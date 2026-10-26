# Case data: NYC weekly MODZCTA rates (suppressed weeks imputed from the borough rate, capped at
# 4 cases) interpolated to days, and Utah daily county counts from the NYT cumulative series.
source("R/setup.R")
library(lubridate)

runs <- select_runs()
NY_WAVES <- runs[state == "NY"]
UT_WAVES <- runs[state == "UT"]

if (nrow(NY_WAVES)) {
  nyc_raw <- fread(CONFIG$nyc_caserate_path)
  modzcta_cols <- grep("^CASERATE_[0-9]+$", names(nyc_raw), value = TRUE)
  nyc_long <- melt(nyc_raw, id.vars = "week_ending", measure.vars = modzcta_cols,
                   variable.name = "modzcta_col", value.name = "case_rate")
  nyc_long[, modzcta := gsub("CASERATE_", "", modzcta_col)][, modzcta_col := NULL]
  nyc_long[, week_ending_date := mdy(week_ending)]
  nyc_long <- nyc_long[!is.na(week_ending_date)]

  cbg_modzcta <- fread(CONFIG$nyc_crosswalk_path)
  cbg_modzcta[, modzcta := as.character(modzcta)]
  saveRDS(cbg_modzcta, file.path(CONFIG$cache_dir, "cbg_modzcta_crosswalk_NY.rds"))
  modzcta_borough <- cbg_modzcta[, .(borough = borough[1]), by = modzcta]

  borough_code <- c(Bronx = "CASERATE_BX", Brooklyn = "CASERATE_BK", Manhattan = "CASERATE_MN",
                    Queens = "CASERATE_QN", `Staten Island` = "CASERATE_SI")
  borough_long <- melt(nyc_raw, id.vars = "week_ending", measure.vars = unname(borough_code),
                       variable.name = "borough_col", value.name = "borough_rate")
  borough_long[, borough := names(borough_code)[match(as.character(borough_col), borough_code)]]
  borough_long[, week_ending_date := mdy(week_ending)]
  borough_long <- borough_long[!is.na(week_ending_date), .(borough, week_ending_date, borough_rate)]

  modzcta_pop <- fread(CONFIG$ny_modzcta_pop_path)[, .(modzcta = as.character(MODIFIED_ZCTA), pop = as.numeric(POP_DENOMINATOR))]
  nyc_long[modzcta_borough, borough := i.borough, on = "modzcta"]
  nyc_long[borough_long, borough_rate := i.borough_rate, on = c("borough", "week_ending_date")]
  nyc_long[modzcta_pop, pop := i.pop, on = "modzcta"]
  nyc_long[, suppressed := is.na(case_rate)]
  nyc_long[suppressed == TRUE & !is.na(borough_rate) & !is.na(pop), case_rate := pmin(borough_rate, 4 / pop * 1e5)]
  nyc_long[, c("borough", "borough_rate", "pop") := NULL]
  nyc_long <- nyc_long[!is.na(case_rate)]

  process_nyc_wave <- function(wave_id, wave_months) {
    wdr <- get_wave_date_range(wave_months)
    out_dir <- get_wave_dir("NY", wave_id)
    dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

    fwrite(nyc_long[week_ending_date >= (wdr$start - 28) & week_ending_date <= (wdr$end + 6),
                    .(modzcta, week_ending_date, case_rate, suppressed)],
           file.path(out_dir, "tables", "case_rates_weekly.csv"))

    relevant <- nyc_long[week_ending_date >= (wdr$start - 6) & week_ending_date <= (wdr$end + 6)]
    midpoints <- relevant[, .(case_rate = mean(case_rate, na.rm = TRUE)), by = .(modzcta, mid_date = week_ending_date - 3L)]
    setorder(midpoints, modzcta, mid_date)

    all_dates <- seq(wdr$start, wdr$end, by = "day")
    daily <- rbindlist(lapply(unique(midpoints$modzcta), function(mz) {
      mz_data <- midpoints[modzcta == mz]
      if (nrow(mz_data) < 2L) return(data.table(modzcta = mz, date = all_dates, case_rate = mz_data$case_rate[1]))
      interp <- approx(as.numeric(mz_data$mid_date), mz_data$case_rate, xout = as.numeric(all_dates), rule = 2)
      data.table(modzcta = mz, date = all_dates, case_rate = pmax(interp$y, 0))
    }))
    daily[, sim_day := as.integer(date - wdr$start) + 1L]
    daily[, source_month := format(date, "%Y-%m")]
    daily[, borough := modzcta_borough[.SD, borough, on = .(modzcta)]]
    daily[is.na(borough), borough := "Unknown"]
    daily[daily[, .(borough_mean = mean(case_rate, na.rm = TRUE)), by = .(date, borough)],
          borough_mean := i.borough_mean, on = .(date, borough)]
    alpha <- CONFIG$spillover_alpha
    daily[, case_rate_spillover := alpha * case_rate + (1 - alpha) * borough_mean]
    daily[is.na(case_rate_spillover), case_rate_spillover := case_rate]
    setorder(daily, date)

    fwrite(daily[, .(geo_id = modzcta, date, sim_day, source_month, case_rate, case_rate_spillover, borough)],
           file.path(out_dir, "tables", "case_rates.csv"))
  }
  for (i in seq_len(nrow(NY_WAVES))) process_nyc_wave(NY_WAVES$wave_id[i], NY_WAVES$wave_months[[i]])
}

if (nrow(UT_WAVES)) {
  ut_years <- unique(unlist(lapply(UT_WAVES$wave_months, substr, 1, 4)))
  nyt_raw <- rbindlist(lapply(ut_years, function(y) fread(CONFIG$nyt_counties_files[[y]])))
  ut_raw <- nyt_raw[state == "Utah" & !is.na(fips)]
  ut_raw[, fips := sprintf("%05d", as.integer(fips))]
  ut_raw[, date := as.Date(date)]
  setorder(ut_raw, fips, date)
  ut_raw[, new_cases := cases - data.table::shift(cases, 1, fill = 0), by = fips]
  ut_raw[new_cases < 0, new_cases := 0]
  # From October 2020 Utah reported some counties in batches; spread batches over the preceding gap
  ut_raw <- redistribute_batch_cases(ut_raw[date >= as.Date("2020-10-01")], id_col = "fips", cases_col = "new_cases") |>
    rbind(ut_raw[date < as.Date("2020-10-01")]) |>
    setorder(fips, date)

  process_ut_wave <- function(wave_id, wave_months, acs_year) {
    wdr <- get_wave_date_range(wave_months)
    out_dir <- get_wave_dir("UT", wave_id)
    dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

    daily <- ut_raw[date >= wdr$start & date <= wdr$end]
    gsf <- st_read(file.path(CONFIG$grids_dir, "grid_100m_UT.gpkg"), quiet = TRUE)
    gdt <- merge(data.table(grid_id = gsf$grid_id, fips_5 = gsf$fips_5),
                 fread(file.path(CONFIG$grids_dir, sprintf("population_UT_%s.csv", acs_year)))[, .(grid_id, population)],
                 by = "grid_id", all.x = TRUE)
    county_pop <- gdt[, .(population = sum(population, na.rm = TRUE)), by = fips_5]

    daily <- merge(daily, county_pop, by.x = "fips", by.y = "fips_5", all.x = TRUE)
    daily[is.na(population) | population < 100, population := 100]
    daily[, case_rate := new_cases / population * 1e5]
    daily[, sim_day := as.integer(date - wdr$start) + 1L]
    daily[, source_month := format(date, "%Y-%m")]
    daily[, region := county]
    daily <- merge(daily, daily[, .(state_mean = mean(case_rate, na.rm = TRUE)), by = date], by = "date", all.x = TRUE)
    alpha <- CONFIG$spillover_alpha
    daily[, case_rate_spillover := alpha * case_rate + (1 - alpha) * state_mean]
    daily[is.na(case_rate_spillover), case_rate_spillover := case_rate]

    fwrite(daily[, .(geo_id = fips, date, sim_day, source_month, new_cases, case_rate, case_rate_spillover, region)],
           file.path(out_dir, "tables", "case_rates.csv"))
  }
  for (i in seq_len(nrow(UT_WAVES))) process_ut_wave(UT_WAVES$wave_id[i], UT_WAVES$wave_months[[i]], UT_WAVES$acs_year[i])
}
