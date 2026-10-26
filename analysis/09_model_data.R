# Unit-day model data: offsets, covariates and counts per unit (100 m cell in NYC, county in Utah) and day,
# the augmented unit graph, and, for NYC, the weekly observations and grid-level pooled prevalence.
source("R/setup.R")

GRID_COVS <- c("ses_contrast", "log_repeat_ratio", "log_contact_density_hourly", "log_daytime_pop",
               "pct_older", "pct_younger", "pct_employed", "pct_unemployed", "building_density")
ACTIVITY_COLS <- c("exposure_import", "pressure_internal", "pressure_external", "internal_activity", "external_activity",
                   "total_activity", "total_visitors")

build_model_data <- function(wave_id, wave_months) {
  out_dir <- get_wave_dir(state, wave_id)
  dir.create(file.path(out_dir, "models"), recursive = TRUE, showWarnings = FALSE)
  dr <- get_wave_date_range(wave_months)
  n_days <- dr$n_days
  gpkg <- file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state))

  # units: visit covariates (0 where a cell has no visits), area covariates and reporting area
  grid_dt <- merge(as.data.table(st_drop_geometry(st_read(gpkg, quiet = TRUE)))[, .(grid_id)],
                   as.data.table(read_feather(file.path(out_dir, "intermediate", "spatial_covariates.feather"))),
                   by = "grid_id", all.x = TRUE)
  for (col in setdiff(names(grid_dt), "grid_id")) grid_dt[is.na(get(col)), (col) := 0]
  demo <- fread(file.path(out_dir, "tables", "demographics.csv"))
  grid_dt <- merge(grid_dt, demo[, .(grid_id, pct_older, pct_younger, pct_employed, pct_unemployed,
                                     building_density, landuse_commercial, landuse_cat)], by = "grid_id", all.x = TRUE)
  geo <- as.data.table(st_drop_geometry(st_read(gpkg, quiet = TRUE)))
  if (state == "NY") {
    xw <- fread(CONFIG$nyc_crosswalk_path)[, .(cbg_geoid = as.character(cbg_geoid), modzcta = as.character(modzcta))]
    geo <- merge(geo[, .(grid_id, cbg_geoid = as.character(cbg_geoid))], xw[!duplicated(cbg_geoid)], by = "cbg_geoid", all.x = TRUE)
    geo[, geo_id := modzcta]
  } else geo[, geo_id := as.character(fips_5)]
  grid_dt <- merge(grid_dt, geo[, .(grid_id, geo_id)], by = "grid_id", all.x = TRUE)[!is.na(geo_id)]

  # population offset; NYC excludes vacant, public and industrial land, and both states drop units with E < 1
  grid_dt <- merge(grid_dt, population_for_wave(state, wave_months)[, .(grid_id, E)], by = "grid_id", all.x = TRUE)
  grid_dt[, E := pmax(fifelse(is.na(E), 0.1, E), 0.1)]
  grid_dt[, log_E := log(E)]
  if (state == "NY") grid_dt <- grid_dt[!landuse_cat %in% c("vacant", "public", "industrial")]
  grid_dt <- grid_dt[E >= 1]
  n_grids <- nrow(grid_dt)
  grid_dt[, grid_idx := .I]

  dates <- dr$start + 0:(n_days - 1)
  time_dt <- data.table(time_id = 1:n_days, sim_day = 1:n_days, date = dates,
                        is_weekend = as.integer(format(dates, "%u") %in% c("6", "7")),
                        is_holiday = as.integer(dates %in% HOLIDAY_DATES))

  gt_dt <- CJ(grid_idx = 1:n_grids, time_id = 1:n_days)
  gt_dt <- merge(gt_dt, grid_dt[, c("grid_id", "grid_idx", "geo_id", GRID_COVS, "landuse_commercial", "E", "log_E"), with = FALSE],
                 by = "grid_idx", all.x = TRUE)
  gt_dt <- merge(gt_dt, time_dt[, .(time_id, sim_day, date, is_weekend, is_holiday)], by = "time_id", all.x = TRUE)

  rho_dt <- fread(file.path(out_dir, "tables", "od_kernel_rho.csv"))
  gt_dt <- merge(gt_dt, rho_dt[, c("grid_id", "sim_day", "rho_norm", "log_rho", "log_exposure_import", ACTIVITY_COLS), with = FALSE],
                 by = c("grid_id", "sim_day"), all.x = TRUE)
  gt_dt[is.na(rho_norm), rho_norm := 1]
  for (col in c("log_rho", "log_exposure_import", ACTIVITY_COLS)) gt_dt[is.na(get(col)), (col) := 0]
  gt_dt[, commercial_x_weekend := landuse_commercial * is_weekend]

  # counts: Utah daily county counts; NYC daily zone counts allocated to cells by population share
  # (the NYC likelihood uses the weekly observations below)
  case_dt <- fread(file.path(out_dir, "tables", "case_rates.csv"))[, geo_id := as.character(geo_id)]
  gt_dt <- merge(gt_dt, case_dt[, intersect(c("geo_id", "sim_day", "case_rate", "case_rate_spillover", "new_cases"), names(case_dt)), with = FALSE],
                 by = c("geo_id", "sim_day"), all.x = TRUE)
  gt_dt[is.na(case_rate), case_rate := 0]
  gt_dt[is.na(case_rate_spillover), case_rate_spillover := 0]
  set.seed(42)
  if (state == "NY") {
    gt_dt[, daily_rate := case_rate / 7]
    gt_dt <- merge(gt_dt, modzcta_population_for_wave(wave_months), by = "geo_id", all.x = TRUE)
    gt_dt[is.na(pop_zone), pop_zone := sum(E), by = .(geo_id, sim_day)]
    gt_dt[, zone_N := as.integer(round(daily_rate[1] * pop_zone[1] / 1e5)), by = .(geo_id, sim_day)]
    gt_dt[, prob_g := E / sum(E), by = .(geo_id, sim_day)]
    gt_dt[is.na(prob_g) | !is.finite(prob_g), prob_g := 0]
    gt_dt[, prob_g := prob_g / sum(prob_g), by = .(geo_id, sim_day)]
    gt_dt[, Y := if (zone_N[1] > 0 && sum(prob_g) > 0) as.integer(rmultinom(1, size = zone_N[1], prob = prob_g)) else 0L,
          by = .(geo_id, sim_day)]
    gt_dt[, c("zone_N", "prob_g", "pop_zone") := NULL]
  } else {
    gt_dt[is.na(new_cases), new_cases := 0L]
    gt_dt[, `:=`(daily_rate = case_rate, Y = as.integer(new_cases))]
  }

  setorder(gt_dt, grid_id, time_id)
  gt_dt[, local_inf_pressure := pressure_internal]
  gt_dt[, log_local_inf_pressure := log1p(local_inf_pressure)]
  gt_dt <- add_contact_prevalence(gt_dt, force = TRUE)

  rho_mean <- mean(gt_dt$rho_norm, na.rm = TRUE)
  if (abs(rho_mean - 1) > 0.05) gt_dt[, rho_norm := rho_norm / rho_mean]
  gt_dt[, rho_norm := pmax(pmin(rho_norm, exp(5)), exp(-5))]
  gt_dt[, `:=`(log_rho = log(rho_norm), offset_total = log_E, idx_spatial = grid_idx, idx_temporal = time_id)]
  gt_dt[, idx_interaction := .I]

  zs <- function(x) if (sd(x, na.rm = TRUE) > 0) scale(x)[, 1] else 0
  for (col in GRID_COVS) {
    gt_dt[is.na(get(col)), (col) := 0]
    gt_dt[, (paste0(col, "_std")) := zs(get(col))]
  }
  gt_dt[, `:=`(local_inf_pressure_std = zs(log_local_inf_pressure), exposure_import_std = zs(log_exposure_import))]

  # NYC: weekly MODZCTA counts allocated to cells by population share (the observation model)
  weekly_obs <- NULL
  if (state == "NY") {
    weekly_obs <- build_weekly_obs(time_dt, modzcta_population_for_wave(wave_months), out_dir)
    cell_E <- gt_dt[, .(E = E[1], n_E = uniqueN(E)), by = .(grid_id, geo_id)]
    stopifnot(all(cell_E$n_E == 1L))
    set.seed(42)
    yw <- merge(cell_E[, .(grid_id, geo_id, E)], weekly_obs$Y_week[, .(geo_id, week_id, Y_week)], by = "geo_id", allow.cartesian = TRUE)
    yw[, Y_cell := if (Y_week[1] > 0 && sum(E) > 0) as.integer(rmultinom(1, Y_week[1], E)) else 0L, by = .(geo_id, week_id)]
    stopifnot(identical(yw[, sum(Y_cell), by = .(geo_id, week_id)]$V1, yw[, Y_week[1], by = .(geo_id, week_id)]$V1))
    weekly_obs$Y_week_cell <- yw[order(week_id, grid_id), .(grid_id, geo_id, week_id, Y_week = Y_cell)]
  }

  # augmented unit graph; NYC prevalence is then pooled over it
  grid_sf <- st_read(gpkg, quiet = TRUE)
  grid_sf <- grid_sf[grid_sf$grid_id %in% unique(gt_dt$grid_id), ]
  nb <- write_unit_graph(grid_sf, gt_dt, file.path(out_dir, "models", "grid_adj_unified.graph"),
                         file.path(out_dir, "tables", "adjacency_crossings.csv"))
  if (state == "NY") {
    xy <- st_coordinates(st_centroid(st_geometry(grid_sf)))
    stopifnot(identical(as.integer(grid_sf$grid_id), unique(gt_dt[order(idx_spatial), grid_id])))
    gt_dt <- pool_prevalence(gt_dt, nb, xy)
  }

  saveRDS(list(gt_dt = gt_dt, n_grids = n_grids, n_days = n_days, time_lookup = time_dt, state = state,
               wave_id = wave_id, wave_months = wave_months, date_range = dr, weekly_obs = weekly_obs),
          file.path(out_dir, "models", "model_data_gt.rds"))
}

for (st in unique(select_runs()$state)) {
  set_transport_globals(st)
  run_waves_by_state(build_model_data, select_runs()[state == st])
}
