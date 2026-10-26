# Dates, identifiers, file paths, population denominators, weekly observations and POI categories.

get_date_range <- function(month_code) {
  parts <- strsplit(month_code, "-")[[1]]
  year <- as.integer(parts[1])
  month <- as.integer(parts[2])
  start_date <- as.Date(sprintf("%d-%02d-01", year, month))
  end_date <- lubridate::ceiling_date(start_date, "month") - 1
  list(start = start_date, end = end_date,
       n_days = as.integer(end_date - start_date) + 1, year = year, month = month)
}

get_wave_date_range <- function(wave_months) {
  ranges <- lapply(wave_months, get_date_range)
  start <- ranges[[1]]$start
  end   <- ranges[[length(ranges)]]$end
  list(start = start, end = end,
       n_days = as.integer(end - start) + 1L,
       month_ranges = ranges)
}

standardize_fips <- function(fips) {
  fips_clean <- trimws(as.character(fips))
  fips_int <- suppressWarnings(as.integer(fips_clean))
  result <- sprintf("%05d", fips_int)
  result[is.na(fips_int)] <- NA_character_
  result
}

standardize_cbg <- function(x) {
  if (inherits(x, "integer64")) {
    x <- format(x, scientific = FALSE)
  } else {
    x <- as.character(x)
  }
  x <- trimws(x)
  sci_idx <- grepl("[eE]", x) & !is.na(x)
  if (any(sci_idx)) {
    x[sci_idx] <- format(as.numeric(x[sci_idx]), scientific = FALSE)
    x[sci_idx] <- trimws(x[sci_idx])
  }
  x <- sub("\\.0*$", "", x)
  short_idx <- !is.na(x) & nchar(x) < 12 & grepl("^[0-9]+$", x)
  if (any(short_idx)) {
    x[short_idx] <- formatC(as.numeric(x[short_idx]), width = 12, format = "d", flag = "0")
  }
  x
}

get_utm_crs <- function(state) {
  if (state == "NY") CONFIG$ny_utm_crs else CONFIG$ut_utm_crs
}

get_metro_fips <- function(state, digits = 5) {
  if (digits == 5) {
    if (state == "NY") CONFIG$nyc_fips_5 else CONFIG$ut_fips_5
  } else {
    if (state == "NY") CONFIG$nyc_fips_3 else CONFIG$ut_fips_3
  }
}

get_cbg_path <- function(state) {
  if (state == "NY") CONFIG$ny_cbg_path else CONFIG$ut_cbg_path
}

get_landuse_path <- function(state, year = NULL) {
  if (state == "NY") CONFIG$ny_landuse_path else file.path(CONFIG$ut_parcels_dir, sprintf("Parcels_LIR_%s.gpkg", year))
}

get_buildings_path <- function(state) {
  if (state == "NY") CONFIG$ny_buildings_path else CONFIG$ut_buildings_path
}

get_state_prefix <- function(state) {
  if (state == "NY") "36" else "49"
}

get_road_paths <- function(state) {
  if (state == "NY") {
    file.path("data", "transport", "newyork_roads2", "tl_2023_36_prisecroads.shp")
  } else {
    list.files(CONFIG$ut_roads_dir, pattern = "\\.shp$", full.names = TRUE)
  }
}

build_poi_area_cache <- function(poi_dt, state, force = FALSE) {
  cache_file <- file.path(CONFIG$cache_dir, sprintf("poi_area_%s.feather", state))

  if (file.exists(cache_file) && !isTRUE(force)) {
    cached <- as.data.table(arrow::read_feather(cache_file))
    if (all(unique(poi_dt$persistent_id) %in% cached$persistent_id))
      return(cached)
  }

  uniq <- unique(poi_dt[!is.na(poi_lon) & !is.na(poi_lat),
                        .(persistent_id, poi_lon, poi_lat)])
  utm_crs <- get_utm_crs(state)
  bld <- sf::st_read(get_buildings_path(state), quiet = TRUE)
  bld <- sf::st_transform(bld, utm_crs)
  bld$.bld_area_sqm <- as.numeric(sf::st_area(bld))
  bld <- bld[, ".bld_area_sqm"]

  pts <- sf::st_as_sf(uniq, coords = c("poi_lon", "poi_lat"), crs = 4326)
  pts <- sf::st_transform(pts, utm_crs)
  joined <- sf::st_join(pts, bld, join = sf::st_within, left = TRUE)

  out <- data.table::data.table(
    persistent_id = joined$persistent_id,
    poi_area_sqm  = joined$.bld_area_sqm)
  out <- unique(out, by = "persistent_id")
  arrow::write_feather(out, cache_file)
  out
}

fit_distance_kernel <- function(distance_m, weights, n_bins = 25) {
  ok <- is.finite(distance_m) & is.finite(weights) &
        distance_m > 0 & weights > 0
  d <- distance_m[ok]; w <- weights[ok]
  if (length(d) < 50) {
    return(list(form = "powerlaw", sigma = 1000, alpha = 2, kappa = NA_real_,
                aic_pl = NA_real_, aic_tpl = NA_real_, aic_gauss = NA_real_,
                n = length(d), note = "insufficient data, using defaults"))
  }
  brks <- exp(seq(log(max(min(d), 50)), log(max(d) + 1), length.out = n_bins + 1))
  bin <- findInterval(d, brks, all.inside = TRUE)
  agg <- tapply(w, bin, sum)
  mid <- (brks[-1] + brks[-length(brks)]) / 2
  y <- as.numeric(agg); x <- mid[as.integer(names(agg))]
  area <- diff(brks)[as.integer(names(agg))]
  nll_pl <- function(par) {
    sg <- exp(par[1]); al <- exp(par[2])
    mu <- (1 + x / sg)^(-al) * area
    mu <- mu * (sum(y) / sum(mu))
    -sum(stats::dpois(round(y), mu, log = TRUE))
  }
  nll_tpl <- function(par) {
    sg <- exp(par[1]); al <- exp(par[2]); kp <- exp(par[3])
    mu <- (1 + x / sg)^(-al) * exp(-x / kp) * area
    mu <- mu * (sum(y) / sum(mu))
    -sum(stats::dpois(round(y), mu, log = TRUE))
  }
  nll_gauss <- function(par) {
    sg <- exp(par[1])
    mu <- exp(-x^2 / (2 * sg^2)) * area
    mu <- mu * (sum(y) / sum(mu))
    -sum(stats::dpois(round(y), mu, log = TRUE))
  }
  fit_pl <- tryCatch(stats::optim(c(log(2000), log(2)), nll_pl,
                                  method = "Nelder-Mead"),
                     error = function(e) NULL)
  fit_tpl <- tryCatch(stats::optim(c(log(2000), log(2), log(20000)), nll_tpl,
                                   method = "Nelder-Mead"),
                      error = function(e) NULL)
  fit_g  <- tryCatch(stats::optim(log(3000), nll_gauss,
                                  method = "Brent",
                                  lower = log(50), upper = log(1e6)),
                     error = function(e) NULL)
  pl_ok  <- !is.null(fit_pl)  && is.finite(fit_pl$value)
  tpl_ok <- !is.null(fit_tpl) && is.finite(fit_tpl$value)
  g_ok   <- !is.null(fit_g)   && is.finite(fit_g$value)
  out <- list(
    sigma       = if (pl_ok)  exp(fit_pl$par[1])  else 1000,
    alpha       = if (pl_ok)  exp(fit_pl$par[2])  else 2,
    sigma_tpl   = if (tpl_ok) exp(fit_tpl$par[1]) else NA_real_,
    alpha_tpl   = if (tpl_ok) exp(fit_tpl$par[2]) else NA_real_,
    kappa       = if (tpl_ok) exp(fit_tpl$par[3]) else NA_real_,
    sigma_gauss = if (g_ok)   exp(fit_g$par)      else NA_real_,
    aic_pl    = if (pl_ok)  2 * 2 + 2 * fit_pl$value  else NA_real_,
    aic_tpl   = if (tpl_ok) 2 * 3 + 2 * fit_tpl$value else NA_real_,
    aic_gauss = if (g_ok)   2 * 1 + 2 * fit_g$value   else NA_real_,
    n = length(d)
  )
  aics <- c(powerlaw = out$aic_pl, tpowerlaw = out$aic_tpl, gauss = out$aic_gauss)
  aics <- aics[is.finite(aics)]
  out$form <- if (length(aics)) names(aics)[which.min(aics)] else "powerlaw"
  out
}

population_for_wave <- function(state, wave_months) {
  yrs <- substr(wave_months, 1, 4)
  w   <- prop.table(table(yrs))

  pop_path <- function(y) file.path(CONFIG$grids_dir,
                                    sprintf("population_%s_%s.csv", state, y))

  available <- list.files(CONFIG$grids_dir,
                          pattern = sprintf("^population_%s_\\d{4}\\.csv$", state))
  avail_yrs <- as.integer(sub(".*_(\\d{4})\\.csv$", "\\1", available))
  if (length(avail_yrs) == 0L)
    stop(sprintf("No population_%s_*.csv files in %s", state, CONFIG$grids_dir))

  resolve_year <- function(y) {
    yi <- as.integer(y)
    if (yi %in% avail_yrs) return(as.character(yi))
    nearest <- avail_yrs[which.min(abs(avail_yrs - yi))]
    warning(sprintf("population_%s_%s.csv missing; using %s", state, y, nearest))
    as.character(nearest)
  }

  parts <- lapply(names(w), function(y) {
    dt <- data.table::fread(pop_path(resolve_year(y)))
    dt[, .(grid_id, population = population * w[[y]])]
  })

  blended <- data.table::rbindlist(parts)[
    , .(population = sum(population)), by = grid_id]
  blended[, E      := pmax(round(population, 2), 0.1)]
  blended[, log_E  := log(E)]
  blended[, population := round(population, 2)]
  blended[]
}

modzcta_population_for_wave <- function(wave_months) {
  yrs <- substr(wave_months, 1, 4)
  w   <- prop.table(table(yrs))
  available <- list.files(CONFIG$grids_dir, pattern = "^population_NY_modzcta_\\d{4}\\.csv$")
  avail_yrs <- as.integer(sub(".*_(\\d{4})\\.csv$", "\\1", available))
  if (!length(avail_yrs)) stop("No population_NY_modzcta_*.csv in ", CONFIG$grids_dir)
  resolve_year <- function(y) { yi <- as.integer(y); if (yi %in% avail_yrs) yi else avail_yrs[which.min(abs(avail_yrs - yi))] }
  parts <- lapply(names(w), function(y) {
    dt <- data.table::fread(file.path(CONFIG$grids_dir, sprintf("population_NY_modzcta_%s.csv", resolve_year(y))))
    dt[, .(grid_id, population = population * w[[y]])]
  })
  blended <- data.table::rbindlist(parts)[, .(pop_zone = sum(population)), by = grid_id]
  zones <- sf::st_drop_geometry(sf::st_read(file.path(CONFIG$grids_dir, "grid_100m_NY_modzcta.gpkg"), quiet = TRUE))
  key <- data.table::as.data.table(zones)[, .(grid_id, geo_id = standardize_fips(rate_geo_id))]
  out <- merge(key, blended, by = "grid_id")[, .(geo_id, pop_zone)]
  out[is.finite(pop_zone) & pop_zone > 0]
}

build_weekly_obs <- function(time_dt, zone_pop, grid_wave_dir) {
  wk <- fread(file.path(grid_wave_dir, "tables", "case_rates_weekly.csv"))
  wk[, `:=`(geo_id = as.character(modzcta), week_ending_date = as.Date(week_ending_date))]
  td <- copy(time_dt)[, date := as.Date(date)]
  td[, week_ending_date := date + (6L - as.integer(format(date, "%w")))]
  wl <- td[, .(week_start = min(date), n_days_in_wave = .N, time_ids = list(time_id)), by = week_ending_date][order(week_ending_date)]
  wl[, week_id := .I]; wl[, observed := n_days_in_wave == 7L]
  td <- merge(td, wl[, .(week_ending_date, week_id, observed)], by = "week_ending_date")
  yw <- merge(wk[week_ending_date %in% wl[observed == TRUE, week_ending_date], .(geo_id, week_ending_date, rate_week = case_rate, suppressed)],
              zone_pop, by = "geo_id")
  yw <- merge(yw, wl[, .(week_ending_date, week_id)], by = "week_ending_date")
  yw[, Y_week := as.integer(round(rate_week * pop_zone / 1e5))]
  list(week_lookup = wl[, .(week_id, week_ending_date, week_start, n_days_in_wave, observed)],
       day_week = td[order(time_id), .(time_id, date, week_id, observed)],
       Y_week = yw[order(week_id, geo_id), .(geo_id, week_id, week_ending_date, rate_week, pop_zone, Y_week, suppressed)])
}

completed_week_rates <- function(case_rates, grid_wave_dir, wave_start) {
  wk <- fread(file.path(grid_wave_dir, "tables", "case_rates_weekly.csv"))
  wk[, `:=`(geo_id = standardize_fips(modzcta), week_ending_date = as.Date(week_ending_date))]
  setkey(wk, geo_id, week_ending_date)
  cr <- copy(case_rates)
  cr[, date := as.Date(wave_start) + sim_day - 1L]
  known <- wk[cr, on = .(geo_id, week_ending_date <= date), mult = "last",
              .(geo_id = i.geo_id, sim_day = i.sim_day, rate_known = x.case_rate)]
  cr <- merge(cr, known, by = c("geo_id", "sim_day"), all.x = TRUE)
  cr[is.na(rate_known), rate_known := 0]
  cr[, case_rate_interp := case_rate]
  cr[, case_rate := rate_known][, c("rate_known", "date") := NULL]
  setkey(cr, geo_id, sim_day)
  cr
}

build_weekly_agg <- function(df, obs, day_week, E_col = "E_inla") {
  df <- data.table::as.data.table(df)
  lat <- data.table::data.table(.col = seq_len(nrow(df)), grid_id = df$grid_id, time_id = df$time_id,
                                E = as.numeric(df[[E_col]]))
  lat <- merge(lat, data.table::as.data.table(day_week)[, .(time_id, week_id)], by = "time_id")
  obs <- data.table::as.data.table(obs)[!is.na(Y_week), .(grid_id, week_id, Y_week)]
  obs[, .row := .I]
  m <- merge(lat, obs, by = c("grid_id", "week_id"))
  missing <- obs[!.row %in% unique(m$.row)]
  if (nrow(missing)) {
    message(sprintf("[build_weekly_agg] %d unit-weeks have no latent rows and are dropped", nrow(missing)))
    obs <- obs[.row %in% unique(m$.row)]; obs[, .row := .I]
    m <- merge(lat, obs, by = c("grid_id", "week_id"))
  }
  m[, E_sum := sum(E), by = .row]
  m[, w := E / E_sum]
  A <- Matrix::sparseMatrix(i = m$.row, j = m$.col, x = m$w, dims = c(nrow(obs), nrow(df)))
  E_obs <- m[, .(E = E_sum[1], n_days = .N), by = .row][order(.row)]
  list(A = A, Y = as.integer(obs$Y_week), E = E_obs$E,
       key = obs[, .(grid_id, week_id, n_days = E_obs$n_days)])
}

redistribute_batch_cases <- function(dt, id_col = "fips", cases_col = "new_cases",
                                     min_gap = 2, threshold_frac = 0.02) {
  setorderv(dt, c(id_col, "date"))
  n_redist <- 0L
  total_moved <- 0

  for (gid in unique(dt[[id_col]])) {
    idx <- which(dt[[id_col]] == gid)
    x   <- dt[[cases_col]][idx]
    n   <- length(x)
    if (n < 3 || sum(x) == 0) next

    med_x  <- median(x[x > 0], na.rm = TRUE)
    thresh <- max(1, threshold_frac * med_x)
    is_low <- x <= thresh

    i <- 1L
    while (i <= n) {
      if (is_low[i]) {
        gap_start <- i
        while (i <= n && is_low[i]) i <- i + 1L
        gap_len <- i - gap_start
        if (gap_len >= min_gap && i <= n && x[i] > thresh) {
          window <- gap_start:i
          window_total <- sum(x[window])
          if (window_total > 0) {
            pre_val  <- if (gap_start > 1) max(x[gap_start - 1], 1) else max(window_total / length(window), 1)
            post_val <- if (i < n) max(x[min(i + 1, n)], 1) else pre_val
            w_len <- length(window)
            weights <- seq(pre_val, post_val, length.out = w_len)
            weights <- weights / sum(weights)
            alloc <- floor(weights * window_total)
            remainder <- as.integer(window_total - sum(alloc))
            if (remainder > 0) {
              frac <- (weights * window_total) - alloc
              top_idx <- order(frac, decreasing = TRUE)[seq_len(remainder)]
              alloc[top_idx] <- alloc[top_idx] + 1L
            }
            x[window] <- as.integer(alloc)
            n_redist <- n_redist + 1L
            total_moved <- total_moved + window_total - alloc[w_len]
          }
          i <- i + 1L
        }
      } else {
        i <- i + 1L
      }
    }
    dt[[cases_col]][idx] <- x
  }
  dt
}

read_inla_graph_as_nb <- function(path) {
  lines <- readLines(path)
  n <- as.integer(lines[1])
  nb <- vector("list", n)
  for (i in seq_len(n)) {
    parts <- as.integer(strsplit(trimws(lines[i + 1]), "\\s+")[[1]])
    nbrs <- if (length(parts) >= 3) parts[3:length(parts)] else integer(0)
    nb[[i]] <- if (length(nbrs) == 0) 0L else nbrs
  }
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_len(n))
  nb
}

HOLIDAY_DATES <- as.Date(c(
  "2020-01-01", "2020-01-20", "2020-02-17", "2020-05-25", "2020-07-04",
  "2020-09-07", "2020-11-26", "2020-12-25",
  "2021-01-01", "2021-01-18", "2021-02-15", "2021-05-31", "2021-07-04",
  "2021-09-06", "2021-11-25", "2021-12-25",
  "2022-01-01", "2022-01-17", "2022-02-21", "2022-05-30", "2022-07-04",
  "2022-09-05", "2022-11-24", "2022-12-25"
))

NAICS_PREFIX_OVERRIDES <- c(
  "5311" = "Commercial & Mixed-Use Buildings",
  "5173" = "Telecoms Retail",
  "5191" = "Libraries & Archives",
  "9261" = "Government Service Offices"
)

NAICS3_CATEGORY_MAP <- c(
  "722" = "Restaurants & Bars", "721" = "Hotels & Lodging",
  "445" = "Grocery & Food Retail", "446" = "Pharmacy & Health Retail",
  "452" = "General Retail", "453" = "General Retail",
  "448" = "Clothing & Accessories",
  "442" = "Home & Furniture Retail", "443" = "Home & Furniture Retail",
  "444" = "Home & Furniture Retail",
  "441" = "Auto Dealers & Repair", "811" = "Auto Dealers & Repair",
  "447" = "Gas Stations",
  "713" = "Fitness & Recreation", "711" = "Fitness & Recreation",
  "812" = "Personal & Consumer Services",
  "611" = "Schools & Education",
  "621" = "Health Care (Ambulatory)", "624" = "Health Care (Ambulatory)",
  "813" = "Religious Organizations",
  "522" = "Banking, Finance & Real Estate",
  "531" = "Banking, Finance & Real Estate",
  "238" = "Industrial & Construction",
  "541" = "Offices & Professional Services",
  "561" = "Offices & Professional Services",
  "712" = "Museums, Arts & Entertainment",
  "512" = "Museums, Arts & Entertainment",
  "488" = "Transportation & Delivery", "492" = "Transportation & Delivery",
  "481" = "Transportation & Delivery"
)

NAICS2_SECTOR_MAP <- c(
  "11" = "Agriculture", "21" = "Mining", "22" = "Utilities",
  "23" = "Construction", "31" = "Manufacturing", "32" = "Manufacturing",
  "33" = "Manufacturing", "42" = "Wholesale Trade", "44" = "Retail Trade",
  "45" = "Retail Trade", "48" = "Transportation", "49" = "Transportation",
  "51" = "Information", "52" = "Finance & Insurance",
  "53" = "Real Estate", "54" = "Professional Services",
  "55" = "Management", "56" = "Administrative Services",
  "61" = "Education", "62" = "Health Care", "71" = "Arts & Recreation",
  "72" = "Accommodation & Food", "81" = "Other Services",
  "92" = "Public Administration"
)

naics_lookup <- function(naics_code, overrides, override_key, map3, default) {
  v  <- unname(overrides[override_key])
  v3 <- unname(map3[substr(as.character(naics_code), 1, 3)])
  v  <- ifelse(is.na(v), v3, v)
  ifelse(is.na(v), default, v)
}

naics_category_label <- function(naics_code) {
  code <- as.character(naics_code)
  sector <- unname(NAICS2_SECTOR_MAP[substr(code, 1, 2)])
  fallback <- ifelse(is.na(sector), "Other", paste0("Other (", sector, ")"))
  naics_lookup(code, NAICS_PREFIX_OVERRIDES, substr(code, 1, 4), NAICS3_CATEGORY_MAP, fallback)
}

BROAD_CATEGORY_MAP <- c(
  "Restaurants & Bars"               = "Food & Beverage",
  "Grocery & Food Retail"            = "Retail & Trade",
  "General Retail"                   = "Retail & Trade",
  "Clothing & Accessories"           = "Retail & Trade",
  "Home & Furniture Retail"          = "Retail & Trade",
  "Other (Retail Trade)"             = "Retail & Trade",
  "Other (Wholesale Trade)"          = "Retail & Trade",
  "Pharmacy & Health Retail"         = "Essential & Services",
  "Health Care (Ambulatory)"         = "Essential & Services",
  "Personal & Consumer Services"     = "Essential & Services",
  "Banking, Finance & Real Estate"   = "Essential & Services",
  "Offices & Professional Services"  = "Essential & Services",
  "Auto Dealers & Repair"            = "Essential & Services",
  "Gas Stations"                     = "Essential & Services",
  "Other (Finance & Insurance)"      = "Essential & Services",
  "Other (Real Estate)"              = "Essential & Services",
  "Other (Administrative Services)"  = "Essential & Services",
  "Other (Utilities)"                = "Essential & Services",
  "Hotels & Lodging"                 = "Recreation & Entertainment",
  "Fitness & Recreation"             = "Recreation & Entertainment",
  "Museums, Arts & Entertainment"    = "Recreation & Entertainment",
  "Schools & Education"              = "Education",
  "Industrial & Construction"        = "Industrial & Construction",
  "Other (Manufacturing)"            = "Industrial & Construction",
  "Other (Construction)"             = "Industrial & Construction",
  "Other (Mining)"                   = "Industrial & Construction",
  "Other (Agriculture)"              = "Industrial & Construction",
  "Transportation & Delivery"        = "Transportation",
  "Other (Transportation)"           = "Transportation",
  "Religious Organizations"          = "Religious Organizations",
  "Commercial & Mixed-Use Buildings" = "Commercial & Mixed-Use Buildings",
  "Telecoms Retail"                  = "Retail & Trade",
  "Libraries & Archives"             = "Education",
  "Government Service Offices"       = "Essential & Services"
)

naics_category_broad <- function(fine_label) {
  b <- unname(BROAD_CATEGORY_MAP[fine_label])
  ifelse(is.na(b), "Other", b)
}

NAICS_INDOOR_OVERRIDES <- c(
  "712190" = 0.5
)

NAICS3_INDOOR_MOD <- c(
  "238" = 0.5,
  "488" = 0.5,
  "492" = 0.5,
  "713" = 0.7,
  "712" = 0.8,
  "811" = 0.7,
  "441" = 0.7,
  "447" = 0.6,
  "481" = 0.7
)

naics_indoor_mod <- function(naics_code)
  naics_lookup(naics_code, NAICS_INDOOR_OVERRIDES, as.character(naics_code), NAICS3_INDOOR_MOD, 1.0)
