# Infection pressure: visits weighted by the recent case rate of the visitor's home area (origin MODZCTA
# or county) and by a contact weight (visits, dwell, indoor amplitude, co-visitation, crowding), split into
# visitors from inside and outside the destination's reporting area, per grid cell (or county) and day.
# Also fits candidate trip-length kernels (Supplementary Table S2) and records category and origin breakdowns.
source("R/setup.R")

# Incubation-period weights over lags of 1-7 days (Utah daily origin rates)
INCUBATION_WEIGHTS <- c(`7` = 0.30, `6` = 0.20, `5` = 0.20, `4` = 0.10, `3` = 0.10, `2` = 0.05, `1` = 0.05)

load_cbg_modzcta_xw <- function(cfg) {
  rds <- file.path(cfg$cache_dir, "cbg_modzcta_crosswalk_NY.rds")
  if (file.exists(rds)) return(readRDS(rds))
  xw <- fread(cfg$nyc_crosswalk_path)
  xw[, cbg_geoid := as.character(cbg_geoid)]
  xw[, modzcta := as.character(modzcta)]
  dir.create(dirname(rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(xw, rds)
  xw
}

build_cbg_lookup <- function(cfg, state) {
  if (state == "NY") {
    xw <- load_cbg_modzcta_xw(cfg)
    unique(data.table(cbg_geoid = standardize_cbg(xw$cbg_geoid),
                      rate_geo_id = standardize_fips(xw$modzcta)), by = "cbg_geoid")
  } else {
    cbg_sf <- st_read(get_cbg_path(state), quiet = TRUE)
    cbg_sf <- cbg_sf[substr(cbg_sf$GEOID, 1, 2) == cfg$ut_state_fips, ]
    data.table(cbg_geoid = standardize_cbg(cbg_sf$GEOID),
               rate_geo_id = substr(standardize_cbg(cbg_sf$GEOID), 1, 5))
  }
}

process_od_kernel <- function(wave_id, wave_months) {
  out_dir <- get_wave_dir(state, wave_id)
  rho_file <- file.path(out_dir, "tables", "od_kernel_rho.csv")
  dr <- get_wave_date_range(wave_months)
  n_days <- dr$n_days
  utm <- get_utm_crs(state)

  case_rates <- fread(file.path(out_dir, "tables", "case_rates.csv"))
  case_rates[, geo_id := standardize_fips(geo_id)]
  if (state == "NY") {
    case_rates <- completed_week_rates(case_rates, out_dir, dr$start)
    INCUBATION_WEIGHTS <- c(`1` = 1)
  }
  setkey(case_rates, geo_id, sim_day)

  cbg_lk <- build_cbg_lookup(CONFIG, state)

  grid_file <- file.path(out_dir, "intermediate", "grid_visits.feather")
  grid_dt <- as.data.table(arrow::read_feather(grid_file))
  if ("rate_geo_id" %in% names(grid_dt)) grid_dt[, rate_geo_id := NULL]
  grid_dt <- unique(grid_dt, by = "grid_id")[visits_total > 0]

  # reporting area of each cell: the county in Utah, the MODZCTA of the cell's CBG in NYC
  gsf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), quiet = TRUE)
  gsf <- gsf[gsf$grid_id %in% grid_dt$grid_id, ]
  gm <- if (state == "UT") data.table(grid_id = gsf$grid_id, rate_geo_id = gsf$fips_5) else
    merge(data.table(grid_id = gsf$grid_id, cbg_geoid = standardize_cbg(gsf$cbg_geoid)),
          unique(cbg_lk, by = "cbg_geoid"), by = "cbg_geoid", all.x = TRUE)[, .(grid_id, rate_geo_id)]
  grid_dt <- merge(grid_dt, unique(gm, by = "grid_id"), by = "grid_id", all.x = TRUE)
  grid_dt[, rate_geo_id := standardize_fips(rate_geo_id)]

  od_file <- file.path(out_dir, "intermediate", "od_home_cbg.feather")
  od_dt <- as.data.table(arrow::read_feather(od_file))
  od_dt[, origin_cbg := standardize_cbg(home_cbg)]
  od_dt[, dest_cbg   := standardize_cbg(poi_cbg)]
  cbg_lk[, cbg_geoid := standardize_cbg(cbg_geoid)]

  od_dt <- merge(od_dt, cbg_lk, by.x = "origin_cbg", by.y = "cbg_geoid", all.x = TRUE)
  setnames(od_dt, "rate_geo_id", "origin_geo_id")
  od_dt <- merge(od_dt, cbg_lk, by.x = "dest_cbg", by.y = "cbg_geoid", all.x = TRUE)
  setnames(od_dt, "rate_geo_id", "dest_geo_id")

  gsf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), quiet = TRUE) |> st_transform(utm)
  osf <- st_as_sf(od_dt[!is.na(poi_lon) & !is.na(poi_lat)], coords = c("poi_lon", "poi_lat"), crs = 4326) |> st_transform(utm)
  poi_xy <- st_coordinates(osf)
  osf$poi_utm_x <- poi_xy[,1]; osf$poi_utm_y <- poi_xy[,2]
  osf <- st_join(osf, gsf[, "grid_id"], left = FALSE)
  od <- as.data.table(st_drop_geometry(osf))
  od <- od[!is.na(origin_geo_id) & !is.na(dest_geo_id)]

  od[, dwell_hours := pmin(median_dwell / 60, 8)]
  od[is.na(dwell_hours) | dwell_hours <= 0, dwell_hours := 0.5]
  od[, activity_weight := visitors * (1 + dwell_hours)]

  od[, activity_weight := activity_weight * naics_indoor_mod(naics_code)]

  cbg_cents <- st_read(get_cbg_path(state), quiet = TRUE) |> st_transform(utm)
  cbg_cents <- cbg_cents[substr(cbg_cents$GEOID, 3, 5) %in% get_metro_fips(state, 3), ]
  cc <- st_coordinates(st_centroid(cbg_cents))
  cbg_xy <- data.table(origin_cbg = standardize_cbg(cbg_cents$GEOID), home_X = cc[,1], home_Y = cc[,2])
  od <- merge(od, cbg_xy, by = "origin_cbg", all.x = TRUE)
  od[, od_distance_m := sqrt((poi_utm_x - home_X)^2 + (poi_utm_y - home_Y)^2)]

  k_fit <- fit_distance_kernel(od$od_distance_m, od$visitors)
  median_trip_m <- {
    ok <- is.finite(od$od_distance_m) & od$od_distance_m > 0 & is.finite(od$visitors) & od$visitors > 0
    d_ok <- od$od_distance_m[ok]; w_ok <- od$visitors[ok]
    ord <- order(d_ok)
    d_ok[ord][which(cumsum(w_ok[ord]) / sum(w_ok) >= 0.5)[1]]
  }
  if (k_fit$form == "tpowerlaw" && is.finite(k_fit$kappa)) {
    sigma_hat <- k_fit$sigma_tpl; alpha_hat <- k_fit$alpha_tpl; kappa_hat <- k_fit$kappa
  } else {
    sigma_hat <- k_fit$sigma;     alpha_hat <- k_fit$alpha;     kappa_hat <- NA_real_
  }
  fwrite(data.table(state = state, wave_id = wave_id,
                    sigma_m = sigma_hat, alpha = alpha_hat, kappa_m = kappa_hat,
                    sigma_gauss_m = k_fit$sigma_gauss,
                    aic_powerlaw = k_fit$aic_pl,
                    aic_tpowerlaw = k_fit$aic_tpl,
                    aic_gauss = k_fit$aic_gauss,
                    chosen_form = k_fit$form, n_pairs = k_fit$n,
                    median_trip_m = median_trip_m),
         file.path(out_dir, "tables", "od_kernel_fit.csv"))

  od[, c("home_X", "home_Y") := NULL]

  # co-visitation: brands whose visitors go on to many other brands the same day or week
  chains <- as.data.table(arrow::read_feather(file.path(out_dir, "intermediate", "trip_chains.feather")))
  if (nrow(chains)) {
    bcs_day  <- chains[scope == "day",  .(chain_deg_day  = sum(covisit_count)), by = source_brand]
    bcs_week <- chains[scope == "week", .(chain_deg_week = sum(covisit_count)), by = source_brand]
    bcs <- merge(bcs_day, bcs_week, by = "source_brand", all = TRUE)
    bcs[is.na(chain_deg_day),  chain_deg_day  := 0]
    bcs[is.na(chain_deg_week), chain_deg_week := 0]
    bcs[, boost_day  := chain_deg_day  / pmax(max(chain_deg_day),  1)]
    bcs[, boost_week := chain_deg_week / pmax(max(chain_deg_week), 1)]
    od <- merge(od, bcs[, .(source_brand, boost_day, boost_week)],
                by.x = "brand", by.y = "source_brand", all.x = TRUE)
    od[is.na(boost_day),  boost_day  := 0]
    od[is.na(boost_week), boost_week := 0]
    od[, activity_weight := activity_weight * (1 + 0.5 * boost_day + 0.25 * boost_week)]
    od[, c("boost_day", "boost_week") := NULL]
  }

  # crowding: cell visit density relative to the median cell, capped at 5x
  vd_grid <- grid_dt[, .(grid_id, avg_visit_density)]
  vd_med <- median(vd_grid$avg_visit_density[vd_grid$avg_visit_density > 0], na.rm = TRUE)
  stopifnot(is.finite(vd_med), vd_med > 0)
  vd_grid[is.na(avg_visit_density) | avg_visit_density <= 0, avg_visit_density := vd_med]
  vd_grid[, crowding_factor := sqrt(pmin(avg_visit_density / vd_med, 5))]
  od <- merge(od, vd_grid[, .(grid_id, crowding_factor)], by = "grid_id", all.x = TRUE)
  od[is.na(crowding_factor), crowding_factor := 1]
  od[, activity_weight := activity_weight * crowding_factor]
  od[, crowding_factor := NULL]

  od[, category_label := naics_category_label(naics_code)]
  od[, category_broad := naics_category_broad(category_label)]

  # day-to-day variation: each cell's daily visits relative to its wave mean, clamped to [0.1, 10]
  gdv <- as.data.table(arrow::read_feather(file.path(out_dir, "intermediate", "grid_daily_visits.feather")))[sim_day >= 1 & sim_day <= n_days]
  gdv[, m := visits / mean(visits), by = grid_id]
  gdv <- gdv[is.finite(m)]
  gdv[, m := pmin(pmax(m, 0.1), 10)]
  m_grids <- sort(unique(gdv$grid_id))
  Mmat <- matrix(1, nrow = length(m_grids), ncol = n_days, dimnames = list(as.character(m_grids), NULL))
  Mmat[cbind(match(gdv$grid_id, m_grids), gdv$sim_day)] <- gdv$m

  # week-to-week origin composition: weekly origin share of each cell's visitors relative to its wave share
  odw_path <- file.path(out_dir, "intermediate", "od_home_cbg_weekly.feather")
  ratio_dt <- NULL
  if (file.exists(odw_path)) {
    odw <- as.data.table(arrow::read_feather(odw_path))
    odw[, origin_cbg := standardize_cbg(home_cbg)]
    odw <- merge(odw, cbg_lk[, .(cbg_geoid, origin_geo_id = rate_geo_id)], by.x = "origin_cbg", by.y = "cbg_geoid", all.x = TRUE)
    osw <- st_as_sf(odw[!is.na(poi_lon) & !is.na(poi_lat)], coords = c("poi_lon", "poi_lat"), crs = 4326) |> st_transform(utm)
    osw <- st_join(osw, gsf[, "grid_id"], left = FALSE)
    odw <- as.data.table(st_drop_geometry(osw))[!is.na(origin_geo_id)]
    odw[, week_start := as.Date(week_start)]
    wk <- odw[, .(v = sum(visitors)), by = .(grid_id, origin_geo_id, week_start)]
    wk[, share_w := v / sum(v), by = .(grid_id, week_start)]
    wave_share <- odw[, .(v = sum(visitors)), by = .(grid_id, origin_geo_id)][, share := v / sum(v), by = grid_id][, .(grid_id, origin_geo_id, share)]
    ratio_dt <- merge(wk, wave_share, by = c("grid_id", "origin_geo_id"))[, .(grid_id, origin_geo_id, week_start, ratio = share_w / share)]
    ratio_dt[!is.finite(ratio), ratio := 0]
    setkey(ratio_dt, grid_id, origin_geo_id, week_start)
    rm(odw, osw, wk, wave_share); gc(verbose = FALSE)
  }
  day_dates <- dr$start + (seq_len(n_days) - 1L)
  day_week  <- day_dates - (as.integer(format(day_dates, "%u")) - 1L)
  rho_list <- vector("list", n_days)
  rho_cat_list <- vector("list", n_days)
  od_origin_list <- vector("list", n_days); origin_prev_list <- vector("list", n_days)
  for (day in seq_len(n_days)) {
    lag_days <- as.integer(names(INCUBATION_WEIGHTS))
    lag_rates <- rbindlist(lapply(seq_along(lag_days), function(k) {
      sd <- day - lag_days[k]
      if (sd < 1) sd <- 1L
      cr <- case_rates[sim_day == sd, .(case_rate = mean(case_rate)), by = .(origin_geo_id = geo_id)]
      cr[, w := INCUBATION_WEIGHTS[k]]
      cr
    }))
    oday <- lag_rates[, .(origin_rate = sum(case_rate * w) / sum(w)), by = origin_geo_id]

    dod <- copy(od)
    if (!is.null(ratio_dt)) {
      rr <- ratio_dt[.(dod$grid_id, dod$origin_geo_id, rep(day_week[day], nrow(dod))), ratio]
      rr[is.na(rr)] <- 0
      dod[, activity_weight := activity_weight * rr]
    }
    mi <- match(as.character(dod$grid_id), rownames(Mmat))
    dod[, activity_weight := activity_weight * ifelse(is.na(mi), 1, Mmat[cbind(mi, rep(day, length(mi)))])]
    dod <- merge(dod, oday, by = "origin_geo_id", all.x = TRUE)
    dod[is.na(origin_rate), origin_rate := 0]
    dod[, is_internal := !is.na(origin_geo_id) & !is.na(dest_geo_id) & origin_geo_id == dest_geo_id]

    dod[, contrib := origin_rate * activity_weight]
    od_origin_list[[day]] <- dod[, .(activity = sum(activity_weight, na.rm = TRUE)),
                                 by = .(dest_grid_id = grid_id, origin_geo_id)][, sim_day := day]
    origin_prev_list[[day]] <- oday[, .(origin_geo_id, sim_day = day, origin_rate)]
    rho_list[[day]] <- dod[, .(
      pressure_internal = sum(contrib[is_internal], na.rm = TRUE),
      pressure_external = sum(contrib[!is_internal], na.rm = TRUE),
      internal_activity = sum(activity_weight[is_internal], na.rm = TRUE),
      external_activity = sum(activity_weight[!is_internal], na.rm = TRUE),
      exposure_import  = sum(contrib[!is_internal], na.rm = TRUE),
      total_activity = sum(activity_weight), total_visitors = sum(visitors)
    ), by = grid_id][, sim_day := day]

    rho_cat_list[[day]] <- dod[, .(
      pressure_internal_c = sum(contrib[is_internal], na.rm = TRUE),
      pressure_external_c = sum(contrib[!is_internal], na.rm = TRUE),
      activity_c          = sum(activity_weight, na.rm = TRUE),
      exposure_import_c   = sum(contrib[!is_internal], na.rm = TRUE)
    ), by = .(grid_id, category_label)][, sim_day := day]
  }
  rho_dt <- rbindlist(rho_list, fill = TRUE)

  rho_cat_dt <- rbindlist(rho_cat_list, fill = TRUE)
  rho_cat_dt[, category_broad := naics_category_broad(category_label)]
  setcolorder(rho_cat_dt, c("grid_id", "sim_day", "category_label", "category_broad"))
  fwrite(rho_cat_dt, file.path(out_dir, "tables", "od_kernel_rho_by_category.csv"))
  rm(rho_cat_list, rho_cat_dt)

  all_gd <- CJ(grid_id = unique(grid_dt$grid_id), sim_day = 1:n_days)
  rho_dt <- merge(all_gd, rho_dt, by = c("grid_id", "sim_day"), all.x = TRUE)
  for (cl in c("exposure_import",
               "pressure_internal", "pressure_external",
               "internal_activity", "external_activity", "total_activity"))
    rho_dt[is.na(get(cl)), (cl) := 0]
  rho_dt[, rho := total_activity]

  gmean <- mean(rho_dt$rho[rho_dt$rho > 0], na.rm = TRUE)
  if (is.na(gmean) || gmean <= 0) gmean <- 1
  rho_dt[, rho_norm := pmax(rho / gmean, 0.01)]
  pmean <- mean(rho_dt$rho_norm)
  if (abs(pmean - 1) > 0.05) rho_dt[, rho_norm := rho_norm / pmean]
  rho_dt[, log_rho := log(rho_norm)]
  rho_dt[, log_exposure_import := log1p(exposure_import)]

  od_origin_dt <- rbindlist(od_origin_list)
  od_origin_dt[, origin_geo_id := as.character(origin_geo_id)]
  arrow::write_feather(od_origin_dt, file.path(out_dir, "tables", "od_activity_by_origin_zone.feather"))
  fwrite(rbindlist(origin_prev_list)[, origin_geo_id := as.character(origin_geo_id)],
         file.path(out_dir, "tables", "origin_prevalence_by_day.csv"))
  fwrite(rho_dt, rho_file)

  gv <- as.data.table(arrow::read_feather(grid_file))
  if ("rate_geo_id" %in% names(gv)) gv[, rate_geo_id := NULL]
  gv <- unique(gv, by = "grid_id")
  gv <- merge(gv, unique(grid_dt[, .(grid_id, rate_geo_id)], by = "grid_id"), by = "grid_id", all.x = TRUE)
  stopifnot("06: grid_visits must have one row per cell" = !anyDuplicated(gv$grid_id))
  arrow::write_feather(gv, grid_file)

  invisible(NULL)
}

invisible(run_waves_by_state(process_od_kernel))
