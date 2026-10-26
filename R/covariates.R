# Covariate construction: infection pressure, residual mobility, lags, binning, screening and POI category shares.

add_contact_prevalence <- function(df, force = FALSE) {
  need <- c("pressure_internal", "exposure_import", "internal_activity", "external_activity")
  if (!all(need %in% names(df))) return(df)
  if (all(c("contact_prev_std", "log_contact_prev", "own_rate_lag", "ext_rate_lag", "own_prev_std", "ext_prev_std", "int_share") %in% names(df)) && !force) return(df)
  was_dt <- data.table::is.data.table(df); df <- data.table::as.data.table(df)
  for (v in need) df[is.na(get(v)), (v) := 0]
  has_zone <- all(c("geo_id", "time_id") %in% names(df))
  if (has_zone) {
    df[, `:=`(zone_internal_activity = sum(internal_activity), zone_external_activity = sum(external_activity),
              zone_pressure_internal = sum(pressure_internal), zone_exposure_import = sum(exposure_import)),
       by = .(geo_id, time_id)]
    df[, contact_prev := {a <- zone_internal_activity + zone_external_activity
                          ifelse(a > 0, (zone_pressure_internal + zone_exposure_import) / a, 0)}]
    df[, own_rate_lag := ifelse(zone_internal_activity > 0, zone_pressure_internal / zone_internal_activity, 0)]
    df[, ext_rate_lag := ifelse(zone_external_activity > 0, zone_exposure_import / zone_external_activity, 0)]
    df[, int_share := {a <- zone_internal_activity + zone_external_activity; ifelse(a > 0, zone_internal_activity / a, 0)}]
    df[, contact_prev_cell := {a <- internal_activity + external_activity
                               ifelse(a > 0, (pressure_internal + exposure_import) / a, NA_real_)}]
    df[, contact_prev_support := "reporting zone"]
  } else {
    df[, contact_prev := {a <- internal_activity + external_activity
                          ifelse(a > 0, (pressure_internal + exposure_import) / a, 0)}]
    df[, own_rate_lag := ifelse(internal_activity > 0, pressure_internal / internal_activity, 0)]
    df[, ext_rate_lag := ifelse(external_activity > 0, exposure_import / external_activity, 0)]
    df[, int_share := {a <- internal_activity + external_activity; ifelse(a > 0, internal_activity / a, 0)}]
    df[, contact_prev_support := "unit"]
  }
  df[, log_contact_prev := log1p(contact_prev)]
  s <- stats::sd(df$log_contact_prev, na.rm = TRUE)
  df[, contact_prev_std := if (is.finite(s) && s > 0) (log_contact_prev - mean(log_contact_prev, na.rm = TRUE)) / s else 0]
  df[, log_own_prev := log1p(own_rate_lag)]
  df[, log_ext_prev := log1p(ext_rate_lag)]
  for (v in c("own", "ext")) {
    lv <- paste0("log_", v, "_prev"); sv <- paste0(v, "_prev_std")
    s_v <- stats::sd(df[[lv]], na.rm = TRUE); m_v <- mean(df[[lv]], na.rm = TRUE)
    data.table::set(df, j = sv, value = if (is.finite(s_v) && s_v > 0) (df[[lv]] - m_v) / s_v else rep(0, nrow(df)))
  }
  if (was_dt) df else as.data.frame(df)
}

PREV_COVS_ALL <- c("contact_prev_std", "own_prev_std", "ext_prev_std",
                   "local_inf_pressure_std", "exposure_import_std")

prevalence_covariates <- function() c("own_prev_std", "ext_prev_std")

CALENDAR_COVS <- c("is_weekend", "is_holiday", "commercial_x_weekend")

POTENTIAL_COVARIATES <- c(
  "ses_contrast_std", "pct_older_std", "pct_younger_std", "pct_employed_std",
  "pct_unemployed_std", "building_density_std", "landuse_commercial",
  "log_contact_density_hourly_std", "log_daytime_pop_std",
  "log_repeat_ratio_std",
  "is_weekend", "is_holiday", "commercial_x_weekend",
  "contact_prev_std"
)

potential_covariates <- function(weekly_mode = FALSE) {
  pc <- setdiff(POTENTIAL_COVARIATES, PREV_COVS_ALL)
  pc <- c(pc, prevalence_covariates())
  if (weekly_mode) pc <- setdiff(pc, CALENDAR_COVS)
  pc
}

select_covariates <- function(df, potential, vif_threshold = 5, force_include = character(0),
                              protect = PREV_COVS_ALL) {
  avail <- c()
  for (v in potential) {
    if (!v %in% names(df)) next
    ok <- if (v %in% c("landuse_commercial")) sum(df[[v]], na.rm = TRUE) > 10
          else if (v %in% c("is_weekend", "is_holiday") || grepl("_x_", v)) TRUE
          else sd(df[[v]], na.rm = TRUE) > 0.01
    if (ok) avail <- c(avail, v)
  }
  force_cols <- intersect(force_include, names(df))
  if (length(force_cols) > 0)
    force_cols <- force_cols[vapply(force_cols, function(v) is.numeric(df[[v]]), logical(1))]
  numeric_avail <- avail[sapply(avail, function(v) is.numeric(df[[v]]) && !v %in% c("is_weekend", "is_holiday"))]
  if (length(numeric_avail) >= 2) {
    all_cols <- c(numeric_avail, force_cols)
    X <- as.data.frame(df[, all_cols, drop = FALSE])
    X <- X[complete.cases(X), , drop = FALSE]
    if (nrow(X) >= 2 * length(all_cols)) {
      repeat {
        n_opt <- length(numeric_avail)
        vifs <- sapply(seq_along(numeric_avail), function(j) {
          r2 <- summary(lm(X[[j]] ~ ., data = X[, -j, drop = FALSE]))$r.squared
          1 / (1 - r2)
        })
        names(vifs) <- numeric_avail
        droppable <- !numeric_avail %in% protect
        if (!any(droppable)) break
        vifs_droppable <- vifs
        vifs_droppable[!droppable] <- -Inf
        worst <- which.max(vifs_droppable)
        if (vifs_droppable[worst] > vif_threshold) {
          avail <- setdiff(avail, numeric_avail[worst])
          numeric_avail <- numeric_avail[-worst]
          X <- X[, -worst, drop = FALSE]
        } else break
      }
    }
  }
  avail
}

EXPOSURE_FLOOR <- 0.05

OFFSET_MIN     <- -4

OFFSET_MAX     <- 20

stabilize_offsets <- function(df) {
  df$offset_total <- pmax(pmin(df$offset_total, OFFSET_MAX), OFFSET_MIN)
  df$E_inla <- exp(df$offset_total)
  efl <- max(quantile(df$E_inla, 0.005, na.rm = TRUE), EXPOSURE_FLOOR)
  df$E_inla <- pmax(df$E_inla, efl)
  df$E_nomob <- df$E_inla
  df
}

add_spatial_orthogonalization <- function(df, grid_sf, nb = NULL) {
  if (is.null(nb))
    nb <- tryCatch(suppressWarnings(spdep::poly2nb(grid_sf, queen = TRUE, snap = 1)), error = function(e) NULL)
  lw <- if (!is.null(nb)) tryCatch(spdep::nb2listw(nb, style = "W", zero.policy = TRUE), error = function(e) NULL) else NULL
  if (is.null(lw)) return(df)
  dt <- as.data.table(df)
  ga <- dt[, .(avg = mean(log_rho, na.rm = TRUE)), by = idx_spatial][order(idx_spatial)]
  lv <- spdep::lag.listw(lw, ga$avg)
  ga[, resid := residuals(lm(avg ~ lv, na.action = na.exclude))]
  ga[, spatial_lag := lv]
  dt <- merge(dt, ga[, .(idx_spatial, avg, resid, spatial_lag)], by = "idx_spatial", all.x = TRUE)
  dt[, log_rho_resid := resid + (log_rho - avg)]
  dt[, mob_spatial_lag_std := scale(spatial_lag)[, 1]]
  dt[, c("avg", "resid", "spatial_lag") := NULL]
  as.data.frame(dt)
}

add_temporal_centering <- function(df) {
  rho_col <- if ("log_rho_resid" %in% names(df)) "log_rho_resid" else "log_rho"
  dt <- as.data.table(df)
  dt[, rho_daily_mean := mean(get(rho_col), na.rm = TRUE), by = idx_temporal]
  if (!"mob_day_mean" %in% names(dt)) dt[, mob_day_mean := rho_daily_mean]
  dt[, (rho_col) := get(rho_col) - rho_daily_mean]
  dt[, rho_daily_mean := NULL]
  as.data.frame(dt)
}

SERIAL_INTERVAL_LAG_WEIGHTS <- c(
  mob_lag3 = 0.4541, mob_lag5 = 0.3409, mob_lag7 = 0.1927,
  mob_lag14 = 0.0119, mob_lag21 = 0.0005
)

add_temporal_lags <- function(df, lags = c(3, 5, 7, 14, 21)) {
  dt <- as.data.table(df); setorder(dt, idx_spatial, idx_temporal)
  rho_col <- if ("log_rho_resid" %in% names(dt)) "log_rho_resid" else "log_rho"
  for (lag_k in lags) {
    lag_col <- paste0("mob_lag", lag_k)
    dt[, (lag_col) := shift(get(rho_col), lag_k, type = "lag"), by = idx_spatial]
    dt[is.na(get(lag_col)), (lag_col) := get(rho_col)[1], by = idx_spatial]
  }
  lag_cols_present <- intersect(names(SERIAL_INTERVAL_LAG_WEIGHTS), names(dt))
  if (length(lag_cols_present) > 0) {
    w <- SERIAL_INTERVAL_LAG_WEIGHTS[lag_cols_present]
    w <- w / sum(w)
    dt[, mob_weighted := Reduce(`+`, Map(function(v, wt) get(v) * wt,
                                          lag_cols_present, w))]
  }
  real_cols <- c(rho_col, paste0("mob_lag", lags))
  dt[, mob_lag_neg_raw := shift(get(rho_col), 7L, type = "lead"), by = idx_spatial]
  dt[, .lead_ok := !is.na(mob_lag_neg_raw)]
  resid_fml <- reformulate(real_cols, response = "mob_lag_neg_raw")
  resid_fit <- lm(resid_fml, data = dt[.lead_ok == TRUE])
  dt[, mob_lag_neg := 0]
  dt[.lead_ok == TRUE, mob_lag_neg := mob_lag_neg_raw - as.numeric(predict(resid_fit, newdata = .SD))]
  dt[, c("mob_lag_neg_raw", ".lead_ok") := NULL]
  as.data.frame(dt)
}

N_RHO_BINS     <- 20

add_rho_bins <- function(df, n_bins = N_RHO_BINS) {
  brks <- seq(min(df$log_rho, na.rm = TRUE) - 0.01,
              max(df$log_rho, na.rm = TRUE) + 0.01, length.out = n_bins + 1)
  df$log_rho_idx <- as.integer(cut(df$log_rho, breaks = brks, include.lowest = TRUE))
  df$log_rho_idx[is.na(df$log_rho_idx)] <- 1L
  df
}

# POI-category shares of each unit-day, computed over the unit's reporting area: of within-area infection
# pressure (local), of incoming pressure (import) and of visit activity (activity). The part not traced
# to any category's visits is "Unallocated".
category_shares <- function(out_dir, gt_dt) {
  zone <- unique(gt_dt[, .(grid_id, geo_id)])
  cat_dt <- merge(fread(file.path(out_dir, "tables", "od_kernel_rho_by_category.csv")), zone, by = "grid_id")
  cat_dt[, category_broad := naics_category_broad(category_label)]
  cat_dt <- cat_dt[, .(local = sum(pressure_internal_c, na.rm = TRUE), import = sum(pressure_external_c, na.rm = TRUE),
                       activity = sum(activity_c, na.rm = TRUE)), by = .(geo_id, time_id = sim_day, category_broad)]
  rho <- merge(fread(file.path(out_dir, "tables", "od_kernel_rho.csv")), zone, by = "grid_id")
  tot <- rho[, .(local = sum(pressure_internal, na.rm = TRUE), import = sum(pressure_external, na.rm = TRUE),
                 activity = sum(total_activity, na.rm = TRUE)), by = .(geo_id, time_id = sim_day)]
  key <- merge(unique(gt_dt[, .(grid_id, time_id, geo_id)]), tot, by = c("geo_id", "time_id"), all.x = TRUE)
  for (v in c("local", "import", "activity")) key[is.na(get(v)), (v) := 0]
  cats <- sort(unique(cat_dt$category_broad))
  share <- function(v) {
    S <- vapply(cats, function(cl) {
      x <- merge(key[, .(geo_id, time_id)], cat_dt[category_broad == cl, c("geo_id", "time_id", v), with = FALSE],
                 by = c("geo_id", "time_id"), all.x = TRUE, sort = FALSE)[[v]]
      x[is.na(x)] <- 0
      ifelse(key[[v]] > 0, x / key[[v]], 0)
    }, numeric(nrow(key)))
    cbind(S, Unallocated = pmax(0, 1 - rowSums(S)))
  }
  list(key = key[, .(grid_id, time_id)], cat_levels = c(cats, "Unallocated"),
       local = share("local"), import = share("import"), activity = share("activity"))
}

# NYC grid-level prevalence: within-area and incoming infection pressure and visit activity are pooled over
# the unit graph with an exponential kernel in network distance (sigma 500 m, cut at 3 km), and their
# ratios replace the reporting-area ratios. Cells with no pooled activity keep the reporting-area value.
pool_prevalence <- function(d, nb, xy, sigma = 500, cutoff = 6 * sigma) {
  G <- max(d$idx_spatial); T <- max(d$idx_temporal)
  stopifnot(length(nb) == G, nrow(xy) == G)
  ii <- rep(seq_along(nb), vapply(nb, function(z) if (identical(z, 0L)) 0L else length(z), 0L))
  jj <- unlist(lapply(nb, function(z) if (identical(z, 0L)) integer(0) else z))
  e <- unique(cbind(pmin(ii, jj), pmax(ii, jj)))
  g <- igraph::graph_from_edgelist(e, directed = FALSE)
  igraph::E(g)$weight <- sqrt((xy[e[, 1], 1] - xy[e[, 2], 1])^2 + (xy[e[, 1], 2] - xy[e[, 2], 2])^2)
  stopifnot(igraph::is_connected(g))

  wi <- wj <- integer(0); wx <- numeric(0)
  for (b in seq(1L, G, by = 2000L)) {
    idx <- b:min(b + 1999L, G)
    D <- igraph::distances(g, v = idx, weights = igraph::E(g)$weight)
    sel <- which(D < cutoff, arr.ind = TRUE)
    wi <- c(wi, idx[sel[, 1]]); wj <- c(wj, sel[, 2]); wx <- c(wx, exp(-D[sel] / sigma))
  }
  W <- Matrix::sparseMatrix(i = wi, j = wj, x = wx, dims = c(G, G))
  pooled <- function(col) {
    M <- matrix(0, G, T)
    M[cbind(d$idx_spatial, d$idx_temporal)] <- d[[col]]
    as.matrix(W %*% M)
  }
  pPI <- pooled("pressure_internal"); pPE <- pooled("pressure_external")
  pAI <- pooled("internal_activity"); pAE <- pooled("external_activity")
  at <- cbind(d$idx_spatial, d$idx_temporal)

  d[, `:=`(z_ai = sum(internal_activity), z_ae = sum(external_activity),
           z_pi = sum(pressure_internal), z_pe = sum(pressure_external)), by = .(geo_id, idx_temporal)]
  d[, own_rate_lag := ifelse(pAI[at] > 0, pPI[at] / pAI[at], ifelse(z_ai > 0, z_pi / z_ai, 0))]
  d[, ext_rate_lag := ifelse(pAE[at] > 0, pPE[at] / pAE[at], ifelse(z_ae > 0, z_pe / z_ae, 0))]
  d[, int_share := { a <- pAI[at] + pAE[at]; ifelse(a > 0, pAI[at] / a, 0) }]
  d[, contact_prev := int_share * own_rate_lag + (1 - int_share) * ext_rate_lag]
  d[, `:=`(log_contact_prev = log1p(contact_prev), log_own_prev = log1p(own_rate_lag), log_ext_prev = log1p(ext_rate_lag))]
  zs <- function(x) { s <- sd(x, na.rm = TRUE); if (is.finite(s) && s > 0) (x - mean(x, na.rm = TRUE)) / s else rep(0, length(x)) }
  d[, `:=`(contact_prev_std = zs(log_contact_prev), own_prev_std = zs(log_own_prev), ext_prev_std = zs(log_ext_prev),
           contact_prev_support = sprintf("graph-pooled exponential %.0f m", sigma))]
  d[, c("z_ai", "z_ae", "z_pi", "z_pe") := NULL]
  d
}
