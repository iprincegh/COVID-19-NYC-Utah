# Rolling-origin temporal cross-validation of the selected model. NYC forecasts the next reporting week
# from each completed week; Utah forecasts the next three days every seven days. Training data after the
# origin are masked and case-derived covariates are frozen at the origin.
source("R/setup.R")
inla.setOption(num.threads = Sys.getenv("INLA_THREADS", "4:9"), inla.mode = "compact", smtp = "taucs")

MIN_TRAIN <- 0.50; MAX_TRAIN <- 0.95; H_UT <- 3L; STEP_UT <- 7L; N_DRAWS <- 200L
SCORE_COLS <- c("RMSE", "MAE", "MdAPE", "CRPS", "Cov50", "Cov80", "Cov95", "logS", "DSS")

ci <- function(d) d[, `:=`(lo = pmax(pred - 1.96 * sd, 0), hi = pred + 1.96 * sd)]

run_temporal_cv <- function(md, arm, warm) {
  fm <- md$fms[[md$best_model]]
  n_days <- md$n_days
  min_train <- ceiling(n_days * MIN_TRAIN)
  max_t0 <- min(n_days - if (is.null(md$agg)) H_UT else 7L, floor(n_days * MAX_TRAIN))
  weekly <- !is.null(md$agg)
  if (weekly) {
    wl <- md$weekly_obs$week_lookup[observed == TRUE]
    tl <- as.data.table(md$time_lookup)
    wl[, t_end := tl$time_id[match(week_ending_date, as.Date(tl$date))]]
    wl <- wl[!is.na(t_end)]
    wl <- wl[t_end >= min_train & t_end <= max_t0 & (week_id + 1L) %in% wl$week_id]
    origins <- wl$t_end; h <- 7L
  } else {
    origins <- seq(min_train, max_t0, by = STEP_UT); h <- H_UT
    units <- unique(as.character(md$df$grid_id))
    eval_row <- if (length(units) <= 100L) as.character(md$df$grid_id) %in% compute_eval_units(md$df, "grid_id", "Y")
                else rep(TRUE, nrow(md$df))
  }
  geo <- if (md$state == "NY") as.character(md$df$geo_id) else NULL

  out <- list()
  for (t0 in origins) {
    if (weekly) {
      w0 <- wl[t_end == t0, week_id]
      holdout <- which(md$agg$key$week_id == w0 + 1L)
      mask <- which(md$agg$key$week_id > w0)
    } else {
      holdout <- which(md$df$idx_temporal > t0 & md$df$idx_temporal <= t0 + h & eval_row)
      mask <- which(md$df$idx_temporal > t0)
    }
    if (!length(holdout)) next
    hm <- holdout_meta(md, holdout)
    hf <- fit_holdout(md, fm, holdout, label = sprintf("%s t0=%d", arm, t0), warmstart_theta = warm,
                      model_name = md$best_model, n_attr_draws = N_DRAWS, mask_idx = mask,
                      cov_fix = if (weekly) NULL else list(mode = "temporal", t0 = t0, h = h), n_pred_draws = N_DRAWS)
    if (!hf$ok) { log_msg(arm, " t0=", t0, " failed"); next }
    if (is.numeric(hf$theta)) warm <- hf$theta
    log_msg(sprintf("%s t0=%d done (%.1f min)", arm, t0, hf$runtime))

    ho <- data.table(time_id = hf$time_id, idx_spatial = hm$idx_spatial, y_true = hf$y_true, mu = hf$mu,
                     sd_fit = hf$sd_fit, pred_var = nb_pred_var(hf$mu, hf$sd_fit, hf$nb_size))
    totals <- function(by, groups) {
      d <- ci(ho[, .(obs = sum(y_true, na.rm = TRUE), pred = sum(mu, na.rm = TRUE),
                     sd = sqrt(sum(pred_var, na.rm = TRUE)), n = .N), by = by])
      dg <- draw_group_totals(hf$y_draws, hf$y_true, groups)
      d[, group := do.call(paste, c(.SD, sep = "|")), .SDcols = by]
      if (!is.null(dg)) d[dg, `:=`(pred_draw = i.pred_draw, pred_med = i.med_draw, lo_draw = i.lo_draw,
                                   hi_draw = i.hi_draw, covered_draw = i.covered_draw), on = "group"]
      d[, `:=`(group = NULL, origin = t0, horizon = h)]
    }
    r <- list(
      scores = cbind(cv_scores(hf$y_true, hf$mu, hf$sd_fit, nb_size = hf$nb_size, family = hf$family %||% "auto"),
                     data.table(origin = t0)),
      fold = data.table(origin = t0, horizon = h, n_oof = length(holdout), runtime_min = round(hf$runtime, 1),
                        nb_size = hf$nb_size %||% NA_real_, family = hf$family %||% NA_character_),
      holdout = hf$holdout_dt[, `:=`(origin = t0, horizon = h)],
      daily = totals("time_id", as.character(hf$time_id))[, lead := time_id - t0],
      unit = totals("idx_spatial", as.character(hm$idx_spatial)))
    if (!is.null(geo)) {
      ho[, geo_id := hm$geo_id]
      r$zone_day <- totals(c("geo_id", "time_id"), paste(hm$geo_id, hf$time_id, sep = "|"))
      r$zone_window <- totals("geo_id", hm$geo_id)
    }
    out[[as.character(t0)]] <- r
    rm(hf); gc()
  }
  bind <- function(x) rbindlist(lapply(out, `[[`, x), fill = TRUE)

  ws <- bind("scores"); daily <- bind("daily")
  daily_agg <- ci(daily[, .(obs = obs[1], pred = mean(pred), sd = sqrt(mean(sd^2))), by = time_id])
  scores <- cbind(data.table(wave = arm, model = md$best_model, horizon = h, n_windows = nrow(ws), N = sum(ws$N)),
                  ws[, lapply(.SD, mean, na.rm = TRUE), .SDcols = intersect(SCORE_COLS, names(ws))])
  scores[, `:=`(Cov95_daily = mean(daily_agg$obs >= daily_agg$lo & daily_agg$obs <= daily_agg$hi),
                daily_RMSE = sqrt(mean((daily_agg$obs - daily_agg$pred)^2)))]

  res <- list(cell_scores = scores, window_scores = ws, fold_diag = bind("fold"), holdout_long = bind("holdout"),
              daily_raw = daily, unit_window_long = bind("unit"))
  if (!is.null(geo)) {
    zd <- bind("zone_day")
    zd <- ci(zd[, .(obs = obs[1], pred = mean(pred), sd = sqrt(mean(sd^2)), pred_med = mean(pred_med, na.rm = TRUE)),
                by = .(geo_id, time_id)])
    zt <- ci(zd[, .(obs = sum(obs), pred = sum(pred), pred_med = sum(pred_med), sd = sqrt(sum(sd^2))), by = geo_id])
    zw <- bind("zone_window")
    res$modz_win <- zw
    res$modz_scores <- data.table(wave = arm, model = md$best_model, horizon = h, n_modzctas = nrow(zt),
      modzcta_day_MdAPE = median(abs(zd$pred_med - zd$obs) / pmax(zd$obs, 1)) * 100,
      modzcta_day_bias = mean(zd$pred_med - zd$obs),
      modzcta_day_Cov95 = mean(zd$obs >= zd$lo & zd$obs <= zd$hi),
      modzcta_total_Corr = cor(zt$obs, zt$pred),
      modzcta_total_MdAPE = median(abs(zt$pred_med - zt$obs) / pmax(zt$obs, 1)) * 100,
      modzcta_total_Cov95 = mean(zt$obs >= zt$lo & zt$obs <= zt$hi),
      modzcta_window_MdAPE = median(abs(zw$pred_med - zw$obs) / pmax(zw$obs, 1)) * 100,
      modzcta_window_Cov95_draw = mean(zw$covered_draw, na.rm = TRUE))
  }
  res
}

arms <- select_runs()$run_id
for (arm in arms) {
  wv <- load_wave_for_cv(arm)
  res <- run_temporal_cv(wv$md, arm, wv$model_results$theta_mode)
  cv_dir <- file.path(wv$out_dir, "validation")
  dir.create(cv_dir, showWarnings = FALSE)
  fwrite(res$cell_scores, file.path(cv_dir, "temporal_cv_scores.csv"))
  if (!is.null(res$modz_scores)) fwrite(res$modz_scores, file.path(cv_dir, "modzcta_temporal_cv_scores.csv"))
  saveRDS(res, file.path(cv_dir, "temporal_cv_results.rds"))
}
