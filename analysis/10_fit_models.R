# Fit the six model variants per wave with INLA and keep the one with the lowest WAIC.
# NYC uses the weekly A-matrix likelihood over latent cell-days; Utah uses daily county counts.
source("R/setup.R")
inla.setOption(num.threads = Sys.getenv("INLA_THREADS", "4:6"), inla.mode = "compact")

FIT_SETTINGS <- list(int_strategy = "ccd", h_val = 0.001, diff_logdens = 4, max_iter = 100, family_rule = "conditional")

fit_wave <- function(wave_id, wave_months) {
  out_dir <- get_wave_dir(state, wave_id)
  md <- readRDS(file.path(out_dir, "models", "model_data_gt.rds"))
  gt_dt <- copy(md$gt_dt)
  mask <- MASK_ARTIFACT_DAYS[[paste0(state, "_", wave_id)]]
  if (!is.null(mask)) gt_dt[gt_dt[mask, on = c("grid_id", "time_id"), which = TRUE], Y := NA_integer_]

  # spatial and temporal indices, unit graph and interaction clusters
  graph_file <- file.path(out_dir, "models", "grid_adj_unified.graph")
  cluster_graph <- file.path(out_dir, "models", "cluster_adj_unified.graph")
  grid_sf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), quiet = TRUE)
  grid_sf <- grid_sf[grid_sf$grid_id %in% unique(gt_dt$grid_id), ]
  nb <- read_inla_graph_as_nb(graph_file)
  stopifnot(length(nb) == nrow(grid_sf))
  df <- gt_dt[order(grid_id, time_id)]
  gi <- data.table(grid_id = sort(unique(df$grid_id)), idx_spatial = seq_len(uniqueN(df$grid_id)))
  df[, idx_spatial := gi$idx_spatial[match(grid_id, gi$grid_id)]]
  df[, idx_temporal := as.integer(factor(time_id))]
  df[, idx_interaction := .I]
  df <- st_clusters(df, grid_sf, nb, graph_file, cluster_graph)
  df[!is.finite(idx_st_cluster), idx_st_cluster := 1L]
  df[, idx_st_cluster := as.integer(idx_st_cluster)]
  if (nrow(grid_sf) > 300) fwrite(unique(df[, .(grid_id, idx_st_cluster)]), file.path(out_dir, "tables", "st_cluster_crosswalk.csv"))
  n_grids <- max(df$idx_spatial)
  df <- as.data.frame(add_rho_bins(stabilize_offsets(as.data.frame(df))))

  weekly <- state == "NY"
  if (weekly) {
    dw <- md$weekly_obs$day_week[, .(time_id, week_id)]
    df$idx_week <- dw$week_id[match(df$time_id, dw$time_id)]
    stopifnot(!anyNA(df$idx_week))
    obs <- md$weekly_obs$Y_week_cell[grid_id %in% unique(df$grid_id), .(grid_id, week_id, Y_week)]
  }
  disp <- if (weekly) var(obs$Y_week) / max(mean(obs$Y_week), 0.01) else var(df$Y, na.rm = TRUE) / max(mean(df$Y, na.rm = TRUE), 0.01)
  zero_frac <- if (weekly) mean(obs$Y_week == 0) else mean(df$Y == 0, na.rm = TRUE)
  covs <- select_covariates(df, potential_covariates(weekly), force_include = "log_rho")

  # residual mobility (spatial, daily-centred, serial-interval lag) and formulas
  df <- df |> add_spatial_orthogonalization(grid_sf, nb = nb) |> add_temporal_centering() |> add_temporal_lags()
  df <- add_st_index(df)
  fms <- build_formulas(covs, graph_file, detect_flags(df), disp, cluster_graph, n_grids, state = state, wave_id = wave_id,
                        st_K = uniqueN(df$idx_st_cluster), st_T = uniqueN(df$idx_temporal))
  if (sd(df$mob_weighted, na.rm = TRUE) > 0.01)
    for (mn in intersect(MOB_MODELS, names(fms))) fms[[mn]] <- update(fms[[mn]], ~ . + mob_weighted)
  agg <- if (weekly) build_weekly_agg(df, obs, dw) else NULL
  if (weekly) stopifnot(all(agg$key$n_days == 7L), all(Matrix::colSums(agg$A > 0) <= 1))

  comp <- list(); outputs <- list(); fits <- list()
  for (mn in names(fms)) {
    E_vec <- get_E(mn, df)
    fit <- do.call(fit_model, c(list(fms[[mn]], df, E_vec, paste0(mn, ": ", MODEL_LABELS[[mn]]), zero_frac, disp,
                                     cpo = FALSE, config = TRUE, n_grids = n_grids, state = state, wave_id = wave_id, agg = agg),
                                FIT_SETTINGS))
    log_msg(sprintf("%s_%s %-12s %s (%.1f min)", state, wave_id, mn, fit$status, fit$runtime))
    r <- fit$result
    gof <- rep(NA_real_, 4)
    if (!is.null(r)) {
      mu <- r$summary.fitted.values$mean
      if (weekly) { mu <- mu[seq_len(nrow(agg$A))]; y <- agg$Y } else y <- df$Y
      ss_tot <- sum((y - mean(y))^2, na.rm = TRUE)
      gof <- c(sqrt(mean((y - mu)^2, na.rm = TRUE)), mean(abs(y - mu), na.rm = TRUE),
               if (ss_tot > 0) 1 - sum((y - mu)^2, na.rm = TRUE) / ss_tot else NA_real_, cor(y, mu, use = "complete.obs"))
    }
    comp[[mn]] <- data.table(Model = mn, Description = MODEL_LABELS[[mn]],
      WAIC = r$waic$waic %||% NA_real_, DIC = r$dic$dic %||% NA_real_, pD = r$dic$p.eff %||% NA_real_,
      RMSE = gof[1], MAE = gof[2], R2 = gof[3], Correlation = gof[4],
      Runtime = fit$runtime, Status = fit$status, Family = fit$family %||% NA_character_,
      FamilyCollapsedFromNB = isTRUE(fit$family_collapsed))
    if (identical(fit$status, "OK")) {
      outputs[[mn]] <- extract_model_outputs(r, df, gi, mn, E_vec, agg = agg)
      fits[[mn]] <- r
    }
    rm(fit, r); gc()
  }

  comp <- rbindlist(comp, fill = TRUE)
  ok <- comp$Status == "OK" & is.finite(comp$WAIC)
  stopifnot("no model converged" = any(ok))
  comp[, c("dWAIC", "dDIC") := NA_real_]
  comp[ok, `:=`(dWAIC = WAIC - min(WAIC), dDIC = DIC - min(DIC))]
  best <- comp$Model[ok][which.min(comp$WAIC[ok])]

  fwrite(comp, file.path(out_dir, "tables", "model_comparison_unified.csv"))
  b <- outputs[[best]]
  for (tbl in c("fixed_effects", "spatial_effects", "temporal_effects", "fitted_values", "fitted_weekly"))
    if (!is.null(b[[tbl]])) fwrite(b[[tbl]], file.path(out_dir, "tables", paste0("unified_", tbl, ".csv")))
  saveRDS(fits[[best]], file.path(out_dir, "models", "best_inla_fit.rds"))
  saveRDS(list(state = state, wave_id = wave_id, comparison = comp, best_model = best, dispersion = disp,
               zero_frac = zero_frac, n_grids = n_grids, model_outputs = outputs, grid_index = gi,
               observation_scale = if (weekly) "cell-week" else "unit-day",
               theta_mode = fits[[best]]$mode$theta, inla_fit_path = file.path("models", "best_inla_fit.rds")),
          file.path(out_dir, "models", "model_results_unified.rds"))
  saveRDS(list(name = best, state = state, wave_id = wave_id, dispersion = disp, zero_frac = zero_frac, outputs = b),
          file.path(out_dir, "models", "best_model.rds"))
  writeLines(best, file.path(out_dir, "models", "best_model.txt"))
}

for (st in unique(select_runs()$state)) {
  set_transport_globals(st)
  run_waves_by_state(fit_wave, select_runs()[state == st])
}
