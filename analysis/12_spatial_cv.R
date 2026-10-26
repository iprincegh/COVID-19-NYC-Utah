# Spatial cross-validation of the selected model on held-out reporting areas: five k-means blocks of
# MODZCTAs in NYC, one county at a time in Utah. Case-derived covariates of held-out areas are rebuilt
# without their own cases. Also pools the out-of-fold PCAtM over folds (case-weighted ratio of sums).
source("R/setup.R")
inla.setOption(num.threads = Sys.getenv("INLA_THREADS", "1:4"), inla.mode = "compact", smtp = "taucs")

N_DRAWS <- 200L
DEGENERATE_DISP <- 1e6

run_spatial_cv <- function(md, folds, arm, warm) {
  fm <- md$fms[[md$best_model]]
  eval_geos <- compute_eval_units(md$df, unit_col = "geo_id", y_col = "Y")
  if (length(eval_geos) && length(eval_geos) < uniqueN(md$df$geo_id) && nrow(folds[area_id %in% eval_geos]))
    folds <- folds[area_id %in% eval_geos]

  detail <- list(); diag <- list(); attr_inputs <- list()
  for (k in seq_len(max(folds$cv_fold))) {
    held <- folds[cv_fold == k, unique(area_id)]
    if (!length(held)) next
    holdout <- if (is.null(md$agg)) which(as.character(md$df$geo_id) %in% held) else
      which(md$agg$key$grid_id %in% unique(md$df$grid_id[as.character(md$df$geo_id) %in% held]))
    if (!length(holdout)) next
    hm <- holdout_meta(md, holdout)
    hf <- fit_holdout(md, fm, holdout, label = sprintf("%s fold %d", arm, k), warmstart_theta = warm,
                      model_name = md$best_model, n_attr_draws = N_DRAWS, compute_cpo = FALSE,
                      cov_fix = list(mode = "spatial"), n_pred_draws = N_DRAWS)
    if (!hf$ok) { log_msg(arm, " fold ", k, " failed"); next }
    if (is.numeric(hf$theta)) warm <- hf$theta
    log_msg(sprintf("%s fold %d done (%.1f min)", arm, k, hf$runtime))

    mu_v <- hf$mu[is.finite(hf$mu)]
    diag[[k]] <- data.table(fold = k, n_oof = length(holdout), n_areas = length(held), nb_size = hf$nb_size %||% NA_real_,
                            family = hf$family %||% NA_character_, runtime_min = round(hf$runtime, 1),
                            pred_disp = if (length(mu_v) > 1) var(mu_v) / max(mean(mu_v), 1e-12) else NA_real_)
    ho <- data.table(geo_id = hm$geo_id, time_id = hf$time_id, y_true = hf$y_true, mu = hf$mu,
                     pred_var = nb_pred_var(hf$mu, hf$sd_fit, hf$nb_size))
    d <- ho[, .(obs = sum(y_true, na.rm = TRUE), pred = sum(mu, na.rm = TRUE), sd = sqrt(sum(pred_var, na.rm = TRUE)),
                n_cells = .N, n_days = uniqueN(time_id)), by = geo_id]
    d[, `:=`(fold = k, lo = pmax(pred - 1.96 * sd, 0), hi = pred + 1.96 * sd)]
    d[, `:=`(covered = obs >= lo & obs <= hi, pct_err = abs(obs - pred) / pmax(obs, 1) * 100)]
    dg <- draw_group_totals(hf$y_draws, hf$y_true, hm$geo_id)
    if (!is.null(dg)) d[dg, `:=`(pred_draw = i.pred_draw, pred_med = i.med_draw, lo_draw = i.lo_draw, hi_draw = i.hi_draw,
                                 covered_draw = i.covered_draw), on = c(geo_id = "group")]
    detail[[k]] <- d
    a <- hf$attribution
    if (!is.null(a)) attr_inputs[[as.character(k)]] <- list(fold = k, mob_covs = a$mob_covs, X_oof = a$X_oof,
                                                            beta_draws = a$beta_draws, mu_ho = a$cell_attr$mu)
    rm(hf); gc()
  }
  detail <- rbindlist(detail, fill = TRUE); diag <- rbindlist(diag)

  # pooled metrics, excluding folds whose predictions are degenerate
  bad <- diag[is.finite(pred_disp) & pred_disp > DEGENERATE_DISP, fold]
  ok <- detail[!fold %in% bad]
  pt <- if ("pred_med" %in% names(ok) && any(is.finite(ok$pred_med))) ok$pred_med else ok$pred
  y <- ok$obs
  pt <- pmin(pt, max(quantile(y, 0.999, na.rm = TRUE), max(y, na.rm = TRUE) * 2, 1))
  pos <- is.finite(y) & y > 0
  r2 <- function(o, p) {
    k <- is.finite(o) & is.finite(p); s <- sum((o[k] - mean(o[k]))^2)
    if (sum(k) >= 2 && is.finite(s) && s > 0) 1 - sum((o[k] - p[k])^2) / s else NA_real_
  }
  metrics <- data.table(arm = arm, model = md$best_model, n_folds_degenerate = length(bad),
    n_areas_scored_ok = nrow(ok), n_areas_mdape_ok = sum(pos),
    area_RMSE_ok = sqrt(mean((y - pt)^2, na.rm = TRUE)), area_MAE_ok = mean(abs(y - pt), na.rm = TRUE),
    area_MdAPE_ok = if (sum(pos) >= 3) median(abs(y[pos] - pt[pos]) / y[pos] * 100) else NA_real_,
    area_Cov95_ok = mean(ok$covered, na.rm = TRUE), area_Cov95_draw_ok = mean(ok$covered_draw, na.rm = TRUE),
    area_Corr_ok = if (nrow(ok) >= 2) suppressWarnings(cor(y, pt)) else NA_real_, area_R2_ok = r2(y, pt))
  list(detail = detail, fold_diag = diag, metrics = metrics, attr_inputs = attr_inputs, pcatm_pooled = pool_pcatm(attr_inputs))
}

# Out-of-fold PCAtM: 1 - sum(mu exp(-eta)) / sum(mu) over all held-out units and folds, per posterior draw,
# where eta is the held-out mobility contrast X_oof beta.
pool_pcatm <- function(ai) {
  ai <- Filter(function(a) is.matrix(a$beta_draws) && NROW(a$X_oof) > 0 && length(a$mu_ho) == NROW(a$X_oof), ai)
  if (!length(ai)) return(NULL)
  M <- max(vapply(ai, function(a) nrow(a$beta_draws), integer(1)))
  num <- numeric(M); den <- 0
  for (a in ai) {
    eta <- as.matrix(a$X_oof) %*% t(a$beta_draws)
    eta <- eta[, if (ncol(eta) == M) seq_len(M) else sample.int(ncol(eta), M, replace = TRUE), drop = FALSE]
    mu <- pmax(a$mu_ho, 1e-12)
    num <- num + (sum(mu) - as.numeric(crossprod(mu, exp(-eta))))
    den <- den + sum(mu)
  }
  p <- num / den
  data.table(n_folds_pooled = length(ai), n_cells_total = sum(vapply(ai, function(a) NROW(a$X_oof), numeric(1))),
             n_draws = M, estimator = "Eq14_case_weighted_ratio_of_sums",
             pcatm_pooled_mean = mean(p), pcatm_pooled_sd = sd(p), pcatm_pooled_q05 = quantile(p, 0.05, names = FALSE),
             pcatm_pooled_q50 = quantile(p, 0.5, names = FALSE), pcatm_pooled_q95 = quantile(p, 0.95, names = FALSE))
}

arms <- select_runs()$run_id
for (arm in arms) {
  wv <- load_wave_for_cv(arm)
  folds <- create_area_folds(wv$md$df, wv$md$grid_sf, wv$state)$fold_dt
  res <- run_spatial_cv(wv$md, folds, arm, wv$model_results$theta_mode)
  cv_dir <- file.path(wv$out_dir, "validation")
  dir.create(cv_dir, showWarnings = FALSE)
  saveRDS(folds, file.path(cv_dir, "area_folds.rds"))
  fwrite(res$detail, file.path(cv_dir, "spatial_cv_area_detail.csv"))
  fwrite(res$fold_diag, file.path(cv_dir, "spatial_cv_fold_diagnostics.csv"))
  fwrite(res$metrics, file.path(cv_dir, "spatial_cv_primary_metrics.csv"))
  fwrite(res$pcatm_pooled, file.path(cv_dir, "pcatm_pooled.csv"))
  saveRDS(res$attr_inputs, file.path(cv_dir, "spatial_cv_attribution_inputs.rds"))
}
