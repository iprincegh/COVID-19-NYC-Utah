# In-sample PCAtM: the share of fitted cases that the mobility-derived terms of the selected model add
# relative to all of them at zero (their average), per wave and per covariate, with the average mobility
# contribution (AMLP) allocated to POI categories by each unit-day's category shares.
source("R/setup.R")

N_POST <- 1000L; N_TERM_DRAWS <- 200L; SEED <- 20261008L
MOB_COVS <- c("own_prev_std", "ext_prev_std", "mob_weighted", "mob_spatial_lag_std", "log_rho_resid")
# category shares used to allocate each term: within-area pressure, incoming pressure, visit activity,
# and the neighbourhood mean of visit activity for the spatial term
SHARE_OF <- c(own_prev_std = "local", ext_prev_std = "import", mob_weighted = "activity",
              log_rho_resid = "activity", mob_spatial_lag_std = "neighbour")

# Joint draws of the requested fixed effects from the configurations stored by config = TRUE, used where
# inla.posterior.sample() aborts (the NYC grid fits). Each configuration is a Gaussian approximation of the
# latent field; the posterior is their mixture weighted by exp(log.posterior). Configurations with
# non-finite variances or means more than 10 posterior sd from summary.fixed are dropped, which on the
# grid leaves the modal configuration only.
beta_draws_from_configs <- function(fit, covs, n) {
  cf <- fit$misc$configs
  k <- match(covs, cf$contents$tag)
  stopifnot(!anyNA(k), all(cf$contents$length[k] == 1L))
  off <- max(cf$contents$start + cf$contents$length - 1L) - nrow(cf$config[[1]]$Q)
  pos <- cf$contents$start[k] - off
  ref_m <- fit$summary.fixed[covs, "mean"]; ref_s <- fit$summary.fixed[covs, "sd"]
  cmean <- function(cc) (cc$improved.mean %||% cc$mean)[pos]
  ok <- vapply(cf$config, function(cc) {
    d <- diag(as.matrix(cc$Qinv[pos, pos, drop = FALSE])); m <- cmean(cc)
    all(is.finite(d)) && all(d > 0) && all(is.finite(m)) && all(abs(m - ref_m) <= 10 * pmax(ref_s, 1e-8))
  }, logical(1))
  stopifnot(any(ok))
  lw <- vapply(cf$config, function(cc) cc$log.posterior, numeric(1))
  w <- ifelse(ok, exp(lw - max(lw[ok])), 0)
  pick <- sample.int(length(w), n, replace = TRUE, prob = w / sum(w))
  B <- matrix(NA_real_, n, length(covs), dimnames = list(NULL, covs))
  for (i in unique(pick)) {
    S <- as.matrix(cf$config[[i]]$Qinv[pos, pos, drop = FALSE]); S <- (S + t(S)) / 2
    e <- eigen(S, symmetric = TRUE)
    A <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
    r <- pick == i
    B[r, ] <- matrix(cmean(cf$config[[i]]), sum(r), length(covs), byrow = TRUE) + matrix(rnorm(sum(r) * length(covs)), sum(r)) %*% t(A)
  }
  B
}

pcatm_wave <- function(wave_id, wave_months) {
  arm_dir <- get_wave_dir(state, wave_id)
  out_dir <- file.path(CONFIG$output_dir, "pcatm", paste0(state, "_", wave_id))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  md <- readRDS(file.path(arm_dir, "models", "model_data_gt.rds"))
  res <- readRDS(file.path(arm_dir, "models", "model_results_unified.rds"))
  best <- res$best_model
  fe <- as.data.table(res$model_outputs[[best]]$fixed_effects)
  fv <- as.data.table(res$model_outputs[[best]]$fitted_values)
  covs <- intersect(MOB_COVS, fe$covariate)
  beta <- setNames(fe$mean, fe$covariate)[covs]; beta_sd <- setNames(fe$sd, fe$covariate)[covs]

  # covariates as in the fit, joined to the fitted cell-day (county-day) counts
  grid_sf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), quiet = TRUE)
  ids <- sort(unique(md$gt_dt$grid_id))
  grid_sf <- grid_sf[grid_sf$grid_id %in% ids, ]
  grid_sf <- grid_sf[order(match(grid_sf$grid_id, ids)), ]
  nb <- read_inla_graph_as_nb(file.path(arm_dir, "models", "grid_adj_unified.graph"))
  stopifnot(length(nb) == length(ids))
  dt <- as.data.frame(md$gt_dt) |> add_spatial_orthogonalization(grid_sf, nb = nb) |> add_temporal_centering() |>
    add_temporal_lags() |> as.data.table()
  d <- merge(dt[, c("grid_id", "time_id", covs), with = FALSE],
             fv[, .(grid_id, time_id, fitted = fitted_mean)], by = c("grid_id", "time_id"))
  X <- as.matrix(d[, ..covs]); X[is.na(X)] <- 0
  w <- d$fitted; sw <- sum(w)

  fit <- readRDS(file.path(arm_dir, res$inla_fit_path))
  set.seed(SEED)
  smp <- tryCatch(INLA::inla.posterior.sample(N_POST, fit, verbose = FALSE), error = function(e) NULL)
  B <- if (is.null(smp)) beta_draws_from_configs(fit, covs, N_POST) else {
    rn <- rownames(smp[[1]]$latent)
    idx <- match(paste0(covs, ":1"), rn)
    if (anyNA(idx)) idx <- match(covs, rn)
    stopifnot(!anyNA(idx))
    t(matrix(vapply(smp, function(s) s$latent[idx, 1], numeric(length(idx))), nrow = length(idx)))
  }
  rm(smp, fit)

  # each term alone
  terms <- rbindlist(lapply(seq_along(covs), function(j) {
    x <- X[, j]; amlp <- sum(w * beta[j] * x) / sw
    pj <- vapply(seq_len(N_TERM_DRAWS), function(i) 1 - sum(w * exp(-x * B[i, j])) / sw, numeric(1))
    aj <- B[seq_len(N_TERM_DRAWS), j] * sum(w * x) / sw
    data.table(component = covs[j], beta = beta[[j]], amlp = amlp, pcatm = 1 - sum(w * exp(-beta[j] * x)) / sw,
               pcatm_loglinear = 1 - exp(-amlp), wave = wave_id,
               pcatm_q025 = quantile(pj, 0.025, names = FALSE), pcatm_q975 = quantile(pj, 0.975, names = FALSE),
               amlp_q025 = quantile(aj, 0.025, names = FALSE), amlp_q975 = quantile(aj, 0.975, names = FALSE),
               n_term_draws = N_TERM_DRAWS)
  }))

  # all terms together
  m <- as.vector(X %*% beta)
  dr <- vapply(seq_len(nrow(B)), function(i) { mi <- as.vector(X %*% B[i, ]); c(1 - sum(w * exp(-mi)) / sw, sum(w * mi) / sw) }, numeric(2))
  q <- function(x, p) quantile(x, p, names = FALSE)
  amlp <- sum(w * m) / sw
  summary <- data.table(state = state, wave = wave_id, model = best, n_rows = nrow(d), covariates = paste(covs, collapse = "+"),
    amlp = amlp, pcatm_loglinear = 1 - exp(-amlp), pcatm = 1 - sum(w * exp(-m)) / sw,
    pcatm_q025 = q(dr[1, ], 0.025), pcatm_q50 = q(dr[1, ], 0.5), pcatm_q975 = q(dr[1, ], 0.975), n_draws = nrow(B),
    amlp_q025 = q(dr[2, ], 0.025), amlp_q50 = q(dr[2, ], 0.5), amlp_q975 = q(dr[2, ], 0.975), amlp_frac_below0 = mean(dr[2, ] < 0))

  # POI-category allocation of each term's AMLP (normal-approximation intervals from the posterior sd of beta)
  sh <- category_shares(arm_dir, md$gt_dt)
  key <- as.data.table(sh$key)
  cats <- sh$cat_levels
  lookup <- function(s) {
    l <- merge(d[, .(grid_id, time_id)], cbind(key, as.data.table(s)), by = c("grid_id", "time_id"), all.x = TRUE)
    for (cl in cats) l[is.na(get(cl)), (cl) := 0]
    l
  }
  act <- lookup(sh$activity)
  shares <- list(local = lookup(sh$local), import = lookup(sh$import), activity = act)
  if ("mob_spatial_lag_std" %in% covs) {
    um <- act[, lapply(.SD, mean), by = grid_id, .SDcols = cats][match(ids, grid_id)]
    M <- as.matrix(um[, ..cats]); M[is.na(M)] <- 0
    NM <- t(vapply(seq_along(ids), function(i) { k <- nb[[i]][nb[[i]] > 0]; if (length(k)) colMeans(M[k, , drop = FALSE]) else M[i, ] },
                   numeric(length(cats))))
    colnames(NM) <- cats
    nl <- merge(d[, .(grid_id, time_id)], data.table(grid_id = ids, NM), by = "grid_id", all.x = TRUE, sort = FALSE)
    shares$neighbour <- nl[match(paste(d$grid_id, d$time_id), paste(grid_id, time_id))]
  }
  categories <- rbindlist(lapply(covs, function(cv) {
    x <- X[, cv]; s <- shares[[SHARE_OF[[cv]]]]
    rbindlist(lapply(cats, function(cl) {
      ws <- sum(w * x * s[[cl]]) / sw
      data.table(category_broad = cl, component = cv, beta = beta[[cv]], amlp = beta[[cv]] * ws,
                 amlp_q025 = beta[[cv]] * ws + qnorm(0.025) * abs(ws) * beta_sd[[cv]],
                 amlp_q975 = beta[[cv]] * ws + qnorm(0.975) * abs(ws) * beta_sd[[cv]])
    }))
  }))
  categories[, wave := wave_id]
  categories[, component_total := sum(amlp), by = component]
  categories[, share_of_component := ifelse(component_total != 0, amlp / component_total, NA_real_)]
  totals <- categories[, .(category_amlp = sum(amlp)), by = category_broad]
  totals[, `:=`(grand_total = sum(categories$amlp), wave = wave_id)]
  totals[, share_of_total := category_amlp / grand_total]

  fwrite(summary, file.path(out_dir, "pcatm_wave_summary.csv"))
  fwrite(terms, file.path(out_dir, "pcatm_covariate_decomposition.csv"))
  fwrite(categories, file.path(out_dir, "pcatm_category_decomposition.csv"))
  fwrite(totals, file.path(out_dir, "pcatm_category_totals.csv"))
}

invisible(run_waves_by_state(pcatm_wave))
