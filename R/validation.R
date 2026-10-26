# Rolling-origin and spatial cross-validation: holdout refits, predictive draws, scores and
# out-of-fold PCAtM. Reference values of the mobility-derived covariates are zero (Section 3.5).

cv_scores <- function(y, mu, sd_fit, nb_size = NULL,
                      family = c("auto", "nbinomial", "poisson", "gaussian")) {
  family <- match.arg(family)
  ok <- !is.na(y) & !is.na(mu) & is.finite(mu)
  y <- y[ok]; mu <- mu[ok]; sd_fit <- sd_fit[ok]
  n <- length(y)
  empty <- data.table(RMSE = NA_real_, MAE = NA_real_,
                       MdAPE = NA_real_, MdAPE_n = 0L, MdAPE_pct_excl = NA_real_,
                       CRPS = NA_real_, Cov50 = NA_real_, Cov80 = NA_real_,
                       Cov95 = NA_real_, Cov95_gauss = NA_real_,
                       logS = NA_real_, DSS = NA_real_,
                       cal_intercept = NA_real_, cal_slope = NA_real_,
                       pit_W1 = NA_real_, N = 0L)
  if (n < 10) return(empty)

  cap <- max(2 * max(y), quantile(y, 0.999))
  mu  <- pmin(mu, cap)
  mu_safe <- pmax(mu, 1e-10)

  rmse <- sqrt(mean((y - mu)^2))
  mae  <- mean(abs(y - mu))

  mdape_ok  <- y > 0
  n_mdape   <- sum(mdape_ok)
  mdape     <- if (n_mdape >= 5) median(abs(y[mdape_ok] - mu[mdape_ok]) / y[mdape_ok]) else NA_real_
  mdape_pct_excl <- 100 * (1 - n_mdape / n)

  has_nb <- !is.null(nb_size) && is.finite(nb_size) && nb_size > 0
  fam <- if (family != "auto") family
         else if (has_nb) "nbinomial"
         else if (length(y) > 0 && all(y >= 0) && all(y == round(y))) "poisson"
         else "gaussian"

  if (fam == "nbinomial" && has_nb) {
    mc_crps <- function() {
      set.seed(42)
      idx <- if (n > 50000) sample.int(n, 50000) else seq_len(n)
      mean(vapply(idx, function(i) {
        m <- mu_safe[i]
        s1 <- rnbinom(200, size = nb_size, mu = m)
        s2 <- rnbinom(200, size = nb_size, mu = m)
        mean(abs(s1 - y[i])) - 0.5 * mean(abs(s1 - s2))
      }, numeric(1)))
    }
    crps <- if (requireNamespace("scoringRules", quietly = TRUE) &&
                requireNamespace("hypergeo", quietly = TRUE)) {
      tryCatch(
        mean(scoringRules::crps_nbinom(y, mu = mu_safe, size = nb_size)),
        error = function(e) mc_crps())
    } else {
      mc_crps()
    }
    cov_pair <- function(p) {
      a <- (1 - p) / 2
      lo <- qnbinom(a,     size = nb_size, mu = mu_safe)
      hi <- qnbinom(1 - a, size = nb_size, mu = mu_safe)
      mean(y >= lo & y <= hi)
    }
    cov50 <- cov_pair(0.50); cov80 <- cov_pair(0.80); cov95 <- cov_pair(0.95)
    logS  <- -mean(dnbinom(y, size = nb_size, mu = mu_safe, log = TRUE))
    pred_sd <- sqrt(mu_safe + mu_safe^2 / nb_size)
    DSS  <- mean(((y - mu_safe) / pred_sd)^2 + 2 * log(pred_sd))
    rpit <- randomized_pit_nb(y, mu_safe, nb_size)
    rpit <- rpit[is.finite(rpit) & rpit >= 0 & rpit <= 1]
    pit_W1 <- if (length(rpit) >= 50) {
      ord <- sort(rpit); m <- length(ord)
      mean(abs(ord - (seq_len(m) - 0.5) / m))
    } else NA_real_
  } else if (fam == "poisson") {
    crps <- if (requireNamespace("scoringRules", quietly = TRUE)) {
      mean(scoringRules::crps_pois(y, lambda = mu_safe))
    } else {
      set.seed(42)
      idx <- if (n > 50000) sample.int(n, 50000) else seq_len(n)
      mean(vapply(idx, function(i) {
        m <- mu_safe[i]
        s1 <- rpois(200, lambda = m); s2 <- rpois(200, lambda = m)
        mean(abs(s1 - y[i])) - 0.5 * mean(abs(s1 - s2))
      }, numeric(1)))
    }
    cov_pair_p <- function(p) {
      a <- (1 - p) / 2
      lo <- qpois(a,     lambda = mu_safe)
      hi <- qpois(1 - a, lambda = mu_safe)
      mean(y >= lo & y <= hi)
    }
    cov50 <- cov_pair_p(0.50); cov80 <- cov_pair_p(0.80); cov95 <- cov_pair_p(0.95)
    logS  <- -mean(dpois(y, lambda = mu_safe, log = TRUE))
    pred_sd <- sqrt(mu_safe)
    DSS  <- mean(((y - mu_safe) / pred_sd)^2 + 2 * log(pred_sd))
    rpit <- randomized_pit_pois(y, mu_safe)
    rpit <- rpit[is.finite(rpit) & rpit >= 0 & rpit <= 1]
    pit_W1 <- if (length(rpit) >= 50) {
      ord <- sort(rpit); m <- length(ord)
      mean(abs(ord - (seq_len(m) - 0.5) / m))
    } else NA_real_
  } else {
    sd_safe <- pmax(sd_fit, 1e-6)
    z <- (y - mu) / sd_safe
    crps <- mean(sd_safe * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)))
    cov_pair_g <- function(p) {
      zq <- qnorm((1 + p) / 2)
      mean(y >= mu - zq * sd_safe & y <= mu + zq * sd_safe)
    }
    cov50 <- cov_pair_g(0.50); cov80 <- cov_pair_g(0.80); cov95 <- cov_pair_g(0.95)
    logS <- -mean(dnorm(y, mean = mu, sd = sd_safe, log = TRUE))
    DSS  <- mean(((y - mu) / sd_safe)^2 + 2 * log(sd_safe))
    pit_W1 <- NA_real_
  }

  sd_safe_g <- pmax(sd_fit, 1e-6)
  cov95_gauss <- mean(y >= mu - 1.96 * sd_safe_g & y <= mu + 1.96 * sd_safe_g)

  cal_int <- NA_real_; cal_slope <- NA_real_
  if (n >= 30 && sd(log(mu_safe)) > 1e-6) {
    fit_cal <- tryCatch(suppressWarnings(
      stats::glm(y ~ log(mu_safe), family = stats::quasipoisson(link = "log"))
    ), error = function(e) NULL)
    if (!is.null(fit_cal) && length(coef(fit_cal)) == 2L &&
        all(is.finite(coef(fit_cal)))) {
      cal_int   <- unname(coef(fit_cal)[1])
      cal_slope <- unname(coef(fit_cal)[2])
    }
  }

  data.table(RMSE = rmse, MAE = mae,
             MdAPE = mdape, MdAPE_n = n_mdape, MdAPE_pct_excl = mdape_pct_excl,
             CRPS = crps,
             Cov50 = cov50, Cov80 = cov80, Cov95 = cov95,
             Cov95_gauss = cov95_gauss,
             logS = logS, DSS = DSS,
             cal_intercept = cal_int, cal_slope = cal_slope,
             pit_W1 = pit_W1, N = n)
}

randomized_pit_pois <- function(y, mu) {
  mu_safe <- pmax(mu, 1e-10)
  F_y    <- ppois(y,     lambda = mu_safe)
  F_y_m1 <- ppois(y - 1, lambda = mu_safe)
  F_y_m1[y == 0] <- 0
  u <- stats::runif(length(y))
  F_y_m1 + u * (F_y - F_y_m1)
}

randomized_pit_nb <- function(y, mu, nb_size) {
  mu_safe <- pmax(mu, 1e-10)
  F_y    <- pnbinom(y,     size = nb_size, mu = mu_safe)
  F_y_m1 <- pnbinom(y - 1, size = nb_size, mu = mu_safe)
  F_y_m1[y == 0] <- 0
  V <- runif(length(y))
  pmin(pmax(F_y_m1 + V * (F_y - F_y_m1), 0), 1)
}

modzcta_to_borough <- function(geo_id) {
  pre <- substr(as.character(geo_id), 1, 3)
  out <- rep("Other", length(geo_id))
  out[pre %in% c("100", "101", "102")]               <- "Manhattan"
  out[pre == "103"]                                   <- "Staten Island"
  out[pre == "104"]                                   <- "Bronx"
  out[pre == "112"]                                   <- "Brooklyn"
  out[pre %in% c("110", "111", "113", "114", "116")] <- "Queens"
  out
}

mask_holdout <- function(df, holdout_idx) {
  df$Y_true <- df$Y
  df$Y[holdout_idx] <- NA_integer_
  df
}

CASE_DERIVED_COVS <- c("contact_prev_std", "own_prev_std", "ext_prev_std", "local_inf_pressure_std", "exposure_import_std")

freeze_case_covariates_forward <- function(df, t0) {
  cols <- intersect(CASE_DERIVED_COVS, names(df))
  if (!length(cols)) return(df)
  ref_rows <- which(df$idx_temporal == t0 + 1L)
  if (!length(ref_rows)) return(df)
  ref_of_unit <- ref_rows[match(df$idx_spatial, df$idx_spatial[ref_rows])]
  later <- which(df$idx_temporal > t0 + 1L & !is.na(ref_of_unit))
  for (cc in cols) df[[cc]][later] <- df[[cc]][ref_of_unit[later]]
  df
}

INCUBATION_WEIGHTS_CV <- c(`7` = 0.30, `6` = 0.20, `5` = 0.20, `4` = 0.10, `3` = 0.10, `2` = 0.05, `1` = 0.05)

replace_own_case_pressure <- function(df, holdout_idx, alpha = NULL) {
  need <- c("case_rate", "case_rate_spillover", "internal_activity", "log_local_inf_pressure",
            "local_inf_pressure_std", "idx_spatial", "idx_temporal")
  if (!all(need %in% names(df)) || !length(holdout_idx)) return(df)
  if (is.null(alpha)) alpha <- tryCatch(CONFIG$spillover_alpha, error = function(e) NULL) %||% 0.7
  parent <- (df$case_rate_spillover - alpha * df$case_rate) / (1 - alpha)
  parent[!is.finite(parent) | parent < 0] <- 0
  key <- data.table::data.table(i = seq_len(nrow(df)), u = df$idx_spatial, t = df$idx_temporal, p = parent)
  lags <- as.integer(names(INCUBATION_WEIGHTS_CV)); w <- unname(INCUBATION_WEIGHTS_CV)
  pw <- numeric(nrow(df))
  for (k in seq_along(lags)) {
    lk <- key[, .(u, t = pmin(pmax(t - lags[k], 1L), max(t)), p)]
    lk <- key[lk[, .(u, t)], on = c("u", "t"), .(pk = p)]
    pw <- pw + w[k] * ifelse(is.finite(lk$pk), lk$pk, 0)
  }
  pw <- pw / sum(w)
  m <- mean(df$log_local_inf_pressure, na.rm = TRUE); s <- sd(df$log_local_inf_pressure, na.rm = TRUE)
  lip_ho <- log1p(df$internal_activity[holdout_idx] * pw[holdout_idx])
  df$local_inf_pressure_std[holdout_idx] <- if (is.finite(s) && s > 0) (lip_ho - m) / s else 0
  df <- rebuild_contact_prev(df, holdout_idx, parent_rate = pw[holdout_idx])
  df
}

rebuild_contact_prev <- function(df, idx, parent_rate) {
  if (!all(c("contact_prev_std", "log_contact_prev", "exposure_import", "internal_activity", "external_activity") %in% names(df))) return(df)
  zone <- all(c("zone_internal_activity", "zone_external_activity", "zone_exposure_import") %in% names(df))
  a_int <- if (zone) df$zone_internal_activity[idx] else df$internal_activity[idx]
  a_ext <- if (zone) df$zone_external_activity[idx] else df$external_activity[idx]
  p_ext <- if (zone) df$zone_exposure_import[idx] else df$exposure_import[idx]
  a <- a_int + a_ext
  cbar <- ifelse(a > 0, (a_int * parent_rate + p_ext) / a, 0)
  m <- mean(df$log_contact_prev, na.rm = TRUE); s <- sd(df$log_contact_prev, na.rm = TRUE)
  df$contact_prev_std[idx] <- if (is.finite(s) && s > 0) (log1p(cbar) - m) / s else 0
  if (all(c("own_prev_std", "log_own_prev") %in% names(df))) {
    mo <- mean(df$log_own_prev, na.rm = TRUE); so <- sd(df$log_own_prev, na.rm = TRUE)
    df$own_prev_std[idx] <- if (is.finite(so) && so > 0) (log1p(parent_rate) - mo) / so else 0
  }
  df
}

load_od_origin_tables <- function(out_dir) {
  fa <- file.path(out_dir, "tables", "od_activity_by_origin_zone.feather")
  fp <- file.path(out_dir, "tables", "origin_prevalence_by_day.csv")
  if (!file.exists(fa) || !file.exists(fp)) return(NULL)
  A <- data.table::as.data.table(arrow::read_feather(fa))
  A[, origin_geo_id := as.character(origin_geo_id)]
  Cp <- data.table::fread(fp)[, origin_geo_id := as.character(origin_geo_id)]
  list(A = A, C = Cp)
}

replace_contact_prev_by_origin <- function(df, md, holdout_geo, parent_of = NULL, pop = NULL) {
  tabs <- md$od_origin
  need <- c("contact_prev_std", "log_contact_prev", "geo_id", "time_id")
  if (is.null(tabs) || !all(need %in% names(df))) return(NULL)
  A <- data.table::copy(tabs$A); Cp <- data.table::copy(tabs$C)
  tl <- data.table::as.data.table(md$time_lookup)
  if (!"sim_day" %in% names(tl)) tl[, sim_day := time_id]
  origins <- unique(Cp$origin_geo_id)
  grp <- if (!is.null(parent_of)) parent_of(origins) else rep("all", length(origins))
  gdt <- data.table::data.table(origin_geo_id = origins, grp = grp)
  Cp <- merge(Cp, gdt, by = "origin_geo_id")
  if (is.null(pop)) pop <- data.table::data.table(origin_geo_id = origins, pop = 1)
  Cp <- merge(Cp, pop, by = "origin_geo_id", all.x = TRUE); Cp[is.na(pop), pop := 1]
  Cp[, held := origin_geo_id %in% holdout_geo]
  par_g <- Cp[held == FALSE, .(parent = sum(origin_rate * pop) / sum(pop)), by = .(grp, sim_day)]
  par_all <- Cp[held == FALSE, .(parent_all = sum(origin_rate * pop) / sum(pop)), by = sim_day]
  Cp <- merge(Cp, par_g, by = c("grp", "sim_day"), all.x = TRUE)
  Cp <- merge(Cp, par_all, by = "sim_day", all.x = TRUE)
  Cp[is.na(parent), parent := parent_all]
  Cp[, rate_cf := fifelse(held, parent, origin_rate)]
  Cp[is.na(rate_cf), rate_cf := 0]
  g2z <- unique(data.table::data.table(dest_grid_id = as.character(df$grid_id), geo_id = as.character(df$geo_id)))
  A[, dest_grid_id := as.character(dest_grid_id)]
  A <- merge(A, g2z, by = "dest_grid_id")
  m <- merge(A, Cp[, .(origin_geo_id, sim_day, rate_cf)], by = c("origin_geo_id", "sim_day"), all.x = TRUE)
  m[is.na(rate_cf), rate_cf := 0]
  cb <- m[, .(cbar = sum(activity * rate_cf) / sum(activity),
              cext = { e <- origin_geo_id != geo_id; if (any(e)) sum(activity[e] * rate_cf[e]) / sum(activity[e]) else 0 }),
          by = .(geo_id, sim_day)]
  cb[!is.finite(cbar), cbar := 0]; cb[!is.finite(cext), cext := 0]
  cb <- merge(cb, Cp[, .(geo_id = origin_geo_id, sim_day, cown = rate_cf)], by = c("geo_id", "sim_day"), all.x = TRUE)
  cb[is.na(cown), cown := 0]
  cb <- merge(cb, tl[, .(sim_day, time_id)], by = "sim_day")
  mu <- mean(df$log_contact_prev, na.rm = TRUE); sd_ <- sd(df$log_contact_prev, na.rm = TRUE)
  key <- data.table::data.table(i = seq_len(nrow(df)), geo_id = as.character(df$geo_id), time_id = df$time_id)
  key <- merge(key, cb[, .(geo_id, time_id, cbar, cext, cown)], by = c("geo_id", "time_id"), all.x = TRUE, sort = FALSE)
  data.table::setorder(key, i)
  ok <- !is.na(key$cbar)
  df$contact_prev_std[key$i[ok]] <- if (is.finite(sd_) && sd_ > 0) (log1p(key$cbar[ok]) - mu) / sd_ else 0
  if (all(c("own_prev_std", "ext_prev_std", "log_own_prev", "log_ext_prev") %in% names(df))) {
    mo <- mean(df$log_own_prev, na.rm = TRUE); so <- sd(df$log_own_prev, na.rm = TRUE)
    me <- mean(df$log_ext_prev, na.rm = TRUE); se <- sd(df$log_ext_prev, na.rm = TRUE)
    df$own_prev_std[key$i[ok]] <- if (is.finite(so) && so > 0) (log1p(key$cown[ok]) - mo) / so else 0
    df$ext_prev_std[key$i[ok]] <- if (is.finite(se) && se > 0) (log1p(key$cext[ok]) - me) / se else 0
  }
  attr(df, "contact_prev_rebuilt") <- sum(ok)
  df
}

compute_eval_units <- function(df,
                                unit_col = "grid_id",
                                y_col = "Y",
                                total_thresh = 50L,
                                pct_zero_thresh = 0.80) {
  if (!unit_col %in% names(df)) return(character(0))
  dt <- as.data.table(df)
  day_col <- intersect(c("sim_day", "time_id", "t", "day"), names(dt))[1]
  if (!is.na(day_col)) {
    per_day <- dt[, .(y_sum = sum(get(y_col), na.rm = TRUE)),
                  by = c(unit_col, day_col)]
    per <- per_day[, .(total = sum(y_sum, na.rm = TRUE),
                       pct_zero = mean(y_sum == 0, na.rm = TRUE)),
                   by = c(unit_col)]
  } else {
    per <- dt[, .(total = sum(get(y_col), na.rm = TRUE),
                  pct_zero = mean(get(y_col) == 0, na.rm = TRUE)),
              by = c(unit_col)]
  }
  ok <- per[total >= total_thresh & pct_zero <= pct_zero_thresh][[unit_col]]
  as.character(ok)
}

nb_pred_var <- function(mu, sd_fit, nb_size) {
  if (!is.null(nb_size) && is.finite(nb_size) && nb_size > 0) {
    mu + mu^2 / nb_size
  } else {
    mu + sd_fit^2
  }
}

ATTR_MOB_PATTERN <- "^(contact_prev_std|own_prev_std|ext_prev_std|local_inf_pressure_std|exposure_import_std|log_rho_resid|mob_weighted|mob_orth_lag0_std|mob_resid_std|mob_spatial_lag_std)$"

agg_design_X <- function(A, df, idx, pattern = ATTR_MOB_PATTERN) {
  cols <- grep(pattern, names(df), value = TRUE)
  X <- as.matrix(as.data.frame(df)[, cols, drop = FALSE]); X[is.na(X)] <- 0
  out <- as.matrix(A[idx, , drop = FALSE] %*% X); colnames(out) <- cols; out
}

compute_attribution_oof <- function(fit_result, df_cv, holdout_idx,
                                    mob_pattern = ATTR_MOB_PATTERN,
                                    n_draws = 200L,
                                    quantiles = c(0.05, 0.50, 0.95),
                                    X_oof = NULL) {
  fe <- fit_result$summary.fixed
  mob_covs <- grep(mob_pattern, rownames(fe), value = TRUE)
  if (length(mob_covs) == 0) return(NULL)

  if (is.null(X_oof)) {
    X_oof <- as.matrix(df_cv[holdout_idx, mob_covs, drop = FALSE])
  } else {
    X_oof <- as.matrix(X_oof[, mob_covs, drop = FALSE])
  }
  X_oof[is.na(X_oof)] <- 0
  mu_hat  <- fit_result$summary.fitted.values$mean[holdout_idx]
  mu_safe <- pmax(mu_hat, 1e-12)

  beta_draws  <- NULL
  approx_flag <- FALSE
  if (n_draws > 0L && requireNamespace("INLA", quietly = TRUE) &&
      !is.null(fit_result$misc$configs)) {
    sel <- setNames(as.list(rep(1L, length(mob_covs))), mob_covs)
    smps <- tryCatch(
      INLA::inla.posterior.sample(as.integer(n_draws), fit_result,
                                  selection = sel, num.threads = 1L,
                                  seed = cv_draw_seed(), verbose = FALSE),
      error = function(e) { warning("inla.posterior.sample failed: ", conditionMessage(e)); NULL })
    if (!is.null(smps)) {
      lat_names <- rownames(smps[[1]]$latent)
      rows <- match(paste0(mob_covs, ":1"), lat_names)
      if (any(is.na(rows))) rows <- match(mob_covs, lat_names)
      if (!any(is.na(rows))) {
        beta_draws <- t(vapply(smps, function(s) s$latent[rows, 1],
                                numeric(length(mob_covs))))
        if (length(mob_covs) == 1L)
          beta_draws <- matrix(beta_draws, ncol = 1L)
        colnames(beta_draws) <- mob_covs
      }
    }
  }
  if (is.null(beta_draws)) {
    bm <- fe[mob_covs, "mean"]; bs <- fe[mob_covs, "sd"]
    M  <- max(as.integer(n_draws), 200L)
    beta_draws <- mapply(function(m, s) rnorm(M, m, s), bm, bs)
    if (!is.matrix(beta_draws)) beta_draws <- matrix(beta_draws, nrow = M)
    colnames(beta_draws) <- mob_covs
    approx_flag <- TRUE
  }
  M <- nrow(beta_draws)

  eta_draws      <- X_oof %*% t(beta_draws)
  mu_nomob_draws <- mu_safe * exp(-eta_draws)
  delta_draws    <- mu_safe - mu_nomob_draws

  row_q <- function(mat, p)
    apply(mat, 1L, function(v) as.numeric(quantile(v, probs = p, na.rm = TRUE)))

  qlo <- quantiles[1]; qmed <- quantiles[2]; qhi <- quantiles[3]
  cell_attr <- data.table::data.table(
    holdout_pos   = seq_along(holdout_idx),
    mu            = mu_hat,
    mu_nomob_mean = rowMeans(mu_nomob_draws),
    mu_nomob_lo   = row_q(mu_nomob_draws, qlo),
    mu_nomob_hi   = row_q(mu_nomob_draws, qhi),
    delta_mean    = rowMeans(delta_draws),
    delta_lo      = row_q(delta_draws, qlo),
    delta_hi      = row_q(delta_draws, qhi))

  total_mu       <- sum(mu_safe)
  total_delta_dr <- colSums(delta_draws)
  pcatm_draws    <- total_delta_dr / max(total_mu, 1e-12)
  pcatm_summary  <- data.table::data.table(
    PCAtM_mean    = mean(pcatm_draws),
    PCAtM_lo      = as.numeric(quantile(pcatm_draws, qlo, na.rm = TRUE)),
    PCAtM_med     = as.numeric(quantile(pcatm_draws, qmed, na.rm = TRUE)),
    PCAtM_hi      = as.numeric(quantile(pcatm_draws, qhi, na.rm = TRUE)),
    attr_cases_mean = mean(total_delta_dr),
    attr_cases_lo = as.numeric(quantile(total_delta_dr, qlo, na.rm = TRUE)),
    attr_cases_hi = as.numeric(quantile(total_delta_dr, qhi, na.rm = TRUE)),
    total_mu      = total_mu,
    n_oof         = length(holdout_idx),
    n_mob_covs    = length(mob_covs),
    n_draws       = M,
    approx_indep  = approx_flag)

  beta_summary <- data.table::data.table(
    covariate = mob_covs,
    mean = fe[mob_covs, "mean"],
    sd   = fe[mob_covs, "sd"],
    q05  = apply(beta_draws, 2L, quantile, probs = qlo,  na.rm = TRUE),
    q50  = apply(beta_draws, 2L, quantile, probs = qmed, na.rm = TRUE),
    q95  = apply(beta_draws, 2L, quantile, probs = qhi,  na.rm = TRUE))

  list(cell_attr     = cell_attr,
       pcatm_summary = pcatm_summary,
       beta_summary  = beta_summary,
       beta_draws    = beta_draws,
       X_oof         = X_oof,
       mob_covs      = mob_covs)
}

predictive_count_draws <- function(fit_result, E_vec, holdout_idx, n_draws = 200L,
                                   family = "poisson", nb_size = NULL,
                                   predictor = "Predictor") {
  if (n_draws <= 0L || is.null(fit_result$misc$configs)) return(NULL)
  n_h <- length(holdout_idx); if (!n_h) return(NULL)
  smps <- tryCatch(
    INLA::inla.posterior.sample(as.integer(n_draws), fit_result,
                                selection = setNames(list(holdout_idx), predictor),
                                num.threads = 1L, seed = cv_draw_seed(), verbose = FALSE),
    error = function(e) { warning("predictive draws: inla.posterior.sample failed: ", conditionMessage(e)); NULL })
  if (is.null(smps)) return(NULL)
  E_h <- E_vec[holdout_idx]
  out <- matrix(NA_integer_, nrow = length(smps), ncol = n_h)
  for (k in seq_along(smps)) {
    eta <- as.numeric(smps[[k]]$latent[seq_len(n_h), 1])
    mu  <- pmin(E_h * exp(eta), 1e9)
    size_k <- nb_size
    hp <- smps[[k]]$hyperpar
    if (!is.null(hp)) { i <- grep("size.*nbinomial", names(hp), ignore.case = TRUE); if (length(i)) size_k <- as.numeric(hp[i[1]]) }
    out[k, ] <- if (identical(family, "nbinomial") && !is.null(size_k) && is.finite(size_k) && size_k > 0)
                  rnbinom(n_h, size = size_k, mu = mu) else rpois(n_h, mu)
  }
  out
}

cv_draw_seed <- function() 42L

draw_group_totals <- function(y_draws, y_true, groups, level = 0.95) {
  if (is.null(y_draws)) return(NULL)
  g <- as.character(groups); a <- (1 - level) / 2
  tot <- rowsum(t(y_draws), g)
  obs <- rowsum(as.numeric(y_true), g)[, 1]
  qs  <- t(apply(tot, 1, quantile, probs = c(a, 0.5, 1 - a), na.rm = TRUE))
  data.table(group = rownames(tot), obs = obs[rownames(tot)], pred_draw = rowMeans(tot),
             med_draw = qs[, 2], lo_draw = qs[, 1], hi_draw = qs[, 3])[, covered_draw := obs >= lo_draw & obs <= hi_draw][]
}

fit_holdout <- function(md, fm, holdout_idx, label, warmstart_theta = NULL,
                        model_name = NULL, n_attr_draws = 200L,
                        compute_cpo = TRUE, mask_idx = NULL, cov_fix = NULL,
                        n_pred_draws = 0L) {
  if (is.null(model_name)) model_name <- md$best_model
  if (is.null(model_name))
    stop("fit_cv_model(): no model_name given and md$best_model is unset; CV must score the model the fit selected, never a default.")
  agg_mode <- !is.null(md$agg)
  agg_cv <- NULL
  spatial_rebuild <- function(df_cv, held_geo, fallback) {
    r <- if (!is.null(md$od_origin)) replace_contact_prev_by_origin(df_cv, md, held_geo, parent_of = md$parent_of, pop = md$origin_pop) else NULL
    if (is.null(r)) fallback(df_cv) else r
  }
  if (agg_mode) {
    df_cv <- md$df
    agg_cv <- md$agg
    agg_cv$Y_true <- agg_cv$Y
    agg_cv$Y[if (!is.null(mask_idx)) mask_idx else holdout_idx] <- NA_integer_
    if (!is.null(cov_fix) && identical(cov_fix$mode, "spatial")) {
      lat_idx <- latent_rows_of(md, holdout_idx)
      held_geo <- unique(as.character(md$df$geo_id[lat_idx]))
      df_cv <- spatial_rebuild(df_cv, held_geo, function(d) replace_own_case_pressure_weekly(d, lat_idx, md, holdout_geo = held_geo))
    }
  } else {
    df_cv <- mask_holdout(md$df, if (!is.null(mask_idx)) mask_idx else holdout_idx)
    if (!is.null(cov_fix)) {
      if (identical(cov_fix$mode, "temporal")) df_cv <- freeze_case_covariates_forward(df_cv, cov_fix$t0)
      if (identical(cov_fix$mode, "spatial")) {
        held_geo <- unique(as.character(md$df$geo_id[holdout_idx]))
        df_cv <- spatial_rebuild(df_cv, held_geo, function(d) replace_own_case_pressure(d, holdout_idx))
      }
    }
  }
  E_vec <- get_E(model_name, df_cv)
  if (isTRUE(md$n_grids <= 300)) warmstart_theta <- NULL
  if (!is.null(cov_fix) && identical(cov_fix$mode, "temporal") || (agg_mode && !is.null(mask_idx) && is.null(cov_fix))) {
    fenv <- environment(fm); st0 <- if (!is.null(fenv)) get0("st_extra", envir = fenv, inherits = FALSE) else NULL
    if (!is.null(st0)) {
      grp <- if (agg_mode) df_cv$idx_week else df_cv$idx_temporal
      grp_i <- as.integer(factor(grp))
      obs_rows <- if (agg_mode) latent_rows_of(md, setdiff(seq_len(nrow(md$agg$A)), mask_idx)) else which(!is.na(df_cv$Y))
      t_obs <- sort(unique(grp_i[obs_rows]))
      st_new <- build_type4_structure(st0$cluster_graph, st0$K, st0$T, t_obs = t_obs)
      fenv2 <- new.env(parent = parent.env(fenv)); for (nm in ls(fenv, all.names = TRUE)) assign(nm, get(nm, fenv), fenv2)
      fenv2$st_extra <- st_new; fm <- fm; environment(fm) <- fenv2
    }
  }
  fit <- fit_model(fm, df_cv, E_vec,
                   model_name = label,
                   zero_frac = md$zero_frac, dispersion = md$dispersion,
                   cpo = compute_cpo, config = (n_attr_draws > 0L || n_pred_draws > 0L), cv_mode = TRUE,
                   n_grids = md$n_grids,
                   warmstart_theta = warmstart_theta,
                   state = md$state, wave_id = md$wave_id,
                   force_family = md$best_model_family, agg = agg_cv)
  if (is.null(fit$result)) return(list(ok = FALSE, runtime = 0))
  mu     <- fit$result$summary.fitted.values$mean[holdout_idx]
  sd_fit <- fit$result$summary.fitted.values$sd[holdout_idx]
  y_true <- if (agg_mode) agg_cv$Y_true[holdout_idx] else df_cv$Y_true[holdout_idx]
  hp <- fit$result$summary.hyperpar
  nb_size_row <- grep("size.*nbinomial", rownames(hp), ignore.case = TRUE)
  nb_size <- if (length(nb_size_row) > 0) hp[nb_size_row[1], "mean"] else NULL

  cpo_train_lcpo     <- NA_real_
  cpo_train_failures <- NA_integer_
  cpo_train_n        <- NA_integer_
  if (!is.null(fit$result$cpo) && length(fit$result$cpo$cpo)) {
    train_mask <- if (agg_mode) !is.na(agg_cv$Y) else !is.na(df_cv$Y)
    cpo_vec  <- fit$result$cpo$cpo[train_mask]
    fail_vec <- fit$result$cpo$failure[train_mask]
    ok_cpo   <- !is.na(cpo_vec) & cpo_vec > 0
    if (!is.null(fail_vec))
      ok_cpo <- ok_cpo & (is.na(fail_vec) | fail_vec == 0)
    if (any(ok_cpo)) {
      cpo_train_lcpo     <- -mean(log(pmax(cpo_vec[ok_cpo], 1e-300)))
      cpo_train_n        <- as.integer(sum(ok_cpo))
      cpo_train_failures <- as.integer(
        sum(!is.na(fail_vec) & fail_vec > 0, na.rm = TRUE))
    }
  }

  fold_family <- fit$result$.args$family %||%
                 (if (!is.null(nb_size)) "nbinomial" else "poisson")

  theta_post <- tryCatch(fit$result$mode$theta, error = function(e) NULL)
  if (is.null(theta_post) || !is.numeric(theta_post)) theta_post <- NULL
  fixed_full <- fit$result$summary.fixed

  attribution <- tryCatch(
    compute_attribution_oof(fit$result, df_cv, holdout_idx,
                            n_draws = n_attr_draws,
                            X_oof = if (agg_mode) agg_design_X(agg_cv$A, df_cv, holdout_idx) else NULL),
    error = function(e) {
      warning("compute_attribution_oof failed: ", conditionMessage(e)); NULL })

  strat_cols <- intersect(c("grid_id", "geo_id", "idx_temporal", "idx_spatial",
                            "area_strat", "urban_class", "wave_phase",
                            "cell_E_q"),
                          names(df_cv))
  holdout_dt <- if (agg_mode) {
    hm <- holdout_meta(md, holdout_idx)
    ud <- unique(as.data.table(df_cv)[, c("grid_id", intersect(strat_cols, names(df_cv))), with = FALSE], by = "grid_id")
    ud[, grid_id := as.character(grid_id)]
    hd <- merge(hm[, .(grid_id, idx_spatial, idx_temporal = time_id, week_id, geo_id)], ud[, setdiff(names(ud), c("idx_spatial", "idx_temporal", "geo_id")), with = FALSE], by = "grid_id", sort = FALSE)
    hd <- hd[match(paste(hm$grid_id, hm$week_id), paste(hd$grid_id, hd$week_id))]
    hd
  } else as.data.table(df_cv[holdout_idx, strat_cols, drop = FALSE])
  holdout_dt[, `:=`(y_true = y_true, mu = mu, sd_fit = sd_fit)]

  y_draws <- if (n_pred_draws > 0L)
    predictive_count_draws(fit$result, if (agg_mode) agg_cv$E else E_vec, holdout_idx, n_draws = n_pred_draws,
                           family = fold_family, nb_size = nb_size,
                           predictor = if (agg_mode) "APredictor" else "Predictor") else NULL

  daily_latent <- NULL
  if (agg_mode && n_pred_draws > 0L && isTRUE(md$n_grids <= 300) && !is.null(fit$result$misc$configs)) {
    lat_idx <- latent_rows_of(md, holdout_idx)
    smp <- tryCatch(INLA::inla.posterior.sample(as.integer(n_pred_draws), fit$result,
                                                selection = list(Predictor = lat_idx), num.threads = 1L,
                                                seed = cv_draw_seed(), verbose = FALSE),
                    error = function(e) NULL)
    if (!is.null(smp)) {
      Ev <- E_vec[lat_idx]; tid <- df_cv$time_id[lat_idx]
      tot <- vapply(smp, function(s) as.numeric(rowsum(Ev * exp(as.numeric(s$latent[seq_along(lat_idx), 1])), tid)), numeric(length(unique(tid))))
      tot <- matrix(tot, ncol = length(smp)); tids <- sort(unique(tid))
      daily_latent <- data.table(time_id = tids, mean = rowMeans(tot), med = apply(tot, 1, median),
                                 lo = apply(tot, 1, quantile, 0.025), hi = apply(tot, 1, quantile, 0.975))
    }
  }
  try({ fit$result$marginals.fitted.values <- NULL
        fit$result$marginals.random        <- NULL
        fit$result$marginals.linear.predictor <- NULL }, silent = TRUE)

  list(ok = TRUE, mu = mu, sd_fit = sd_fit, y_true = y_true,
       nb_size = nb_size, family = fold_family,
       runtime = fit$runtime,
       time_id = if (agg_mode) holdout_meta(md, holdout_idx)$time_id else df_cv$idx_temporal[holdout_idx],
       fixed   = fixed_full,
       fixed_full = fixed_full,
       theta   = theta_post,
       cpo_train_lcpo     = cpo_train_lcpo,
       cpo_train_failures = cpo_train_failures,
       cpo_train_n        = cpo_train_n,
       attribution = attribution,
       holdout_dt = holdout_dt,
       y_draws = y_draws,
       daily_latent = daily_latent)
}

rebuild_model_data <- function(gt_dt, out_dir, wave_id = NULL, weekly_obs = NULL) {
  weekly_mode <- !is.null(weekly_obs)
  grid_sf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), quiet = TRUE)
  grid_sf <- grid_sf[grid_sf$grid_id %in% unique(gt_dt$grid_id), ]

  aug <- build_augmented_adjacency(grid_sf, bridge_links, WORMHOLES, gt_dt = gt_dt)
  nb_aug <- aug$nb

  comp_check <- safe_n_comp_nb(nb_aug)
  if (comp_check$nc > 1) {
    cs <- table(comp_check$comp.id)
    keep <- which(comp_check$comp.id == as.integer(names(which.max(cs))))
    grid_sf <- grid_sf[keep, ]
    gt_dt <- gt_dt[grid_id %in% grid_sf$grid_id]
    nb_keep <- nb_aug[keep]
    old_to_new <- setNames(seq_along(keep), keep)
    nb_aug <- lapply(nb_keep, function(nbrs)
      as.integer(old_to_new[as.character(nbrs[nbrs %in% keep])]))
    attr(nb_aug, "class") <- "nb"
    attr(nb_aug, "region.id") <- as.character(seq_along(keep))
    nb_aug <- sanitize_nb(nb_aug)
  }

  dir.create(file.path(out_dir, "validation"), recursive = TRUE, showWarnings = FALSE)
  .gtag <- Sys.getpid()
  graph_file <- file.path(out_dir, "validation", sprintf("cv_grid_adj_unified_%s.graph", .gtag))
  safe_nb2INLA(graph_file, nb_aug)
  n_grids <- length(unique(grid_sf$grid_id))

  df <- as.data.table(gt_dt)[order(grid_id, time_id)]
  gi <- data.table(grid_id = sort(unique(df$grid_id)),
                   idx_spatial = seq_len(uniqueN(df$grid_id)))
  if ("idx_spatial" %in% names(df)) df[, idx_spatial := NULL]
  df <- merge(df, gi, by = "grid_id", all.x = TRUE)
  df[, idx_temporal := as.integer(factor(time_id))]
  df[, idx_interaction := .I]

  cluster_graph <- file.path(out_dir, "validation", sprintf("cv_cluster_adj_unified_%s.graph", .gtag))
  df <- st_clusters(df, grid_sf, nb_aug, graph_file, cluster_graph)
  st_res <- list(df = df, cluster_graph = cluster_graph)
  df <- as.data.frame(add_rho_bins(stabilize_offsets(df)))
  df <- df |>
    add_spatial_orthogonalization(grid_sf, nb = nb_aug) |>
    add_temporal_centering() |>
    add_temporal_lags()

  obs <- NULL; dw <- NULL
  if (weekly_mode) {
    dw <- weekly_obs$day_week[, .(time_id, week_id)]
    df$idx_week <- dw$week_id[match(df$time_id, dw$time_id)]
    stopifnot("every day must map to a reporting week" = !anyNA(df$idx_week))
    obs <- weekly_obs$Y_week_cell[grid_id %in% unique(df$grid_id), .(grid_id, week_id, Y_week)]
  }
  potential <- potential_covariates(weekly_mode)
  df <- add_contact_prevalence(df)
  avail <- select_covariates(df, potential,
                             force_include = c("log_rho"))
  flags <- detect_flags(df)
  disp  <- if (weekly_mode) var(obs$Y_week) / max(mean(obs$Y_week), 0.01)
           else var(df$Y, na.rm = TRUE) / max(mean(df$Y, na.rm = TRUE), 0.01)
  zero_frac <- if (weekly_mode) mean(obs$Y_week == 0) else mean(df$Y == 0, na.rm = TRUE)
  df <- add_st_index(df)
  fms <- build_formulas(avail, graph_file, flags, disp,
                        st_res$cluster_graph, n_grids,
                        state = state, wave_id = wave_id,
                        st_K = uniqueN(df$idx_st_cluster), st_T = uniqueN(df$idx_temporal))

  best_lags <- intersect("mob_weighted", names(df))
  m1_extras <- best_lags
  m1_extras <- m1_extras[sapply(m1_extras, function(v) sd(df[[v]], na.rm = TRUE) > 0.01)]
  add_extras <- function(flist) {
    if (length(m1_extras) == 0) return(flist)
    mob_models <- intersect(MOB_MODELS, names(flist))
    for (mn in mob_models)
      flist[[mn]] <- update(flist[[mn]],
                            as.formula(paste("~ . +", paste(m1_extras, collapse = " + "))))
    flist
  }
  fms <- add_extras(fms)

  agg <- NULL
  if (weekly_mode) {
    stopifnot("df rows must be unique unit-days" = !anyDuplicated(as.data.table(df)[, .(grid_id, time_id)]))
    agg <- build_weekly_agg(df, obs, dw)
    stopifnot("aggregation rows must be complete 7-day weeks" = all(agg$key$n_days == 7L))
  }

  list(df = df, grid_sf = grid_sf, graph_file = graph_file,
       fms = fms, dispersion = disp, zero_frac = zero_frac,
       n_grids = n_grids, n_days = max(df$idx_temporal),
       grid_index = gi, nb = nb_aug, cluster_graph = st_res$cluster_graph,
       state = state, wave_id = wave_id,
       agg = agg, weekly_obs = weekly_obs)
}

holdout_meta <- function(md, idx) {
  if (is.null(md$agg)) {
    return(data.table(grid_id = as.character(md$df$grid_id[idx]), idx_spatial = md$df$idx_spatial[idx],
                      time_id = md$df$idx_temporal[idx],
                      geo_id = if ("geo_id" %in% names(md$df)) as.character(md$df$geo_id[idx]) else NA_character_))
  }
  key <- md$agg$key[idx]
  u <- unique(as.data.table(md$df)[, .(grid_id, idx_spatial, geo_id = if ("geo_id" %in% names(md$df)) as.character(geo_id) else NA_character_)])
  wl <- md$weekly_obs$week_lookup; tl <- md$time_lookup
  wend <- if (!is.null(tl)) tl$time_id[match(wl$week_ending_date, as.Date(tl$date))] else NA_integer_
  m <- merge(key, u, by = "grid_id", sort = FALSE)
  m <- m[match(paste(key$grid_id, key$week_id), paste(m$grid_id, m$week_id))]
  data.table(grid_id = as.character(m$grid_id), idx_spatial = m$idx_spatial,
             time_id = wend[match(m$week_id, wl$week_id)], week_id = m$week_id, geo_id = m$geo_id)
}

latent_rows_of <- function(md, idx) {
  if (is.null(md$agg)) return(idx)
  A <- md$agg$A[idx, , drop = FALSE]
  which(Matrix::colSums(A != 0) > 0)
}

replace_own_case_pressure_weekly <- function(df, lat_idx, md, holdout_geo) {
  need <- c("internal_activity", "log_local_inf_pressure", "local_inf_pressure_std", "geo_id", "time_id")
  if (!all(need %in% names(df)) || !length(lat_idx) || is.null(md$weekly_rates)) return(df)
  df <- data.table::as.data.table(df)
  wr <- md$weekly_rates[!geo_id %in% holdout_geo]
  par_b <- wr[, .(parent = sum(case_rate * pop) / sum(pop)), by = .(borough, week_ending_date)]
  par_c <- wr[, .(parent_city = sum(case_rate * pop) / sum(pop)), by = week_ending_date]
  zb <- unique(md$weekly_rates[, .(geo_id, borough)])
  tl <- as.data.table(md$time_lookup)[, .(time_id, date = as.Date(date))]
  rows <- data.table(i = lat_idx, geo_id = as.character(df$geo_id[lat_idx]), time_id = df$time_id[lat_idx])
  rows <- merge(rows, tl, by = "time_id", sort = FALSE)
  rows[, date_known := date - 1L]
  rows[zb, borough := i.borough, on = "geo_id"]
  setkey(par_b, borough, week_ending_date); setkey(par_c, week_ending_date)
  pb <- par_b[rows, on = .(borough, week_ending_date <= date_known), mult = "last", .(i = i.i, parent = x.parent)]
  pc <- par_c[rows, on = .(week_ending_date <= date_known), mult = "last", .(i = i.i, parent_city = x.parent_city)]
  rows[pb, parent := i.parent, on = "i"]; rows[pc, parent_city := i.parent_city, on = "i"]
  rows[is.na(parent), parent := parent_city]; rows[is.na(parent), parent := 0]
  m <- mean(df$log_local_inf_pressure, na.rm = TRUE); s <- sd(df$log_local_inf_pressure, na.rm = TRUE)
  lip <- log1p(df$internal_activity[rows$i] * rows$parent)
  df$local_inf_pressure_std[rows$i] <- if (is.finite(s) && s > 0) (lip - m) / s else 0
  df <- rebuild_contact_prev(df, rows$i, parent_rate = rows$parent)
  df
}

create_area_folds <- function(df, grid_sf, st, k_ny = 5L) {
  df_dt <- as.data.table(df)
  y_by_area <- df_dt[, .(n_cases = sum(Y, na.rm = TRUE)),
                     by = .(area_id = as.character(geo_id))]
  y_by_area <- y_by_area[!area_id %in% c("unknown", "99999")]
  if (st == "NY") {
    gm <- unique(data.table(grid_id = as.character(df$grid_id),
                            area_id = as.character(df$geo_id)))
    gm <- gm[area_id != "unknown" & area_id != "99999"]
    cents <- suppressWarnings(st_centroid(grid_sf))
    cents_dt <- data.table(grid_id = as.character(cents$grid_id),
                           x = st_coordinates(cents)[, 1],
                           y = st_coordinates(cents)[, 2])
    gm <- merge(gm, cents_dt, by = "grid_id")
    area_cents <- gm[, .(x = mean(x), y = mean(y), n_grids = .N), by = area_id]
    k_use <- min(k_ny, nrow(area_cents))
    set.seed(42)
    km <- kmeans(area_cents[, .(x, y)], centers = k_use, nstart = 25)
    area_cents[, cv_fold := km$cluster]
    fold_dt <- area_cents[, .(area_id, cv_fold)]
    n_folds <- max(fold_dt$cv_fold)
  } else {
    county_lkp <- data.table(grid_id = as.character(grid_sf$grid_id),
                             county_name = grid_sf$county_name)
    gm <- unique(data.table(grid_id = as.character(df$grid_id),
                            area_id = as.character(df$geo_id)))
    gm <- merge(gm, county_lkp, by = "grid_id", all.x = TRUE)
    gm_valid <- gm[!is.na(county_name)]
    counties_sorted <- sort(unique(gm_valid$county_name))
    county_fold <- data.table(county_name = counties_sorted,
                              cv_fold = seq_along(counties_sorted))
    gm_valid <- merge(gm_valid, county_fold, by = "county_name")
    fold_dt <- unique(gm_valid[, .(area_id, cv_fold)])

    n_folds <- max(fold_dt$cv_fold, na.rm = TRUE)
  }
  fold_dt <- merge(fold_dt, y_by_area, by = "area_id", all.x = TRUE)
  fold_dt[is.na(n_cases), n_cases := 0]
  fold_case_totals <- fold_dt[, .(fold_n_cases = sum(n_cases), n_areas = .N), by = cv_fold]
  if (nrow(fold_case_totals) > 1) {
    imb <- max(fold_case_totals$fold_n_cases) /
           max(1, min(fold_case_totals$fold_n_cases))
    if (imb > 5)
      warning(sprintf(
        "create_area_folds (%s): fold case-total imbalance is %.1fx (max/min); folds are spatially blocked, not case-balanced -- see fold_case_totals attribute.",
        st, imb))
  }
  attr(fold_dt, "fold_case_totals") <- fold_case_totals

  grid_area <- unique(data.table(grid_id = as.character(df$grid_id),
                                 area_id = as.character(df$geo_id)))
  grid_area <- grid_area[area_id != "unknown" & area_id != "99999"]
  list(fold_dt = fold_dt, grid_area = grid_area, n_folds = n_folds)
}

get_nb_size <- function(model_results, wave_name = NULL, model_name = NULL) {
  if (is.null(model_name)) model_name <- model_results$best_model
  stopifnot("the fit records no selected model" = !is.null(model_name))
  hp <- model_results$model_outputs[[model_name]]$hyperparameters
  if (is.null(hp)) return(NULL)
  nb_row <- hp[grepl("size.*nbinomial", hp$parameter, ignore.case = TRUE), ]
  if (nrow(nb_row) > 0) nb_row$mean[1] else NULL
}

model_family_of <- function(model_results, model_name) {
  hp <- model_results$model_outputs[[model_name]]$hyperparameters
  if (!is.null(hp) && any(grepl("size.*nbinomial", hp$parameter, ignore.case = TRUE))) "nbinomial" else "poisson"
}

# Rebuild the selected model's data and formulas for cross-validation, exactly as fitted.
load_wave_for_cv <- function(wave_name) {
  st <- sub("_.*", "", wave_name)
  wave_id <- sub("^[A-Z]+_(wave[0-9]+).*$", "\\1", wave_name)
  out_dir <- file.path(CONFIG$output_dir, wave_name)
  set_transport_globals(st)

  model_data    <- readRDS(file.path(out_dir, "models", "model_data_gt.rds"))
  model_results <- readRDS(file.path(out_dir, "models", "model_results_unified.rds"))

  gt_clean <- copy(model_data$gt_dt)
  mask_tbl <- MASK_ARTIFACT_DAYS[[paste0(st, "_", wave_id)]]
  if (!is.null(mask_tbl)) {
    idx <- gt_clean[mask_tbl, on = c("grid_id", "time_id"), which = TRUE]
    gt_clean[idx, Y := NA_integer_]
  }
  weekly <- identical(st, "NY")
  stopifnot("CV must use the likelihood of the fit" =
              identical(model_results$observation_scale, if (weekly) "cell-week" else "unit-day"))
  weekly_obs <- if (weekly) model_data$weekly_obs else NULL
  md <- rebuild_model_data(gt_clean, out_dir, wave_id = wave_id, weekly_obs = weekly_obs)
  md$time_lookup <- model_data$time_lookup

  if (weekly) {
    wk <- fread(file.path(out_dir, "tables", "case_rates_weekly.csv"))
    wk[, `:=`(geo_id = standardize_fips(modzcta), week_ending_date = as.Date(week_ending_date))]
    pop <- unique(weekly_obs$Y_week[, .(geo_id = as.character(geo_id), pop = pop_zone)])
    wk <- merge(wk[, .(geo_id, week_ending_date, case_rate)], pop, by = "geo_id")
    wk[, borough := modzcta_to_borough(geo_id)]
    md$weekly_rates <- wk
  }

  md$od_origin <- load_od_origin_tables(out_dir)
  if (!is.null(md$od_origin)) {
    md$parent_of <- if (weekly) function(g) modzcta_to_borough(g) else function(g) rep("state", length(g))
    md$origin_pop <- if (weekly) unique(as.data.table(weekly_obs$Y_week)[, .(origin_geo_id = as.character(geo_id), pop = pop_zone)])
                     else as.data.table(md$df)[, .(pop = E[1]), by = .(origin_geo_id = as.character(geo_id))]
  }

  best_model <- model_results$best_model
  stopifnot("the fit's selected model must be among the rebuilt formulas" = best_model %in% names(md$fms))
  md$best_model <- best_model
  md$best_model_family <- model_family_of(model_results, best_model)

  list(md = md, model_results = model_results, nb_size = get_nb_size(model_results),
       state = st, wave_id = wave_id, out_dir = out_dir, best_model = best_model)
}
