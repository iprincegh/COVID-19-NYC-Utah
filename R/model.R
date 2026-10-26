# Model variants, priors, INLA fitting and output extraction.

N_ST_CLUSTERS <- 40L

hyper_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                   phi  = list(prior = "pc",      param = c(0.5, 0.5)))

# County-days set to missing as reporting artefacts (Utah Wave 2).
MASK_ARTIFACT_DAYS <- list(
  UT_wave2 = data.table::data.table(grid_id = c(20L, 24L, 27L), time_id = c(45L, 75L, 42L)))

MODEL_LABELS <- c(
  Base          = "Baseline (BYM2+RW2, no mobility)",
  Temp          = "Temporal-only (RW2, no spatial RE)",
  LinMob        = "Linear mobility (BYM2+RW2+rho)",
  NLMob         = "Nonlinear mobility (BYM2+RW1+rw1(rho))",
  LinMob_ST_IV  = "LinMob + Knorr-Held Type IV (ICAR x RW1)",
  NLMob_ST_IV   = "NLMob + Knorr-Held Type IV (ICAR x RW1)"
)

MOB_MODELS <- c("Temp", "LinMob", "NLMob", "LinMob_ST_IV", "NLMob_ST_IV")

rename_quants <- function(dt) {
  nms <- names(dt)
  nms[nms == "0.025quant"] <- "q025"; nms[nms == "0.5quant"] <- "q50"; nms[nms == "0.975quant"] <- "q975"
  setnames(dt, names(dt), nms)
}

get_E <- function(mn, df) {
  if (mn %in% c("Base", "Temp")) df$E_nomob
  else                           df$E_inla
}

detect_flags <- function(df) {
  has <- function(x) x %in% names(df) && isTRUE(sd(df[[x]], na.rm = TRUE) > 0.01)
  list(rho = "log_rho" %in% names(df),
       exposure_import = has("exposure_import_std"), spatial_lag = has("mob_spatial_lag_std"))
}

is_tight_wave <- function(state, wave_id)
  isTRUE(state == "NY") && isTRUE(wave_id == "wave3")

build_type4_structure <- function(cluster_graph, K, T, t_obs = NULL) {
  g <- INLA::inla.read.graph(cluster_graph)
  stopifnot(g$n == K)
  W <- Matrix::sparseMatrix(i = rep(seq_len(K), g$nnbs), j = unlist(g$nbs), x = 1, dims = c(K, K))
  R_s <- Matrix::Diagonal(K, Matrix::rowSums(W)) - W
  R_s <- INLA::inla.scale.model(R_s, constr = list(A = matrix(1, 1, K), e = 0))
  D1 <- Matrix::bandSparse(T - 1, T, k = 0:1, diagonals = list(rep(-1, T - 1), rep(1, T - 1)))
  R_t <- Matrix::t(D1) %*% D1
  R_t <- INLA::inla.scale.model(R_t, constr = list(A = matrix(1, 1, T), e = 0))
  Q <- Matrix::kronecker(R_t, R_s)
  ind_t <- rep(1, T); if (!is.null(t_obs)) { ind_t[] <- 0; ind_t[t_obs] <- 1 }
  A_k <- Matrix::kronecker(matrix(ind_t, 1, T), Matrix::Diagonal(K))
  A_t <- Matrix::kronecker(Matrix::Diagonal(T), matrix(1, 1, K))
  A <- as.matrix(rbind(A_k, A_t[-T, , drop = FALSE]))
  list(Q = Q, A = A, e = rep(0, nrow(A)), rankdef = K + T - 1, K = K, T = T, cluster_graph = cluster_graph, t_obs = t_obs)
}

add_st_index <- function(df, weekly_mode = FALSE) {
  if (!"idx_st_cluster" %in% names(df)) return(df)
  K <- length(unique(df$idx_st_cluster)); grp <- if (weekly_mode) df$idx_week else df$idx_temporal
  df$idx_st <- as.integer((as.integer(factor(grp)) - 1L) * K + df$idx_st_cluster)
  df
}

build_formulas <- function(covs, graph, flags, dispersion = 1.1, cluster_graph = NULL, n_grids = Inf,
                           state = NULL, wave_id = NULL, st_group = "idx_temporal",
                           st_K = NULL, st_T = NULL) {
  fx <- paste(covs, collapse = " + ")
  small_domain <- is.finite(n_grids) && n_grids < 50

  tight <- is_tight_wave(state, wave_id)

  if (small_domain) {
    bym2 <- sprintf("f(idx_spatial,model='bym2',graph='%s',scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(2,0.01)),phi=list(prior='pc',param=c(0.5,0.8))))", graph)
  } else if (tight) {
    bym2 <- sprintf("f(idx_spatial,model='bym2',graph='%s',scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(0.3,0.01)),phi=list(prior='pc',param=c(0.5,0.5))))", graph)
  } else {
    bym2 <- sprintf("f(idx_spatial,model='bym2',graph='%s',scale.model=TRUE,hyper=hyper_bym2)", graph)
  }
  rw_u <- if (dispersion < 1.5) 1 else if (dispersion < 3) 0.75 else 0.5
  if (tight) rw_u <- min(rw_u, 0.3)
  rho_u <- if (tight) 0.3 else 0.5
  rw2 <- sprintf("f(idx_temporal,model='rw2',constr=TRUE,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))))", rw_u)
  rw1_t <- sprintf("f(idx_temporal,model='rw1',constr=TRUE,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))))", rw_u)
  rw1_rho <- sprintf("f(log_rho_idx,model='rw1',constr=TRUE,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))))", rho_u)
  rv <- if (flags$spatial_lag) "log_rho_resid" else "log_rho"

  fx_base <- paste(setdiff(covs, PREV_COVS_ALL), collapse = " + ")
  fms <- list(Base = as.formula(paste(c("Y ~ 1", fx_base, bym2, rw2), collapse = " + ")))
  if (flags$rho) {
    lm <- c("Y ~ 1", fx, bym2, rw2, rv)
    if (flags$spatial_lag) lm <- c(lm, "mob_spatial_lag_std")
    fms$LinMob <- as.formula(paste(lm, collapse = " + "))
    nl <- c("Y ~ 1", fx, bym2, rw1_t, rw1_rho)
    if (flags$spatial_lag) nl <- c(nl, "mob_spatial_lag_std")
    fms$NLMob <- as.formula(paste(nl, collapse = " + "))
    if (!is.null(cluster_graph)) {
      st_u <- if (small_domain) 0.5
              else if (dispersion > 5) 0.05
              else if (dispersion > 3) 0.08
              else 0.15
      if (tight) st_u <- 0.3
      st_diag <- if (small_domain) 0.5 else 1e-1
      st_extra <- NULL
      if (!is.null(st_K) && !is.null(st_T) && st_K > 1) {
        st_extra <- build_type4_structure(cluster_graph, st_K, st_T)
      }
      st_iv <- if (is.null(st_extra)) {
        sprintf("f(idx_st_cluster,model='besag',graph='%s',scale.model=TRUE,group=%s,control.group=list(model='rw1',hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01)))),hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))),diagonal=%s)", cluster_graph, st_group, st_u, st_u, st_diag)
      } else {
        sprintf("f(idx_st,model='generic0',Cmatrix=st_extra$Q,rankdef=%d,constr=FALSE,extraconstr=list(A=st_extra$A,e=st_extra$e),hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))),diagonal=%s)", st_extra$rankdef, st_u, st_diag)
      }

      lm_st4 <- c("Y ~ 1", fx, bym2, rw2, rv, st_iv)
      if (flags$spatial_lag)  lm_st4 <- c(lm_st4, "mob_spatial_lag_std")
      fms$LinMob_ST_IV <- as.formula(paste(lm_st4, collapse = " + "))

      st4_terms <- c("Y ~ 1", fx, bym2, rw1_t, rw1_rho, st_iv)
      if (flags$spatial_lag) st4_terms <- c(st4_terms, "mob_spatial_lag_std")
      fms$NLMob_ST_IV <- as.formula(paste(st4_terms, collapse = " + "))
    }
    tm <- c("Y ~ 1", fx, rw2, rv)
    if (flags$spatial_lag) tm <- c(tm, "mob_spatial_lag_std")
    fms$Temp <- as.formula(paste(tm, collapse = " + "))
  }

  ord <- intersect(names(MODEL_LABELS), names(fms))
  fms <- fms[ord]
  fenv <- new.env(parent = globalenv()); fenv$st_extra <- if (exists("st_extra", inherits = FALSE)) st_extra else NULL
  for (nm in names(fms)) environment(fms[[nm]]) <- fenv
  fms
}

nb_hyper_from <- function(dispersion, tight = FALSE) {
  u <- if (dispersion < 1.5) 1
       else if (dispersion < 3) 0.7
       else if (dispersion < 5) 0.5
       else if (dispersion < 100) 0.3
       else 0.1
  alpha <- if (dispersion < 5) 0.01 else 0.05
  if (dispersion >= 100) alpha <- 0.1
  if (tight) u <- min(u, 0.3)
  c(u, alpha)
}

fit_model <- function(formula, df, E_vec, model_name = NULL,
                      zero_frac = 0, dispersion = 1.1,
                      cpo = FALSE, config = TRUE, cv_mode = FALSE,
                      n_grids = Inf, warmstart_theta = NULL,
                      state = NULL, wave_id = NULL, force_family = NULL,
                      fixed_theta = NULL, agg = NULL,
                      int_strategy = NULL, h_val = NULL, diff_logdens = NULL,
                      max_iter = 25, family_rule = "marginal") {
  t0 <- Sys.time()
  if (!is.null(E_vec)) {
    stopifnot("E_vec length must equal nrow(df)" = length(E_vec) == nrow(df))
    df$.inla_E <- as.numeric(E_vec)
  } else {
    df$.inla_E <- 1
  }
  if (!is.null(agg)) {
    stopifnot("agg$A must have nrow(df) columns" = ncol(agg$A) == nrow(df),
              "agg$Y / agg$E must match nrow(agg$A)" = length(agg$Y) == nrow(agg$A) && length(agg$E) == nrow(agg$A))
    .inla_data <- c(as.list(df)[setdiff(names(df), c("Y", ".inla_E"))],
                    list(Y = as.numeric(agg$Y), .inla_E = as.numeric(agg$E)))
    .ctrl_pred <- list(A = agg$A, compute = TRUE, link = 1)
    n_obs_rows <- nrow(agg$A)
  } else {
    .inla_data <- df
    .ctrl_pred <- list(compute = TRUE, link = 1)
    n_obs_rows <- nrow(df)
  }
  is_st <- !is.null(model_name) && grepl("ST", model_name)
  small_domain <- is.finite(n_grids) && n_grids < 50
  tight <- is_tight_wave(state, wave_id)
  strat <- if (cv_mode) "gaussian"
           else if (small_domain || zero_frac > 0.90 || is_st) "gaussian"
           else "simplified.laplace"
  diag  <- if (small_domain) 0.5
           else if (is_st) 1e-1 else if (zero_frac > 0.90) 1e-2 else 1e-3
  nb_hyper <- nb_hyper_from(dispersion, tight = tight)
  use_poisson <- if (!is.null(force_family)) identical(force_family, "poisson")
                 else isTRUE(dispersion < 1.1)
  fam <- if (use_poisson) "poisson" else "nbinomial"
  int_strat <- int_strategy %||% "eb"
  dld <- diff_logdens %||% if (cv_mode) 2 else 3
  cv_h_val     <- 0.015
  cv_quantiles <- c(0.025, 0.975)

  primary_threads <- INLA::inla.getOption("num.threads") %||% "4:1"

  ctrl_mode <- if (!is.null(fixed_theta) && length(fixed_theta))
                 list(theta = fixed_theta, restart = FALSE)
               else if (!is.null(warmstart_theta) && length(warmstart_theta))
                 list(theta = warmstart_theta, restart = TRUE)
               else list()
  if (!is.null(fixed_theta) && length(fixed_theta)) config <- TRUE

  last_inla_msg <- character(0)
  run_inla <- function(safe_mode = FALSE, diag_override = NULL, strat_override = NULL) {
    d <- if (!is.null(diag_override)) diag_override
         else if (safe_mode) max(diag, 5e-2) else diag
    s <- if (!is.null(strat_override)) strat_override else strat
    h <- h_val %||% if (cv_mode) cv_h_val else 0.005
    cc <- list(dic = !cv_mode, waic = !cv_mode,
               cpo = cpo, config = config,
               return.marginals.predictor = FALSE)
    if (cv_mode) cc$quantiles <- cv_quantiles
    tryCatch({
      msgs <- capture.output({
        .res <- INLA::inla(formula, family = fam, data = .inla_data, E = .inla_E,
             control.inla = list(strategy = s, int.strategy = int_strat,
                                 diff.logdens = dld,
                                 diagonal = d,
                                 h = h, cmin = 0,
                                 tolerance.step = 1e-2,
                                 max.iter = max_iter,
                                 control.vb = list(enable = FALSE)),
             control.compute = cc,
             control.family = if (use_poisson) list() else
               list(hyper = list(theta = list(prior = "pc.prec", param = nb_hyper))),
             control.fixed = list(mean = 0, prec = 0.3),
             control.predictor = .ctrl_pred,
             control.mode = ctrl_mode,
             safe = safe_mode,
             num.threads = primary_threads, verbose = FALSE)
      }, type = "message")
      if (length(msgs) > 0) {
        keep <- grep("(failed|error|crash|lambda|abort|cannot|invalid|NaN|Inf|singular|max.iter|Illegal)",
                     msgs, ignore.case = TRUE, value = TRUE)
        if (length(keep)) last_inla_msg <<- c(last_inla_msg,
                                              paste(head(keep, 4), collapse = " | "))
      }
      if (is.null(.res$summary.fitted.values)) {
        last_inla_msg <<- c(last_inla_msg, "summary.fitted.values=NULL")
        return(NULL)
      }
      .res
    }, error = function(e) {
      last_inla_msg <<- c(last_inla_msg, paste0("R-error: ", conditionMessage(e)))
      NULL
    })
  }
  result <- run_inla(safe_mode = FALSE)
  if (is.null(result) && length(ctrl_mode)) {
    gc(full = TRUE)
    if (cv_mode) last_inla_msg <- c(last_inla_msg, "<retry: drop warmstart>")
    ctrl_mode <- list(restart = TRUE)
    warmstart_theta <- NULL
    result <- run_inla(safe_mode = FALSE)
  }
  if (is.null(result)) {
    gc(full = TRUE)
    result <- run_inla(safe_mode = TRUE, diag_override = if (small_domain) 0.8 else max(diag, 5e-2))
  }
  if (is.null(result)) {
    gc(full = TRUE)
    result <- run_inla(safe_mode = TRUE, diag_override = if (small_domain) 1.0 else 0.1, strat_override = "gaussian")
  }
  if (is.null(result) && cv_mode) {
    gc(full = TRUE)
    .saved_fam <- fam
    .saved_use_p <- use_poisson
    fam <- "poisson"
    use_poisson <- TRUE
    result <- run_inla(safe_mode = TRUE,
                       diag_override = if (small_domain) 1.0 else 0.1,
                       strat_override = "gaussian")
    if (!is.null(result)) {
      last_inla_msg <- c(last_inla_msg, "<rescued by Poisson fallback>")
    } else {
      fam <- .saved_fam
      use_poisson <- .saved_use_p
    }
  }
  if (cv_mode && is.null(result) && length(last_inla_msg)) {
    message("[fit_model CV-mode FAIL] ", paste(last_inla_msg, collapse = " ;; "))
  }
  if (is.null(result) && !cv_mode) {
    gc(full = TRUE)
    r3_diag <- if (small_domain) 1.0 else 0.1
    r3_int <- int_strategy %||% "eb"
    r3_h   <- h_val %||% 0.005
    r3_dld <- diff_logdens %||% 3
    tryCatch({
      invisible(capture.output({
        .res <- INLA::inla(formula, family = fam, data = .inla_data, E = .inla_E,
             control.inla = list(strategy = "gaussian", int.strategy = r3_int,
                                 diff.logdens = r3_dld,
                                 diagonal = r3_diag, h = r3_h, cmin = 0,
                                 tolerance.step = 1e-2, max.iter = max_iter,
                                 control.vb = list(enable = FALSE)),
             control.compute = list(dic = TRUE, waic = TRUE, cpo = cpo, config = config),
             control.family = if (use_poisson) list() else
               list(hyper = list(theta = list(prior = "pc.prec", param = nb_hyper))),
             control.fixed = list(mean = 0, prec = 0.3),
             control.predictor = .ctrl_pred,
             safe = TRUE,
             num.threads = "2:1", verbose = FALSE)
      }, type = "message"))
      if (!is.null(.res$summary.fitted.values)) result <- .res
    }, error = function(e) NULL)
  }
  if (is.null(result) && !cv_mode) {
    gc(full = TRUE)
    r4_diag <- if (small_domain) 1.0 else 0.2
    r4_strat <- if (small_domain) "gaussian" else "simplified.laplace"
    r4_int <- int_strategy %||% "eb"
    r4_h   <- h_val %||% 0.005
    r4_dld <- diff_logdens %||% 3
    tryCatch({
      invisible(capture.output({
        .res <- INLA::inla(formula, family = fam, data = .inla_data, E = .inla_E,
             control.inla = list(strategy = r4_strat, int.strategy = r4_int,
                                 diff.logdens = r4_dld,
                                 diagonal = r4_diag, h = r4_h, cmin = 0,
                                 tolerance.step = 1e-2, max.iter = max_iter,
                                 control.vb = list(enable = FALSE)),
             control.compute = list(dic = TRUE, waic = TRUE, cpo = cpo, config = config),
             control.family = if (use_poisson) list() else
               list(hyper = list(theta = list(prior = "pc.prec", param = nb_hyper))),
             control.fixed = list(mean = 0, prec = 0.3),
             control.predictor = .ctrl_pred,
             control.mode = list(restart = TRUE),
             safe = TRUE,
             num.threads = "1:1", verbose = FALSE)
      }, type = "message"))
      if (!is.null(.res$summary.fitted.values)) result <- .res
    }, error = function(e) NULL)
  }

  family_collapsed <- FALSE
  if (identical(family_rule, "conditional") && is.null(force_family) && !use_poisson) {
    .collapse_at <- 1e4
    if (!is.null(result)) {
      hp <- result$summary.hyperpar
      sz_row <- grep("size.*nbinomial", rownames(hp), ignore.case = TRUE)
      sz_med <- if (length(sz_row)) hp[sz_row[1], "0.5quant"] else NA_real_
      if (is.finite(sz_med) && sz_med > .collapse_at) {
        gc(full = TRUE)
        .nb_result <- result
        fam <- "poisson"; use_poisson <- TRUE
        p_res <- run_inla(safe_mode = FALSE)
        if (is.null(p_res)) p_res <- run_inla(safe_mode = TRUE,
                                              diag_override = if (small_domain) 0.8 else max(diag, 5e-2))
        if (!is.null(p_res)) {
          result <- p_res; family_collapsed <- TRUE
        } else {
          fam <- "nbinomial"; use_poisson <- FALSE; result <- .nb_result
        }
        rm(.nb_result)
      }
    } else {
      gc(full = TRUE)
      fam <- "poisson"; use_poisson <- TRUE
      result <- run_inla(safe_mode = FALSE)
      if (is.null(result)) result <- run_inla(safe_mode = TRUE,
                                              diag_override = if (small_domain) 0.8 else max(diag, 5e-2),
                                              strat_override = "gaussian")
      if (!is.null(result)) family_collapsed <- TRUE
      else { fam <- "nbinomial"; use_poisson <- FALSE }
    }
  }

  runtime <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(result)) {
    if (!cv_mode) message("[fit_model FAIL] ", model_name %||% "", " :: ",
                          if (length(last_inla_msg)) paste(last_inla_msg, collapse = " ;; ")
                          else "(no INLA/R message captured)")
    return(list(result = NULL, runtime = runtime, status = "FAILED",
                convergence = NULL, family = fam,
                family_collapsed = family_collapsed,
                fail_msg = paste(last_inla_msg, collapse = " ;; ")))
  }

  fitted_rescaled <- FALSE
  {
    fv <- result$summary.fitted.values
    .Y_obs <- if (is.null(agg)) df$Y else agg$Y
    .E_obs <- if (is.null(agg)) E_vec else agg$E
    y_ok <- !is.na(.Y_obs)
    if (!is.null(fv) && nrow(fv) >= n_obs_rows && length(.E_obs) == n_obs_rows && any(y_ok)) {
      .obs <- seq_len(n_obs_rows)
      ratio <- sum(fv$mean[.obs][y_ok], na.rm = TRUE) / max(sum(.Y_obs[y_ok], na.rm = TRUE), 1)
      if (is.finite(ratio) && ratio < 1e-2) {
        for (cc in intersect(c("mean", "sd", "0.025quant", "0.5quant", "0.975quant", "mode"),
                             names(fv)))
          fv[[cc]][.obs] <- fv[[cc]][.obs] * .E_obs
        result$summary.fitted.values <- fv
        fitted_rescaled <- TRUE
        message("[fit_model] summary.fitted.values were on the rate scale (INLA ",
                as.character(utils::packageVersion("INLA")), "); rescaled by E to counts.")
      }
    }
  }

  if (cv_mode) return(list(result = result, runtime = runtime, family = fam,
                           family_collapsed = family_collapsed,
                           fitted_rescaled = fitted_rescaled))

  waic_val <- result$waic$waic
  pD <- result$dic$p.eff
  max_fit <- max(result$summary.fitted.values$mean[seq_len(n_obs_rows)], na.rm = TRUE)
  waic_bad <- !is.null(waic_val) && !is.na(waic_val) &&
              (is.infinite(waic_val) || waic_val > 1e7)
  pD_bad   <- !is.null(pD) && !is.na(pD) && pD < 0
  status <- if (waic_bad || pD_bad) "UNSTABLE"
            else if (max_fit > 1e6) "IMPLAUSIBLE" else "OK"
  if (status == "OK" && cpo && is.null(result$cpo)) {
    result <- tryCatch(INLA::inla.cpo(result, verbose = FALSE),
                       error = function(e) result)
  }
  if (status == "OK" && cpo && !is.null(result$cpo$failure)) {
    n_fail <- sum(result$cpo$failure > 0, na.rm = TRUE)
    if (n_fail > 0 && n_fail <= 5000) {
      result <- tryCatch(INLA::inla.cpo(result, verbose = FALSE),
                         error = function(e) result)
    }
  }

  list(result = if (status == "OK") result else NULL, runtime = runtime, status = status,
       convergence = list(
         mlik_int    = result$mlik["log marginal-likelihood (integration)", 1],
         mlik_gauss  = result$mlik["log marginal-likelihood (Gaussian)", 1],
         mode_status = result$mode$mode.status,
         n_hyperpar  = nrow(result$summary.hyperpar),
         cpu_total_s = sum(result$cpu.used, na.rm = TRUE),
         max_fitted  = max_fit, waic = waic_val, pD = pD
       ),
       family = fam, family_collapsed = family_collapsed,
       fitted_rescaled = fitted_rescaled)
}

extract_model_outputs <- function(result, df, grid_index, model_name, E_vec = df$E_inla, agg = NULL) {
  if (is.null(result)) return(NULL)
  out <- list()
  n_obs <- if (is.null(agg)) nrow(df) else nrow(agg$A)
  lat_rows <- if (is.null(agg)) seq_len(nrow(df)) else n_obs + seq_len(nrow(df))

  if (!is.null(result$summary.fixed)) {
    fx <- data.table::as.data.table(result$summary.fixed, keep.rownames = "covariate")
    rename_quants(fx)
    fx[, `:=`(RR = exp(mean), RR_lower = exp(q025), RR_upper = exp(q975),
              significant = (q025 > 0) | (q975 < 0))]
    out$fixed_effects <- fx
  }

  if ("idx_spatial" %in% names(result$summary.random)) {
    sp <- data.table::as.data.table(result$summary.random$idx_spatial)
    names(sp)[1] <- "idx_spatial"; rename_quants(sp)
    sp <- merge(sp, grid_index, by = "idx_spatial", all.x = TRUE)
    sp[, RR_spatial := exp(mean)]
    out$spatial_effects <- sp
    out$n_areas <- length(unique(stats::na.omit(grid_index$idx_spatial)))
    out$spatial_model <- if (nrow(sp) >= 2L * out$n_areas) "bym2" else "besag"
  }

  if ("idx_temporal" %in% names(result$summary.random)) {
    tp <- data.table::as.data.table(result$summary.random$idx_temporal)
    names(tp)[1] <- "idx_temporal"; rename_quants(tp)
    tp[, RR_temporal := exp(mean)]
    out$temporal_effects <- tp[order(idx_temporal)]
  }

  if (!is.null(agg) && !is.null(result$summary.linear.predictor)) {
    lp <- data.table::as.data.table(result$summary.linear.predictor)[lat_rows]
    m <- lp$mean; s2 <- lp$sd^2
    ft <- data.table::data.table(
      fitted_mean = E_vec * exp(m + s2 / 2),
      fitted_sd   = E_vec * sqrt((exp(s2) - 1) * exp(2 * m + s2)),
      `fitted_0.025quant` = E_vec * exp(lp[["0.025quant"]]),
      `fitted_0.5quant`   = E_vec * exp(m),
      `fitted_0.975quant` = E_vec * exp(lp[["0.975quant"]]))
    ft[, `:=`(Y = df$Y, Y_observed = FALSE, E = E_vec, grid_id = df$grid_id, time_id = df$time_id,
              pearson_resid = (df$Y - fitted_mean) / sqrt(pmax(fitted_mean, 1e-6)),
              RR = exp(m), RR_lower = exp(lp[["0.025quant"]]), RR_upper = exp(lp[["0.975quant"]]),
              RR_exceed1   = 1 - pnorm(0,        mean = m, sd = lp$sd),
              RR_exceed1.5 = 1 - pnorm(log(1.5), mean = m, sd = lp$sd),
              RR_exceed2   = 1 - pnorm(log(2),   mean = m, sd = lp$sd),
              CI_width = exp(lp[["0.975quant"]]) - exp(lp[["0.025quant"]]))]
    out$fitted_values <- ft
    out$daily_outcome_note <- "Y on fitted_values is the interpolated daily series (diagnostic only); the likelihood used fitted_weekly"
    fw <- data.table::as.data.table(result$summary.fitted.values)[seq_len(n_obs)]
    data.table::setnames(fw, names(fw), paste0("fitted_", names(fw)))
    lpo <- data.table::as.data.table(result$summary.linear.predictor)[seq_len(n_obs)]
    fw[, `:=`(Y = agg$Y, E = agg$E, grid_id = agg$key$grid_id, week_id = agg$key$week_id,
              n_days = agg$key$n_days,
              pearson_resid = (agg$Y - fitted_mean) / sqrt(pmax(fitted_mean, 1e-6)),
              RR = exp(lpo$mean), RR_lower = exp(lpo[["0.025quant"]]), RR_upper = exp(lpo[["0.975quant"]]))]
    out$fitted_weekly <- fw
    out$observation_scale <- "unit-week"
  } else if (!is.null(result$summary.fitted.values)) {
    ft <- data.table::as.data.table(result$summary.fitted.values)
    data.table::setnames(ft, names(ft), paste0("fitted_", names(ft)))
    ft[, `:=`(Y = df$Y, E = E_vec, grid_id = df$grid_id, time_id = df$time_id,
              pearson_resid = (df$Y - fitted_mean) / sqrt(pmax(fitted_mean, 1e-6)))]
    if (!is.null(result$summary.linear.predictor)) {
      lp <- data.table::as.data.table(result$summary.linear.predictor)
      ft[, `:=`(RR = exp(lp$mean), RR_lower = exp(lp[["0.025quant"]]),
                RR_upper = exp(lp[["0.975quant"]]),
                RR_exceed1   = 1 - pnorm(0,          mean = lp$mean, sd = lp$sd),
                RR_exceed1.5 = 1 - pnorm(log(1.5),   mean = lp$mean, sd = lp$sd),
                RR_exceed2   = 1 - pnorm(log(2),     mean = lp$mean, sd = lp$sd),
                CI_width = exp(lp[["0.975quant"]]) - exp(lp[["0.025quant"]]))]
    }
    out$fitted_values <- ft
    out$observation_scale <- "unit-day"
  }

  if (!is.null(result$summary.hyperpar))
    out$hyperparameters <- data.table::as.data.table(result$summary.hyperpar, keep.rownames = "parameter")
  if ("log_rho_idx" %in% names(result$summary.random))
    out$rw1_curve <- data.table::as.data.table(result$summary.random$log_rho_idx)

  if (!is.null(result$cpo$cpo)) {
    out$cpo <- if (is.null(agg)) data.table::data.table(
      grid_id = df$grid_id, time_id = df$time_id,
      cpo = result$cpo$cpo, pit = result$cpo$pit, failure = result$cpo$failure)
    else data.table::data.table(
      grid_id = agg$key$grid_id, week_id = agg$key$week_id,
      cpo = result$cpo$cpo, pit = result$cpo$pit, failure = result$cpo$failure)
    out$lcpo <- -mean(log(pmax(result$cpo$cpo, 1e-300)), na.rm = TRUE)
  }

  if ("idx_st_cluster" %in% names(result$summary.random)) {
    st <- data.table::as.data.table(result$summary.random$idx_st_cluster)
    names(st)[1] <- "idx_st"; rename_quants(st)
    st[, RR_st := exp(mean)]
    out$st_interaction <- st
  } else if ("idx_st" %in% names(result$summary.random)) {
    st <- data.table::as.data.table(result$summary.random$idx_st)
    K <- length(unique(df$idx_st_cluster))
    st[, `:=`(idx_st = ((ID - 1L) %% K) + 1L, idx_time = ((ID - 1L) %/% K) + 1L)]
    st[, ID := NULL]; rename_quants(st)
    st[, RR_st := exp(mean)]
    out$st_interaction <- st
  }

  out
}

# rgeoda::skater() can crash the R session (GeoDaCenter/rgeoda#29), so it runs in a child process (R/skater.R).
run_skater_isolated <- function(nb_aug, co_df, n_grids, nc, grid_id, max_attempts = 4) {
  inp <- tempfile(fileext = ".rds"); out <- tempfile(fileext = ".rds"); log <- tempfile(fileext = ".log")
  on.exit(unlink(c(inp, out, log)), add = TRUE)
  attr(nb_aug, "region.id") <- as.character(seq_len(n_grids))
  saveRDS(list(nb_aug = nb_aug, co_df = co_df, n_grids = n_grids, nc = nc, grid_id = grid_id), inp)
  for (attempt in seq_len(max_attempts)) {
    status <- system2("Rscript", shQuote(c(file.path("R", "skater.R"), inp, out)), stdout = log, stderr = log)
    if (identical(status, 0L) && file.exists(out)) return(readRDS(out))
    log_msg(sprintf("SKATER attempt %d/%d failed:\n%s", attempt, max_attempts, paste(tail(readLines(log, warn = FALSE), 15), collapse = "\n")))
  }
  stop("rgeoda::skater() failed after ", max_attempts, " attempts")
}

# Spatial clusters for the Type IV interaction: SKATER on the unit graph, about 100 units per cluster
# and at most N_ST_CLUSTERS, for domains above 300 units; smaller domains use the unit graph itself.
# Writes the cluster graph and returns df with idx_st_cluster.
st_clusters <- function(df, grid_sf, nb, graph_file, cluster_graph) {
  n <- nrow(grid_sf)
  if (n <= 300) {
    file.copy(graph_file, cluster_graph, overwrite = TRUE)
    return(df[, idx_st_cluster := idx_spatial])
  }
  nc <- max(1L, min(N_ST_CLUSTERS, round(n / 100)))
  co <- as.data.frame(st_coordinates(suppressWarnings(st_centroid(grid_sf))))
  names(co) <- c("X", "Y")
  cl <- run_skater_isolated(nb, co, n, nc, grid_sf$grid_id)
  if ("idx_st_cluster" %in% names(df)) df[, idx_st_cluster := NULL]
  df <- merge(df, cl, by = "grid_id", all.x = TRUE)
  edges <- rbindlist(lapply(seq_along(nb), function(i) {
    ni <- nb[[i]][nb[[i]] > 0L]
    if (length(ni)) data.table(a = cl$idx_st_cluster[i], b = cl$idx_st_cluster[ni])
  }))
  edges <- if (nrow(edges)) unique(edges)[a != b] else data.table(a = integer(), b = integer())
  nb_cl <- lapply(seq_len(nc), function(i) {
    nbs <- sort(unique(edges[a == i, b]))
    if (length(nbs)) as.integer(nbs) else 0L
  })
  class(nb_cl) <- "nb"
  attr(nb_cl, "region.id") <- as.character(seq_len(nc))
  safe_nb2INLA(cluster_graph, nb_cl)
  df
}
