# Negative control: refit the selected model with same-day residual mobility seven days ahead, after
# removing the part predicted by current and past mobility. Future mobility cannot cause current cases,
# so a clear association points to confounding.
source("R/setup.R")
inla.setOption(num.threads = Sys.getenv("INLA_THREADS", "4:6"), inla.mode = "compact")

arms <- select_runs()$run_id
dir.create(file.path(CONFIG$output_dir, "negative_control"), showWarnings = FALSE)
for (arm in arms) {
  wv <- load_wave_for_cv(arm)
  md <- wv$md
  fm <- update(md$fms[[wv$best_model]], ~ . + mob_lag_neg)
  fit <- fit_model(fm, md$df, get_E(wv$best_model, md$df), model_name = paste0(arm, "_placebo"),
                   zero_frac = md$zero_frac, dispersion = md$dispersion, config = FALSE, n_grids = md$n_grids,
                   state = wv$state, wave_id = wv$wave_id, force_family = md$best_model_family, agg = md$agg)
  stopifnot(!is.null(fit$result))
  fe <- as.data.table(fit$result$summary.fixed, keep.rownames = "covariate")
  setnames(fe, c("0.025quant", "0.5quant", "0.975quant"), c("q025", "q50", "q975"), skip_absent = TRUE)
  fe[, `:=`(wave = arm, family = md$best_model_family)]
  fwrite(fe, file.path(CONFIG$output_dir, "negative_control", paste0(arm, "_fixed_effects.csv")))
}
