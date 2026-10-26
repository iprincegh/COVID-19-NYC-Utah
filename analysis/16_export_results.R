# Collect the results behind the journal figures and tables into results/, which make_figures.R reads.
# Everything written here is a model output or derives from public case, population and boundary data;
# nothing derives from the Advan visit records except fitted model quantities.
source("R/setup.R")
res_dir <- "results"
dir.create(res_dir, showWarnings = FALSE)
# Names in results/ follow the paper: models as in Table 2, mobility-derived covariates by their symbols
# (rho_inf, rho_imp, z_resid, L, Wz_bar) and the other covariates and hyperparameters by their table names.
MODEL_NAME <- c(Base = "Base", Temp = "Non-spatial", LinMob = "Linear-mobility", NLMob = "Nonlinear-mobility",
                LinMob_ST_IV = "Linear + interaction", NLMob_ST_IV = "Nonlinear + interaction")
COV_NAME <- c(own_prev_std = "rho_inf", ext_prev_std = "rho_imp", log_rho_resid = "z_resid", mob_weighted = "L",
              mob_spatial_lag_std = "Wz_bar", `(Intercept)` = "intercept", building_density_std = "building_density",
              landuse_commercial = "commercial_land_use", commercial_x_weekend = "commercial_land_use_x_weekend",
              log_repeat_ratio_std = "repeat_visit_ratio", log_contact_density_hourly_std = "contact_density",
              log_daytime_pop_std = "daytime_population_density", ses_contrast_std = "ses_contrast",
              pct_younger_std = "younger_share", pct_older_std = "older_share", pct_unemployed_std = "unemployment_rate",
              pct_employed_std = "employment_rate", is_weekend = "weekend", is_holiday = "holiday")
HYPER_NAME <- c(`Precision for idx_spatial` = "BYM2 precision", `Phi for idx_spatial` = "BYM2 mixing (phi)",
                `Precision for idx_temporal` = "Temporal precision", `Precision for idx_st` = "Interaction precision",
                `Precision for idx_st_cluster` = "Interaction precision",
                `Precision for log_rho_idx` = "Nonlinear activity effect precision",
                `size for the nbinomial observations (1/overdispersion)` = "Negative binomial size")
readable <- function(x) {
  x <- as.data.table(x)
  for (m in list(list("model", MODEL_NAME), list("covariate", COV_NAME), list("term", COV_NAME), list("parameter", HYPER_NAME)))
    if (m[[1]] %in% names(x)) x[, (m[[1]]) := fifelse(get(m[[1]]) %in% names(m[[2]]), unname(m[[2]][get(m[[1]])]), get(m[[1]]))]
  x
}
out <- function(x, name) fwrite(readable(x), file.path(res_dir, name))
arm_of <- function(x, arm) cbind(data.table(setting = sub("_.*", "", arm), wave = as.integer(sub(".*wave", "", arm))), x)
per_arm <- function(f) rbindlist(lapply(WAVE_RUNS$run_id, function(arm) arm_of(f(arm, get_wave_dir(sub("_.*", "", arm), sub(".*_", "", arm))), arm)), fill = TRUE)
best <- function(dir) readRDS(file.path(dir, "models", "best_model.rds"))
# Utah counts as fitted: the masked reporting artefacts are missing
counts <- function(arm, md) {
  d <- as.data.table(md$gt_dt)
  m <- MASK_ARTIFACT_DAYS[[arm]]
  if (!is.null(m)) d[d[m, on = c("grid_id", "time_id"), which = TRUE], Y := NA_integer_]
  d
}

# model selection, fixed effects and hyperparameters
out(per_arm(function(arm, dir) fread(file.path(dir, "tables", "model_comparison_unified.csv"))[Status == "OK",
  .(model = Model, WAIC, DIC, pD, RMSE, Correlation, Family)]), "model_comparison.csv")
bm <- lapply(setNames(nm = WAVE_RUNS$run_id), function(arm) best(get_wave_dir(sub("_.*", "", arm), sub(".*_", "", arm)))$outputs[
  c("fixed_effects", "hyperparameters", "spatial_effects", "temporal_effects", "st_interaction", "rw1_curve")])
models <- per_arm(function(arm, dir) { b <- best(dir); data.table(model = b$name, dispersion = b$dispersion) })
out(models, "selected_models.csv")
stack <- function(part, cols) rbindlist(lapply(names(bm), function(arm) if (!is.null(bm[[arm]][[part]]))
  arm_of(as.data.table(bm[[arm]][[part]])[, intersect(cols, names(bm[[arm]][[part]])), with = FALSE], arm)), fill = TRUE)
out(stack("fixed_effects", c("covariate", "mean", "sd", "q025", "q975", "RR", "RR_lower", "RR_upper")), "fixed_effects.csv")
out(stack("hyperparameters", c("parameter", "mean", "sd", "0.025quant", "0.5quant", "0.975quant")), "hyperparameters.csv")
out(stack("temporal_effects", c("idx_temporal", "mean", "sd", "q025", "q975")), "temporal_effects.csv")
out(stack("st_interaction", c("idx_st", "mean", "q025", "q975")), "st_interaction.csv.gz")
sp <- rbindlist(lapply(names(bm), function(arm) {
  s <- as.data.table(bm[[arm]]$spatial_effects)
  arm_of(s[idx_spatial <= nrow(s) %/% 2L, .(grid_id, mean, sd, q025, q975)], arm)   # first half: the BYM2 effect
}))
out(sp, "spatial_effects.csv.gz")

# nonlinear activity effect, at the centres of the 20 equal-width bins of each wave's activity range
out(rbindlist(lapply(names(bm), function(arm) {
  cv <- bm[[arm]]$rw1_curve
  if (is.null(cv)) return(NULL)
  z <- readRDS(file.path(get_wave_dir(sub("_.*", "", arm), sub(".*_", "", arm)), "models", "model_data_gt.rds"))$gt_dt$log_rho
  brks <- seq(min(z, na.rm = TRUE) - 0.01, max(z, na.rm = TRUE) + 0.01, length.out = N_RHO_BINS + 1)
  arm_of(as.data.table(cv)[, .(bin = ID, z = (brks[ID] + brks[ID + 1]) / 2, mean, q025 = `0.025quant`, q975 = `0.975quant`)], arm)
})), "activity_effect.csv")

# cases: weekly totals per setting and wave, NYC allocation check and the cells of one MODZCTA
weekly <- rbindlist(lapply(WAVE_RUNS$run_id, function(arm) {
  st <- sub("_.*", "", arm); dir <- get_wave_dir(st, sub(".*_", "", arm))
  md <- readRDS(file.path(dir, "models", "model_data_gt.rds"))
  if (st == "NY") {
    wl <- md$weekly_obs$week_lookup[observed == TRUE]
    y <- md$weekly_obs$Y_week_cell[, .(cases = sum(Y_week)), by = week_id]
    d <- merge(y, wl[, .(week_id, day = as.numeric(week_start - md$date_range$start) + 4, n_days = 7L)], by = "week_id")
  } else {
    d <- counts(arm, md)[, .(cases = sum(Y, na.rm = TRUE)), by = time_id][, week_id := (time_id - 1L) %/% 7L + 1L]
    d <- d[, .(cases = sum(cases), n_days = .N, day = mean(time_id)), by = week_id][n_days >= 4L]
  }
  arm_of(d[, .(week = week_id, day, cases_per_week = 7 * cases / n_days)], arm)
}))
out(weekly, "weekly_cases.csv")
out(per_arm(function(arm, dir) {
  cr <- fread(file.path(dir, "tables", "case_rates.csv"))[, .(rate = sum(case_rate, na.rm = TRUE)), by = date]
  cr[, .(day = as.integer(as.Date(date) - min(as.Date(date))) + 1L, rate)][order(day)]
}), "case_rate_daily_total.csv")
out(rbindlist(lapply(grep("^UT", WAVE_RUNS$run_id, value = TRUE), function(arm)
  arm_of(counts(arm, readRDS(file.path(get_wave_dir("UT", sub(".*_", "", arm)), "models", "model_data_gt.rds")))[,
    .(cases = sum(Y)), by = .(day = time_id)], arm))), "utah_daily_cases.csv")

ny1 <- readRDS(file.path(get_wave_dir("NY", "wave1"), "models", "model_data_gt.rds"))
rep <- fread(file.path(get_wave_dir("NY", "wave1"), "tables", "case_rates_weekly.csv"), colClasses = list(character = "modzcta"))
rt <- merge(ny1$weekly_obs$Y_week_cell[, .(Y_sum = sum(Y_week)), by = .(geo_id, week_id)],
            ny1$weekly_obs$Y_week[, .(geo_id, week_id, week_ending_date, pop_zone)], by = c("geo_id", "week_id"))
rt <- merge(rt, rep[, .(geo_id = modzcta, week_ending_date = as.Date(week_ending_date), rate_reported = case_rate, suppressed)],
            by = c("geo_id", "week_ending_date"))[suppressed == FALSE & rate_reported > 0]
out(rt[, .(geo_id, week_ending_date, rate_reported, rate_allocated = Y_sum / pop_zone * 1e5)], "allocation_check.csv")
# cells fitted in any NYC wave, with the Wave 1 population weight and allocated cases (Figure 3)
cells <- unique(rbindlist(lapply(c("wave2", "wave3"), function(w)
  unique(readRDS(file.path(get_wave_dir("NY", w), "models", "model_data_gt.rds"))$gt_dt[, .(grid_id, geo_id)]))))
w1 <- merge(unique(as.data.table(ny1$gt_dt)[, .(grid_id, geo_id, E)]),
            ny1$weekly_obs$Y_week_cell[, .(cases_wave1 = sum(Y_week)), by = grid_id], by = "grid_id", all.x = TRUE)
cells <- rbind(w1, cells[!grid_id %in% w1$grid_id], fill = TRUE)
xy <- st_read(file.path(CONFIG$grids_dir, "grid_100m_NY.gpkg"), quiet = TRUE)
xy <- data.table(grid_id = xy$grid_id, round(st_coordinates(st_centroid(st_geometry(xy))), 1))
out(merge(cells, xy, by = "grid_id")[, .(grid_id, x = X, y = Y, geo_id, E, cases_wave1)], "nyc_grid_cells.csv.gz")
rm(ny1); gc()

# public boundaries used by the maps (UTM 18N for NYC, UTM 12N for Utah); MODZCTA 99999 is parks and airports
mz <- st_read(CONFIG$ny_modzcta_path, quiet = TRUE) |> st_transform(get_utm_crs("NY"))
st_write(st_sf(geo_id = as.character(mz$modzcta), geometry = st_geometry(mz)), file.path(res_dir, "boundaries.gpkg"),
         layer = "modzcta", delete_dsn = TRUE, quiet = TRUE)
ut <- st_read(file.path(CONFIG$grids_dir, "grid_100m_UT.gpkg"), quiet = TRUE)
st_write(st_sf(grid_id = ut$grid_id, geo_id = ut$fips_5, name = ut$county_name, geometry = st_geometry(ut)),
         file.path(res_dir, "boundaries.gpkg"), layer = "utah_counties", quiet = TRUE)

# temporal cross-validation
tcv <- lapply(setNames(nm = WAVE_RUNS$run_id), function(arm)
  readRDS(file.path(get_wave_dir(sub("_.*", "", arm), sub(".*_", "", arm)), "validation", "temporal_cv_results.rds")))
tstack <- function(f) rbindlist(lapply(names(tcv), function(arm) { x <- f(tcv[[arm]], arm); if (!is.null(x)) arm_of(x, arm) }), fill = TRUE)
out(tstack(function(r, arm) {
  s <- as.data.table(r$cell_scores)[1, .(model, MdAPE, CRPS, Cov95)]
  if (!is.null(r$modz_scores)) s <- cbind(s, as.data.table(r$modz_scores)[, .(modzcta_total_MdAPE, modzcta_total_Corr, modzcta_total_Cov95)])
  s
}), "temporal_cv_scores.csv")
out(tstack(function(r, arm) as.data.table(r$fold_diag)[, .(origin, family, nb_size)]), "temporal_cv_folds.csv")
out(tstack(function(r, arm) as.data.table(r$holdout_long)[, .(origin, y_true, mu)]), "temporal_cv_holdout.csv.gz")
out(tstack(function(r, arm) as.data.table(r$daily_raw)[, .(origin, time_id, obs, pred_med, lo_draw, hi_draw)]), "temporal_cv_daily.csv")
out(tstack(function(r, arm) if (startsWith(arm, "NY")) as.data.table(r$modz_win)[, .(origin, unit = geo_id, obs, pred_med)]
           else as.data.table(r$unit_window_long)[, .(origin, unit = as.character(idx_spatial), obs, pred_med)]),
    "temporal_cv_windows.csv")

# observed counts for the naive benchmark: NYC MODZCTA-weeks and Utah county-days
out(per_arm(function(arm, dir) {
  md <- readRDS(file.path(dir, "models", "model_data_gt.rds"))
  if (startsWith(arm, "NY")) {
    wl <- md$weekly_obs$week_lookup
    merge(md$weekly_obs$Y_week[, .(unit = geo_id, week_id, y = Y_week)],
          wl[, .(week_id, last_day = as.integer(week_ending_date - md$date_range$start) + 1L)], by = "week_id")
  } else counts(arm, md)[, .(unit = as.character(idx_spatial), last_day = time_id, y = Y)]
}), "observed_counts.csv.gz")

# spatial cross-validation and out-of-fold PCAtM
cv_file <- function(name) per_arm(function(arm, dir) fread(file.path(dir, "validation", name), colClasses = list(character = "geo_id")))
out(cv_file("spatial_cv_area_detail.csv")[, .(setting, wave, geo_id, fold, obs, pred_med, lo_draw, hi_draw, covered_draw)], "spatial_cv_areas.csv")
out(cv_file("spatial_cv_primary_metrics.csv")[, .(setting, wave, n_areas = n_areas_scored_ok, MdAPE = area_MdAPE_ok, R2 = area_R2_ok,
                                                  Cov95 = area_Cov95_draw_ok)], "spatial_cv_metrics.csv")
out(cv_file("pcatm_pooled.csv")[, .(setting, wave, pcatm = pcatm_pooled_mean, q05 = pcatm_pooled_q05, q95 = pcatm_pooled_q95)], "pcatm_holdout.csv")

# in-sample PCAtM and the negative control
pc <- function(name) per_arm(function(arm, dir) fread(file.path(CONFIG$output_dir, "pcatm", arm, name))[, wave := NULL])
out(pc("pcatm_wave_summary.csv")[, .(setting, wave, model, pcatm, pcatm_q025, pcatm_q975, amlp)], "pcatm_summary.csv")
out(pc("pcatm_covariate_decomposition.csv")[, .(setting, wave, term = component, pcatm, pcatm_q025, pcatm_q975, amlp)], "pcatm_terms.csv")
out(pc("pcatm_category_decomposition.csv")[, .(setting, wave, term = component, category = category_broad, amlp, term_total = component_total)],
    "pcatm_categories.csv")
out(per_arm(function(arm, dir) fread(file.path(CONFIG$output_dir, "negative_control", paste0(arm, "_fixed_effects.csv")))[
  covariate == "mob_lag_neg", .(mean, q025, q975, family)]), "negative_control.csv")
