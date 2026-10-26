# Descriptive diagnostics quoted in the Results: share of contact-weighted visits made from within the
# destination's reporting area, between- and within-MODZCTA variation of covariates and of the NYC
# spatial effect, and how far NYC cells on a MODZCTA border are pulled towards the neighbouring MODZCTA.
source("R/setup.R")
out <- file.path(CONFIG$output_dir, "diagnostics")
dir.create(out, showWarnings = FALSE)

# within-area share of visit activity, per arm
fwrite(rbindlist(lapply(WAVE_RUNS$run_id, function(arm) {
  d <- fread(file.path(CONFIG$output_dir, arm, "tables", "od_kernel_rho.csv"), select = c("internal_activity", "total_activity"))
  data.table(arm = arm, within_area_share = sum(d$internal_activity) / sum(d$total_activity))
})), file.path(out, "within_area_visit_share.csv"))

# percentage of variance within groups
within_pct <- function(x, g) {
  ok <- is.finite(x); x <- x[ok]; g <- g[ok]
  100 * mean((x - ave(x, g))^2) * length(x) / (length(x) - 1) / var(x)
}
spatial_effect <- function(arm) {
  sp <- fread(file.path(CONFIG$output_dir, arm, "tables", "unified_spatial_effects.csv"))
  sp[idx_spatial <= nrow(sp) %/% 2L, .(grid_id, b = mean)]   # first half: the combined BYM2 effect
}

nyc <- rbindlist(lapply(c("wave1", "wave2", "wave3"), function(w) {
  arm <- paste0("NY_", w)
  gt <- as.data.table(readRDS(file.path(get_wave_dir("NY", w), "models", "model_data_gt.rds"))$gt_dt)[geo_id != "99999"]
  cells <- gt[, .(geo_id = geo_id[1], z = mean(log_rho)), by = grid_id]
  cells <- merge(cells, spatial_effect(arm), by = "grid_id")
  gt[, `:=`(own_w = own_prev_std - ave(own_prev_std, geo_id, time_id), ext_w = ext_prev_std - ave(ext_prev_std, geo_id, time_id))]

  # border cells: within-MODZCTA deviation of b against the gap to neighbouring MODZCTAs
  g <- INLA::inla.read.graph(file.path(get_wave_dir("NY", w), "models", "grid_adj_unified.graph"))
  cz <- merge(unique(gt[, .(grid_id, geo_id, idx = idx_spatial)]), cells[, .(grid_id, b)], by = "grid_id")
  setkey(cz, idx)
  zone <- cz[J(seq_len(g$n)), geo_id]; b <- cz[J(seq_len(g$n)), b]
  zmean <- tapply(b, zone, mean, na.rm = TRUE)
  e <- rbindlist(lapply(seq_len(g$n), function(i) {
    other <- g$nbs[[i]][!is.na(zone[g$nbs[[i]]]) & zone[g$nbs[[i]]] != zone[i]]
    if (length(other) && !is.na(zone[i]))
      data.table(bw = b[i] - zmean[[zone[i]]], pull = mean(zmean[zone[other]]) - zmean[[zone[i]]])
  }))

  data.table(arm = arm, cells = nrow(cells), zones = uniqueN(cells$geo_id),
             unit_days_with_visits = mean(gt$total_activity > 0),
             within_pct_spatial_effect = within_pct(cells$b, cells$geo_id),
             within_pct_activity = within_pct(cells$z, cells$geo_id),
             within_pct_local_pressure = 100 * var(gt$own_w) / var(gt$own_prev_std),
             within_pct_imported_exposure = 100 * var(gt$ext_w) / var(gt$ext_prev_std),
             border_pull_slope = coef(lm(bw ~ pull, data = e))[[2]])
}))
fwrite(nyc, file.path(out, "nyc_within_modzcta_variation.csv"))
