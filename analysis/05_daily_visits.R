# Daily total POI visits and weekend indicator per wave (used for Figure 2).
source("R/setup.R")

process_daily <- function(wave_id, wave_months) {
  out_dir <- get_wave_dir(state, wave_id)
  dates <- unique(fread(file.path(out_dir, "tables", "case_rates.csv"))[, .(date = as.Date(date))])
  setorder(dates, date)
  daily <- fread(file.path(out_dir, "intermediate", "daily_stats.csv"))[, .(date = as.Date(date), total_visits, median_dwell)]
  covs <- merge(dates, daily, by = "date", all.x = TRUE)
  covs[, sim_day := as.integer(date - min(date)) + 1L]
  covs[, is_weekend := as.integer(weekdays(date) %in% c("Saturday", "Sunday"))]
  fwrite(covs[, .(date, sim_day, total_visits, median_dwell, is_weekend)],
         file.path(out_dir, "tables", "temporal_covariates.csv"))
}

invisible(run_waves_by_state(process_daily))
