# Cell-level visit covariates: SES contrast (iOS share of devices), repeat-visit ratio,
# hourly contact density and daytime visitor population.
source("R/setup.R")

process_covariates <- function(wave_id, wave_months) {
  out_dir <- get_wave_dir(state, wave_id)
  utm <- get_utm_crs(state)
  grid_dt <- as.data.table(read_feather(file.path(out_dir, "intermediate", "grid_visits.feather")))
  grid_dt[is.na(avg_repeat_ratio), avg_repeat_ratio := 1]
  grid_dt[is.na(avg_contact_density_hourly), avg_contact_density_hourly := 0]
  grid_dt[, `:=`(ses_contrast = visits_ios / pmax(visits_ios + visits_android, 1) - 0.5,
                 log_repeat_ratio = log1p(avg_repeat_ratio - 1),
                 log_contact_density_hourly = log1p(avg_contact_density_hourly))]
  grid_dt[is.na(ses_contrast), ses_contrast := 0]

  daytime <- as.data.table(read_feather(file.path(out_dir, "intermediate", "daytime_visitors.feather")))
  daytime <- st_as_sf(daytime[!is.na(poi_lon) & !is.na(poi_lat)], coords = c("poi_lon", "poi_lat"), crs = 4326) |>
    st_transform(utm)
  gsf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), quiet = TRUE) |> st_transform(utm)
  daytime <- as.data.table(st_drop_geometry(st_join(daytime, gsf[, "grid_id"], left = FALSE)))[,
    .(daytime_visitors = sum(daytime_visitors, na.rm = TRUE)), by = grid_id]
  grid_dt <- merge(grid_dt, daytime, by = "grid_id", all.x = TRUE)
  grid_dt[is.na(daytime_visitors), daytime_visitors := 0]
  grid_dt[, log_daytime_pop := log1p(daytime_visitors)]

  write_feather(grid_dt[, .(grid_id, ses_contrast, log_repeat_ratio, log_contact_density_hourly, log_daytime_pop)],
                file.path(out_dir, "intermediate", "spatial_covariates.feather"))
}

invisible(run_waves_by_state(process_covariates))
