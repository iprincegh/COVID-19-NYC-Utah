# Spatial units: the NYC 100 m grid, Utah counties, and the NYC MODZCTA zones used for case allocation.
source("R/setup.R")

for (state in CONFIG$states) {
  utm_crs <- get_utm_crs(state)
  fips_3  <- get_metro_fips(state, 3)

  if (state == "UT") {
    counties_sf <- st_transform(st_read(CONFIG$ut_counties_path, quiet = TRUE), utm_crs)
    counties_sf$COUNTYFP10 <- sprintf("%03d", as.integer(counties_sf$FIPS))
    counties_sf <- counties_sf[counties_sf$COUNTYFP10 %in% fips_3, ]
    counties_sf <- counties_sf[order(counties_sf$COUNTYFP10), ]
    grid_sf <- st_sf(
      grid_id     = seq_len(nrow(counties_sf)),
      county_fips = counties_sf$COUNTYFP10,
      fips_5      = paste0(get_state_prefix(state), counties_sf$COUNTYFP10),
      cbg_geoid   = NA_character_,
      county_name = counties_sf$NAME,
      geometry    = st_geometry(counties_sf))
    cc <- st_coordinates(st_centroid(grid_sf))
    grid_sf$cent_x <- cc[, 1]
    grid_sf$cent_y <- cc[, 2]
  } else {
    cbg_sf <- st_read(get_cbg_path(state), quiet = TRUE)
    metro_cbg <- st_transform(cbg_sf[substr(cbg_sf$GEOID, 3, 5) %in% fips_3, ], utm_crs)

    grid_sf <- st_sf(geometry = st_make_grid(metro_cbg, cellsize = CONFIG$grid_size_m, what = "polygons"))
    coords <- st_coordinates(st_centroid(grid_sf))
    grid_sf$cent_x <- coords[, 1]
    grid_sf$cent_y <- coords[, 2]
    grid_sf <- grid_sf[st_intersects(grid_sf, st_union(metro_cbg), sparse = FALSE)[, 1], ]
    grid_sf$grid_id <- seq_len(nrow(grid_sf))

    cbg_assign <- st_within(st_centroid(grid_sf), metro_cbg, sparse = TRUE)
    grid_sf$cbg_geoid   <- vapply(cbg_assign, function(x) if (length(x)) metro_cbg$GEOID[x[1]] else NA_character_, "")
    grid_sf$county_fips <- substr(grid_sf$cbg_geoid, 3, 5)
    grid_sf$fips_5      <- paste0(get_state_prefix(state), grid_sf$county_fips)
  }
  st_write(grid_sf, file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)), delete_dsn = TRUE, quiet = TRUE)
}

# MODZCTA zones (the 99999 placeholder is not a reporting area)
modzcta_sf <- st_read(CONFIG$ny_modzcta_path, quiet = TRUE)
modzcta_sf <- modzcta_sf[modzcta_sf$modzcta != "99999", ]
stopifnot(nrow(modzcta_sf) == 177)
modzcta_sf <- st_transform(modzcta_sf, get_utm_crs("NY"))
modzcta_sf <- modzcta_sf[order(modzcta_sf$modzcta), ]
zones <- st_sf(
  grid_id     = seq_len(nrow(modzcta_sf)),
  modzcta     = standardize_fips(modzcta_sf$modzcta),
  rate_geo_id = standardize_fips(modzcta_sf$modzcta),
  label       = modzcta_sf$label,
  geometry    = st_geometry(modzcta_sf))
st_write(zones, file.path(CONFIG$grids_dir, "grid_100m_NY_modzcta.gpkg"), delete_dsn = TRUE, quiet = TRUE)
