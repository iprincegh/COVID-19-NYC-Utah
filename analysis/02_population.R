# Population offsets: WorldPop 100 m sums for NYC grid cells, ACS county totals for Utah,
# and the NYC Health MODZCTA denominators used to convert weekly rates to counts.
source("R/setup.R")
library(terra)
library(exactextractr)

study_years <- unique(WAVE_RUNS$acs_year)

# NYC grid cells: sum of the year-matched WorldPop constrained 100 m raster
grid_ny <- st_read(file.path(CONFIG$grids_dir, "grid_100m_NY.gpkg"), quiet = TRUE)
for (year in study_years) {
  r <- terra::rast(file.path(CONFIG$pop_raster_dir, sprintf("usa_pop_%s_CN_100m_R2025A_v1.tif", year)))
  v <- exactextractr::exact_extract(r, st_transform(grid_ny, terra::crs(r)), "sum", progress = FALSE)
  v[is.na(v)] <- 0
  pop <- round(v, 2)
  fwrite(data.table(grid_id = grid_ny$grid_id, population = pop, E = pmax(pop, 0.1), log_E = log(pmax(pop, 0.1))),
         file.path(CONFIG$grids_dir, sprintf("population_NY_%s.csv", year)))
}

# Utah counties: ACS 5-year total population (S0101_C01_001E), tracts summed to counties
grid_ut <- st_read(file.path(CONFIG$grids_dir, "grid_100m_UT.gpkg"), quiet = TRUE)
key <- data.table(grid_id = grid_ut$grid_id, fips_5 = as.character(grid_ut$fips_5))
for (year in study_years) {
  a <- fread(CONFIG$age_pop_paths[[year]])
  a <- a[as.integer(state) == 49L & !is.na(suppressWarnings(as.numeric(S0101_C01_001E)))]
  a[, fips_5 := paste0(sprintf("%02d", as.integer(state)), sprintf("%03d", as.integer(county)))]
  cty <- a[, .(population = sum(as.numeric(S0101_C01_001E))), by = fips_5]
  dt <- merge(key, cty, by = "fips_5", all.x = TRUE)[order(grid_id)]
  stopifnot("every Utah county needs an ACS total" = !anyNA(dt$population))
  fwrite(dt[, .(grid_id, population = round(population, 2), E = pmax(population, 0.1), log_E = log(pmax(population, 0.1)))],
         file.path(CONFIG$grids_dir, sprintf("population_UT_%s.csv", year)))
}

# NYC MODZCTA: NYC Health population denominators (data-by-modzcta.csv)
zones <- st_read(file.path(CONFIG$grids_dir, "grid_100m_NY_modzcta.gpkg"), quiet = TRUE)
pop_raw <- fread(CONFIG$ny_modzcta_pop_path)[, modzcta := standardize_fips(MODIFIED_ZCTA)]
merged <- merge(as.data.table(st_drop_geometry(zones))[, .(grid_id, modzcta)],
                pop_raw[, .(modzcta, POP_DENOMINATOR)], by = "modzcta", all.x = TRUE)
stopifnot(!anyNA(merged$POP_DENOMINATOR), nrow(merged) == nrow(zones))
setorder(merged, grid_id)
pop_dt <- merged[, .(grid_id, population = POP_DENOMINATOR, E = pmax(POP_DENOMINATOR, 0.1), log_E = log(pmax(POP_DENOMINATOR, 0.1)))]
for (year in study_years)
  fwrite(pop_dt, file.path(CONFIG$grids_dir, sprintf("population_NY_modzcta_%s.csv", year)))
