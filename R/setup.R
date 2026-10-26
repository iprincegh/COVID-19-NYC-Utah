# Packages, configuration and run helpers shared by every analysis script.
# Run all scripts from the repository root, e.g. Rscript analysis/10_fit_models.R

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(spdep)
  library(INLA)
  library(arrow)
  library(Matrix)
})

CONFIG <- list(
  data_dir    = "data",
  advan_dir   = file.path("data", "advan_weekly_plus"),
  output_dir  = "outputs",
  grids_dir   = file.path("outputs", "grids"),
  cache_dir   = file.path("outputs", "cache"),

  ny_cbg_path      = file.path("data", "cbg_shapefiles", "newyork_shp", "tl_2024_36_bg.shp"),
  ut_cbg_path      = file.path("data", "cbg_shapefiles", "utah_shp", "tl_2024_49_bg.shp"),
  ut_counties_path = file.path("data", "ut_counties", "CensusCounties2010.shp"),
  ny_modzcta_path  = file.path("data", "MODZCTA", "geo_export_5392c927-ee1e-4b9c-9787-2cad9456c83b.shp"),
  ny_modzcta_pop_path = file.path("data", "MODZCTA", "data-by-modzcta.csv"),

  nyc_caserate_path  = file.path("data", "us_CovidCases", "caseRate", "NYCcaserate-by-modzcta.csv"),
  nyc_crosswalk_path = file.path("data", "us_CovidCases", "caseRate", "NYCcbg_to_modzcta.csv"),
  nyt_counties_files = list(
    "2020" = file.path("data", "us_CovidCases", "us-counties-2020.csv"),
    "2021" = file.path("data", "us_CovidCases", "us-counties-2021.csv"),
    "2022" = file.path("data", "us_CovidCases", "us-counties-2022.csv")),

  age_pop_paths = list(
    "2020" = file.path("data", "demographics", "AgePop", "ACSST5Y2020.S0101-Data.csv"),
    "2021" = file.path("data", "demographics", "AgePop", "ACSST5Y2021.S0101-Data.csv"),
    "2022" = file.path("data", "demographics", "AgePop", "ACSST5Y2022.S0101-Data.csv")),
  income_empl_paths = list(
    "2020" = file.path("data", "demographics", "IncomeEmpl", "ACSDP5Y2020.DP03-Data.csv"),
    "2021" = file.path("data", "demographics", "IncomeEmpl", "ACSDP5Y2021.DP03-Data.csv"),
    "2022" = file.path("data", "demographics", "IncomeEmpl", "ACSDP5Y2022.DP03-Data.csv")),

  ny_landuse_path    = file.path("data", "landuse", "nyc_mappluto_25v4_shp", "MapPLUTO.shp"),
  ny_buildings_path  = file.path("data", "building_footprints", "NewYorkBuildings.geojson"),
  ut_parcels_dir     = file.path("data", "landuse", "utah_parcels_lir"),
  ut_buildings_path  = file.path("data", "building_footprints", "UtahBuildings.geojson"),
  pop_raster_dir     = file.path("data", "population"),

  ny_roads_shp        = file.path("data", "transport", "newyork_roads2", "tl_2023_36_prisecroads.shp"),
  ut_roads_dir        = file.path("data", "transport", "utah_roads_statewide"),
  bridge_csv          = file.path("data", "transport", "Bridge_Ratings_20260213.csv"),
  subway_stations_csv = file.path("data", "transport", "MTA_Subway_Stations_20260217.csv"),
  subway_routes_shp   = file.path("data", "transport", "nyu_subway", "routes_nyc_subway_may2016.shp"),

  grid_size_m = 100,
  ny_utm_crs  = 32618,
  ut_utm_crs  = 32612,

  nyc_fips_5 = c("36061", "36047", "36081", "36005", "36085"),
  nyc_fips_3 = c("061", "047", "081", "005", "085"),
  ut_state_fips = "49",
  ut_fips_5 = sprintf("49%03d", seq(1, 57, by = 2)),
  ut_fips_3 = sprintf("%03d", seq(1, 57, by = 2)),

  waves = list(
    wave1 = list(id = "wave1", months = c("2020-03", "2020-04", "2020-05"),
                 label = "Wave 1: Ancestral (Mar-May 2020)", acs_year = "2020"),
    wave2 = list(id = "wave2", months = c("2020-12", "2021-01", "2021-02"),
                 label = "Wave 2: Alpha (Dec 2020-Feb 2021)", acs_year = "2021"),
    wave3 = list(id = "wave3", months = c("2021-12", "2022-01", "2022-02"),
                 label = "Wave 3: Omicron (Dec 2021-Feb 2022)", acs_year = "2022")),

  spillover_alpha = 0.7,
  states = c("NY", "UT")
)

WAVE_RUNS <- rbindlist(lapply(CONFIG$states, function(s) rbindlist(lapply(CONFIG$waves, function(w)
  data.table(state = s, wave_id = w$id, wave_label = w$label,
             wave_months = list(w$months), acs_year = w$acs_year)))))
WAVE_RUNS[, run_id := paste0(state, "_", wave_id)]

# Optional command-line selection: Rscript analysis/<script>.R [NY|UT] [wave1|wave2|wave3]
select_runs <- function(args = commandArgs(trailingOnly = TRUE)) {
  runs <- WAVE_RUNS
  if (length(args) >= 1) runs <- runs[state %in% strsplit(args[1], ",")[[1]]]
  if (length(args) >= 2) runs <- runs[wave_id %in% strsplit(args[2], ",")[[1]]]
  stopifnot("no state/wave matches the arguments" = nrow(runs) > 0)
  runs
}

# Call fn(wave_id, wave_months) for each selected run, with `state` set in the caller.
run_waves_by_state <- function(fn, wave_runs = select_runs()) {
  results <- list()
  for (st in unique(wave_runs$state)) {
    st_waves <- wave_runs[state == st]
    assign("state", st, envir = parent.frame())
    for (i in seq_len(nrow(st_waves)))
      results[[st_waves$run_id[i]]] <- fn(st_waves$wave_id[i], st_waves$wave_months[[i]])
  }
  invisible(results)
}

get_wave_dir <- function(state, wave_id) file.path(CONFIG$output_dir, paste0(state, "_", wave_id))

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)

invisible(lapply(c(CONFIG$output_dir, CONFIG$grids_dir, CONFIG$cache_dir),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

for (f in c("data.R", "covariates.R", "adjacency.R", "model.R", "validation.R"))
  source(file.path("R", f))
