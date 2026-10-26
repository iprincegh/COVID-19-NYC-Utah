# Area covariates per cell: ACS age and employment shares (tract in NYC, county in Utah),
# residential density (MapPLUTO units per ha in NYC, building footprints per ha in Utah) and land use.
source("R/setup.R")

SF_TO_M2 <- 0.092903

read_acs_age <- function(path) {
  dt <- fread(path, header = TRUE)[grepl("^[0-9]+US", GEO_ID)]
  out <- data.table(tract_geoid = sub(".*US", "", dt$GEO_ID), total_pop = as.numeric(dt$S0101_C01_001E),
                    pct_older = as.numeric(dt$S0101_C01_030E) / pmax(as.numeric(dt$S0101_C01_001E), 1),
                    pct_younger = as.numeric(dt$S0101_C01_022E) / pmax(as.numeric(dt$S0101_C01_001E), 1))
  out <- out[nchar(tract_geoid) == 11]
  out[total_pop < 10, c("pct_older", "pct_younger") := NA_real_]
  out
}

read_acs_employment <- function(path) {
  dt <- fread(path, header = TRUE)[grepl("^[0-9]+US", GEO_ID)]
  out <- data.table(tract_geoid = sub(".*US", "", dt$GEO_ID), pop_16plus = as.numeric(dt$DP03_0001E),
                    pct_employed = as.numeric(dt$DP03_0004E) / pmax(as.numeric(dt$DP03_0001E), 1),
                    pct_unemployed = as.numeric(dt$DP03_0005E) / pmax(as.numeric(dt$DP03_0003E), 1))
  out <- out[nchar(tract_geoid) == 11]
  out[pop_16plus < 10, c("pct_employed", "pct_unemployed") := NA_real_]
  out
}

# NYC: residential units per ha, parcels split across cells by area
building_density_pluto <- function(grid_sf) {
  pluto <- st_read(CONFIG$ny_landuse_path, quiet = TRUE) |> st_transform(st_crs(grid_sf)) |> st_make_valid()
  for (cc in c("UnitsRes", "ResArea", "BldgArea")) {
    pluto[[cc]] <- as.numeric(pluto[[cc]])
    pluto[[cc]][is.na(pluto[[cc]])] <- 0
  }
  pluto <- pluto[(pluto$UnitsRes + pluto$ResArea + pluto$BldgArea) > 0, ]
  pluto$parcel_area_m2 <- as.numeric(st_area(pluto))
  pluto <- pluto[pluto$parcel_area_m2 > 0, ]
  inter <- suppressWarnings(st_intersection(st_set_precision(grid_sf, 1e3)[, "grid_id"],
                                            st_set_precision(pluto, 1e3)[, c("UnitsRes", "parcel_area_m2")]))
  inter <- suppressWarnings(st_collection_extract(inter, "POLYGON"))
  inter <- inter[!st_is_empty(inter), ]
  inter$frac <- as.numeric(st_area(inter)) / inter$parcel_area_m2
  agg <- as.data.table(st_drop_geometry(inter))[, .(units_res = sum(UnitsRes * frac, na.rm = TRUE)), by = grid_id]
  out <- merge(data.table(grid_id = grid_sf$grid_id, area_ha = as.numeric(st_area(grid_sf)) / 1e4), agg,
               by = "grid_id", all.x = TRUE)
  out[is.na(units_res), units_res := 0]
  out[, .(grid_id, building_density = units_res / area_ha)]
}

# Utah: building footprints (by centroid) per ha
building_density_footprints <- function(grid_sf) {
  centroids <- suppressWarnings(st_centroid(st_read(CONFIG$ut_buildings_path, quiet = TRUE) |> st_transform(st_crs(grid_sf))))
  data.table(grid_id = grid_sf$grid_id,
             building_density = lengths(st_intersects(grid_sf, centroids)) / (as.numeric(st_area(grid_sf)) / 1e4))
}

classify_lir_prop_class <- function(x) {
  s <- toupper(trimws(as.character(x)))
  cat <- rep("other", length(s))
  cat[grepl("MIXED", s)] <- "mixed"
  set_if <- function(hit, label) cat[hit & cat == "other"] <<- label
  set_if((s == "I") | grepl("INDUSTRIAL", s), "industrial")
  set_if((s == "C") | grepl("COMMERCIAL|RETAIL|OFFICE", s), "commercial")
  set_if((s == "R") | s %in% c("MH", "PUD", "CONDO") |
           grepl("RESIDENT|SINGLE FAMILY|SNGL FAM|MULTI[- ]?FAM|APARTM|TOWNHOUSE|MOBILE HOME|CONDO|PUD|^111|^115|^116|^119|^BR0|^LR0|^RES", s),
         "residential")
  set_if(grepl("AGRICULT|GREENBELT|FAA|FARM|GRAZ", s), "agricultural")
  set_if(grepl("VACANT|UNDEVELOP|^LAND$|LAND VACANT|LAND RES|LAND COM", s), "vacant")
  set_if(grepl("EXEMPT|TAX EXEMPT|GOVERN|RECREATION|RELIGIOUS|CHARITABLE", s), "exempt")
  cat[is.na(s) | s == "" | s == "DEFAULT" | s == "UNKNOWN"] <- "other"
  cat
}

# NYC: MapPLUTO land use at the cell centroid. Utah: dominant developed land use by parcel acreage.
landuse_categories <- function(state, grid_sf, year) {
  if (state == "NY") {
    lu <- st_read(CONFIG$ny_landuse_path, quiet = TRUE) |> st_transform(st_crs(grid_sf))
    lu_col <- grep("LandUse|landuse|LANDUSE", names(lu), value = TRUE)[1]
    j <- st_join(suppressWarnings(st_centroid(grid_sf)), lu[, lu_col], left = TRUE)
    out <- data.table(grid_id = j$grid_id, landuse = j[[lu_col]])
    out[, landuse_cat := fcase(landuse %in% c("01", "02", "03"), "residential", landuse %in% c("04", "05"), "commercial",
                               landuse %in% c("06", "07"), "industrial", landuse %in% c("08", "09"), "mixed",
                               landuse == "10", "public", landuse == "11", "vacant", default = "other")]
    return(out[, .(grid_id, landuse_cat)])
  }
  parcels <- st_read(get_landuse_path(state, year), quiet = TRUE) |> st_transform(st_crs(grid_sf))
  pdt <- as.data.table(st_drop_geometry(parcels))
  pdt[, landuse_cat := classify_lir_prop_class(PROP_CLASS)]
  pdt[, acres := suppressWarnings(as.numeric(PARCEL_ACRES))]
  pdt[is.na(acres) | acres < 0, acres := 0]
  pdt[, grid_id := st_join(suppressWarnings(st_centroid(parcels))[, c()], grid_sf[, "grid_id"], join = st_within, left = TRUE)$grid_id]
  agg <- pdt[!is.na(grid_id), .(acres = sum(acres, na.rm = TRUE)), by = .(grid_id, landuse_cat)]
  dev <- agg[landuse_cat %in% c("residential", "commercial", "industrial", "mixed")]
  dom <- dev[dev[, .I[which.max(acres)], by = grid_id]$V1, .(grid_id, landuse_cat)]
  out <- merge(data.table(grid_id = grid_sf$grid_id), dom, by = "grid_id", all.x = TRUE)
  out[is.na(landuse_cat), landuse_cat := "other"]
  out
}

runs <- select_runs()
age <- lapply(CONFIG$age_pop_paths, read_acs_age)
empl <- lapply(CONFIG$income_empl_paths, read_acs_employment)

for (st in unique(runs$state)) {
  grid_sf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", st)), quiet = TRUE) |> st_transform(get_utm_crs(st))
  grid_dt <- as.data.table(st_drop_geometry(grid_sf))[, .(grid_id, fips_5,
    tract_geoid = if (st == "NY") substr(as.character(cbg_geoid), 1, 11) else NA_character_)]
  bld <- if (st == "NY") building_density_pluto(grid_sf) else building_density_footprints(grid_sf)
  lu_by_year <- list()

  for (i in which(runs$state == st)) {
    year <- runs$acs_year[i]
    demo <- copy(grid_dt)
    if (st == "UT") {
      a <- age[[year]][, .(pct_older = weighted.mean(pct_older, total_pop, na.rm = TRUE),
                           pct_younger = weighted.mean(pct_younger, total_pop, na.rm = TRUE)), by = .(fips_5 = substr(tract_geoid, 1, 5))]
      e <- empl[[year]][, .(pct_employed = weighted.mean(pct_employed, pop_16plus, na.rm = TRUE),
                            pct_unemployed = weighted.mean(pct_unemployed, pop_16plus, na.rm = TRUE)), by = .(fips_5 = substr(tract_geoid, 1, 5))]
      demo <- merge(merge(demo, a, by = "fips_5", all.x = TRUE), e, by = "fips_5", all.x = TRUE)
    } else {
      demo <- merge(demo, age[[year]][, .(tract_geoid, pct_older, pct_younger)], by = "tract_geoid", all.x = TRUE)
      demo <- merge(demo, empl[[year]][, .(tract_geoid, pct_employed, pct_unemployed)], by = "tract_geoid", all.x = TRUE)
    }
    demo <- merge(demo, bld, by = "grid_id", all.x = TRUE)
    lu_key <- if (st == "NY") "NY" else year
    if (is.null(lu_by_year[[lu_key]])) lu_by_year[[lu_key]] <- landuse_categories(st, grid_sf, year)
    demo <- merge(demo, lu_by_year[[lu_key]], by = "grid_id", all.x = TRUE)

    for (col in c("pct_older", "pct_younger", "pct_employed", "pct_unemployed", "building_density")) {
      med <- median(demo[[col]], na.rm = TRUE)
      demo[is.na(get(col)), (col) := if (is.na(med)) 0 else med]
    }
    demo[is.na(landuse_cat), landuse_cat := "other"]
    demo[, `:=`(landuse_residential = as.integer(landuse_cat == "residential"),
                landuse_commercial = as.integer(landuse_cat == "commercial"),
                landuse_industrial = as.integer(landuse_cat == "industrial"))]
    fwrite(demo, file.path(get_wave_dir(st, runs$wave_id[i]), "tables", "demographics.csv"))
  }
}
