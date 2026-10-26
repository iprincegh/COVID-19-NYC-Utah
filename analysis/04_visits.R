# Advan Weekly Patterns Plus: POI visits within each wave, aggregated to grid cells (NYC) or counties (Utah).
# Writes per-cell visit summaries, daily visits, visitor home-CBG flows and brand co-visitation edges.
# Requires the licensed Advan files in data/advan_weekly_plus/ (not distributed).
source("R/setup.R")
library(DBI)
library(duckdb)
library(lubridate)
library(jsonlite)

naics_lookup <- data.table(
  sub_category = c(
    "Religious Organizations", "Drinking Places (Alcoholic Beverages)",
    "Full-Service Restaurants", "Cafeterias, Grill Buffets, and Buffets",
    "Fitness and Recreational Sports Centers", "Personal Care Services",
    "Barber Shops", "Beauty Salons", "Nail Salons",
    "Offices of Dentists", "Offices of Physicians",
    "Offices of Other Health Practitioners",
    "Child Day Care Services", "Elementary and Secondary Schools",
    "Colleges, Universities, and Professional Schools",
    "Museums, Historical Sites, and Similar Institutions",
    "Bowling Centers", "All Other Amusement and Recreation Industries",
    "Limited-Service Restaurants", "Snack and Nonalcoholic Beverage Bars",
    "Specialty Food Stores", "Meat Markets",
    "Other Schools and Instruction", "Clothing Stores",
    "Shoe Stores", "Jewelry Stores", "Sporting Goods Stores",
    "Drycleaning and Laundry Services",
    "Supermarkets and Other Grocery (except Convenience) Stores",
    "Convenience Stores", "Pharmacies and Drug Stores",
    "Cosmetics, Beauty Supplies, and Perfume Stores",
    "Depository Credit Intermediation",
    "Electronics and Appliance Stores",
    "Home Furnishings Stores", "Furniture Stores",
    "Automobile Dealers", "Tire Dealers",
    "Investigation and Security Services",
    "General Automotive Repair", "Car Washes",
    "Home Centers", "Paint and Wallpaper Stores",
    "Bus and Other Motor Vehicle Transit Systems",
    "Gasoline Stations with Convenience Stores", "Other Gasoline Stations",
    "Building Contractors",
    "Manufacturing", "Nature Parks and Other Similar Institutions",
    "Amusement and Theme Parks", "Golf Courses and Country Clubs",
    "Couriers and Express Delivery Services"),
  category_risk = c(
    1.8, 1.7, 1.5, 1.5, 1.5, 1.5, 1.5, 1.4, 1.4,
    1.4, 1.3, 1.3, 1.3, 1.2, 1.1, 1.2, 1.3, 1.2,
    1.2, 1.1, 1.2, 1.1, 1.2, 1.1, 1.1, 1.1, 1.1,
    1.1, 1.0, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.9, 0.9, 0.9, 0.8, 0.9, 0.9, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7),
  enclosed_mod = c(
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 0.5, 1.0, 0.5, 0.5, 1.0, 1.0, 0.0, 0.5,
    0.5, 0.0, 1.0, 0.0, 0.5, 0.0, 0.0),
  effective_hours = c(
    8, 14, 12, 12, 16, 10, 10, 10, 10,
    9, 9, 9, 11, 10, 10, 8, 12, 10,
    16, 12, 11, 11, 10, 11, 11, 11, 11,
    10, 16, 24, 14, 11, 10, 11, 11, 11,
    11, 10, 10, 10, 12, 12, 11, 16, 24,
    24, 10, 10, 14, 10, 12, 10))

top_category_risk <- data.table(
  top_category = c(
    "Religious Organizations", "Drinking Places (Alcoholic Beverages)",
    "Restaurants and Other Eating Places", "Fitness and Recreational Sports Centers",
    "Personal Care Services", "Offices of Dentists",
    "Offices of Physicians", "Offices of Other Health Practitioners",
    "Child Day Care Services", "Elementary and Secondary Schools",
    "Museums, Historical Sites", "Other Amusement and Recreation",
    "Limited-Service Restaurants", "Specialty Food Stores",
    "Other Schools and Instruction", "Clothing Stores",
    "Sporting Goods Stores", "Jewelry Stores",
    "Drycleaning and Laundry Services", "Grocery Stores", "Supermarkets",
    "Health and Personal Care Stores", "Depository Credit Intermediation",
    "Electronics and Appliance Stores", "Home Furnishings Stores",
    "Furniture Stores", "Office-based services", "Automobile Dealers",
    "Investigation and Security Services", "Automotive Repair and Maintenance",
    "Building Material Dealers", "Urban Transit Systems",
    "Gasoline Stations", "Building Contractors",
    "Manufacturing", "Nature Parks",
    "Amusement and Theme Parks", "Golf Courses and Country Clubs",
    "Couriers and Express Delivery"),
  category_risk_top = c(
    1.8, 1.7, 1.5, 1.5, 1.5, 1.4, 1.3, 1.3,
    1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1,
    1.1, 1.1, 1.1, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.8,
    0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7))

parse_json_safe <- function(x) tryCatch(fromJSON(x), error = function(e) NULL)

merge_count_jsons <- function(json_vec) {
  json_vec <- json_vec[!is.na(json_vec) & nchar(json_vec) > 2]
  if (length(json_vec) == 0) return(NA_character_)
  agg <- list()
  for (j in json_vec) {
    p <- parse_json_safe(j)
    if (is.null(p) || length(p) == 0) next
    nms <- names(p)
    if (is.null(nms)) next
    vals <- suppressWarnings(as.numeric(unlist(p)))
    keep <- !is.na(vals)
    if (!any(keep)) next
    for (i in which(keep)) {
      k <- nms[i]
      if (is.null(k) || is.na(k) || nchar(k) == 0) next
      agg[[k]] <- (if (is.null(agg[[k]])) 0 else agg[[k]]) + vals[i]
    }
  }
  if (length(agg) == 0) return(NA_character_)
  toJSON(agg, auto_unbox = TRUE)
}

sum_array_jsons <- function(json_vec, n) {
  acc <- numeric(n)
  any_ok <- FALSE
  for (j in json_vec) {
    if (is.na(j) || nchar(j) < 3) next
    v <- suppressWarnings(as.numeric(parse_json_safe(j)))
    if (length(v) == 0 || all(is.na(v))) next
    v <- v[1:min(length(v), n)]
    v[is.na(v)] <- 0
    if (length(v) < n) v <- c(v, rep(0, n - length(v)))
    acc <- acc + v
    any_ok <- TRUE
  }
  if (!any_ok) return(NA_character_)
  toJSON(acc, auto_unbox = FALSE)
}

collapse_vbeh_to_24 <- function(j) {
  if (is.na(j) || nchar(j) < 3) return(NA_character_)
  v <- suppressWarnings(as.numeric(parse_json_safe(j)))
  if (length(v) < 24) return(NA_character_)
  v[is.na(v)] <- 0
  if (length(v) >= 168) {
    m <- matrix(v[1:168], nrow = 7, byrow = TRUE)
    h24 <- colSums(m)
  } else {
    h24 <- v[1:24]
  }
  toJSON(h24, auto_unbox = FALSE)
}

parse_wd_frac <- function(j) {
  v <- tryCatch(as.numeric(fromJSON(j)), error = function(e) NULL)
  if (is.null(v) || length(v) < 7 || sum(v, na.rm = TRUE) == 0) return(NA_real_)
  v <- v[1:7]; v[is.na(v)] <- 0
  tot <- sum(v); if (tot == 0) return(NA_real_)
  sum(v[1:5]) / tot
}

parse_hourly <- function(j) {
  h <- tryCatch(as.numeric(fromJSON(j)), error = function(e) NULL)
  if (is.null(h) || length(h) < 24 || sum(h, na.rm = TRUE) == 0)
    return(list(pk = NA_integer_, hc = NA_real_))
  h <- h[1:24]; h[is.na(h)] <- 0; props <- h / sum(h)
  list(pk = which.max(h) - 1L, hc = sum(props^2))
}

parse_daily_cv <- function(j) {
  v <- tryCatch(as.numeric(fromJSON(j)), error = function(e) NULL)
  if (is.null(v) || length(v) < 2 || sum(v, na.rm = TRUE) == 0)
    return(list(cv = NA_real_, wd_dom = NA_real_))
  v[is.na(v)] <- 0
  props <- v / sum(v)
  list(cv = sd(v) / mean(v), wd_dom = max(props))
}

parse_bucketed_dwell <- function(j) {
  d <- parse_json_safe(j)
  if (is.null(d) || length(d) == 0) return(NA_real_)
  midpoints <- c(2.5, 12.5, 40, 150, 360)
  counts <- suppressWarnings(as.numeric(unlist(d)))
  if (length(counts) > 5) counts <- counts[1:5]
  counts[is.na(counts)] <- 0
  if (sum(counts) > 0)
    sum(midpoints[1:length(counts)] * counts) / sum(counts)
  else NA_real_
}

parse_device_ios_frac <- function(j) {
  if (is.na(j) || nchar(j) < 3) return(NA_real_)
  p <- parse_json_safe(j)
  if (is.null(p) || length(p) == 0) return(NA_real_)
  nms <- tolower(names(p))
  vals <- suppressWarnings(as.numeric(unlist(p)))
  vals[is.na(vals)] <- 0
  ios <- sum(vals[nms == "ios"])
  android <- sum(vals[nms == "android"])
  tot <- ios + android
  if (tot <= 0) return(NA_real_)
  ios / tot
}

parse_visitor_home_cbgs <- function(poi_home_dt, state) {
  poi_home_dt <- poi_home_dt[!is.na(visitor_home_cbgs_json) & nchar(visitor_home_cbgs_json) > 2]
  if (nrow(poi_home_dt) == 0) return(data.table())
  state_fips_2 <- substr(get_metro_fips(state, 5)[1], 1, 2)
  chunks <- split(seq_len(nrow(poi_home_dt)), ceiling(seq_len(nrow(poi_home_dt)) / 5000))
  rbindlist(lapply(chunks, function(idx) {
    out <- poi_home_dt[idx][, {
      parsed <- parse_json_safe(.BY$visitor_home_cbgs_json)
      if (is.null(parsed) || length(parsed) == 0)
        data.table(home_cbg = character(0), visitors = integer(0))
      else
        data.table(home_cbg = names(parsed),
                   visitors = as.integer(unlist(parsed)))
    }, by = .(poi_cbg, poi_lon, poi_lat, visits_total, median_dwell, brand,
              naics_code, visitor_home_cbgs_json)]
    out[, visitor_home_cbgs_json := NULL]
    out[!is.na(home_cbg) & substr(home_cbg, 1, 2) == state_fips_2 & visitors > 0]
  }))
}

parse_covisitation <- function(od_raw, json_col) {
  rows <- od_raw[!is.na(get(json_col)) & nchar(get(json_col)) > 2 & nchar(brand) > 0]
  if (nrow(rows) == 0) return(data.table())
  rbindlist(lapply(seq_len(nrow(rows)), function(i) {
    parsed <- parse_json_safe(rows[i][[json_col]])
    if (is.null(parsed) || length(parsed) == 0) return(NULL)
    data.table(source_brand = rows[i]$brand,
               target_brand = names(parsed),
               covisit_count = suppressWarnings(as.integer(unlist(parsed))))
  }))
}

parse_daytime_visitors <- function(j) {
  d <- parse_json_safe(j)
  if (is.null(d) || length(d) == 0) return(0L)
  v <- suppressWarnings(as.integer(unlist(d)))
  sum(v, na.rm = TRUE)
}

aggregate_weekly_to_poi <- function(weekly_dt) {
  weekly_dt[, w := pmax(visits_total_w, 0)]
  weekly_dt[, ios_frac := vapply(device_type_json, parse_device_ios_frac,
                                 numeric(1))]
  agg <- weekly_dt[, .(
    poi_lon         = poi_lon[1],
    poi_lat         = poi_lat[1],
    top_category    = top_category[1],
    sub_category    = sub_category[1],
    brand           = brand[1],
    naics_code      = naics_code[1],
    poi_cbg         = poi_cbg[1],
    region          = region[1],

    visits_total    = sum(visits_total_w, na.rm = TRUE),
    visitors_unique = sum(visitors_unique_w, na.rm = TRUE),
    median_dwell    = if (sum(w, na.rm = TRUE) > 0)
                        weighted.mean(median_dwell_w, w, na.rm = TRUE)
                      else NA_real_,
    distance_from_home = if (sum(w, na.rm = TRUE) > 0)
                           weighted.mean(distance_from_home_w, w, na.rm = TRUE)
                         else NA_real_,
    ios_frac        = if (sum(w, na.rm = TRUE) > 0)
                        weighted.mean(ios_frac, w, na.rm = TRUE)
                      else NA_real_,

    visits_by_day_json       = sum_array_jsons(visits_by_day_json_w, 7),
    popularity_by_hour_json  = {
      h24 <- vapply(visits_by_each_hour_json_w, collapse_vbeh_to_24,
                    character(1))
      sum_array_jsons(h24, 24)
    },
    visitor_home_cbgs_json   = merge_count_jsons(visitor_home_cbgs_json_w),
    visitor_daytime_cbgs_json = merge_count_jsons(visitor_daytime_cbgs_json_w),
    bucketed_dwell_times_json = merge_count_jsons(bucketed_dwell_times_json_w),
    related_same_day_brand_json  = merge_count_jsons(related_same_day_brand_json_w),
    related_same_week_brand_json = merge_count_jsons(related_same_week_brand_json_w)
  ), by = .(persistent_id)]

  agg[, wd_frac := vapply(visits_by_day_json, parse_wd_frac, numeric(1))]
  agg
}

lookup_parse <- function(dt, json_col, parse_fn, out_cols) {
  valid <- !is.na(dt[[json_col]]) & nchar(dt[[json_col]]) > 2
  if (!any(valid)) { dt[, (json_col) := NULL]; return(dt) }
  uniq <- unique(dt[[json_col]][valid])
  lk <- rbindlist(lapply(uniq, function(j) {
    res <- parse_fn(j)
    if (is.list(res)) c(setNames(list(j), json_col), res)
    else c(setNames(list(j), json_col), setNames(list(res), out_cols))
  }), fill = TRUE)
  dt <- merge(dt, lk, by = json_col, all.x = TRUE)
  dt[, (json_col) := NULL]
  dt
}

process_pois <- function(poi_dt, state) {
  poi_dt <- merge(poi_dt, naics_lookup, by = "sub_category", all.x = TRUE)
  poi_dt <- merge(poi_dt, top_category_risk, by = "top_category", all.x = TRUE)
  poi_dt[is.na(category_risk),    category_risk    := category_risk_top]
  poi_dt[is.na(category_risk),    category_risk    := 1.0]
  poi_dt[is.na(enclosed_mod),     enclosed_mod     := 1.0]
  poi_dt[is.na(effective_hours),  effective_hours  := 12]
  poi_dt[, category_risk_top := NULL]

  area_dt <- build_poi_area_cache(poi_dt, state)
  poi_dt[area_dt, poi_area_sqm := i.poi_area_sqm, on = "persistent_id"]
  match_rate <- mean(!is.na(poi_dt$poi_area_sqm))
  message(sprintf("  POI footprint match rate: %.1f%% (%s/%s)",
                  100 * match_rate,
                  fmt(sum(!is.na(poi_dt$poi_area_sqm))),
                  fmt(nrow(poi_dt))))
  if (match_rate < 0.80)
    warning(sprintf("Low POI footprint coverage (%.1f%%); crowding factor will lean on median imputation.",
                    100 * match_rate))

  poi_dt[, dwell_hours := pmin(median_dwell / 60, 8)]
  poi_dt[is.na(dwell_hours) | dwell_hours <= 0, dwell_hours := 0.5]
  poi_dt[, activity_per_visit := (1 + dwell_hours) * category_risk * enclosed_mod]

  poi_dt[, visits_norm := as.numeric(visits_total)]
  poi_dt[is.na(visits_norm) | visits_norm <= 0, visits_norm := 0]
  poi_dt[, activity_total := visits_norm * activity_per_visit]
  poi_dt[is.na(distance_from_home) | distance_from_home < 0,
         distance_from_home := NA_real_]

  poi_dt[, visit_density := fifelse(
    is.na(poi_area_sqm) | poi_area_sqm <= 0, NA_real_,
    visits_norm / poi_area_sqm)]
  poi_dt[, contact_density_hourly := fifelse(
    is.na(poi_area_sqm) | poi_area_sqm <= 0 | is.na(effective_hours) |
      effective_hours <= 0, NA_real_,
    visits_norm / poi_area_sqm / effective_hours)]

  poi_dt[is.na(wd_frac), wd_frac := 5/7]
  poi_dt[, visits_weekday := round(visits_norm * wd_frac)]
  poi_dt[, visits_weekend := visits_norm - visits_weekday][, wd_frac := NULL]

  poi_dt[, c("peak_hour", "hourly_concentration") := .(NA_integer_, NA_real_)]
  poi_dt <- lookup_parse(poi_dt, "popularity_by_hour_json", parse_hourly,
                         c("pk", "hc"))
  poi_dt[!is.na(pk), peak_hour := pk]
  poi_dt[!is.na(hc), hourly_concentration := hc]
  poi_dt[, c("pk", "hc") := NULL]

  poi_dt[, c("daily_visit_cv", "weekday_dominance") := .(NA_real_, NA_real_)]
  poi_dt <- lookup_parse(poi_dt, "visits_by_day_json", parse_daily_cv,
                         c("cv", "wd_dom"))
  poi_dt[!is.na(cv), daily_visit_cv := cv][, cv := NULL]
  poi_dt[!is.na(wd_dom), weekday_dominance := wd_dom][, wd_dom := NULL]

  poi_dt[, repeat_ratio := visits_total / pmax(visitors_unique, 1L)]

  poi_dt[, bdwell := vapply(bucketed_dwell_times_json,
                            parse_bucketed_dwell, numeric(1))]
  poi_dt[!is.na(bdwell), median_dwell := bdwell]
  poi_dt[, bdwell := NULL]

  poi_dt[is.na(ios_frac), ios_frac := 0.5]
  poi_dt[, visits_ios     := visits_norm * ios_frac]
  poi_dt[, visits_android := visits_norm * (1 - ios_frac)]
  poi_dt[, activity_ios     := activity_total * ios_frac]
  poi_dt[, activity_android := activity_total * (1 - ios_frac)]
  poi_dt[, activity_weekday := visits_weekday * activity_per_visit]
  poi_dt[, activity_weekend := visits_weekend * activity_per_visit]
  poi_dt
}

aggregate_to_grid <- function(dt) {
  dt[, .(
    n_pois = .N,
    visits_total       = sum(visits_total, na.rm = TRUE),
    visits_norm        = sum(visits_norm, na.rm = TRUE),
    visits_weekday     = sum(visits_weekday, na.rm = TRUE),
    visits_weekend     = sum(visits_weekend, na.rm = TRUE),
    activity_total     = sum(activity_total, na.rm = TRUE),
    activity_weekday   = sum(activity_weekday, na.rm = TRUE),
    activity_weekend   = sum(activity_weekend, na.rm = TRUE),
    visits_ios         = sum(visits_ios, na.rm = TRUE),
    visits_android     = sum(visits_android, na.rm = TRUE),
    activity_ios       = sum(activity_ios, na.rm = TRUE),
    activity_android   = sum(activity_android, na.rm = TRUE),
    avg_dwell                  = weighted.mean(dwell_hours, visits_norm, na.rm = TRUE),
    avg_category_risk          = weighted.mean(category_risk, visits_norm, na.rm = TRUE),
    enclosed_frac              = weighted.mean(enclosed_mod, visits_norm, na.rm = TRUE),
    avg_distance_from_home     = weighted.mean(distance_from_home, visits_norm, na.rm = TRUE),
    avg_poi_area               = mean(poi_area_sqm, na.rm = TRUE),
    avg_visit_density          = weighted.mean(visit_density, visits_norm, na.rm = TRUE),
    avg_peak_hour              = weighted.mean(peak_hour, visits_norm, na.rm = TRUE),
    avg_hourly_conc            = weighted.mean(hourly_concentration, visits_norm, na.rm = TRUE),
    avg_daily_cv               = weighted.mean(daily_visit_cv, visits_norm, na.rm = TRUE),
    avg_weekday_dominance      = weighted.mean(weekday_dominance, visits_norm, na.rm = TRUE),
    avg_repeat_ratio           = weighted.mean(repeat_ratio, visits_norm, na.rm = TRUE),
    avg_effective_hours        = weighted.mean(effective_hours, visits_norm, na.rm = TRUE),
    avg_contact_density_hourly = weighted.mean(contact_density_hourly, visits_norm, na.rm = TRUE)
  ), by = grid_id]
}

select_files_for_range <- function(advan_files, wave_start, wave_end) {
  bn <- basename(advan_files)
  pk <- substr(bn, 1, 10)
  ok <- suppressWarnings(!is.na(as.Date(pk, format = "%Y-%m-%d")))
  pk <- as.Date(pk[ok], format = "%Y-%m-%d")
  files <- advan_files[ok]
  keep <- (pk + 6) >= wave_start & pk <= wave_end
  files[keep]
}

process_wave <- function(wave_id, wave_months) {
  out_dir <- get_wave_dir(state, wave_id)
  dir.create(file.path(out_dir, "intermediate"), recursive = TRUE,
             showWarnings = FALSE)

  fips_5  <- get_metro_fips(state, 5)
  utm_crs <- get_utm_crs(state)
  wdr     <- get_wave_date_range(wave_months)
  total_n_days <- wdr$n_days

  grid_sf   <- st_read(file.path(CONFIG$grids_dir,
                                 sprintf("grid_100m_%s.gpkg", state)),
                       quiet = TRUE)
  grid_sf   <- st_transform(grid_sf, utm_crs)
  grid_join <- grid_sf[, "grid_id"]

  advan_files <- list.files(CONFIG$advan_dir, pattern = "\\.parquet$",
                            full.names = TRUE)
  wave_files <- select_files_for_range(advan_files, wdr$start, wdr$end)
  message(sprintf("[%s/%s] %d/%d weekly files in [%s, %s]",
                  state, wave_id, length(wave_files), length(advan_files),
                  wdr$start, wdr$end))
  if (length(wave_files) == 0)
    stop("No weekly files matched the wave window. Check filename format.")

  con <- dbConnect(duckdb(), dbdir = ":memory:")
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  fsql <- paste0("'", fips_5, "'", collapse = ", ")
  region_filter <- if (state == "UT") "REGION = 'UT'"
                   else "REGION = 'NY'"

  files_glob <- paste0("'", paste(wave_files, collapse = "','"), "'")
  weekly_query <- sprintf("
    SELECT
      CAST(PERSISTENT_ID_STORE AS VARCHAR)              AS persistent_id,
      CAST(PERSISTENT_ID AS VARCHAR)                    AS company_id,
      CAST(LONGITUDE AS DOUBLE)                         AS poi_lon,
      CAST(LATITUDE  AS DOUBLE)                         AS poi_lat,
      CAST(POI_CBG AS VARCHAR)                          AS poi_cbg,
      CAST(REGION AS VARCHAR)                           AS region,
      TOP_CATEGORY                                      AS top_category,
      COALESCE(SUB_CATEGORY,'')                         AS sub_category,
      COALESCE(BRAND,'')                                AS brand,
      COALESCE(NAICS_CODE,'')                           AS naics_code,

      DATE_RANGE_START                                  AS week_start,
      CAST(VISIT_COUNTS  AS DOUBLE)                     AS visits_total_w,
      CAST(VISITOR_COUNTS AS DOUBLE)                    AS visitors_unique_w,
      CAST(MEDIAN_DWELL  AS DOUBLE)                     AS median_dwell_w,
      CAST(DISTANCE_FROM_HOME AS DOUBLE)                AS distance_from_home_w,
      CAST(OPEN_DATE  AS VARCHAR)                       AS open_date,
      CAST(CLOSE_DATE AS VARCHAR)                       AS close_date,
      IS_DISTRIBUTOR                                    AS is_distributor,

      COALESCE(CAST(VISITS_BY_DAY        AS VARCHAR),'') AS visits_by_day_json_w,
      COALESCE(CAST(VISITS_BY_EACH_HOUR  AS VARCHAR),'') AS visits_by_each_hour_json_w,
      COALESCE(CAST(VISITOR_HOME_CBGS    AS VARCHAR),'') AS visitor_home_cbgs_json_w,
      COALESCE(CAST(VISITOR_DAYTIME_CBGS AS VARCHAR),'') AS visitor_daytime_cbgs_json_w,
      COALESCE(CAST(BUCKETED_DWELL_TIMES AS VARCHAR),'') AS bucketed_dwell_times_json_w,
      COALESCE(CAST(RELATED_SAME_DAY_BRAND  AS VARCHAR),'') AS related_same_day_brand_json_w,
      COALESCE(CAST(RELATED_SAME_WEEK_BRAND AS VARCHAR),'') AS related_same_week_brand_json_w,
      COALESCE(CAST(DEVICE_TYPE          AS VARCHAR),'') AS device_type_json
    FROM read_parquet([%s])
    WHERE %s
      AND SUBSTR(CAST(POI_CBG AS VARCHAR),1,5) IN (%s)
      AND VISIT_COUNTS > 0
      AND LONGITUDE IS NOT NULL AND LATITUDE IS NOT NULL
      AND PERSISTENT_ID_STORE IS NOT NULL",
    files_glob, region_filter, fsql)

  weekly_dt <- as.data.table(dbGetQuery(con, weekly_query))
  message(sprintf("  read %s weekly POI rows", fmt(nrow(weekly_dt))))
  if (nrow(weekly_dt) == 0) {
    return(list(state = state, wave_id = wave_id, n_pois = 0,
                n_days = total_n_days))
  }

  n_before <- nrow(weekly_dt)
  weekly_dt <- weekly_dt[is.na(is_distributor) | is_distributor != 1]
  if (nrow(weekly_dt) < n_before)
    message(sprintf("  dropped %s distributor rows",
                    fmt(n_before - nrow(weekly_dt))))

  weekly_dt[, open_d  := suppressWarnings(as.Date(open_date))]
  weekly_dt[, close_d := suppressWarnings(as.Date(close_date))]
  weekly_dt[is.na(open_d),  open_d  := as.Date("1970-01-01")]
  weekly_dt[is.na(close_d), close_d := as.Date("2038-01-01")]
  weekly_dt[, week_end_d := as.Date(week_start) + 6L]
  n_before <- nrow(weekly_dt)
  weekly_dt <- weekly_dt[open_d <= week_end_d & close_d >= as.Date(week_start)]
  if (nrow(weekly_dt) < n_before)
    message(sprintf("  dropped %s rows outside POI active window",
                    fmt(n_before - nrow(weekly_dt))))
  weekly_dt[, c("open_d", "close_d", "week_end_d",
                "open_date", "close_date", "is_distributor") := NULL]
  if (nrow(weekly_dt) == 0) {
    return(list(state = state, wave_id = wave_id, n_pois = 0,
                n_days = total_n_days))
  }

  weekly_dt[, week_start_d := as.Date(week_start)]
  vbyday <- lapply(weekly_dt$visits_by_day_json_w, function(j) {
    v <- suppressWarnings(as.numeric(parse_json_safe(j)))
    if (length(v) >= 7) { v <- v[1:7]; v[is.na(v)] <- 0; v }
    else rep(0, 7)
  })
  vmat <- do.call(rbind, vbyday)
  n <- nrow(weekly_dt)
  start_int <- as.integer(wdr$start)
  end_int   <- as.integer(wdr$end)
  daily_dt <- data.table(
    day_int = rep(as.integer(weekly_dt$week_start_d), each = 7) +
              rep(0L:6L, times = n),
    vd      = as.numeric(t(vmat)),
    dwell   = rep(weekly_dt$median_dwell_w, each = 7)
  )
  daily_dt[, wt := vd]
  daily_dt <- daily_dt[day_int >= start_int & day_int <= end_int]

  ds_all <- daily_dt[, .(
    total_visits = sum(vd, na.rm = TRUE),
    median_dwell = if (sum(wt, na.rm = TRUE) > 0)
                     weighted.mean(dwell, wt, na.rm = TRUE) else NA_real_
  ), by = day_int]
  setorder(ds_all, day_int)
  ds_all[, date := as.Date(day_int, origin = "1970-01-01")]
  ds_all[, sim_day := day_int - start_int + 1L]
  ds_all[, source_month := format(date, "%Y-%m")]
  ds_all[, is_weekend := as.integer(lubridate::wday(date, week_start = 1) >= 6)]
  ds_all[, day_int := NULL]
  setcolorder(ds_all, c("sim_day", "date", "source_month", "is_weekend",
                        "total_visits", "median_dwell"))
  poi_daily <- data.table(
    persistent_id = rep(weekly_dt$persistent_id, each = 7),
    day_int = rep(as.integer(weekly_dt$week_start_d), each = 7) + rep(0L:6L, times = n),
    vd = as.numeric(t(vmat)))
  poi_daily <- poi_daily[day_int >= start_int & day_int <= end_int & is.finite(vd)]
  rm(daily_dt, vbyday, vmat); weekly_dt[, week_start_d := NULL]
  gc(verbose = FALSE)

  od_w_src <- weekly_dt[!is.na(visitor_home_cbgs_json_w) & nchar(visitor_home_cbgs_json_w) > 2,
                        .(poi_cbg, poi_lon, poi_lat, week_start,
                          visits_total = visits_total_w, median_dwell = median_dwell_w,
                          brand, naics_code, visitor_home_cbgs_json = visitor_home_cbgs_json_w)]
  if (nrow(od_w_src) > 0) {
    od_w <- rbindlist(lapply(split(od_w_src, od_w_src$week_start), function(ws) {
      x <- parse_visitor_home_cbgs(ws, state); if (nrow(x)) x[, week_start := as.Date(ws$week_start[1])]; x
    }), fill = TRUE)
    if (nrow(od_w) > 0) {
      od_w <- od_w[, .(visitors = sum(visitors)), by = .(poi_cbg, poi_lon, poi_lat, week_start, home_cbg)]
      write_feather(od_w, file.path(out_dir, "intermediate", "od_home_cbg_weekly.feather"))
    }
    rm(od_w_src); gc(verbose = FALSE)
  }

  poi_dt <- aggregate_weekly_to_poi(weekly_dt)
  rm(weekly_dt); gc(verbose = FALSE)
  message(sprintf("  %s unique POIs after wave aggregation", fmt(nrow(poi_dt))))

  poi_dt <- process_pois(poi_dt, state)

  poi_sf <- st_join(
    st_transform(st_as_sf(poi_dt[!is.na(poi_lon) & !is.na(poi_lat)],
                          coords = c("poi_lon", "poi_lat"), crs = 4326),
                 utm_crs),
    grid_join, left = FALSE)
  if (nrow(poi_sf) > 0) {
    grid_dt <- aggregate_to_grid(as.data.table(st_drop_geometry(poi_sf)))
    grid_dt[, `:=`(
      log_visits         = log1p(visits_norm),
      log_activity       = log1p(activity_total),
      log_ios_visits     = log1p(visits_ios),
      log_android_visits = log1p(visits_android),
      ios_frac           = visits_ios / pmax(visits_ios + visits_android, 1)
    )]
    grid_dt[, ses_contrast := ios_frac - 0.5]
    write_feather(grid_dt, file.path(out_dir, "intermediate",
                                     "grid_visits.feather"))
  }

  if (nrow(poi_sf) > 0) {
    pg <- as.data.table(st_drop_geometry(poi_sf))[, .(persistent_id, grid_id)]
    pg <- unique(pg, by = "persistent_id")
    gd <- merge(poi_daily, pg, by = "persistent_id")[, .(visits = sum(vd, na.rm = TRUE), n_pois = uniqueN(persistent_id)),
                                                     by = .(grid_id, day_int)]
    gd[, `:=`(date = as.Date(day_int, origin = "1970-01-01"), sim_day = day_int - start_int + 1L)]
    gd[, day_int := NULL]
    setcolorder(gd, c("grid_id", "sim_day", "date", "visits", "n_pois"))
    write_feather(gd, file.path(out_dir, "intermediate", "grid_daily_visits.feather"))
  }
  rm(poi_daily)
  fwrite(ds_all, file.path(out_dir, "intermediate", "daily_stats.csv"))
  od_src <- poi_dt[!is.na(visitor_home_cbgs_json) &
                   nchar(visitor_home_cbgs_json) > 2,
                   .(poi_cbg, poi_lon, poi_lat, visits_total, median_dwell,
                     brand, naics_code, visitor_home_cbgs_json)]
  od_long <- parse_visitor_home_cbgs(od_src, state)
  if (nrow(od_long) > 0) {
    od_long[, `:=`(origin_fips = substr(home_cbg, 1, 5),
                   dest_fips   = substr(poi_cbg, 1, 5))]
    write_feather(od_long, file.path(out_dir, "intermediate",
                                     "od_home_cbg.feather"))
  }

  cov_src <- poi_dt[, .(
    brand = brand,
    covisit_json      = related_same_day_brand_json,
    covisit_week_json = related_same_week_brand_json
  )]
  cv_day  <- parse_covisitation(cov_src, "covisit_json")
  if (nrow(cv_day) > 0) cv_day[, scope := "day"]
  cv_week <- parse_covisitation(cov_src, "covisit_week_json")
  if (nrow(cv_week) > 0) cv_week[, scope := "week"]
  cv <- rbindlist(list(cv_day, cv_week), fill = TRUE)
  if (nrow(cv) > 0) {
    cv_agg <- cv[, .(covisit_count = sum(covisit_count), n_pois = .N),
                 by = .(source_brand, target_brand, scope)][order(-covisit_count)]
    write_feather(cv_agg, file.path(out_dir, "intermediate",
                                    "trip_chains.feather"))
  }

  dt_src <- poi_dt[!is.na(visitor_daytime_cbgs_json) &
                   nchar(visitor_daytime_cbgs_json) > 2]
  if (nrow(dt_src) > 0) {
    daytime_agg <- dt_src[, .(
      poi_lon = poi_lon[1], poi_lat = poi_lat[1],
      daytime_visitors = parse_daytime_visitors(visitor_daytime_cbgs_json)
    ), by = poi_cbg]
    write_feather(daytime_agg, file.path(out_dir, "intermediate",
                                         "daytime_visitors.feather"))
  }

  list(state = state, wave_id = wave_id, n_pois = nrow(poi_dt),
       n_days = total_n_days)
}

invisible(run_waves_by_state(process_wave))
