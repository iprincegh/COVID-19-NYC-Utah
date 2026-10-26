# Connected neighbour graph for the 100 m grid and Utah counties: queen contiguity plus
# transport links (bridges, tunnels, ferries, subway, roads) and nearest-neighbour joins.

suppressPackageStartupMessages({ library(RANN); library(igraph) })

SNAP_QUEEN         <- 10
SNAP_TIERS         <- c(300, 500, 750, 1000, 1500)
SNAP_MAX           <- 2000
OD_CORR_THRESHOLD  <- 0.6
OD_MIN_ACTIVITY    <- 5.0
OD_MIN_DAYS        <- 10
OD_MIN_ROWS        <- 1000
OD_MAX_FRAG_SIZE   <- 10
OD_MAX_FRAG_DEGREE <- 1.5
MAX_WORMHOLES      <- 10

NYC_WORMHOLES <- list(
  list(name="LGA-Midtown",           from=c(-73.8740,40.7769), to=c(-73.9857,40.7484), priority="critical"),
  list(name="LGA-EastElmhurst",      from=c(-73.8740,40.7769), to=c(-73.8700,40.7610), priority="critical"),
  list(name="JFK-DowntownBK",         from=c(-73.7781,40.6413), to=c(-73.9857,40.6892), priority="critical"),
  list(name="JFK-Jamaica",            from=c(-73.7781,40.6413), to=c(-73.7930,40.7028), priority="critical"),
  list(name="JFK-HowardBeach",        from=c(-73.7781,40.6413), to=c(-73.8300,40.6590), priority="critical"),
  list(name="StGeorge-Whitehall",     from=c(-74.0766,40.6437), to=c(-74.0130,40.7013), priority="critical", crossing=TRUE),
  list(name="VerrazzanoNarrows-SI-BK", from=c(-74.0539,40.6058), to=c(-74.0330,40.6070), priority="critical", crossing=TRUE),
  list(name="Greenpt-E34thFerry",     from=c(-73.9615,40.7320), to=c(-73.9718,40.7440), priority="critical", crossing=TRUE),
  list(name="BKCruiseTerminal-DTBK",  from=c(-74.0168,40.6808), to=c(-73.9857,40.6892), priority="manual"),
  list(name="BKCruiseTerminal-BayR",  from=c(-74.0168,40.6808), to=c(-74.0302,40.6346), priority="manual"),
  list(name="SunnysideYd-PennSta",    from=c(-73.9215,40.7425), to=c(-73.9937,40.7506), priority="manual"),
  list(name="AtlanticYd-Atlantic",    from=c(-73.9762,40.6846), to=c(-73.9776,40.6862), priority="manual"),
  list(name="GrandCentral-Harlem125", from=c(-73.9772,40.7527), to=c(-73.9390,40.8043), priority="manual"),
  list(name="Fulton-Barclays",        from=c(-74.0071,40.7092), to=c(-73.9762,40.6846), priority="manual"),
  list(name="CrossBronx-GWBridge",    from=c(-73.8800,40.8380), to=c(-73.9520,40.8510), priority="manual"),
  list(name="BQE-HughCarey",          from=c(-74.0020,40.6870), to=c(-74.0145,40.6895), priority="manual")
)

NYC_SUBWAY_CROSSINGS <- list(
  list(name="L-14thSt-tunnel",        from=c(-73.9815,40.7308), to=c(-73.9566,40.7172)),
  list(name="EM-53rdSt-tunnel",       from=c(-73.9691,40.7576), to=c(-73.9455,40.7472)),
  list(name="NRW-60thSt-tunnel",      from=c(-73.9675,40.7626), to=c(-73.9402,40.7508)),
  list(name="7-Steinway-tunnel",      from=c(-73.9766,40.7513), to=c(-73.9538,40.7427)),
  list(name="23-ClarkSt-tunnel",      from=c(-74.0091,40.7069), to=c(-73.9930,40.6975)),
  list(name="R-MontagueSt-tunnel",    from=c(-74.0130,40.7031), to=c(-73.9920,40.6941)),
  list(name="AC-CranberrySt-tunnel",  from=c(-74.0077,40.7101), to=c(-73.9906,40.6992)),
  list(name="F-RutgersSt-tunnel",     from=c(-73.9902,40.7139), to=c(-73.9866,40.7014)),
  list(name="45-JoralemonSt-tunnel",  from=c(-74.0141,40.7049), to=c(-73.9902,40.6929)),
  list(name="456-LexAv-HarlemRiver",  from=c(-73.9375,40.8043), to=c(-73.9299,40.8132)),
  list(name="BD-Concourse-HarlemRiver", from=c(-73.9384,40.8302), to=c(-73.9256,40.8279)),
  list(name="25-HarlemRiver",         from=c(-73.9477,40.8140), to=c(-73.9177,40.8161))
)

ADJ_AUG_VERSION <- 3L

CROSS_SAME_MAX  <- 2500

CROSS_OTHER_MAX <- 2500

CROSSING_EXCLUDE_PATTERN <- "George Washington|Lincoln Tunnel|Holland Tunnel"

unit_borough <- function(grid_sf) {
  xw_path <- tryCatch(CONFIG$nyc_crosswalk_path, error = function(e) NULL)
  if (is.null(xw_path) || !file.exists(xw_path)) return(NULL)
  xw <- data.table::fread(xw_path)
  if ("cbg_geoid" %in% names(grid_sf)) {
    key <- as.character(grid_sf$cbg_geoid)
    map <- xw[, .(borough = borough[1]), by = .(k = as.character(cbg_geoid))]
  } else if ("rate_geo_id" %in% names(grid_sf) || "modzcta" %in% names(grid_sf)) {
    key <- as.character(if ("rate_geo_id" %in% names(grid_sf)) grid_sf$rate_geo_id else grid_sf$modzcta)
    map <- xw[, .(borough = borough[1]), by = .(k = as.character(modzcta))]
  } else return(NULL)
  b <- map$borough[match(key, map$k)]
  if (all(is.na(b))) return(NULL)
  b
}

NYC_EXTRA_CROSSINGS <- list(
  list(name="MacombsDam-Bridge",       kind="bridge", from=c(-73.9380,40.8290), to=c(-73.9290,40.8285)),
  list(name="AlexanderHamilton-Bridge", kind="bridge", from=c(-73.9310,40.8465), to=c(-73.9210,40.8450)),
  list(name="Washington-Bridge-181st", kind="bridge", from=c(-73.9300,40.8475), to=c(-73.9200,40.8470)),
  list(name="Pulaski-Bridge",          kind="bridge", from=c(-73.9530,40.7350), to=c(-73.9520,40.7415)),
  list(name="Kosciuszko-Bridge",       kind="bridge", from=c(-73.9330,40.7255), to=c(-73.9280,40.7315)),
  list(name="GreenpointAve-Bridge",    kind="bridge", from=c(-73.9400,40.7325), to=c(-73.9380,40.7360)),
  list(name="MetropolitanAve-Bridge",  kind="bridge", from=c(-73.9270,40.7135), to=c(-73.9245,40.7150)),
  list(name="GrandSt-Bridge",          kind="bridge", from=c(-73.9230,40.7150), to=c(-73.9205,40.7170)),
  list(name="G-Greenpoint-tube",       kind="subway", from=c(-73.9540,40.7313), to=c(-73.9495,40.7440)),
  list(name="F-63rdSt-tunnel",         kind="subway", from=c(-73.9662,40.7645), to=c(-73.9425,40.7540)),
  list(name="RooseveltIsland-Bridge",  kind="bridge", from=c(-73.9500,40.7615), to=c(-73.9385,40.7625)),
  list(name="RooseveltIsland-Tram",    kind="bridge", from=c(-73.9508,40.7574), to=c(-73.9640,40.7615), same_boro=TRUE),
  list(name="NorthChannel-Addabbo-Bridge", kind="bridge", from=c(-73.8390,40.6420), to=c(-73.8250,40.6210), same_boro=TRUE),
  list(name="CrossBay-Veterans-Bridge", kind="bridge", from=c(-73.8200,40.6010), to=c(-73.8175,40.5875), same_boro=TRUE),
  list(name="A-JamaicaBay-north",      kind="subway", from=c(-73.8305,40.6605), to=c(-73.8160,40.6085), same_boro=TRUE),
  list(name="A-JamaicaBay-south",      kind="subway", from=c(-73.8160,40.6085), to=c(-73.8135,40.5880), same_boro=TRUE),
  list(name="Rockaway-Ferry",          kind="ferry",  from=c(-74.0080,40.7015), to=c(-73.8330,40.5825))
)

UT_WORMHOLES <- list(
  list(name="SLCairport-DowntownSLC", from=c(-111.9791,40.7899), to=c(-111.8910,40.7608), priority="critical"),
  list(name="SLCairport-WestValley",  from=c(-111.9791,40.7899), to=c(-111.9390,40.6916), priority="critical"),
  list(name="I15-Weber-Davis",        from=c(-111.9730,41.2230), to=c(-111.9020,40.9880), priority="critical"),
  list(name="I15-Utah-Juab",          from=c(-111.6580,40.2330), to=c(-111.7480,39.7310), priority="critical"),
  list(name="I15-Iron-Washington",    from=c(-113.0610,37.6770), to=c(-113.5080,37.1080), priority="critical"),
  list(name="I80-SL-Tooele",          from=c(-111.8910,40.7608), to=c(-112.3260,40.6220), priority="critical"),
  list(name="I80-Summit-Wasatch",     from=c(-111.4980,40.7420), to=c(-111.5010,40.6540), priority="critical"),
  list(name="I70-Sevier-Emery",       from=c(-111.8240,38.7830), to=c(-111.2540,38.9960), priority="manual"),
  list(name="US89-Sanpete-Sevier",    from=c(-111.5790,39.3640), to=c(-111.8240,38.7830), priority="manual"),
  list(name="US6-Carbon-Duchesne",    from=c(-110.9530,39.6070), to=c(-110.4030,40.1630), priority="manual"),
  list(name="I84-Weber-BoxElder",     from=c(-111.9730,41.2230), to=c(-112.0150,41.7390), priority="manual")
)

safe_n_comp_nb <- function(nb_obj) {
  n <- length(nb_obj)
  el <- vector("list", n)
  for (i in seq_len(n)) {
    nbrs <- nb_obj[[i]]
    nbrs <- nbrs[nbrs > 0L & nbrs <= n]
    if (length(nbrs) > 0) el[[i]] <- data.frame(from = i, to = nbrs)
  }
  edges <- do.call(rbind, el)
  g <- make_empty_graph(n = n, directed = FALSE)
  if (!is.null(edges) && nrow(edges) > 0)
    g <- add_edges(g, as.vector(t(edges[, 1:2])))
  mem <- components(g)
  list(nc = mem$no, comp.id = as.integer(mem$membership))
}

nb_set_edge <- function(nb, i, j) {
  ni <- nb[[i]]; ni <- ni[ni > 0L]
  nj <- nb[[j]]; nj <- nj[nj > 0L]
  nb[[i]] <- sort(unique(c(ni, as.integer(j))))
  nb[[j]] <- sort(unique(c(nj, as.integer(i))))
  nb
}

sanitize_nb <- function(nb_obj) {
  n <- length(nb_obj)
  for (i in seq_len(n)) {
    nbrs <- nb_obj[[i]]
    nbrs <- as.integer(nbrs[!is.na(nbrs)])
    nbrs <- sort(unique(nbrs[nbrs != i & nbrs > 0L & nbrs <= n]))
    nb_obj[[i]] <- if (length(nbrs) == 0) 0L else nbrs
  }
  attr(nb_obj, "class") <- "nb"
  attr(nb_obj, "region.id") <- as.character(1:n)
  nb_obj
}

safe_nb2INLA <- function(graph_file, nb_obj) {
  nb2INLA(graph_file, sanitize_nb(nb_obj))
}

apply_wormholes <- function(nb, wh_list, coords_proj, layer_name) {
  n_wh <- 0
  if (length(wh_list) == 0) return(list(nb = nb, n = 0))
  comp <- safe_n_comp_nb(nb)
  for (wh in wh_list) {
    if (n_wh >= MAX_WORMHOLES) break
    i1 <- which.min(sqrt((coords_proj[,1] - wh$from[1])^2 + (coords_proj[,2] - wh$from[2])^2))
    i2 <- which.min(sqrt((coords_proj[,1] - wh$to[1])^2   + (coords_proj[,2] - wh$to[2])^2))
    if (i1 != i2 && comp$comp.id[i1] != comp$comp.id[i2] && !(i2 %in% nb[[i1]])) {
      nb <- nb_set_edge(nb, i1, i2)
      n_wh <- n_wh + 1
      comp <- safe_n_comp_nb(nb)
    }
  }
  list(nb = nb, n = n_wh)
}

build_augmented_adjacency <- function(grid_sf, bridge_links, wormholes, gt_dt = NULL) {
  n_grid <- nrow(grid_sf)
  grid_ids <- grid_sf$grid_id
  utm_crs <- st_crs(grid_sf)

  centroids_proj <- st_centroid(grid_sf)
  coords_proj <- st_coordinates(centroids_proj)

  wh_pts <- do.call(rbind, lapply(wormholes, function(w) rbind(w$from, w$to)))
  wh_sf <- st_as_sf(data.frame(x = wh_pts[,1], y = wh_pts[,2]),
                    coords = c("x", "y"), crs = 4326) |> st_transform(utm_crs)
  wh_utm <- st_coordinates(wh_sf)
  for (i in seq_along(wormholes)) {
    wormholes[[i]]$from <- wh_utm[2*i - 1, 1:2]
    wormholes[[i]]$to   <- wh_utm[2*i, 1:2]
  }

  if (nrow(bridge_links) > 0) {
    br_sf <- st_as_sf(bridge_links[!is.na(lat) & !is.na(lon)],
                      coords = c("lon", "lat"), crs = 4326) |> st_transform(utm_crs)
    br_utm <- st_coordinates(br_sf)
    bridge_links[!is.na(lat) & !is.na(lon), c("utm_x", "utm_y") := .(br_utm[,1], br_utm[,2])]
  }

  suppressWarnings({ nb <- poly2nb(grid_sf, queen = TRUE, snap = SNAP_QUEEN) })
  comp0 <- safe_n_comp_nb(nb)

  wh_critical <- Filter(function(w) w$priority == "critical", wormholes)
  wh_manual   <- Filter(function(w) w$priority != "critical", wormholes)
  res0 <- apply_wormholes(nb, wh_critical, coords_proj, "L0")
  nb <- res0$nb

  n_bridge <- 0
  if (nrow(bridge_links) > 0) {
    for (i in 1:nrow(bridge_links)) {
      bx <- bridge_links$utm_x[i]; by <- bridge_links$utm_y[i]
      if (is.na(bx) || is.na(by)) next
      dists <- sqrt((coords_proj[,1] - bx)^2 + (coords_proj[,2] - by)^2)
      nearest_idx <- order(dists)[1:min(20, n_grid)]
      cur <- safe_n_comp_nb(nb)
      nc <- cur$comp.id[nearest_idx]
      for (j in 1:(length(nearest_idx) - 1)) {
        linked <- FALSE
        for (k in (j + 1):min(length(nearest_idx), j + 10)) {
          if (nc[j] != nc[k]) {
            ia <- nearest_idx[j]; ib <- nearest_idx[k]
            if (!(ib %in% nb[[ia]])) {
              nb <- nb_set_edge(nb, ia, ib); n_bridge <- n_bridge + 1
            }
            linked <- TRUE; break
          }
        }
        if (linked) break
      }
    }
  }
  n_cross <- 0
  cross_log <- list()
  grp <- unit_borough(grid_sf)
  if (!is.null(grp)) {
    unit_label <- if ("rate_geo_id" %in% names(grid_sf)) as.character(grid_sf$rate_geo_id) else as.character(grid_ids)
    pts <- list()
    if (nrow(bridge_links) > 0 && "BORO" %in% names(bridge_links)) {
      for (i in seq_len(nrow(bridge_links))) {
        if (is.na(bridge_links$lat[i]) || is.na(bridge_links$lon[i]) || is.na(bridge_links$BORO[i])) next
        if (grepl(CROSSING_EXCLUDE_PATTERN, as.character(bridge_links[[1]][i]))) next
        pts[[length(pts) + 1]] <- list(name = as.character(bridge_links[[1]][i]), kind = "bridge",
          role = "single", boro = bridge_links$BORO[i], lon = bridge_links$lon[i], lat = bridge_links$lat[i])
      }
    }
    wh_cross <- Filter(function(w) isTRUE(w$crossing), wormholes)
    for (x in c(lapply(NYC_SUBWAY_CROSSINGS, function(x) c(x, list(kind = "subway"))),
                NYC_EXTRA_CROSSINGS,
                lapply(wh_cross, function(w) c(w, list(kind = "wormhole"))))) {
      pts[[length(pts) + 1]] <- list(name = x$name, kind = x$kind, role = "from", boro = NA,
        lon = x$from[1], lat = x$from[2], utm = identical(x$kind, "wormhole"), same_boro = isTRUE(x$same_boro))
      pts[[length(pts) + 1]] <- list(name = x$name, kind = x$kind, role = "to", boro = NA,
        lon = x$to[1], lat = x$to[2], utm = identical(x$kind, "wormhole"), same_boro = isTRUE(x$same_boro))
    }
    if (length(pts)) {
      is_utm <- vapply(pts, function(x) isTRUE(x$utm), logical(1))
      xy <- do.call(rbind, lapply(pts, function(x) c(x$lon, x$lat)))
      p_ll <- st_as_sf(data.frame(x = xy[!is_utm, 1], y = xy[!is_utm, 2]), coords = c("x", "y"), crs = 4326) |> st_transform(utm_crs)
      geom <- vector("list", length(pts))
      geom[!is_utm] <- st_geometry(p_ll)
      if (any(is_utm)) geom[is_utm] <- st_geometry(st_as_sf(data.frame(x = xy[is_utm, 1], y = xy[is_utm, 2]), coords = c("x", "y"), crs = utm_crs))
      p_sf <- st_sf(geometry = st_sfc(geom, crs = utm_crs))
      D <- matrix(as.numeric(st_distance(p_sf, grid_sf)), nrow = length(pts))
      has_grp <- !is.na(grp)
      nearest_unit <- function(k, mask, cap) {
        idx <- which(mask); if (!length(idx)) return(NA_integer_)
        j <- idx[which.min(D[k, idx])]; if (D[k, j] > cap) NA_integer_ else j
      }
      add_cross <- function(nb, ia, ib, name, kind, same_boro_ok = FALSE) {
        if (is.na(ia) || is.na(ib) || ia == ib) return(nb)
        if (grp[ia] == grp[ib] && !same_boro_ok) return(nb)
        dd <- sqrt((coords_proj[ia,1] - coords_proj[ib,1])^2 + (coords_proj[ia,2] - coords_proj[ib,2])^2)
        cross_log[[length(cross_log) + 1]] <<- data.frame(name = name, kind = kind,
          from = grp[ia], to = grp[ib], unit_from = unit_label[ia], unit_to = unit_label[ib],
          centroid_dist_m = round(dd), stringsAsFactors = FALSE)
        if (!(ib %in% nb[[ia]])) { nb <- nb_set_edge(nb, ia, ib); n_cross <<- n_cross + 1 }
        nb
      }
      k <- 1
      while (k <= length(pts)) {
        x <- pts[[k]]
        if (identical(x$role, "single")) {
          ia <- nearest_unit(k, has_grp & grp == x$boro, CROSS_SAME_MAX)
          ib <- nearest_unit(k, has_grp & grp != x$boro, CROSS_OTHER_MAX)
          nb <- add_cross(nb, ia, ib, x$name, x$kind); k <- k + 1
        } else {
          ia <- nearest_unit(k, has_grp, CROSS_SAME_MAX)
          ib <- nearest_unit(k + 1, has_grp, CROSS_SAME_MAX)
          if (!is.na(ia) && !is.na(ib) && grp[ia] == grp[ib] && !isTRUE(x$same_boro))
            cross_log[[length(cross_log) + 1]] <- data.frame(name = x$name, kind = x$kind,
              from = grp[ia], to = grp[ib], unit_from = unit_label[ia], unit_to = unit_label[ib],
              centroid_dist_m = NA_integer_, stringsAsFactors = FALSE)
          nb <- add_cross(nb, ia, ib, x$name, x$kind, same_boro_ok = isTRUE(x$same_boro)); k <- k + 2
        }
      }
    }
  }

  n_snap <- 0
  for (ti in seq_along(SNAP_TIERS)) {
    snap_dist <- min(SNAP_TIERS[ti], SNAP_MAX)
    cur <- safe_n_comp_nb(nb)
    if (cur$nc <= 1) break
    cs <- table(cur$comp.id)
    lg_id <- as.integer(names(which.max(cs)))
    main_g <- which(cur$comp.id == lg_id)
    small_ids <- as.integer(names(cs[cs < 50 & names(cs) != as.character(lg_id)]))
    if (length(small_ids) == 0) break
    n_this <- 0
    for (cid in small_ids) {
      frag <- which(cur$comp.id == cid)
      frag_s <- frag[seq(1, length(frag), length.out = min(10, length(frag)))]
      nn <- nn2(coords_proj[main_g, , drop = FALSE], coords_proj[frag_s, , drop = FALSE], k = 1)
      bi <- which.min(nn$nn.dists)
      if (nn$nn.dists[bi] <= snap_dist) {
        ia <- frag_s[bi]; ib <- main_g[nn$nn.idx[bi]]
        if (!(ib %in% nb[[ia]])) { nb <- nb_set_edge(nb, ia, ib); n_snap <- n_snap + 1; n_this <- n_this + 1 }
      }
    }
    if (n_this == 0 && ti > 2) break
  }

  n_rail <- 0
  subway_loaded <- FALSE
  if (file.exists(SUBWAY_ROUTES_SHP)) {
    tryCatch({
      subway <- st_transform(st_read(SUBWAY_ROUTES_SHP, quiet = TRUE), st_crs(grid_sf))
      hits <- st_intersects(grid_sf, st_buffer(subway, 75), sparse = TRUE)
      cur <- safe_n_comp_nb(nb)
      for (li in 1:nrow(subway)) {
        touching <- which(lengths(lapply(hits, function(x) x[x == li])) > 0)
        if (length(touching) < 2) next
        lc <- cur$comp.id[touching]; uc <- unique(lc)
        if (length(uc) <= 1) next
        for (ci in 2:length(uc)) {
          g1 <- touching[which(lc == uc[1])[1]]
          g2 <- touching[which(lc == uc[ci])[1]]
          if (!(g2 %in% nb[[g1]])) { nb <- nb_set_edge(nb, g1, g2); n_rail <- n_rail + 1 }
        }
        if (li %% 50 == 0) cur <- safe_n_comp_nb(nb)
      }
      subway_loaded <- TRUE
    }, error = function(e) NULL)
  }
  if (!subway_loaded && exists("SUBWAY_STATIONS_CSV") && file.exists(SUBWAY_STATIONS_CSV)) {
    tryCatch({
      mta <- fread(SUBWAY_STATIONS_CSV)
      mta <- mta[!is.na(`GTFS Latitude`) & !is.na(`GTFS Longitude`)]
      mta_sf <- st_transform(st_as_sf(mta, coords = c("GTFS Longitude", "GTFS Latitude"), crs = 4326), st_crs(grid_sf))
      mta_coords <- st_coordinates(mta_sf)
      cur <- safe_n_comp_nb(nb)
      for (ln in unique(mta$Line)) {
        ln_idx <- which(mta$Line == ln)
        if (length(ln_idx) < 2) next
        grids_on_line <- vapply(ln_idx, function(si) {
          which.min(sqrt((coords_proj[,1] - mta_coords[si,1])^2 + (coords_proj[,2] - mta_coords[si,2])^2))
        }, integer(1))
        for (si in 1:(length(grids_on_line) - 1)) {
          ia <- grids_on_line[si]; ib <- grids_on_line[si + 1]
          if (ia != ib && cur$comp.id[ia] != cur$comp.id[ib]) {
            nb <- nb_set_edge(nb, ia, ib); n_rail <- n_rail + 1
          }
        }
      }
    }, error = function(e) NULL)
  }

  n_road <- 0
  road_files <- if (file.exists(ROAD_SHP)) ROAD_SHP else get_road_paths(state)
  for (rf in road_files) {
    if (!file.exists(rf)) next
    tryCatch({
      roads_sf <- st_transform(st_read(rf, quiet = TRUE), st_crs(grid_sf))
      mtfcc_col <- intersect(c("MTFCC", "mtfcc"), names(roads_sf))[1]
      if (!is.na(mtfcc_col)) roads_major <- roads_sf[roads_sf[[mtfcc_col]] %in% TIGER_BRIDGE_TYPES, ]
      else roads_major <- roads_sf
      comp_pre <- safe_n_comp_nb(nb)
      if (comp_pre$nc > 1 && nrow(roads_major) > 0) {
        road_hits <- st_intersects(roads_major, grid_sf)
        for (r in seq_along(road_hits)) {
          tidx <- road_hits[[r]]
          if (length(tidx) < 2) next
          tcomps <- comp_pre$comp.id[tidx]; ucomps <- unique(tcomps)
          if (length(ucomps) < 2) next
          for (ci in 1:(length(ucomps) - 1)) {
            g_c1 <- tidx[tcomps == ucomps[ci]]; g_c2 <- tidx[tcomps == ucomps[ci + 1]]
            best_d <- Inf; best_i <- NA; best_j <- NA
            for (gi in g_c1[1:min(5, length(g_c1))]) {
              for (gj in g_c2[1:min(5, length(g_c2))]) {
                d <- sqrt((coords_proj[gi,1] - coords_proj[gj,1])^2 + (coords_proj[gi,2] - coords_proj[gj,2])^2)
                if (d < best_d) { best_d <- d; best_i <- gi; best_j <- gj }
              }
            }
            if (!is.na(best_i) && best_i != best_j && !(best_j %in% nb[[best_i]])) {
              nb <- nb_set_edge(nb, best_i, best_j); n_road <- n_road + 1
              comp_pre$comp.id[comp_pre$comp.id == ucomps[ci + 1]] <- ucomps[ci]
            }
          }
        }
      }
    }, error = function(e) NULL)
  }

  n_od <- 0
  comp2 <- safe_n_comp_nb(nb)
  if (!is.null(gt_dt) && nrow(gt_dt) > OD_MIN_ROWS && comp2$nc > 1) {
    cs <- table(comp2$comp.id)
    lg_id <- as.integer(names(which.max(cs)))
    main_idx <- which(comp2$comp.id == lg_id)
    main_gids <- grid_ids[main_idx]

    isolate_grids <- integer(0)
    for (fc in unique(comp2$comp.id[comp2$comp.id != lg_id])) {
      fc_idx <- which(comp2$comp.id == fc)
      fc_deg <- mean(sapply(fc_idx, function(i) length(nb[[i]][nb[[i]] > 0L])))
      if (length(fc_idx) <= OD_MAX_FRAG_SIZE || fc_deg < OD_MAX_FRAG_DEGREE)
        isolate_grids <- c(isolate_grids, fc_idx)
    }
    isolate_gids <- grid_ids[isolate_grids]

    rho_col <- if ("log_rho" %in% names(gt_dt)) "log_rho" else "Y"
    act_col <- if ("rho_norm" %in% names(gt_dt)) "rho_norm" else "Y"
    gstats <- gt_dt[, .(total_act = sum(get(act_col), na.rm = TRUE),
                         n_days = sum(get(act_col) > 0, na.rm = TRUE)), by = grid_id]
    elig <- gstats[grid_id %in% isolate_gids & total_act >= OD_MIN_ACTIVITY & n_days >= OD_MIN_DAYS, grid_id]

    if (length(elig) > 0) {
      fp <- dcast(gt_dt[grid_id %in% elig], grid_id ~ time_id, value.var = rho_col, fill = 0, fun.aggregate = mean)
      fg <- fp$grid_id; fm <- as.matrix(fp[, -1, with = FALSE])
      ms <- gstats[grid_id %in% main_gids][order(-total_act)]
      ms_ids <- ms$grid_id[1:min(300, nrow(ms))]
      mp <- dcast(gt_dt[grid_id %in% ms_ids], grid_id ~ time_id, value.var = rho_col, fill = 0, fun.aggregate = mean)
      mg <- mp$grid_id; mm <- as.matrix(mp[, -1, with = FALSE])
      for (fi in 1:nrow(fm)) {
        fgid <- fg[fi]; fidx <- which(grid_ids == fgid)
        if (length(fidx) == 0) next
        fv <- fm[fi, ]; if (sd(fv, na.rm = TRUE) < 1e-10) next
        cors <- apply(mm, 1, function(mv) if (sd(mv, na.rm = TRUE) < 1e-10) 0 else cor(fv, mv, use = "complete.obs"))
        bi <- which.max(cors)
        if (!is.na(cors[bi]) && cors[bi] > OD_CORR_THRESHOLD) {
          midx <- which(grid_ids == mg[bi])
          if (length(midx) > 0 && !(midx %in% nb[[fidx]])) { nb <- nb_set_edge(nb, fidx, midx); n_od <- n_od + 1 }
        }
      }
    }
  }
  n_wh_l3 <- 0
  if (length(wh_manual) > 0 && safe_n_comp_nb(nb)$nc > 1) {
    res3 <- apply_wormholes(nb, wh_manual, coords_proj, "L3")
    nb <- res3$nb; n_wh_l3 <- res3$n
  }
  n_abs <- 0
  for (pass in 1:10) {
    sc <- safe_n_comp_nb(nb)
    if (sc$nc <= 1) break
    cs <- table(sc$comp.id)
    lg_id <- as.integer(names(which.max(cs)))
    main_g <- which(sc$comp.id == lg_id)
    n_this <- 0
    for (cid in unique(sc$comp.id[sc$comp.id != lg_id])) {
      frag <- which(sc$comp.id == cid)
      frag_s <- frag[seq(1, length(frag), length.out = min(5, length(frag)))]
      nn <- nn2(coords_proj[main_g, , drop = FALSE], coords_proj[frag_s, , drop = FALSE], k = 1)
      bi <- which.min(nn$nn.dists)
      ia <- frag_s[bi]; ib <- main_g[nn$nn.idx[bi]]
      if (!(ib %in% nb[[ia]])) { nb <- nb_set_edge(nb, ia, ib); n_abs <- n_abs + 1; n_this <- n_this + 1 }
    }
    if (n_this == 0) break
  }

  n_isl <- 0
  zeros <- which(card(nb) == 0)
  if (length(zeros) > 0) {
    others <- setdiff(1:n_grid, zeros)
    for (iso in zeros) {
      nn <- nn2(coords_proj[others, , drop = FALSE], matrix(coords_proj[iso, ], nrow = 1), k = 1)
      nb <- nb_set_edge(nb, iso, others[nn$nn.idx[1]]); n_isl <- n_isl + 1
    }
  }
  comp_f <- safe_n_comp_nb(nb)
  list(nb = nb, comp = comp_f, comp_before = comp0$nc, comp_after = comp_f$nc,
       n_wh_l0 = res0$n, n_bridge = n_bridge, n_cross = n_cross,
       crossings = if (length(cross_log)) do.call(rbind, cross_log) else NULL,
       n_snap = n_snap, n_rail = n_rail,
       n_road = n_road, n_od = n_od, n_wormhole = res0$n + n_wh_l3,
       n_absorbed = n_abs, n_island_fix = n_isl, aug_version = ADJ_AUG_VERSION)
}

# Transport layers used to connect the neighbour graph; set per state before building it.
TIGER_BRIDGE_TYPES <- c("S1100", "S1200")
set_transport_globals <- function(st) {
  if (st == "NY") {
    SUBWAY_ROUTES_SHP   <<- CONFIG$subway_routes_shp
    SUBWAY_STATIONS_CSV <<- CONFIG$subway_stations_csv
    ROAD_SHP            <<- CONFIG$ny_roads_shp
    if (file.exists(CONFIG$bridge_csv)) {
      bridges <- fread(CONFIG$bridge_csv, na.strings = c("", "NA"))
      setnames(bridges, c("X - COORD  (LAT)", "Y - COORD   (LON)"), c("lat", "lon"), skip_absent = TRUE)
      bridge_links <<- bridges[nchar(BORO) >= 2 & !is.na(lat) & !is.na(lon)]
    } else bridge_links <<- data.table(lat = numeric(), lon = numeric())
  } else {
    SUBWAY_ROUTES_SHP   <<- ""
    SUBWAY_STATIONS_CSV <<- ""
    ROAD_SHP            <<- ""
    bridge_links <<- data.table(lat = numeric(), lon = numeric())
  }
  WORMHOLES <<- if (st == "NY") NYC_WORMHOLES else UT_WORMHOLES
  state <<- st
}

# Write the augmented unit graph for INLA and return it as read back from the file; BYM2 needs one connected component.
write_unit_graph <- function(grid_sf, gt_dt, graph_file, crossings_file = NULL) {
  aug <- build_augmented_adjacency(grid_sf, bridge_links, WORMHOLES, gt_dt = gt_dt)
  nc <- safe_n_comp_nb(aug$nb)$nc
  if (nc > 1) stop(sprintf("the unit graph for %s has %d components", graph_file, nc))
  safe_nb2INLA(graph_file, aug$nb)
  if (!is.null(crossings_file) && !is.null(aug$crossings)) fwrite(aug$crossings, crossings_file)
  read_inla_graph_as_nb(graph_file)
}
