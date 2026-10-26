# Journal figures (figures/) and tables (tables/) from the shared results in results/.
# Run from the repository root: Rscript make_figures.R
suppressPackageStartupMessages({
  library(data.table); library(sf); library(ggplot2); library(patchwork)
  library(scales); library(scico); library(surveillance)
})

res <- function(f, ...) fread(file.path("results", f), ...)
for (d in c("figures", "tables")) dir.create(d, showWarnings = FALSE)
save_fig <- function(p, name, w, h, dpi = 600) ggsave(file.path("figures", name), p, width = w, height = h, dpi = dpi, bg = "white")
save_tab <- function(x, name) fwrite(x, file.path("tables", name))

WAVE_COLS <- c("Wave 1" = "#1B9E77", "Wave 2" = "#D95F02", "Wave 3" = "#7570B3")
wave_f <- function(w) factor(paste("Wave", w), names(WAVE_COLS))
SETTING <- c(NY = "NYC (100 m grid)", UT = "Utah (counties)")
SITE <- c(NY = "New York City", UT = "Utah")
POP <- c(NY = 8804190, UT = 3271616)   # 2020 Census resident population
MODEL <- c(Base = "Base", Temp = "Non-spatial", LinMob = "Linear-mobility", NLMob = "Nonlinear-mobility",
           LinMob_ST_IV = "Linear + interaction", NLMob_ST_IV = "Nonlinear + interaction")
TERM <- c(own_prev_std = "Local infection pressure", ext_prev_std = "Imported exposure",
          log_rho_resid = "Same-day residual mobility", mob_weighted = "Lagged residual mobility",
          mob_spatial_lag_std = "Spatial mobility")
COVARIATE <- c(TERM, log_repeat_ratio_std = "Repeat-visit ratio", building_density_std = "Building density",
               log_contact_density_hourly_std = "Contact density (hourly)", log_daytime_pop_std = "Daytime population density",
               landuse_commercial = "Commercial land use", commercial_x_weekend = "Commercial land use x weekend",
               ses_contrast_std = "SES contrast", pct_younger_std = "Younger-population share",
               pct_older_std = "Older-population share", pct_unemployed_std = "Unemployment rate",
               pct_employed_std = "Employment rate", is_weekend = "Weekend", is_holiday = "Holiday")
borough <- function(g) {
  p <- substr(g, 1, 3)
  fcase(p %in% c("100", "101", "102"), "Manhattan", p == "103", "Staten Island", p == "104", "Bronx",
        p == "112", "Brooklyn", p %in% c("110", "111", "113", "114", "116"), "Queens", default = "Other")
}
theme_map <- function(...) theme_void(base_size = 7) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(b = 3)),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4),
        legend.text = element_text(size = 6), legend.title = element_text(size = 7), ...)
theme_panels <- function(...) theme_bw(base_size = 16) +
  theme(strip.text = element_text(face = "bold", size = 15), axis.text = element_text(size = 12),
        axis.title = element_text(size = 16), panel.grid.minor = element_line(linewidth = 0.2, colour = "grey92"),
        panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"), ...)
# Diverging relative-risk scale on log(RR): clipped symmetrically at the 1st/99th percentiles (all values
# below 500 units), with ticks at constant steps in log(RR).
rr_scale <- function(rr, step = NULL, barheight = 2.6) {
  q <- if (length(rr) < 500) range(rr) else quantile(rr, c(0.01, 0.99))
  clip <- max(abs(log(q)))
  if (is.null(step)) step <- if (clip < 0.5) 0.1 else if (clip < 1) 0.25 else 0.5
  n <- ceiling(clip / step)
  scale_fill_gradientn(colours = scico(256, palette = "vik"), name = "RR", limits = c(-n, n) * step, oob = squish,
    breaks = seq(-n, n) * step, labels = function(x) signif(exp(x), 2),
    guide = guide_colourbar(ticks.colour = "white", ticks.linewidth = 0.4, frame.colour = NA, draw.ulim = FALSE,
                            draw.llim = FALSE, barheight = unit(barheight, "cm"), barwidth = unit(0.32, "cm")))
}

modzcta <- st_read("results/boundaries.gpkg", "modzcta", quiet = TRUE)
counties <- st_read("results/boundaries.gpkg", "utah_counties", quiet = TRUE)
cells <- res("nyc_grid_cells.csv.gz", colClasses = list(character = "geo_id"))
sp <- res("spatial_effects.csv.gz")[, log_rr := mean]
te <- res("temporal_effects.csv")
hp <- res("hyperparameters.csv")
selected <- res("selected_models.csv")

# ---- Figure 1: study areas ----
modzcta$borough <- borough(modzcta$geo_id)
counties$class <- ifelse(toupper(counties$name) %in% c("SALT LAKE", "UTAH", "DAVIS", "WEBER", "CACHE"), "Urban", "Rural")
area_map <- function(d, fill, pal, title, border, lw) ggplot(d) + geom_sf(aes(fill = .data[[fill]]), colour = border, linewidth = lw) +
  scale_fill_manual(values = pal, name = NULL) + labs(title = title) + coord_sf(datum = NA) +
  theme_map(legend.position = "bottom", legend.key.size = unit(0.55, "cm"))
save_fig((area_map(modzcta, "borough", c(Manhattan = "#496A9A", Brooklyn = "#6F91B5", Queens = "#6F9B9C", Bronx = "#A8C2CD",
                                         `Staten Island` = "#ADA5C2"), "New York City", "white", 0.25) |
         area_map(counties, "class", c(Urban = "#287B91", Rural = "#DCE6EB"), "Utah", "grey40", 0.3)) + plot_layout(widths = c(1, 1)),
         "fig01_study_areas.png", 7.2, 4.11)

# ---- Figure 2: weekly cases (and POI visits where the licensed visit data are present) ----
wk <- res("weekly_cases.csv")[, y := cases_per_week / POP[setting] * 1e5]
visits <- rbindlist(lapply(c("NY", "UT"), function(st) rbindlist(lapply(1:3, function(w) {
  f <- file.path("outputs", sprintf("%s_wave%d", st, w), "tables", "temporal_covariates.csv")
  if (!file.exists(f)) return(NULL)
  tc <- fread(f)[, date := as.Date(date)]
  tc[, week := if (st == "NY") match(date + (6L - as.POSIXlt(date)$wday), sort(unique(date + (6L - as.POSIXlt(date)$wday))))
               else (sim_day - 1L) %/% 7L + 1L]
  tc[, .(visits = 7 * sum(total_visits) / .N), by = week][, `:=`(setting = st, wave = w)]
}))))
if (nrow(visits)) wk <- merge(wk, visits, by = c("setting", "wave", "week"), all.x = TRUE)
panel_ts <- function(st, what, x_title) {
  d <- wk[setting == st][, v := if (what == "cases") y else visits / POP[st]][v > 0]
  p <- ggplot(d, aes(day, v, colour = wave_f(wave))) + geom_line(linewidth = 0.6) + geom_point(size = 0.9) +
    scale_colour_manual(values = WAVE_COLS, name = NULL) + scale_x_continuous(breaks = seq(0, 90, 15), limits = c(0, 92)) +
    labs(x = if (x_title) "Day of wave" else NULL, title = if (what == "cases") SITE[[st]] else NULL,
         y = if (st == "UT") NULL else if (what == "cases") "Reported cases per 100,000\nper week (log scale)" else "POI visits per resident\nper week") +
    theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 9, hjust = 0.5))
  if (what == "cases") p + scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c("0.1", "1", "10", "100", "1,000"), limits = c(0.04, 4000))
  else p + scale_y_continuous(limits = c(0, 50), expand = expansion(mult = c(0, 0.02)))
}
has_visits <- "visits" %in% names(wk)
p <- panel_ts("NY", "cases", !has_visits) | panel_ts("UT", "cases", !has_visits)
if (has_visits) p <- p / (panel_ts("NY", "visits", TRUE) | panel_ts("UT", "visits", TRUE)) + plot_layout(heights = c(1.15, 1))
save_fig(p + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "fig02_cases_visits.png", 7.2, if (has_visits) 5.2 else 2.9)

# ---- Figure 3: allocation of MODZCTA counts to cells (East Harlem, Wave 1) and the allocation check ----
zone <- cells[geo_id == "10029"][is.na(cases_wave1), cases_wave1 := 0L]
lim <- function(v) { r <- range(v) + c(-50, 50); r + c(-1, 1) * 0.03 * diff(r) }
theme_alloc <- theme_map(legend.position = "right", legend.margin = margin(l = 10), legend.key.height = unit(0.45, "cm")) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(b = 4)))
zone_map <- function(fill, title, scale) ggplot() + geom_tile(data = zone, aes(x, y, fill = .data[[fill]]), width = 100, height = 100) +
  geom_sf(data = modzcta[modzcta$geo_id == "10029", ], fill = NA, colour = "grey20", linewidth = 0.5) + scale +
  labs(title = title) + coord_sf(datum = NA, xlim = lim(zone$x), ylim = lim(zone$y), expand = FALSE) + theme_alloc
p_in <- zone_map("E", "Population weight (input)", scale_fill_viridis_c(option = "mako", direction = -1, trans = "sqrt",
                 name = expression(E[g]), breaks = c(2, 75, 250, 525, 900)))
p_out <- zone_map("cases_wave1", "Allocated cases (output)", scale_fill_viridis_c(option = "rocket", direction = -1, trans = "log1p",
                  name = expression(Sigma * Y[g]), breaks = unique(round(expm1(seq(0, log1p(max(zone$cases_wave1)), length.out = 5))))))
brk <- c(10, 30, 100, 300, 1000, 3000)
p_chk <- ggplot(res("allocation_check.csv"), aes(rate_reported, rate_allocated)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#d7191c", linewidth = 0.5) +
  geom_point(alpha = 0.25, size = 0.6, colour = "#2c7bb6") +
  scale_x_log10(labels = comma, breaks = brk) + scale_y_log10(labels = comma, breaks = brk) +
  labs(title = "Allocation check, all MODZCTA-weeks", x = "Reported weekly rate per 100,000", y = "Rate from allocated cases") +
  theme_minimal(base_size = 7) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8, margin = margin(b = 4)), axis.text = element_text(size = 6),
        axis.title.x = element_text(size = 7, margin = margin(t = 4)), axis.title.y = element_text(size = 7, margin = margin(r = 4)),
        panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.4))
save_fig((p_in | p_out | p_chk) + plot_layout(widths = c(1, 1, 1.15)) & theme(plot.margin = margin(6, 12, 6, 12)),
         "fig03_case_allocation.png", 7.2, 2.46)

# ---- Figures 4-5: spatial relative risk, one colour scale for both settings ----
rr_fill <- rr_scale(exp(sp$mean))
ny_maps <- lapply(1:3, function(w) ggplot() + geom_sf(data = modzcta, fill = NA, colour = "grey60", linewidth = 0.25) +
  geom_tile(data = merge(cells, sp[setting == "NY" & wave == w], by = "grid_id"), aes(x, y, fill = log_rr), width = 100, height = 100) +
  rr_fill + coord_sf(datum = NA) + labs(title = paste("Wave", w)) + theme_map(legend.position = if (w == 3) "right" else "none"))
save_fig(wrap_plots(ny_maps, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "right"),
         "fig04_spatial_rr_nyc.png", 7.2, 5.5 * 7.2 / 14)
ut_maps <- lapply(1:3, function(w) {
  d <- merge(counties, sp[setting == "UT" & wave == w], by = "grid_id")
  d$alpha <- ifelse(d$q025 > 0 | d$q975 < 0, 1, 0.2)
  ggplot() + geom_sf(data = counties, fill = "grey95", colour = "grey60", linewidth = 0.3, alpha = 0.3) +
    geom_sf(data = d, aes(fill = log_rr, alpha = alpha), colour = "grey40", linewidth = 0.3) + rr_fill + scale_alpha_identity() +
    coord_sf(datum = NA) + labs(title = paste("Wave", w)) + theme_map(legend.position = if (w == 3) "right" else "none")
})
save_fig(wrap_plots(ut_maps, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "right"),
         "fig05_spatial_rr_utah.png", 7.2, 6 * 7.2 / 14)

# ---- Figure 6 and Supplementary Figure S2: Type IV interaction ----
st <- res("st_interaction.csv.gz")
interaction <- function(s, w) {
  d <- st[setting == s & wave == w]
  K <- uniqueN(d$idx_st)
  stopifnot(nrow(d) / K == te[setting == s & wave == w, uniqueN(idx_temporal)])   # grouped by day
  d[, `:=`(cluster = idx_st, day = (seq_len(.N) - 1L) %/% K + 1L)]
  d[, demeaned := mean - mean(mean), by = cluster]
  prec <- hp[setting == s & wave == w & grepl("idx_st_cluster|Precision for idx_st$", parameter), mean][1]
  degenerate <- prec > 1e5
  ord <- if (degenerate) d[, .(v = mean(mean)), by = cluster][order(-v), cluster]
         else d[, .(pk = day[which.max(demeaned)]), by = cluster][order(pk), cluster]
  d[, rank := match(cluster, rev(ord))]
  list(d = d, K = K, degenerate = degenerate)
}
ia <- unlist(lapply(c("NY", "UT"), function(s) lapply(1:3, function(w) c(interaction(s, w), setting = s, wave = w))), recursive = FALSE)
vik_mid <- function(q, name, breaks) scale_fill_gradientn(colours = scico(256, palette = "vik"), rescaler = ~ rescale_mid(.x, mid = 0),
  limits = c(-q, q), breaks = breaks[abs(breaks) <= q + 1e-9], labels = label_number(accuracy = 0.1), name = name)
x <- ia[[3]]   # NYC Wave 3
v <- if (x$degenerate) x$d$mean else x$d$demeaned
q <- max(quantile(abs(v), 0.99), 1e-6)
p6 <- ggplot(x$d[, val := pmin(pmax(v, -q), q)], aes(day, factor(rank), fill = val)) + geom_tile(width = 1) +
  vik_mid(q, if (x$degenerate) "ST effect (log-RR)" else "Deviation from\ncluster mean", pretty(c(-q, q), n = 5)) +
  scale_x_continuous(breaks = seq(0, 90, 15), expand = c(0, 0)) + labs(x = "Day of wave", y = "Spatial cluster") +
  theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.key.width = unit(0.9, "cm"), legend.key.height = unit(0.25, "cm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 7, margin = margin(r = 12)), axis.title = element_text(size = 7),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 6))
save_fig(p6, "fig06_interaction_nyc_wave3.png", 7.2, (if (x$K > 60) 8 else 7) * 7.2 / 12)
q <- max(quantile(abs(unlist(lapply(ia, function(x) x$d$demeaned))), 0.99), 1e-6)
all_ia <- rbindlist(lapply(ia, function(x) x$d[, .(day, rank, val = pmin(pmax(demeaned, -q), q), wave = wave_f(x$wave),
  site = if (x$setting == "NY") sprintf("NYC (%d clusters)", x$K) else "Utah (counties)")]))
pS2 <- ggplot(all_ia, aes(day, rank, fill = val)) + geom_tile(width = 1) +
  facet_grid(site ~ wave, scales = "free_y", space = "free_y") +
  vik_mid(q, "Deviation from the unit mean (log RR)", pretty(c(-q, q), n = 5)) +
  scale_x_continuous(breaks = seq(0, 90, 15), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Day of wave", y = "Spatial units, ordered by peak day") + theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        strip.text = element_text(face = "bold", size = 11), legend.position = "bottom", legend.key.width = unit(1.6, "cm"),
        legend.key.height = unit(0.3, "cm"), legend.title = element_text(size = 10, margin = margin(r = 10)), panel.spacing = unit(0.6, "lines"))
save_fig(pS2, "figS2_interaction_all_waves.png", 11, 6.5, dpi = 400)

# ---- Figure 7: temporal relative risk, with the week around the peak of reported cases shaded ----
trr <- te[idx_temporal > 1][!(setting == "NY" & wave == 1 & idx_temporal < 11)]   # NYC reporting starts 11 March 2020
peaks <- res("case_rate_daily_total.csv")[order(setting, wave, day)][, roll := frollmean(rate, 7, align = "center"), by = .(setting, wave)][!is.na(roll), .(peak = day[which.max(roll)]), by = .(setting, wave)]
facets <- function(d) d[, `:=`(state = factor(SETTING[setting], SETTING), wave_f = wave_f(wave))]
p7 <- ggplot(facets(trr), aes(idx_temporal, exp(mean), colour = wave_f, fill = wave_f)) +
  geom_rect(data = facets(peaks), inherit.aes = FALSE, aes(xmin = peak - 3, xmax = peak + 3, fill = wave_f), ymin = -Inf, ymax = Inf, alpha = 0.14) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_ribbon(aes(ymin = exp(q025), ymax = exp(q975)), alpha = 0.16, colour = NA) + geom_line(linewidth = 0.8) +
  facet_grid(state ~ wave_f, scales = "free") +
  scale_colour_manual(values = unname(WAVE_COLS), guide = "none") + scale_fill_manual(values = unname(WAVE_COLS), guide = "none") +
  scale_x_continuous(breaks = seq(0, 90, 15), expand = c(0.02, 0)) +
  scale_y_log10(breaks = c(0.03, 0.1, 0.3, 1, 3, 10, 30, 100), labels = function(x) format(x, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)) +
  labs(x = "Day of wave", y = "Temporal relative risk (log scale)") + theme_panels(legend.position = "none")
save_fig(p7, "fig07_temporal_rr.png", 15, 1.5 + 4.2 * 2)

# ---- Figure 8: nonlinear activity effect ----
ae <- res("activity_effect.csv")[, wave_l := wave_f(wave)]
site_row <- function(st, x_title) {
  p <- ggplot(ae[setting == st], aes(z, mean, colour = wave_l, fill = wave_l)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.25, colour = NA) + geom_line(linewidth = 0.9) + geom_point(size = 1.3) +
    facet_wrap(vars(wave_l), nrow = 1) +
    scale_colour_manual(values = WAVE_COLS, guide = "none") + scale_fill_manual(values = WAVE_COLS, guide = "none") +
    labs(x = if (x_title) expression("Visit activity, " * italic(z)[italic(g) * "," * italic(t)]) else NULL, y = "Effect on log-RR") +
    theme_bw(base_size = 16) +
    theme(strip.text = element_text(face = "bold", size = 15), strip.background = element_rect(fill = "grey94", colour = NA),
          panel.grid.minor = element_blank(), panel.spacing = unit(1.1, "lines"), axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))
  (p | wrap_elements(full = grid::textGrob(SITE[[st]], rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 14)))) +
    plot_layout(widths = c(1, 0.055))
}
save_fig(site_row("NY", FALSE) / site_row("UT", TRUE), "fig08_dose_response.png", 11.5, 6.4)

# ---- Figure 9 and Supplementary Data 1: AMLP by POI category ----
pc <- res("pcatm_categories.csv")
save_tab(pc[, .(setting = fifelse(setting == "NY", "NYC", "Utah"), wave, term = TERM[term], category,
                amlp_contribution = signif(amlp, 6), term_total_amlp = signif(term_total, 6))][order(setting, wave, term, category)],
         "supplementary_data_1_amlp_by_category.csv")
ct <- pc[, .(amlp = sum(amlp)), by = .(setting, wave, category)][, pct := 100 * amlp / sum(amlp), by = .(setting, wave)]
ord <- ct[!category %in% c("Other", "Unallocated"), mean(pct), by = category][order(V1), category]
ct[, `:=`(category = factor(category, c("Unallocated", "Other", ord)), area = factor(SITE[setting], SITE), wave_l = wave_f(wave),
          lab = sprintf("%.1f%%", ifelse(abs(pct) < 0.05, 0, pct)))]
ct[, key := ifelse(category %in% c("Other", "Unallocated"), "grey", as.character(wave_l))]
p9 <- ggplot(ct, aes(pct, category, fill = key)) + geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) + geom_col(width = 0.7) +
  geom_text(aes(x = ifelse(pct >= 0, pct + 1, pct - 1), label = lab, hjust = ifelse(pct >= 0, 0, 1)), size = 2.6) +
  facet_grid(area ~ wave_l) + scale_fill_manual(values = c(WAVE_COLS, grey = "#9E9E9E"), guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 60, 20), expand = expansion(mult = c(0.22, 0.28))) +
  labs(x = "Share of the total mobility contribution", y = NULL) + theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), strip.text = element_text(face = "bold", size = 10.5),
        strip.background = element_rect(fill = "grey88", colour = NA), axis.text.y = element_text(size = 9))
save_fig(p9, "fig09_pcatm_categories.png", 10.5, 6.8)

# ---- Supplementary Figure S1: 100 m spatial relative risk in four MODZCTAs ----
zones <- c("10029", "11372", "10459", "11225")
zc <- merge(cells[geo_id %in% zones, .(grid_id, x, y, geo_id)], sp[setting == "NY"], by = "grid_id")
clip_s1 <- max(abs(log(quantile(exp(zc$mean), c(0.01, 0.99)))))
fill_s1 <- rr_scale(exp(zc$mean), step = if (clip_s1 < 1) 0.25 else 0.5, barheight = 3.2)
panels <- unlist(lapply(zones, function(z) lapply(1:3, function(w) {
  d <- zc[geo_id == z & wave == w]; b <- modzcta[modzcta$geo_id == z, ]; bb <- st_bbox(b)
  xr <- range(c(d$x - 50, d$x + 50, bb[c("xmin", "xmax")])); yr <- range(c(d$y - 50, d$y + 50, bb[c("ymin", "ymax")]))
  ggplot() + geom_sf(data = b, fill = NA, colour = "grey15", linewidth = 0.8) +
    geom_tile(data = d, aes(x, y, fill = log_rr), width = 100, height = 100) + fill_s1 +
    coord_sf(datum = NA, xlim = xr + c(-1, 1) * 0.04 * diff(xr), ylim = yr + c(-1, 1) * 0.04 * diff(yr), expand = FALSE) +
    labs(title = if (z == zones[1]) paste("Wave", w) else NULL, y = if (w == 1) sprintf("%s (%s)", borough(z), z) else NULL) +
    theme_map(legend.position = if (z == tail(zones, 1) && w == 3) "right" else "none",
              axis.title.y = if (w == 1) element_text(size = 7, face = "bold", angle = 90, margin = margin(r = 8)) else element_blank())
})), recursive = FALSE)
save_fig(wrap_plots(panels, ncol = 3, byrow = TRUE) + plot_layout(guides = "collect") & theme(legend.position = "right"),
         "figS1_spatial_rr_100m_zones.png", 7.2, 8.64)

# ---- Supplementary Figure S3: PIT histograms and calibration tests of the temporal forecasts ----
ho <- merge(res("temporal_cv_holdout.csv.gz"), res("temporal_cv_folds.csv"), by = c("setting", "wave", "origin"))
ho <- ho[is.finite(mu) & is.finite(y_true) & !is.na(family)][, x := as.integer(round(y_true))]
ho[, size := if (family[1] == "nbinomial") nb_size else NA_real_, by = .(setting, wave, origin)]
calib <- lapply(split(ho, by = c("setting", "wave")), function(h) {
  J <- 10L; mass <- numeric(J)
  for (fam in unique(h$family)) {
    g <- h[family == fam]
    p <- if (fam == "nbinomial") pit(g$x, pdistr = "pnbinom", mu = g$mu, size = g$size, J = J, plot = NULL)
         else pit(g$x, pdistr = "ppois", lambda = g$mu, J = J, plot = NULL)
    mass <- mass + p$density / J * nrow(g)
  }
  nb <- h$family == "nbinomial"
  ct <- if (all(nb)) calibrationTest(h$x, mu = h$mu, size = h$size, which = "dss")
        else if (!any(nb)) calibrationTest(h$x, mu = h$mu, which = "dss")
  sc <- rbindlist(lapply(split(h, by = "origin"), function(g) {
    s <- if (g$family[1] == "nbinomial") scores(g$x, mu = g$mu, size = g$size, which = c("logs", "rps"), individual = TRUE)
         else scores(g$x, mu = g$mu, which = c("logs", "rps"), individual = TRUE)
    data.table(origin = g$origin[1], logS = mean(s[, "logs"]), RPS = mean(s[, "rps"]))
  }))
  lab <- h[1, .(site = factor(fifelse(setting == "NY", "NYC", "Utah"), c("NYC", "Utah")), wave = wave_f(wave))]
  list(pit = cbind(lab, data.table(mid = (seq_len(J) - 0.5) / J, density = mass / nrow(h) * J)),
       test = cbind(lab, data.table(n = nrow(h), z = if (is.null(ct)) NA_real_ else unname(ct$statistic),
                                    p = if (is.null(ct)) NA_real_ else ct$p.value, below_half = sum(mass[1:5]) / nrow(h))),
       scores = cbind(lab, sc[order(origin)]))
})
pick <- function(k) rbindlist(lapply(calib, `[[`, k))
save_tab(pick("test"), "pit_calibration_tests.csv")
save_tab(pick("scores"), "forecast_scores_by_origin.csv")
pS3 <- ggplot(pick("pit"), aes(mid, density, fill = wave)) + geom_col(width = 0.1, colour = "white", linewidth = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30", linewidth = 0.4) + facet_grid(site ~ wave) +
  scale_fill_manual(values = WAVE_COLS, guide = "none") + scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) + labs(x = "PIT (non-randomised, 10 bins)", y = "Relative frequency") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), strip.text = element_text(face = "bold"))
save_fig(pS3, "figS3_pit_calibration.png", 7.2, 4.4)

# ---- Supplementary Figure S4: temporal forecasts against observations, with the naive forecast ----
LAB_S4 <- c(NY = "NYC (citywide, weekly)", UT = "Utah (statewide, daily)")
obs_ts <- rbind(wk[setting == "NY", .(setting, wave, day = day + 3, obs = cases_per_week)],   # NYC: last day of each reporting week
             res("utah_daily_cases.csv")[, .(setting, wave, day, obs = cases)])
fc <- res("temporal_cv_daily.csv")[, day := time_id]
fc[, naive := obs_ts[setting == .BY$setting & wave == .BY$wave][match(origin, day), obs], by = .(setting, wave)]
obs_ts <- merge(obs_ts, fc[, .(start = min(origin)), by = .(setting, wave)], by = c("setting", "wave"))[day >= start]
segs <- fc[, .(x0 = min(day) - ifelse(setting == "NY", 2, 0.4), x1 = max(day) + ifelse(setting == "NY", 2, 0.4), y = naive[1]),
           by = .(setting, wave, origin)]
for (d in list(obs_ts, fc, segs)) d[, `:=`(set_l = factor(LAB_S4[setting], LAB_S4), wave_l = wave_f(wave))]
pS4 <- ggplot() + geom_line(data = obs_ts, aes(day, pmax(obs, 1)), colour = "grey30", linewidth = 0.5) +
  geom_point(data = obs_ts, aes(day, pmax(obs, 1)), colour = "grey30", size = 0.7) +
  geom_ribbon(data = fc[setting == "UT"], aes(day, ymin = pmax(lo_draw, 1), ymax = hi_draw, fill = wave_l, group = origin), alpha = 0.3) +
  geom_linerange(data = fc[setting == "NY"], aes(day, ymin = pmax(lo_draw, 1), ymax = hi_draw, colour = wave_l), linewidth = 2.4, alpha = 0.35) +
  geom_point(data = fc, aes(day, pmax(pred_med, 1), colour = wave_l), size = 1.8) +
  geom_segment(data = segs, aes(x = x0, xend = x1, y = y, yend = y, linetype = "Naive forecast"), colour = "black", linewidth = 0.6) +
  facet_grid(set_l ~ wave_l, scales = "free_y") +
  scale_colour_manual(values = WAVE_COLS, guide = "none") + scale_fill_manual(values = WAVE_COLS, guide = "none") +
  scale_linetype_manual(values = c("Naive forecast" = "22"), name = NULL) +
  scale_y_log10(labels = label_comma(), breaks = breaks_log(n = 5), expand = expansion(mult = c(0.05, 0.18))) +
  scale_x_continuous(breaks = seq(45, 90, 15), limits = c(44, 92)) + labs(x = "Day of wave", y = "Cases (log scale)") +
  theme_panels(legend.position = "bottom", legend.text = element_text(size = 13))
save_fig(pS4, "figS4_temporal_cv.png", 15, 10, dpi = 400)

# Naive benchmark: last observed MODZCTA week (NYC) or county day carried forward (Utah, three days)
win <- res("temporal_cv_windows.csv", colClasses = list(character = "unit"))
oc <- res("observed_counts.csv.gz", colClasses = list(character = "unit"))
win[oc, naive := fifelse(setting == "NY", i.y, 3 * fcoalesce(as.numeric(i.y), 0)), on = .(setting, wave, unit, origin = last_day)]
bench <- win[, .(mae_model = mean(abs(obs - pred_med)), mae_naive = mean(abs(obs - naive)), ratio_model = sum(pred_med) / sum(obs),
                 ratio_naive = sum(naive) / sum(obs), n = .N), by = .(setting, wave, origin)][, rel_mae := mae_model / mae_naive]
save_tab(bench, "forecast_benchmark_by_origin.csv")
save_tab(win[, .(rel_mae = mean(abs(obs - pred_med)) / mean(abs(obs - naive))), by = .(setting, wave)], "forecast_benchmark_rel_mae.csv")

# ---- Supplementary Figure S5: held-out area totals in spatial cross-validation ----
ar <- res("spatial_cv_areas.csv", colClasses = list(character = "geo_id"))[obs > 0 & is.finite(pred_med)]
ar[is.na(covered_draw), covered_draw := obs >= lo_draw & obs <= hi_draw]
ar[, wave_l := wave_f(wave)]
row_cv <- function(st, top) {
  d <- ar[setting == st][, row := c(NY = "NYC (MODZCTAs)", UT = "Utah (counties)")[[st]]]
  l <- range(c(d$obs, pmax(d$pred_med, 1), pmax(d$lo_draw, 1), d$hi_draw))
  ggplot(d, aes(obs, pmax(pred_med, 1), colour = wave_l)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.5) +
    geom_linerange(aes(ymin = pmax(lo_draw, 1), ymax = hi_draw), alpha = 0.3, linewidth = 0.5) +
    geom_point(aes(shape = covered_draw), size = 2, stroke = 0.8) + facet_grid(row ~ wave_l) +
    scale_x_log10(labels = label_comma(), limits = l) + scale_y_log10(labels = label_comma(), limits = l) +
    scale_colour_manual(values = WAVE_COLS, guide = "none") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Inside 95% interval", `FALSE` = "Outside"), name = NULL) +
    labs(x = "Observed held-out total (log scale)", y = "Predicted total, posterior median (log scale)") +
    theme_panels(strip.text.x = if (top) element_text(face = "bold", size = 15) else element_blank(),
                 strip.background.x = if (top) element_rect() else element_blank(), legend.position = "bottom", legend.text = element_text(size = 13))
}
save_fig((row_cv("NY", TRUE) / row_cv("UT", FALSE)) + plot_layout(guides = "collect", axis_titles = "collect") & theme(legend.position = "bottom"),
         "figS5_spatial_cv.png", 15, 10.5, dpi = 400)

# ---- Tables ----
fmt_ci <- function(m, lo, hi, d = 1) {
  z <- function(x) ifelse(round(x, d) == 0, 0, x)   # no "-0.0"
  sprintf(paste0("%.", d, "f (%.", d, "f, %.", d, "f)"), z(m), z(lo), z(hi))
}
sw <- function(d) d[, setting := SETTING[setting]][]

mc <- res("model_comparison.csv")[, dWAIC := WAIC - min(WAIC), by = .(setting, wave)][order(setting, wave, dWAIC)]
save_tab(sw(mc[, .(setting, wave, model = paste0(MODEL[model], ifelse(dWAIC == 0, "*", "")), dWAIC = round(dWAIC, 2), RMSE = round(RMSE, 3))]),
         "table2_model_comparison.csv")

fe <- res("fixed_effects.csv")[covariate != "(Intercept)"]
fe[, `:=`(covariate = factor(COVARIATE[covariate], unique(COVARIATE)), RR = fmt_ci(RR, RR_lower, RR_upper, 2))]
nl <- selected[startsWith(model, "NLMob"), .(setting, wave, covariate = factor(TERM[["log_rho_resid"]], unique(COVARIATE)), RR = "nonlinear")]
fe <- rbind(fe[, .(setting, wave, covariate, RR)], nl)
for (s in c("NY", "UT"))
  save_tab(dcast(fe[setting == s], covariate ~ paste("Wave", wave), value.var = "RR", fill = "--"),
           sprintf("table%d_fixed_effects_%s.csv", if (s == "NY") 3 else 4, if (s == "NY") "nyc" else "utah"))

re <- merge(sp[, .(units = .N, sd_b = sd(mean), rr_b_min = exp(min(mean)), rr_b_max = exp(max(mean)), var_b = var(mean)), by = .(setting, wave)],
            te[, .(sd_g = sd(mean), rr_g_min = exp(min(mean)), rr_g_max = exp(max(mean)), var_g = var(mean)), by = .(setting, wave)])
save_tab(sw(re[, .(setting, wave, units, sd_b = round(sd_b, 3), spatial_rr = sprintf("%.2f-%.2f", rr_b_min, rr_b_max), sd_g = round(sd_g, 2),
                   temporal_rr = sprintf("%.2f-%.2f", rr_g_min, rr_g_max), spatial_share_pct = round(100 * var_b / (var_b + var_g)))]),
         "table5_random_effects.csv")

pt <- res("pcatm_terms.csv")
tab6 <- merge(merge(res("pcatm_summary.csv"), pt[term == "ext_prev_std", .(setting, wave, imp = pcatm, imp_lo = pcatm_q025, imp_hi = pcatm_q975)]),
              res("pcatm_holdout.csv")[, .(setting, wave, ho = pcatm, q05, q95)], by = c("setting", "wave"))
save_tab(sw(tab6[, .(setting, wave, in_sample_95 = fmt_ci(100 * pcatm, 100 * pcatm_q025, 100 * pcatm_q975),
                     imported_exposure_only_95 = fmt_ci(100 * imp, 100 * imp_lo, 100 * imp_hi),
                     spatial_holdout_90 = fmt_ci(100 * ho, 100 * q05, 100 * q95))]), "table6_pcatm.csv")

nb <- hp[grepl("^size for the nbinomial", parameter), .(setting, wave, nb_size = `0.5quant`)]
fam <- merge(merge(selected, res("model_comparison.csv")[, .(setting, wave, model, Family)], by = c("setting", "wave", "model")), nb, by = c("setting", "wave"), all.x = TRUE)
save_tab(sw(fam[, .(setting, wave, family = ifelse(Family == "nbinomial", "NegBin", "Poisson"), dispersion_ratio = round(dispersion, 2),
                    nb_size = ifelse(nb_size >= 10, round(nb_size), round(nb_size, 1)))]), "tableS6_family.csv")

nc <- res("negative_control.csv")
save_tab(sw(nc[, .(setting, wave, rr_95 = sprintf("%.3f [%.3f, %.3f]", exp(mean), exp(q025), exp(q975)), spans_1 = q025 < 0 & q975 > 0)]),
         "tableS7_negative_control.csv")

s8 <- merge(pt, res("pcatm_summary.csv")[, .(setting, wave, total_amlp = amlp)], by = c("setting", "wave"))
s8 <- s8[order(setting, wave, match(term, c("own_prev_std", "ext_prev_std", "mob_spatial_lag_std", "mob_weighted", "log_rho_resid")))]
save_tab(sw(rbind(s8[, .(setting, wave, covariate = TERM[term], pcatm_95 = fmt_ci(100 * pcatm, 100 * pcatm_q025, 100 * pcatm_q975),
                         share_of_amlp = round(100 * amlp / total_amlp, 1))],
                  tab6[, .(setting, wave, covariate = "Total", pcatm_95 = fmt_ci(100 * pcatm, 100 * pcatm_q025, 100 * pcatm_q975), share_of_amlp = 100)])[
                    order(setting, wave, covariate == "Total")]), "tableS8_pcatm_terms.csv")

tcv <- res("temporal_cv_scores.csv")
save_tab(tcv[setting == "NY", .(wave, model = MODEL[model], cell_MdAPE = round(100 * MdAPE, 1), modzcta_total_MdAPE = round(modzcta_total_MdAPE, 1),
                                corr = round(modzcta_total_Corr, 3), cov95 = round(100 * modzcta_total_Cov95, 1))], "tableS9_temporal_cv_nyc.csv")
save_tab(tcv[setting == "UT", .(wave, model = MODEL[model], MdAPE = round(100 * MdAPE, 1), CRPS = round(CRPS, 2), cov95 = round(100 * Cov95, 1))],
         "tableS10_temporal_cv_utah.csv")
save_tab(sw(res("spatial_cv_metrics.csv")[, .(setting, wave, areas = n_areas, MdAPE = round(MdAPE, 1), cov95 = round(100 * Cov95))]),
         "tableS11_spatial_cv.csv")
