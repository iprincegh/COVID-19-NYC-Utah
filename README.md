# Bayesian Spatiotemporal Attribution of COVID-19 Risk to Human Mobility Across Spatial Resolutions

Code and results for the paper by Prince Oppong Boakye, Poshan Niraula and Edzer Pebesma (Institute for Geoinformatics, University of Münster).

The analysis models COVID-19 cases in three waves (March–May 2020, December 2020–February 2021 and December 2021–February 2022) at two resolutions. New York City (NYC) uses 100 m grid cells, with weekly cases for Modified ZIP Code Tabulation Areas (MODZCTAs) allocated to cells by population. Utah uses daily county counts. Bayesian spatiotemporal models fitted with INLA combine a BYM2 spatial random effect, a random-walk temporal random effect, a Type IV space-time interaction and covariates built from point-of-interest (POI) visit data. The Proportion of Cases Attributable to Mobility (PCAtM) is estimated in-sample and under spatial holdouts.

## Reproducing the figures and tables

`results/` holds the model outputs behind every figure and table. To redraw them, run from the repository root:

```sh
Rscript make_figures.R
```

This writes `figures/` and `tables/` (CSV) in under a minute. The licensed visit data are not included, so Figure 2 shows reported cases only. If the visit stages have been run (see below), Figure 2 also shows POI visits.

To draw the journal figures in the browser without installing anything, open the [`binder` branch](https://github.com/iprincegh/COVID-19-NYC-Utah/tree/binder) in Binder and run `source("make_figures.R")` in the R console:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/iprincegh/COVID-19-NYC-Utah/binder?urlpath=rstudio)

| Output | Content |
|---|---|
| `fig01`–`fig09` | main-text Figures 1–9 |
| `figS1`–`figS5` | Supplementary Figures S1–S5 |
| `table2`–`table6` | main-text Tables 2–6 |
| `tableS6`–`tableS11` | Supplementary Tables S6–S11 |
| `supplementary_data_1_amlp_by_category.csv` | Supplementary Data 1 |
| `pit_calibration_tests.csv`, `forecast_*.csv` | calibration tests and the naive-forecast benchmark |

Names in `results/` follow the paper: models carry their Table 2 names, and the mobility-derived covariates are named after their symbols.

| Name | Covariate |
|---|---|
| `rho_inf` | Local infection pressure, $\rho^{\mathrm{inf}}_{g,t}$ |
| `rho_imp` | Imported exposure, $\rho^{\mathrm{imp}}_{g,t}$ |
| `z_resid` | Same-day residual mobility, $z^{\mathrm{resid}}_{g,t}$ |
| `L` | Lagged residual mobility, $L_{g,t}$ |
| `Wz_bar` | Spatial mobility, $[\mathbf{W}\bar{\mathbf{z}}]_g$ |

## Running the full analysis

The analysis is a sequence of scripts in `analysis/`, each run from the repository root. Each reads the outputs of the previous ones in `outputs/`. Stages that loop over waves take an optional setting and wave, for example `Rscript analysis/10_fit_models.R NY wave1`.

| Stage | Script | Output |
|---|---|---|
| 01 | `01_spatial_units.R` | NYC 100 m grid, MODZCTAs and Utah counties |
| 02 | `02_population.R` | WorldPop cell populations, ACS county totals, MODZCTA denominators |
| 03 | `03_cases.R` | NYC weekly MODZCTA rates and daily series, Utah daily county counts |
| 04 | `04_visits.R` | POI visits per cell and day, visitor home-CBG flows, co-visitation (needs Advan data) |
| 05 | `05_daily_visits.R` | daily visit totals for Figure 2 |
| 06 | `06_infection_pressure.R` | infection pressure and visit activity per unit-day, trip-length kernels (Supplementary Table S2) |
| 07 | `07_spatial_covariates.R` | cell-level visit covariates |
| 08 | `08_demographics.R` | ACS shares, building density, land use |
| 09 | `09_model_data.R` | unit-day model data, NYC weekly observations, the neighbour graph, and local infection pressure and imported exposure pooled to nearby grid cells |
| 10 | `10_fit_models.R` | six model variants per wave and the WAIC-selected model |
| 11 | `11_temporal_cv.R` | rolling-origin forecasts |
| 12 | `12_spatial_cv.R` | spatial holdouts and PCAtM under spatial holdouts |
| 13 | `13_pcatm.R` | in-sample PCAtM, per covariate and by POI category |
| 14 | `14_negative_control.R` | negative-control refit |
| 15 | `15_diagnostics.R` | share of visits from within the destination's reporting area, and the within-MODZCTA variation reported in the Results |
| 16 | `16_export_results.R` | `results/` for `make_figures.R` |

`R/` contains the shared functions. `setup.R` holds the configuration and file paths, `data.R` input handling and population offsets, `covariates.R` covariate construction, `adjacency.R` the neighbour graph, `model.R` the INLA models and `validation.R` cross-validation.

### Computing requirements

Utah runs on a laptop in minutes. NYC has about 54,000 cells and five million cell-days, so its stages need a large-memory machine. Stages 09–14 used 24 cores and 250–700 GB of memory per wave, each finishing within 24 hours. The number of INLA threads is set with the environment variable `INLA_THREADS` (for example `4:6`). The reported results used INLA 26.8.7. R 4.4 or later is required for the native `%||%`. The required packages are:

`data.table`, `R.utils`, `sf`, `spdep`, `INLA`, `arrow`, `Matrix`, `igraph`, `RANN`, `rgeoda`, `lubridate`, `terra`, `exactextractr`, `DBI`, `duckdb`, `jsonlite`; for `make_figures.R` also `ggplot2`, `patchwork`, `scales`, `scico` and `surveillance`.

### Randomness

The case allocation to cells (stage 09) and the k-means spatial folds (stage 12) use `set.seed(42)`. Cross-validation draws use INLA seed 42 and are reproducible. In-sample PCAtM (stage 13) sets `set.seed(20261008)`. In the NYC grid fits, only the modal configuration stored with the fit is valid, so `inla.posterior.sample()` aborts and stage 13 draws the coefficients from that configuration. These draws are reproducible, and their intervals exclude hyperparameter uncertainty. Utah uses `inla.posterior.sample()`, whose draws depend on INLA's own random state, so its PCAtM intervals vary slightly between runs.

## Data

Input data go in `data/`, at the paths set in `R/setup.R`. All sources except the visit data are public:

| Data | Source | Path |
|---|---|---|
| NYC weekly case rates by MODZCTA, MODZCTA population | NYC DOHMH [coronavirus-data](https://github.com/nychealth/coronavirus-data): `caserate-by-modzcta.csv`, `data-by-modzcta.csv` | `data/us_CovidCases/caseRate/NYCcaserate-by-modzcta.csv`, `data/MODZCTA/data-by-modzcta.csv` |
| MODZCTA boundaries | NYC Open Data, Modified Zip Code Tabulation Areas | `data/MODZCTA/` |
| Block group to MODZCTA crosswalk | NYC DOHMH [coronavirus-data](https://github.com/nychealth/coronavirus-data) (`cbg_geoid`, `modzcta`, `borough`) | `data/us_CovidCases/caseRate/NYCcbg_to_modzcta.csv` |
| Utah county cases | [New York Times COVID-19 data](https://github.com/nytimes/covid-19-data) (`us-counties-2020/21/22.csv`) | `data/us_CovidCases/` |
| Population | WorldPop constrained 100 m, R2025A, 2020–2022 | `data/population/` |
| Age and employment | US Census ACS 5-year Tables S0101 and DP03, 2020–2022 | `data/demographics/` |
| Block groups | US Census TIGER/Line 2024 | `data/cbg_shapefiles/` |
| Utah counties | Utah Geospatial Resource Center | `data/ut_counties/` |
| Land use, residential units | NYC MapPLUTO 25v4; Utah parcels (LIR), UGRC | `data/landuse/` |
| Building footprints | Microsoft US Building Footprints | `data/building_footprints/` |
| Bridges, subway, roads (NYC neighbour graph) | NYC DOT Bridge Ratings (NYC Open Data); MTA Subway Stations (Open NY); 2016 New York City Subway Routes (NYU Spatial Data Repository); US Census TIGER/Line 2023 primary and secondary roads | `data/transport/` |
| POI visits | Advan Weekly Patterns Plus (commercial licence) | `data/advan_weekly_plus/` |

The Advan data cannot be redistributed and are available from [Advan](https://advanresearch.com) under a commercial licence. `results/` contains no Advan records: it holds model outputs, case counts and boundaries.

## Licence

Code: MIT licence (`LICENSE`). Each data source keeps its own licence.

## Contact

Prince Oppong Boakye, prince.ob@uni-muenster.de
