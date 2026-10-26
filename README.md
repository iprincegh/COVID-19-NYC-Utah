# Journal figures: Bayesian Spatiotemporal Attribution of COVID-19 Risk to Human Mobility Across Spatial Resolutions

This branch draws the figures of the paper by Prince Oppong Boakye, Poshan Niraula and Edzer Pebesma (Institute for Geoinformatics, University of Münster) from the model outputs in `results/`. The full analysis code is on the [`main` branch](https://github.com/iprincegh/COVID-19-NYC-Utah).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/iprincegh/COVID-19-NYC-Utah/binder?urlpath=rstudio)

Click the badge to open RStudio in the browser, then run in the R console:

```r
source("make_figures.R")
```

The figures are written to `figures/`: `fig01`–`fig09` are main-text Figures 1–9 and `figS1`–`figS5` are Supplementary Figures S1–S5. The script also writes the table values to `tables/`. The first launch builds the environment and can take several minutes. The licensed visit data are not included, so Figure 2 shows reported cases only.

| File | Content |
|---|---|
| `make_figures.R` | draws every figure from `results/` |
| `results/` | model outputs, case counts and boundaries |
| `runtime.txt`, `install.R`, `apt.txt` | R 4.4, the R packages and the system libraries Binder installs |

Names in `results/` follow the paper: models carry their Table 2 names, and the mobility-derived covariates are named after their symbols.

| Name | Covariate |
|---|---|
| `rho_inf` | Local infection pressure, $\rho^{\mathrm{inf}}_{g,t}$ |
| `rho_imp` | Imported exposure, $\rho^{\mathrm{imp}}_{g,t}$ |
| `z_resid` | Same-day residual mobility, $z^{\mathrm{resid}}_{g,t}$ |
| `L` | Lagged residual mobility, $L_{g,t}$ |
| `Wz_bar` | Spatial mobility, $[\mathbf{W}\bar{\mathbf{z}}]_g$ |

## Licence

Code: MIT licence (`LICENSE`).

## Contact

Prince Oppong Boakye, prince.ob@uni-muenster.de
