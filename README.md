# Journal figures: Bayesian Spatiotemporal Attribution of COVID-19 Risk to Human Mobility Across Spatial Resolutions

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/iprincegh/COVID-19-NYC-Utah/binder?urlpath=rstudio) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![R ≥ 4.4](https://img.shields.io/badge/R-%E2%89%A5%204.4-276DC3.svg)](https://www.r-project.org/)

This branch draws the figures of the paper by Prince Oppong Boakye, Poshan Niraula and Edzer Pebesma (Institute for Geoinformatics, University of Münster) from the model outputs in `results/`. The full analysis code is on the [`main` branch](https://github.com/iprincegh/COVID-19-NYC-Utah).

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

## Contact

Prince Oppong Boakye, prince.ob@uni-muenster.de
