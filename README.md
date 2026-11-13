# qqrtools: Tools for Quantile-on-Quantile Regression

<!-- badges: start -->
[![R-CMD-check](https://github.com/Baffour536/qqrtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Baffour536/qqrtools/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`qqrtools` estimates **quantile-on-quantile regression (QQR)** surfaces in
the sense of Sim and Zhou (2015). Standard quantile regression (QR) tells you
how a regressor *X* shifts the τ-quantile of an outcome *Y*. QQR asks a finer
question: how does the **θ-quantile of X** affect the **τ-quantile of Y**?
The answer is a two-dimensional surface β(θ, τ) that can reveal
asymmetric, tail-specific or state-dependent relationships which a single
QR slope averages away.

## What the package estimates

For each grid point x<sup>θ</sup> (by default the empirical θ-quantile of
*X*) and each τ, `qqr()` solves the kernel-weighted local linear quantile
regression

min over (α, β, γ) of Σ<sub>t</sub> ρ<sub>τ</sub>( Y<sub>t</sub> − α − β (X<sub>t</sub> − x<sup>θ</sup>) − γ′W<sub>t</sub> ) · K(d<sub>t</sub> / h),

where ρ<sub>τ</sub> is the check function, *K* is the Gaussian kernel and
*W* are optional controls. The local slope β̂(θ, τ) is the QQR estimate.
The kernel distance is either d<sub>t</sub> = X<sub>t</sub> − x<sup>θ</sup>
(`kernel_space = "x"`, default; bandwidth in units of *X*) or
d<sub>t</sub> = F<sub>n</sub>(X<sub>t</sub>) − θ
(`kernel_space = "cdf"`, the exact Sim–Zhou weighting; bandwidth on the
probability scale, Sim and Zhou use h = 0.05). Each local problem is solved
by `quantreg::rq()`.

## Features

- **Estimation** – `qqr()` returns coefficient and p-value surfaces for the
  intercept, the quantile regressor and every control (numeric or factor),
  with the grid, its θ levels, Kish effective local sample sizes and a
  record of failed local fits.
- **Bandwidth rules** – `"silverman"` (default), `"robust"`, `"iqr"`,
  `"cv_fast"` (leave-one-out check-loss cross-validation), `"sj"`
  (Sheather–Jones), `"yj"` (Ruppert–Sheather–Wand, needs *KernSmooth*),
  `"npcv"` and `"sc"` (need *np*), or a number. The chosen rule and a
  sentence for your methods section are stored in the fit.
- **Comparison with QR** – `compare.qqr()` averages the surface over θ for
  each τ and compares it with `quantreg::rq()` (Sim and Zhou's robustness
  check); `plot.compare.qqr()`, `excel.compare.qqr()`, `word.compare.qqr()`.
- **Graphics** – `surface3d()`, `panel_surface()`, `plot()` (base
  graphics), `Int.plot()` (*rgl*), and `qqr_heatmap()`, `qqr_contour()`,
  `qqr_plotly()` (*plotly*) with optional significance overlays.
- **Data diagnostics** – `qqr_cor_heatmap()` summarises the bivariate
  quantile dependence of *X* and *Y* (Spearman, Kendall, Somers' D, Pearson,
  Blomqvist's β, co-exceedance with an independence benchmark).
- **Tables and export** – `qqr_surface()` (wide or long), `rowave.qqr()`,
  `excel.qqr()`, `word.qqr()` (optionally with significance stars).

## Installation

`qqrtools` is not on CRAN. Install it from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Baffour536/qqrtools")
```

R ≥ 4.1 is required (tested on the current, previous and development
versions of R on Windows, macOS and Linux). Only *quantreg* is a hard dependency; plotting and
export functions use *plot3D*, *rgl*, *plotly*, *ggplot2*, *openxlsx*,
*officer* and *flextable* if they are installed.

## Quick start

```r
library(qqrtools)

set.seed(42)
n   <- 400
x   <- rnorm(n)
w   <- rnorm(n)
## the effect of x is larger in the upper tail of y when x is high
y   <- 0.5 + x + 0.3 * w + (1 + 0.5 * pmax(x, 0)) * rnorm(n)
dat <- data.frame(y = y, x = x, w = w)

fit <- qqr(y ~ x + w, data = dat, qvar = "x",
           taus = seq(0.1, 0.9, by = 0.1),   # quantiles of y
           m    = 9)                         # grid at the 0.05..0.95 quantiles of x
fit
summary(fit)

## QQR slope surface: rows = theta (quantiles of x), columns = tau
round(fit$qcoef, 2)

## Robustness check against standard quantile regression
compare.qqr(fit, digits = 3)

## Graphics
surface3d(fit, palette = "RdBu", add_zero_plane = TRUE)
panel_surface(fit, coef = c("x", "w"))
if (requireNamespace("plotly", quietly = TRUE)) {
  qqr_heatmap(fit, sig_overlay = "stars")
}

## Tables and export
head(qqr_surface(fit, coef = "x", digits = 3))
if (requireNamespace("openxlsx", quietly = TRUE)) {
  excel.qqr(fit, file = file.path(tempdir(), "qqr_surfaces.xlsx"), digits = 4)
}
```

### The Sim–Zhou weighting

```r
fit_sz <- qqr(y ~ x + w, data = dat, qvar = "x",
              taus = seq(0.1, 0.9, by = 0.1), m = 9,
              kernel_space = "cdf", h = 0.05)   # warns: manual bandwidth
```

## Interpreting the output

- `fit$qcoef[i, j]` is β̂(θ<sub>i</sub>, τ<sub>j</sub>): the local effect of
  *X* on the τ<sub>j</sub>-quantile of *Y* when *X* is near its
  θ<sub>i</sub>-quantile (`fit$thetas`, `fit$grid`).
- The intercept surface is the fitted τ-quantile of *Y* at
  *X* = x<sup>θ</sup> and *W* = 0; control surfaces are local effects of
  the controls.
- P-values are **pointwise** results from `quantreg::summary.rq()` for each
  weighted local problem. They ignore bandwidth selection, smoothing bias and
  the many cells tested, so treat them as descriptive.
- Results depend on the bandwidth. Automatic rules make the choice
  reproducible, but they come from density or mean-regression theory and
  are not optimal for QQR slopes. Rules of the form h ∝ n<sup>−1/5</sup>,
  such as Silverman's, can give noisy slope surfaces. Report *h* and check
  how sensitive the results are to it (see the methodology vignette).

## Documentation

- `vignette("qqrtools")` – getting started (runs on your machine).
- `vignette("qqr-methodology")` – methodology, bandwidths, inference,
  visualisation and dependence diagnostics.
- Help pages: `?qqr`, `?compare.qqr`, `?surface3d`, `?qqr_heatmap`, …

## Development and testing

The test suite (`tests/testthat`) checks every cell of the surface against
direct weighted `quantreg::rq()` fits, the large-bandwidth limit against
standard QR, exact recovery of noiseless models, bandwidth formulas, the
dependence measures against brute-force definitions, input validation and
all exported functions. Run it with `devtools::test()`. GitHub Actions runs
`R CMD check` on Windows, macOS and Linux. Validation reports live in
[`audit/`](audit/).

## Citation

```r
citation("qqrtools")
```

Please also cite Sim and Zhou (2015) for the QQR method.

## References

Koenker, R. (2005) *Quantile Regression*. Cambridge University Press.
doi:10.1017/CBO9780511754098

Koenker, R. and Bassett, G. (1978) Regression quantiles. *Econometrica*,
46(1), 33–50. doi:10.2307/1913643

Ruppert, D., Sheather, S. J. and Wand, M. P. (1995) An effective bandwidth
selector for local least squares regression. *Journal of the American
Statistical Association*, 90(432), 1257–1270.
doi:10.1080/01621459.1995.10476630

Shahzad, F., Shahzad, U., Fareed, Z., Iqbal, N., Hashmi, S. H. and
Ahmad, F. (2020) Asymmetric nexus between temperature and COVID-19 in the top
ten affected provinces of China: A current application of
quantile-on-quantile approach. *Science of the Total Environment*, 736,
139115. doi:10.1016/j.scitotenv.2020.139115

Sheather, S. J. and Jones, M. C. (1991) A reliable data-based bandwidth
selection method for kernel density estimation. *Journal of the Royal
Statistical Society Series B*, 53(3), 683–690.
doi:10.1111/j.2517-6161.1991.tb01857.x

Silverman, B. W. (1986) *Density Estimation for Statistics and Data
Analysis*. London: Chapman & Hall.

Sim, N. and Zhou, H. (2015) Oil prices, US stock return, and the dependence
between their quantiles. *Journal of Banking & Finance*, 55, 1–8.
doi:10.1016/j.jbankfin.2015.01.013

Yu, K. and Jones, M. C. (1998) Local linear quantile regression. *Journal
of the American Statistical Association*, 93(441), 228–237.
doi:10.1080/01621459.1998.10474104

## Authors

Collins Baffour Kyei (maintainer), William Godfred Cantah and Peterson
Owusu Junior (Department of Economic Studies and Department of Data Science
and Economic Policy, University of Cape Coast, Ghana), and Dominic Owusu
Abeyie. ORCID iDs are listed in `DESCRIPTION`.

## Licence

MIT © 2025–2026 Collins Baffour Kyei et al. See [LICENSE.md](LICENSE.md).
