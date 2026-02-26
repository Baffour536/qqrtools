# qqr: Quantile-on-Quantile Regression in R

<!-- badges: start -->
<!-- badges: end -->

## Overview

`qqr` implements **quantile-on-quantile regression (QQR)** as introduced
by Sim and Zhou (2015). Where standard quantile regression (QR) models the
conditional quantiles of a response as a function of predictors, QQR goes a
step further: it models how the effect of a regressor on the outcome varies
jointly across the *distribution of the outcome* **and** the *distribution
of the regressor*. The result is a two-dimensional coefficient surface
indexed by both quantile dimensions.

This richer structure is particularly useful in applications where the
relationship between two variables is asymmetric, non-linear, or
concentrated in specific regions of the joint distribution — for example,
tail dependence in financial markets, regime-dependent transmission in
macroeconomics, or heterogeneous treatment effects in policy evaluation.

Recent applied uses of QQR include studies of oil price spillovers
(Sim and Zhou, 2015), COVID-19 and climate interactions
(Shahzad et al., 2020), and carbon-green bond market linkages (Wu, 2023).

---

## Key features

- **Full QQR coefficient surfaces** — fits local linear quantile regressions
  over a grid of regressor values and a vector of outcome quantile levels,
  storing coefficient and p-value surfaces for every coefficient in the model.

- **Automatic bandwidth selection** — eight principled, data-driven methods
  remove researcher discretion over the most critical tuning parameter and
  promote reproducibility:

## Choosing a bandwidth method

| Situation | Recommended method | Rationale |
|-----------|-------------------|-----------|
| Normal or near-normal regressor | `"silverman"` | AMISE-optimal under normality; fast |
| Financial returns, exchange rates, commodity prices | `"robust"` | IQR correction guards against heavy tails |
| Distribution better described by quartiles | `"iqr"` | Scale matched to interquartile spread |
| Fast, assumption-free selection (no extra packages) | `"cv_fast"` | LOO quantile loss; no distributional assumption |
| Best density-based selector, base R only | `"sj"` | Sheather-Jones; best available density selector |
| Regression data, best structural match to QQR | `"yj"` | Uses both X and Y; regression MISE-optimal |
| Maximum data-driven rigour, mixed-type data | `"npcv"` | np full LS cross-validation |
| Varying-coefficient structural analogue | `"sc"` | np smooth-coefficient CV; mirrors QQR structure |
| Replicating a specific prior result | Numeric `h` with documentation | Manual, triggers warning |

The `"yj"` selector (Ruppert-Sheather-Wand) is particularly noteworthy
because it is the only method that uses both the quantile regressor `x`
and the outcome `y` together to minimise regression MISE -- making it the
most theoretically appropriate selector for QQR, which is fundamentally a
regression problem.


  Manual numeric `h` is accepted but triggers a reproducibility warning.

- **Publication-quality 3D plots** — two complementary surface visualisations
  using only base R graphics (no extra dependencies):
  - `surface3d()` — a smooth coloured 3D surface, equivalent to MATLAB's
    `surf()`, with colour legend and optional zero-reference plane.
  - `panel_surface()` — a ribbon panel plot with one coloured flat panel per
    coefficient, matching the style used in several published QQR studies.
  - `plot.qqr()` and `Int.plot()` — static and interactive (via `rgl`) plots.

- **Word export with significance stars** — `word.qqr(which = "coef_stars")`
  produces a single combined table with coefficients and stars appended
  directly (`***`, `**`, `*`), plus an automatic footnote. Unavailable
  p-values are marked `[NA]` rather than left blank.

- **QR vs. QQR comparison** — `compare.qqr()`, `plot.compare.qqr()`,
  `excel.compare.qqr()`, and `word.compare.qqr()` collapse the QQR surface
  to a univariate function of `tau` and compare it with standard QR
  coefficients — a robustness check recommended in the applied literature.

- **Tidy surface extraction and export** — `qqr_surface()`, `rowave.qqr()`,
  `excel.qqr()`, and `word.qqr()` for downstream reporting and analysis.

---

## Installation


---

## Quick start

```r
library(qqr)

set.seed(42)
n   <- 300
x   <- rnorm(n)
w   <- rnorm(n)
y   <- 0.5 + x * qnorm(runif(n)) + 0.3 * w
dat <- data.frame(y = y, x = x, w = w)

## --- 1. Fit with automatic bandwidth (Silverman -- default) ---
fit <- qqr(y ~ x + w, data = dat, qvar = "x",
           taus    = seq(0.05, 0.95, by = 0.05),
           m       = 19,
           h       = "silverman",   # default; also try "robust", "sj", "yj"
           se_type = "nid")

## --- 1b. Other automatic selectors (same interface) ---
# fit_sj <- qqr(y ~ x + w, data = dat, qvar = "x", h = "sj")
# fit_yj <- qqr(y ~ x + w, data = dat, qvar = "x", h = "yj")  # needs KernSmooth
# fit_cv <- qqr(y ~ x + w, data = dat, qvar = "x", h = "cv_fast")

print(fit)
summary(fit)

## --- 2. 3D surface plot ---
surface3d(fit, coef = "x", palette = "viridis")

## --- 3. Panel ribbon plot ---
panel_surface(fit)

## --- 4. Export to Word with significance stars ---
word.qqr(fit, coef = "x", which = "coef_stars", digits = 3,
         file = "QQR_results_stars.docx")

## --- 5. Compare with standard QR ---
cmp <- compare.qqr(fit, coef = "x")
plot.compare.qqr(fit, coef = "x")

## --- 6. Export all surfaces to Excel ---
excel.qqr(fit, which = "both", digits = 4,
          file = "QQR_surfaces.xlsx")
```

---

## Bandwidth selection — why it matters

The bandwidth `h` is the single most influential tuning parameter in QQR. A
bandwidth that is too small produces a noisy, erratic surface; one that is
too large smooths over genuine local heterogeneity. Because researchers can,
intentionally or not, choose `h` to obtain more favourable results, free
bandwidth selection is a form of specification searching that undermines the
credibility of QQR findings.

`qqr` makes automatic selection the default and stores the method and full
formula in the fitted object so it can be reported transparently:

```r
fit <- qqr(y ~ x, data = dat, qvar = "x", h = "robust")

## The summary prints the selection details for the methods section:
summary(fit)
#> Selection method    :  Robust rule of thumb
#> Details             :  h = 0.9 * min(sd(x), IQR(x)/1.34) * n^(-1/5)
#>                        sd(x) = 0.9874,  IQR(x) = 1.3021,  n = 300,  h = 0.2814
```

Manual `h` is still allowed for reproducibility of existing work, but
triggers a warning:

```r
fit2 <- qqr(y ~ x, data = dat, qvar = "x", h = 0.05)
#> Warning: Bandwidth h was set manually to 0.05. Manual bandwidth selection
#> may affect results and reduce reproducibility. For automatic selection,
#> use h = 'silverman' (default), 'robust', 'iqr', 'cv_fast', 'sj', 'yj', 'npcv', or 'sc'.
```


---

## Main functions

| Function | Description |
|----------|-------------|
| `qqr()` | Fit a QQR model |
| `print.qqr()` | Compact model summary |
| `summary.qqr()` | Detailed surface statistics and bandwidth diagnostics |
| `surface3d()` | Smooth 3D coefficient surface (base R) |
| `panel_surface()` | Ribbon panel 3D plot (base R) |
| `plot.qqr()` | Static 3D plot via `plot3D` or base R |
| `Int.plot()` | Interactive 3D plot via `rgl` |
| `qqr_surface()` | Extract surfaces as wide or long data frames |
| `rowave.qqr()` | Row-wise averages across outcome quantiles |
| `compare.qqr()` | Compare row-mean QQR with standard QR |
| `plot.compare.qqr()` | `ggplot2` comparison plot |
| `excel.qqr()` | Export surfaces to Excel |
| `word.qqr()` | Export surfaces to Word (with optional stars) |
| `excel.compare.qqr()` | Export QR vs. QQR comparison to Excel |
| `word.compare.qqr()` | Export QR vs. QQR comparison to Word |

---

## References

Fan, J. and Gijbels, I. (1996) *Local Polynomial Modelling and Its
Applications*. London: Chapman & Hall.

Koenker, R. (2005) *Quantile Regression*. Cambridge: Cambridge University
Press. doi:10.1017/CBO9780511754098

Koenker, R. and Bassett, G. (1978) Regression quantiles. *Econometrica*,
46(1), 33--50. doi:10.2307/1913643

Shahzad, F., Shahzad, U., Fareed, Z., Iqbal, N., Hashmi, S.H. and Ahmad, F.
(2020) Asymmetric nexus between temperature and COVID-19 in the top ten
affected provinces of China: A current application of quantile-on-quantile
approach. *Science of the Total Environment*, 736, 139115.
doi:10.1016/j.scitotenv.2020.139115

Silverman, B.W. (1986) *Density Estimation for Statistics and Data Analysis*.
London: Chapman & Hall.

Sim, N. and Zhou, H. (2015) Oil prices, US stock return and the dependence
between their quantiles. *Journal of Banking & Finance*, 55, 1--8.
doi:10.1016/j.jbankfin.2015.02.030

Wu, D. (2023) The linkage between carbon market and green bond market:
Evidence from a quantile-on-quantile regression approach. *Sustainability*,
15(13), 10634. doi:10.3390/su151310634

Hayfield, T. and Racine, J.S. (2008) Nonparametric econometrics: The np
package. *Journal of Statistical Software*, 27(5), 1--32.

Jones, M.C., Marron, J.S. and Sheather, S.J. (1996) A brief survey of
bandwidth selection for density estimation. *Journal of the American
Statistical Association*, 91(433), 401--407.

Li, Q. and Racine, J.S. (2010) Smooth varying-coefficient estimation and
inference for qualitative and quantitative data. *Econometric Theory*,
26(6), 1--31.

Ruppert, D., Sheather, S.J. and Wand, M.P. (1995) An effective bandwidth
selector for local least squares regression. *Journal of the American
Statistical Association*, 90(432), 1257--1270.

Scott, D.W. (1992) *Multivariate Density Estimation: Theory, Practice, and
Visualization*. New York: Wiley.

Sheather, S.J. and Jones, M.C. (1991) A reliable data-based bandwidth
selection method for kernel density estimation. *Journal of the Royal
Statistical Society Series B*, 53(3), 683--690.

Yu, K. and Jones, M.C. (1998) Local linear quantile regression.
*Journal of the American Statistical Association*, 93(441), 228--237.
doi:10.1080/01621459.1998.10473680

---

## Authors

Collins Baffour Kyei, William Godfred Cantah, Peterson Owusu Junior  
Department of Economic Studies and Department of Data Science and Economic
Policy, University of Cape Coast, Ghana.

Corresponding author: [collinskb536@gmail.com](mailto:collinskb536@gmail.com)
