# qqrtools 0.3.0

First release under the name **qqrtools** (previously **qqr**). This
version consolidates the local development of the visualisation suite with
a full review of the estimator, and corrects several defects. Results
obtained with qqr 0.2.0 should be re-estimated.

## Migration from qqr

* Replace `library(qqr)` with `library(qqrtools)`. Function names, the
  `"qqr"` class and argument names are unchanged.
* To reproduce the qqr 0.2.0 grid, use `qqr(..., grid = "range")`.
* `compare.qqr()` results change (see below); recompute them.

## Breaking changes (methodological corrections)

* **Grid at quantiles of the regressor.** The documentation, plots and Sim
  and Zhou (2015) define the QQR surface at quantiles of the regressor,
  but qqr 0.2.0 used `m` equally spaced points over its range and labelled
  them 0.05-0.95 in the 3D plots. `qqr()` now places the grid at the
  empirical quantiles `quantile(x, thetas)` (new arguments `grid =
  "quantile"` (default) and `thetas`), stores `thetas`, and all plots use
  them. `grid = "range"` keeps the old grid, with correctly computed
  `thetas` (the empirical CDF at the grid points).
* **`compare.qqr()` now performs the Sim-Zhou comparison.** It previously
  averaged each surface over tau and then paired grid point j with tau_j,
  which is not a meaningful comparison. It now averages the surface over
  theta for each tau and compares it with `quantreg::rq()` on the same
  design. The output gains a `difference` column.
* **`h = "cv_fast"` cross-validates y, not x.** The old selector predicted
  each `x_i` by a kernel-weighted median of the other `x` values, ignoring
  the outcome, and therefore almost always returned the smallest candidate.
  It now minimises the leave-one-out check loss of a kernel-weighted
  conditional median of `y` given `x`.
* `effective_n` now holds the Kish effective sample size
  `(sum w)^2 / sum w^2` instead of the raw sum of kernel weights.

## New features

* `kernel_space = "cdf"` applies the kernel to `F_n(X_t) - theta`, the
  exact weighting of Sim and Zhou (2015); the bandwidth is then on the
  probability scale.
* Factor controls (expanded with treatment contrasts), transformed
  regressors such as `qvar = "log(x)"`, transformed responses and models
  without an intercept are supported. Previously, factor controls and
  transformed regressors silently produced all-`NA` surfaces.
* The fitted object stores `formula`, `thetas`, `grid_type`,
  `kernel_space`, `n_failed`, `y` and `design`.
* `qqr_cor_heatmap()` returns the independence benchmark for
  co-exceedance (`$independence`).
* Visualisation suite (from the development version): `qqr_heatmap()`,
  `qqr_contour()`, `qqr_plotly()`, `qqr_cor_heatmap()`.

## Bug fixes

* Input validation: `taus`/`thetas` must be strictly increasing values in
  (0, 1); `m` a whole number >= 2; `h` a positive number or a known rule;
  informative errors for empty data, non-numeric or constant regressors,
  infinite values, interactions with the quantile regressor and too few
  observations.
* Local fits that cannot be estimated (too few effective observations,
  rank-deficient local design, solver errors) are counted in `n_failed`
  and reported with one warning, instead of failing silently.
* Somers' D(Y|X) used the pairs tied on X in its denominator (that is
  D(X|Y)); it now uses pairs tied on Y only, as defined by Somers (1962).
* `qqr_cor_heatmap()`: the first bin break is now the minimum (it was the
  `tau_1 / 2` quantile); co-exceedance is computed at the tau-quantiles
  over the full sample, consistent with its independence benchmark;
  constant cells return `NA`; the independence line now marks where the
  observed co-exceedance equals the benchmark.
* plotly colour scales: diverging scales were reversed (negative values
  shown red). Scales are now built explicitly (low = blue for "RdBu");
  every `grDevices::hcl.pals()` name, e.g. "plasma", is supported (plotly
  silently ignored unknown names before).
* `qqr_heatmap()` accepts `coef = "intercept"`; all plotting functions
  validate `clim` and give clear errors for missing p-value surfaces.
* `qqr_surface(tau_subset = )` matched tau values with exact floating-point
  equality, so values such as 0.1 from `seq(0.05, 0.95, by = 0.05)` were
  silently dropped; matching now uses a tolerance and unmatched values
  give an error.
* `word.qqr(which = "coef_stars")` marks missing p-values with `[NA]` (as
  documented) and formats coefficients with fixed decimals.
* `print.qqr()`/`print.summary.qqr()` label all eight bandwidth methods.
* A duplicate internal `.qqr_match_coef()` definition was removed.
* Corrected references: the DOI of Sim and Zhou (2015) is
  10.1016/j.jbankfin.2015.01.013; Yu and Jones (1998) is
  10.1080/01621459.1998.10474104; Li and Racine (2004) appeared in
  *Statistica Sinica* 14, 485-512; Li and Racine (2010) is *Econometric
  Theory* 26(6), 1607-1637. An unverifiable reference to "Li, Racine and
  Wooldridge (2002)" was replaced by Li, Huang, Li and Fu (2002), and a
  README reference to "Wu (2023)" with an incorrect title was removed.

## Dependencies and infrastructure

* `dplyr` and `tidyr` are no longer used (long tables are built in base
  R). Imports are now graphics, grDevices, quantreg and stats only.
* All examples run (file-writing parts write to `tempdir()`).
* New test files with reference calculations (every surface cell against a
  direct weighted `quantreg::rq()` fit, large-bandwidth limit equals QR,
  exact recovery of noiseless models, closed-form bandwidths, dependence
  measures against brute-force definitions); line coverage about 94%.
* New "Getting started" vignette that is evaluated when the package is
  built; GitHub Actions `R CMD check` on Windows, macOS and Linux.
* Validation and audit reports in `audit/`.

## Errata for the 0.2.0 notes below

* The 0.2.0 notes describe `h = "cv"`, which never existed under that
  name in 0.2.0 (it was `"cv_fast"`), and refer to "four" selectors; there
  were eight. Removing `@keywords internal` from the package help page does
  not affect exports.

# qqr 0.2.0

## Breaking changes
* The default value of `h` in `qqr()` has changed from `0.05` to
  `"silverman"`. Existing code that relies on `h = 0.05` will receive a
  deprecation warning and should be updated to pass `h = 0.05` explicitly
  or switch to an automatic method.

* `h = "cv"` has been **renamed to `h = "cv_fast"`** to distinguish it
  from the new full np cross-validation selector (`"npcv"`). Code using
  `h = "cv"` must be updated to `h = "cv_fast"`. The selector itself is
  unchanged; only the name has changed.
  
* Internal coefficient names now use the actual variable name for the
  quantile regressor (e.g. `"XOF_USD"`) instead of the generic alias
  `"z"`. Objects fitted with earlier versions of the package should be
  refitted.

## New features

### Bandwidth selection (`qqr()`)

* **Eight automatic bandwidth selectors** have been added via the `h` argument. 
All eight methods are accessed via the `h` argument of `qqr()`. * The fitted object now stores `h_method` (character label) and
`h_details` (full description suitable for a paper's methods section) alongside the numeric `bandwidth`.
Pass a character string to select:

  - `h = "silverman"` (new default): Silverman's rule of thumb,
    `h = 1.06 * sd(x) * n^(-1/5)`. Optimal under normality; fast and
    easy to report.
  - `h = "robust"`: robust rule of thumb,
    `h = 0.9 * min(sd(x), IQR(x)/1.34) * n^(-1/5)`. Preferred for
    heavy-tailed or skewed regressors, as common in financial and
    macroeconomic series.
  - `h = "iqr"`: IQR-based rule, `h = IQR(x) * n^(-1/5)`.
  - `h = "cv"`: leave-one-out cross-validation over a 15-point grid
    centred on the Silverman bandwidth. The most defensible choice for
    publication; computationally heavier than rule-of-thumb methods.
  - **`h = "sj"` -- Sheather-Jones solve-the-equation plug-in.**
    Uses base R `bw.SJ(x)` to solve iteratively for the bandwidth that
    satisfies

      h = [ R(K) / (n * psi4_hat(g(h))) ]^{1/5}

  where `psi4_hat` is a pilot estimate of `integral [f''(x)]^2 dx`
  obtained with a preliminary bandwidth `g(h)`. Widely regarded as the
  best density-based selector available in base R. Requires no
  additional packages. Falls back to the Silverman rule with a warning
  if the iteration fails to converge.

  *Reference: Sheather, S.J. and Jones, M.C. (1991) A reliable
  data-based bandwidth selection method for kernel density estimation.
  Journal of the Royal Statistical Society Series B, 53(3), 683--690.*

  *Additional survey: Jones, M.C., Marron, J.S. and Sheather, S.J.
  (1996) A brief survey of bandwidth selection for density estimation.
  Journal of the American Statistical Association, 91(433), 401--407.*

* **`h = "yj"` -- Ruppert-Sheather-Wand direct plug-in for local
  linear regression.**
  Uses `KernSmooth::dpill(x, y)` to minimise the MISE of a local
  linear regression of the outcome on the quantile regressor:

      h = [ sigma2_eps * integral K^2 /
            (n * f_X(x) * [m''(x)]^2 * mu2^2) ]^{1/5}

  where `m''(x)` is estimated via pilot fits using Mallows Cp block
  sizes. This is the **only selector that uses both `x` and `y`** and
  is tuned to regression (not density) loss, making it the most
  structurally appropriate automatic selector for QQR. Requires
  `KernSmooth`, a recommended package that ships with standard R
  installations.

  *Reference: Ruppert, D., Sheather, S.J. and Wand, M.P. (1995) An
  effective bandwidth selector for local least squares regression.
  Journal of the American Statistical Association, 90(432), 1257--1270.*

* **`h = "npcv"` -- np full least-squares cross-validation.**
  Uses `np::npregbw(xdat = x, ydat = y, regtype = "ll")` to minimise
  the leave-one-out squared-error loss:

      h_CV = argmin_h  (1/n) sum_t [y_t - m_hat_{-t}(x_t; h)]^2

  where `m_hat_{-t}` is the leave-one-out local-linear fit at
  observation `t`. Supports mixed continuous and categorical data via
  generalised product kernels. The most assumption-free selector
  available; also the most computationally intensive. Requires the
  `np` package (listed in `Suggests`).

  *References: Hayfield, T. and Racine, J.S. (2008) Nonparametric
  econometrics: The np package. Journal of Statistical Software,
  27(5), 1--32.*

  *Li, Q. and Racine, J.S. (2004) Cross-validated local linear
  nonparametric regression. Biometrika, 91(1), 639--643.*

* **`h = "sc"` -- np smooth-coefficient cross-validation.**
  Uses `np::npscoefbw(xdat = x, ydat = y)` to select the bandwidth
  for the smooth varying-coefficient model

      Y_i = W_i' gamma(Z_i) + u_i

  This is the closest structural analogue to QQR in least-squares
  form: QQR estimates `beta(tau, z_i)` -- a coefficient that varies
  with `z_i` -- which is structurally identical to a varying-
  coefficient model, except within a quantile rather than mean
  regression framework. Requires the `np` package (listed in
  `Suggests`).

  *References: Li, Q. and Racine, J.S. (2010) Smooth varying-
  coefficient estimation and inference for qualitative and
  quantitative data. Econometric Theory, 26(6), 1--31.*

  *Hayfield, T. and Racine, J.S. (2008) Journal of Statistical
  Software, 27(5), 1--32.*

### Updated `.qqr_select_h()` internal function

* The internal bandwidth helper now accepts an optional `y` argument
  (the outcome vector). It is passed automatically by `qqr()` when
  `h` is one of the regression-based selectors (`"yj"`, `"npcv"`,
  `"sc"`).

* All methods return full bibliographic references inside `h_details`,
  so `cat(summary(fit)$h_details)` produces a paragraph ready to paste
  into a paper's methods section.

* Graceful fallback: `"sj"`, `"yj"`, `"npcv"`, and `"sc"` catch errors
  and fall back to the Silverman rule with an informative warning,
  rather than stopping the entire `qqr()` call.

### Updated manual-bandwidth warning

The reproducibility warning for numeric `h` now lists all eight
automatic methods.

### Printing and summaries

* `print.qqr()` has been redesigned to display bandwidth, selection
  method and a manual-h warning in a formatted block.

* `print.summary.qqr()` now reports a complete bandwidth diagnostics
  section including the formula used, parameter values, and the selected
  `h`. The `summary.qqr()` object gains fields `bandwidth`, `h_method`,
  `h_details`, `kernel`, `se_type`, and `n`.

### Three-dimensional surface plots

* **`surface3d()`** (generic) and **`surface3d.qqr()`** (S3 method): a
  new publication-quality 3D surface plot of the raw QQR coefficient
  matrix using `persp()`. The x-axis spans quantiles of the quantile
  regressor, the y-axis spans outcome quantile levels, and the z-axis
  shows coefficient values. Equivalent to MATLAB's `surf()` applied to
  the raw QQR output. Key arguments: `coef`, `what`, `palette`, `theta`,
  `phi`, `clim`, `add_zero_plane`, `alpha`, `legend`.

* **`panel_surface()`** (generic) and **`panel_surface.qqr()`** (S3
  method): a ribbon panel 3D plot with one coloured flat panel per
  coefficient. Panels are arranged side-by-side along the variable axis
  with quantile levels on the other axis, reproducing the style used in
  several published QQR studies.

### Word export with significance stars

* `word.qqr()` gains a new `which = "coef_stars"` option that produces a
  single combined table with coefficients and significance stars appended
  directly to each cell. Stars follow the convention `*** p < 0.01`,
  `** p < 0.05`, `* p < 0.10`. Cells where p-values could not be
  computed (e.g. insufficient local observations) are marked `[NA]`
  rather than left blank. A formatted footnote is appended automatically.

## Bug fixes

* Fixed two blank lines inside roxygen documentation blocks in
  `qqr-methods.R` that prevented `@export` tags from being parsed
  correctly, causing `compare.qqr`, `plot.compare.qqr`, and
  `Int.plot.qqr` to be missing from the installed package's namespace.

* Removed `@keywords internal` from the `_PACKAGE` sentinel in
  `qqr-package.R`, which was suppressing package-wide exports.

* All S3 methods now carry explicit `@rawNamespace export()` directives
  so they are callable by full name (e.g. `compare.qqr(fit)`) as well as
  through dispatch.

* Replaced `graphics::trans3d` (not formally exported by the `graphics`
  package) with an internal helper `.qqr_trans3d()` to eliminate an
  `R CMD CHECK` warning.

* Removed non-ASCII characters (em-dashes, dagger symbol) from R source
  files to comply with CRAN's portability requirements.

* `dplyr` and `tidyr` moved from `Imports` to `Suggests` since they are
  only used in `qqr_surface()` when `long = TRUE`.

## Internal changes

* `.qqr_select_h()`: new internal function implementing all four
  automatic bandwidth selectors.

* `.qqr_trans3d()`: internal projection helper replacing
  `graphics::trans3d`.

* `qqr-core.R` restructured for clarity: bandwidth resolution is now a
  distinct step (Step 3) between formula parsing and grid construction.

### Updated DESCRIPTION

* `KernSmooth` and `np` added to `Suggests`.
* Version bumped to `0.2.0`.
