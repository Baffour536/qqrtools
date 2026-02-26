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
