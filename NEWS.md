# qqrtools 0.4.0

## Bandwidth and bias

* New default bandwidth rule `h = "slope"`:
  h = 2 * min(sd, IQR/1.34) * n^(-1/7). It uses the rate that is optimal
  for local slopes (Fan and Gijbels, 1996); the constant was calibrated by
  simulation (`tools/validation/bandwidth_calibration.R`). The previous
  default, `"silverman"`, is a density bandwidth that under-smooths slopes
  and remains available.
* New rule `h = "cv_slope"`: the cross-validated bandwidth of `"cv_fast"`
  rescaled to the slope rate.
* New argument `degree` in `qqr()`: `2` (default) fits a local quadratic
  model, which removes the first-order bias of the local linear slope at
  the low and high quantiles of the regressor; `1` gives the local linear
  estimator of Sim and Zhou (2015).
* Results obtained with version 0.3.0 defaults are reproduced with
  `h = "silverman", degree = 1, se_type = "nid"`.

## Inference

* The default standard-error type is now `se_type = "ker"` instead of
  `"nid"`. In null simulations with the default bandwidth and degree, 5%
  tests rejected in 3.5-4.4% of cells with `"ker"`, about 5% with `"boot"`
  and 7-9% with `"nid"` (`tools/validation/size_default.R`). `"boot"` is
  recommended for final results.

## Robustness tools

* `qqr_sensitivity()` refits a model over a range of bandwidths and
  summarises the stability of the surface (with `print()` and `plot()`
  methods).
* `qqr_padjust()` adjusts the p-value surfaces for the number of cells
  tested (`stats::p.adjust()`; Benjamini-Hochberg by default).

# qqrtools 0.3.0

First public release.

## Estimation

* `qqr()` fits quantile-on-quantile regression (Sim and Zhou, 2015) by
  Gaussian-kernel weighted local linear quantile regression, returning
  coefficient and p-value surfaces for the intercept, the quantile
  regressor and every control.
* Grid at the empirical quantiles of the regressor (`grid = "quantile"`,
  levels `thetas`) or equally spaced over its range (`grid = "range"`).
* Kernel weights on the regressor scale (`kernel_space = "x"`) or on the
  empirical distribution function scale of Sim and Zhou
  (`kernel_space = "cdf"`).
* Numeric and factor controls, transformed regressors (e.g.
  `qvar = "log(x)"`), transformed responses and models without an
  intercept.
* Pointwise p-values from `quantreg::summary.rq()` (`se_type = "nid"`,
  `"ker"`, `"boot"` or `"none"`).
* Kish effective local sample sizes (`effective_n`) and a count of local
  fits that could not be estimated (`n_failed`), reported with a warning.

## Bandwidth selection

* Automatic rules via `h`: `"silverman"` (default), `"robust"`, `"iqr"`,
  `"cv_fast"` (leave-one-out check-loss cross-validation), `"sj"`
  (Sheather-Jones), `"yj"` (Ruppert-Sheather-Wand, needs 'KernSmooth'),
  `"npcv"` and `"sc"` (need 'np'), or a single positive number.
* The rule, the bandwidth and a reporting sentence are stored in the fit
  (`h_method`, `bandwidth`, `h_details`).

## Comparison with standard quantile regression

* `compare.qqr()` averages each surface over the regressor quantiles for
  every tau and compares it with `quantreg::rq()`;
  `plot.compare.qqr()`, `excel.compare.qqr()` and `word.compare.qqr()`
  plot and export the comparison.

## Graphics

* Base graphics: `surface3d()`, `panel_surface()`, `plot()`.
* Interactive: `Int.plot()` ('rgl'); `qqr_heatmap()`, `qqr_contour()` and
  `qqr_plotly()` ('plotly') with significance overlays and zero-centred
  diverging colour scales.

## Data diagnostics

* `qqr_cor_heatmap()` summarises the quantile dependence between the
  response and the regressor with Spearman, Kendall, Somers' D, Pearson,
  Blomqvist's beta and co-exceedance (with an independence benchmark).

## Tables and export

* `qqr_surface()` (wide or long tables), `rowave.qqr()`, `excel.qqr()` and
  `word.qqr()` (optionally with significance stars).

## Documentation and testing

* Vignettes: "Getting started with qqrtools" and a methodology vignette.
* Test suite with reference calculations against `quantreg` (about 94%
  line coverage) and a Monte Carlo study in `tools/validation`.
