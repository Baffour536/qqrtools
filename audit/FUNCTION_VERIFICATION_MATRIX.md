# Function verification matrix – qqrtools 0.3.0

Every exported function and every internal function on a numerical or
data-handling path is listed. "Tests" names the test files in
`tests/testthat/` that exercise the function (`basic` = `test-qqr-basic.R`,
`math` = `test-math-reference.R`, `val` = `test-validation.R`,
`plots` = `test-plots-dependence.R`, `mx` = `test-methods-export.R`).
"Mathematical validation" describes the independent reference against which
results are compared. Status: **Verified** = reviewed, corrected where
needed, and tested against an independent reference; **Verified (run)** =
reviewed and tested for correct structure/behaviour (graphics and file
output, where there is no numerical reference); **Partial** = see note.

## Exported functions

| Function | Purpose | Tests | Mathematical validation | Status |
|---|---|---|---|---|
| `qqr()` | Fit QQR surfaces by kernel-weighted local linear quantile regression | basic, math, val, mx | Every cell equals a direct `quantreg::rq(weights=)` fit (coefficients and `nid` p-values, tol 1e-8); Sim–Zhou CDF weights recomputed by hand; huge-h limit equals global QR; noiseless linear model recovered exactly; location-scale slope recovered; Monte Carlo study (audit/ECONOMETRIC_VALIDATION_REPORT.md) | Verified |
| `print.qqr()` | Console summary of a fit | basic, mx, val | Output content checked (grid type, kernel space, bandwidth label) | Verified (run) |
| `summary.qqr()` | Global slope-surface summaries | basic, mx | `row_means`/`col_means` equal `rowMeans`/`colMeans` of `qcoef` | Verified |
| `print.summary.qqr()` | Print a `summary.qqr` | mx | Output content checked | Verified (run) |
| `qqr_surface()` | Wide/long tables of surfaces | basic, mx | Long values equal `as.vector(t(surface))`; tau matching with tolerance; subsets | Verified |
| `rowave.qqr()` | Averages over tau for each grid point | basic, mx | Equals `rowMeans` of the surface | Verified |
| `compare()` | S3 generic | math | Dispatch to `compare.qqr` | Verified (run) |
| `compare.qqr()` | Theta-averaged QQR vs standard QR | basic, math, val, mx | `qr` equals `coef(rq())` on the same design (tol 1e-8); `qqr_mean` equals `colMeans` of the surface; equality of both in the huge-h limit; factor controls | Verified |
| `plot.compare.qqr()` | ggplot of the comparison | basic, mx | Plot data equal comparison table | Verified (run) |
| `excel.compare.qqr()` | Comparison to Excel | basic, mx | Sheet names checked | Verified (run) |
| `word.compare.qqr()` | Comparison to Word | basic | File written | Verified (run) |
| `excel.qqr()` | Surfaces to Excel | basic, mx | Sheet contents read back equal the rounded surface | Verified |
| `word.qqr()` | Surfaces to Word (optionally with stars) | basic, mx | Text read back: `[NA]` marker, fixed-decimal coefficients | Verified (run) |
| `plot.qqr()` | Static 3D plot (plot3D or base) | basic, plots | Axis = stored thetas | Verified (run) |
| `Int.plot()`, `Int.plot.qqr()` | Interactive 3D plot (rgl) | basic, plots | – (rgl null device) | Verified (run) |
| `surface3d()`, `surface3d.qqr()` | Coloured 3D surface (base) | plots, val | Axis = stored thetas; argument validation | Verified (run) |
| `panel_surface()`, `panel_surface.qqr()` | Side-by-side 3D surfaces with shared scale | plots | Shared limits; argument validation | Verified (run) |
| `qqr_heatmap()` | plotly heatmap with significance overlay | plots | Colour scale orientation checked (negative = blue) | Verified (run) |
| `qqr_contour()` | plotly contour plot | plots | as above | Verified (run) |
| `qqr_plotly()` | plotly surface/heatmap/contour hub | plots | as above | Verified (run) |
| `qqr_cor_heatmap()` | Quantile dependence heatmap | plots | Cell Spearman recomputed by hand; bin breaks; co-exceedance and its independence benchmark recomputed; binary-mode Pearson recomputed | Verified |
| `print.qqr_cor_heatmap()` | Print dependence summary | plots | Output checked | Verified (run) |

## Internal functions

| Function | Purpose | Tests | Mathematical validation | Status |
|---|---|---|---|---|
| `.qqr_select_h()` | Eight bandwidth rules | math, val, mx | Silverman, robust, IQR equal closed forms; `sj` equals `bw.SJ`; `yj` equals `KernSmooth::dpill`; `cv_fast` equals a brute-force LOO check-loss minimiser; error and fallback branches | Verified (`npcv`, `sc`: run only, results are those of *np*) |
| `.qqr_local_fit()` | One weighted local `rq()` fit with p-values | math, mx | Equal to direct `rq()`; failure paths (rank deficiency, Kish ESS < p+1, solver error) | Verified |
| `.qqr_kish()` | Kish effective sample size | math, mx | Equals (Σw)²/Σw² | Verified |
| `.qqr_kernel()` | Gaussian kernel | math, mx | `dnorm` | Verified |
| `.qqr_ecdf_ranks()` | Empirical CDF F_n(X_t) | math | Equals `ecdf(x)(x)` | Verified |
| `.qqr_check_probs()` | Validate taus/thetas | val | Error cases | Verified |
| `.qqr_match_coef()`, `.qqr_user_coefs()`, `.qqr_resolve_coefs()` | Map user names to surfaces | basic, val, plots, mx | Intercept aliases, old `"z"` alias, errors | Verified |
| `.qqr_thetas()` | Grid quantile levels (also for 0.2.0 objects) | val, plots | Equals ECDF at grid for old objects | Verified |
| `.qqr_get_surface()` | Extract a surface; error on empty p-values | plots | Error message for `se_type = "none"` | Verified |
| `.qqr_design()` | Design matrix for `compare.qqr` (also old objects) | val | Same QR as stored design | Verified |
| `.qqr_somers_d()` | Somers' D(Y\|X) | plots | Brute-force definition with ties; equals Kendall tau-a without ties | Verified (corrected) |
| `.qqr_blomqvist()` | Blomqvist's beta | plots | Sign definition | Verified |
| `.qqr_coexceedance()` | Joint exceedance probability | plots | Empirical mean | Verified |
| `.qqr_cell_stat()` | Dispatch of dependence measures | plots | via the above | Verified |
| `.qqr_colorscale()`, `.qqr_rdbu()`, `.qqr_diverging()` | plotly colour scales | plots | Low end blue, high end red; hcl palettes | Verified |
| `.qqr_zlim()` | Colour limits | plots | Symmetric limits, degenerate range | Verified |
| `.qqr_trans3d()` | Projection for zero planes | plots | Used by drawn zero planes (visual) | Verified (run) |
| `.qqr_h_label()` | Labels for print methods | mx | All methods | Verified |

## Removed

| Function | Reason |
|---|---|
| duplicate `.qqr_match_coef()` in `qqr-methods.R` | Silently overrode the version in `qqr-core.R` (collation order); one definition kept |
| `.qqr_lprq_point()` | Replaced by `.qqr_local_fit()` (formula-safe names, failure accounting) |
