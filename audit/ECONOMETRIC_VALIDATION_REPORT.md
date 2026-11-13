# Econometric validation report – qqrtools 0.4.0

This report states what `qqr()` estimates, checks it against the cited
method and against independent computations, and reports a Monte Carlo
study. All numbers come from code in this repository
(`tests/testthat/test-math-reference.R`,
`tests/testthat/test-bandwidth-degree.R`, `tools/validation/simulation.R`,
`tools/validation/bandwidth_calibration.R`, `tools/validation/size_default.R`);
they can be reproduced with R ≥ 4.1 and the package installed.

## 1. Method actually implemented

For a response *Y*, a designated regressor *X* (`qvar`) and controls *W*
(numeric or factor dummies), `qqr()` computes, for each grid point
x<sup>θ</sup> and each τ, the solution of

  min<sub>α,β,γ</sub> Σ<sub>t</sub> ρ<sub>τ</sub>(Y<sub>t</sub> − α − β(X<sub>t</sub> − x<sup>θ</sup>) − γ′W<sub>t</sub>) K(d<sub>t</sub>/h),

with ρ<sub>τ</sub>(u) = u(τ − 1{u < 0}), Gaussian *K*, and

* d<sub>t</sub> = X<sub>t</sub> − x<sup>θ</sup> (`kernel_space = "x"`, default), or
* d<sub>t</sub> = F<sub>n</sub>(X<sub>t</sub>) − θ (`kernel_space = "cdf"`).

Grid: x<sup>θ</sup> = Q̂<sub>X</sub>(θ) (type-7 sample quantile) at
`thetas` (default `seq(0.05, 0.95, length.out = m)`), or `m` equally
spaced points over [min X, max X] (`grid = "range"`).
With `degree = 2` (default) the term δ(X<sub>t</sub> − x<sup>θ</sup>)² is
added to the local model (local quadratic); δ̂ is not stored. With
`degree = 1` the problem is as displayed (local linear).
The weighted problem is solved exactly by `quantreg::rq(method = "br")`.
P-values are those of `quantreg::summary.rq(se = se_type)` for the weighted
problem.

### Correspondence with Sim and Zhou (2015)

Sim and Zhou estimate β<sub>1</sub>(θ, τ) from a local linear quantile
regression of *Y* on X<sub>t</sub> − X<sup>θ</sup> with weights
K((F<sub>n</sub>(X<sub>t</sub>) − θ)/h), at θ on a grid, with h = 0.05,
and compare the θ-average of β̂<sub>1</sub>(θ, τ) with the QR slope.

| Element | Sim and Zhou (2015) | qqrtools |
|---|---|---|
| Location of local fits | sample θ-quantiles of X | θ-quantiles (default); equally spaced range optional |
| Kernel distance | F<sub>n</sub>(X<sub>t</sub>) − θ | X<sub>t</sub> − x (default) or F<sub>n</sub>(X<sub>t</sub>) − θ |
| Local regressor | X<sub>t</sub> − X<sup>θ</sup> | same |
| Controls | not in the original model (bivariate) | linear, with local coefficients |
| QR comparison | average over θ for each τ | average over θ for each τ (`compare.qqr()`) |

With `grid = "quantile"`, `kernel_space = "cdf"` and `h = 0.05` (and no
controls) `qqr()` implements the Sim–Zhou estimator. The default
`kernel_space = "x"` is the local linear quantile regression of Yu and
Jones (1998) evaluated at the quantiles of *X*. When *X* is uniform the
two weightings coincide (F<sub>n</sub>(x) ≈ x); the Monte Carlo study
confirms that both settings then give the same estimates to three decimals
(§3, compare the `x` and `cdf` rows for h = 0.1, 0.2).

## 2. Deterministic reference checks

All checks pass (see `TESTING_REPORT.md`).

| Check | Reference | Result |
|---|---|---|
| Each surface cell (intercept, slope, control) and its `nid` p-value | direct `quantreg::rq(y ~ I(x − g) + w, weights = dnorm((x − g)/h))` and `summary(…, se = "nid")` | equal to 1e-8 for all 4 × 3 cells × 3 coefficients |
| CDF kernel | weights `dnorm((ecdf(x)(x) − θ)/h)` computed by hand | equal to 1e-8 |
| Large-bandwidth limit (h = 10<sup>6</sup>) | global `rq()`; local intercept = QR intercept + slope · x<sup>θ</sup> | equal to 1e-6; `compare.qqr()` difference ≈ 0 |
| Noiseless linear model *y* = 2 − 3x + 0.5w | exact | every cell equals −3 and 0.5 to 1e-8 |
| Location-scale model, n = 4000, h = 0.3 | β = 2 + 0.5 Φ<sup>−1</sup>(τ) | max cell error < 0.35, θ-average error < 0.25 |
| `compare.qqr()` | `colMeans` of the surface; `coef(rq())` | exact |
| Bandwidth rules | closed forms; `bw.SJ`; `KernSmooth::dpill`; brute-force LOO check-loss | exact |
| Kish effective sample size | (Σw)²/Σw² | exact |
| Transformed regressor | `qvar = "log(ex)"` with ex = e<sup>x</sup> vs `qvar = "x"` | identical surfaces |
| Somers' D, Blomqvist's β, co-exceedance, cell correlations | brute-force definitions | exact |

## 3. Monte Carlo study of the local linear estimator

This section documents the local linear estimator (`degree = 1`) with
density-type and manual bandwidths. It motivates the defaults introduced
in version 0.4.0, which are evaluated in section 3a.

Script `tools/validation/simulation.R`; 200 replications per setting; seed
base 20261113; θ = τ ∈ {0.1, 0.25, 0.5, 0.75, 0.9}; X ~ U(0, 1),
ε ~ N(0, 1), no controls. Full table: `tools/validation/simulation_summary.md`.

* **A (location-scale):** Y = 1 + 2X + (0.5 + X)ε, true β(θ, τ) = 2 + Φ<sup>−1</sup>(τ) (constant in θ; the local linear model is exact).
* **B (curved):** Y = 1 + 2X² + (0.5 + X)ε, true β(θ, τ) = 4θ + Φ<sup>−1</sup>(τ).
* **C (null):** Y = 1 + ε, β = 0; used for the size of the pointwise tests.

Average over the 25 cells (`kernel_space = "x"`):

| Design | n | h | mean \|bias\| | mean SD | mean RMSE | RMSE of θ-average |
|---|---|---|---|---|---|---|
| A | 250 | Silverman (≈0.10) | 0.129 | 1.395 | 1.400 | 0.624 |
| A | 1000 | Silverman (≈0.08) | 0.054 | 0.974 | 0.974 | 0.386 |
| A | 1000 | 0.1 | 0.036 | 0.707 | 0.706 | 0.316 |
| A | 1000 | 0.2 | 0.014 | 0.309 | 0.308 | 0.181 |
| A | 1000 | 0.3 | 0.010 | 0.203 | 0.203 | 0.155 |
| B | 1000 | Silverman (≈0.08) | 0.086 | 0.978 | 0.982 | 0.388 |
| B | 1000 | 0.1 | 0.110 | 0.708 | 0.721 | 0.317 |
| B | 1000 | 0.2 | 0.404 | 0.315 | 0.521 | 0.196 |
| B | 1000 | 0.3 | 0.686 | 0.211 | 0.728 | 0.198 |

Bias by cell, design B, n = 1000, h = 0.2 (rows θ, columns τ):

| θ \ τ | 0.1 | 0.25 | 0.5 | 0.75 | 0.9 |
|---|---|---|---|---|---|
| 0.1 | +0.58 | +0.57 | +0.58 | +0.56 | +0.58 |
| 0.25 | +0.24 | +0.25 | +0.24 | +0.23 | +0.24 |
| 0.5 | −0.09 | −0.10 | −0.12 | −0.14 | −0.13 |
| 0.75 | −0.36 | −0.40 | −0.41 | −0.43 | −0.45 |
| 0.9 | −0.65 | −0.65 | −0.69 | −0.70 | −0.74 |

Size of the pointwise 5 % test of β(θ, τ) = 0 (design C, pooled over 25
cells; per-cell rates in the summary file):

| n | h | rejection rate |
|---|---|---|
| 250 | Silverman (≈0.10) | 0.079 |
| 1000 | Silverman (≈0.08) | 0.072 |
| 250 | 0.2 | 0.062 |
| 1000 | 0.2 | 0.051 |

### Interpretation

1. **Consistency of the implementation.** In design A, where the local
   linear model is exact, the bias is small at every bandwidth and shrinks
   with n; the SD falls with n and with h. This is the behaviour expected
   of a correctly implemented local linear quantile estimator.
2. **Bias–variance trade-off.** In design B the slope varies with θ. With
   a wide bandwidth the slope estimates at the edges of the θ-grid are
   strongly biased (+0.58 at θ = 0.1, −0.70 at θ = 0.9 for h = 0.2): the
   kernel window is cut off by the boundary of the support and averages
   the slope over an asymmetric region. Interior cells (θ = 0.5) are
   nearly unbiased. This is the known boundary behaviour of local linear
   derivative estimates, not a coding error.
3. **Density-oriented bandwidths give noisy surfaces.** The density-oriented
   rules (Silverman ≈ 0.08–0.10 here, robust even smaller) give cell SDs
   around 1.0–1.4 at n = 250–1000 for slopes between about 0.7 and 3.3. The
   variance of a local slope scales like 1/(n h³), so rules of order
   n<sup>−1/5</sup> designed for density estimation under-smooth slopes. The
   θ-averaged slope used by `compare.qqr()` is much more precise.
   **Practical consequence:** single-cell estimates at a density-type bandwidth
   should not be interpreted in isolation in samples of a few hundred
   observations; report *h*, show results for larger *h*, and focus on
   patterns that are stable across bandwidths.
4. **Pointwise p-values** are close to nominal for h = 0.2 (5.1–6.2 %) and
   moderately oversized for the small Silverman bandwidth (7.2–7.9 %), with
   the largest distortions in the corner cells (up to 13 % at θ = τ = 0.1,
   n = 250). They are not adjusted for testing 25 or more cells; with 361
   cells (19 × 19) many "significant" cells are expected under the null.
5. **CDF vs X weighting.** With X uniform the two kernel spaces gave
   identical summaries (e.g. A, n = 1000, h = 0.2: RMSE 0.308 vs 0.309),
   as theory implies. With h = 0.05 on the CDF scale (Sim and Zhou's
   choice) the cell RMSE at n = 1000 is 1.68, i.e. that bandwidth is very
   small for slope estimation at this sample size.

## 3a. Calibration of the default bandwidth, degree and standard errors (version 0.4.0)

**Bandwidth and degree.** Script `tools/validation/bandwidth_calibration.R`;
200 replications; θ = τ ∈ {0.1, 0.25, 0.5, 0.75, 0.9}; n = 250 and 1000.
Designs (Y = m(X) + s(X)ε, ε ~ N(0, 1)):

* **D1:** X ~ U(0, 1), m = 1 + 2x², s = 0.5 + x (slope linear in x);
* **D2:** X ~ N(0, 1), m = sin(1.5x), s = 0.5 + 0.25|x| (rapidly varying slope);
* **D3:** X ~ χ²(3), m = 2(1 − e^(−0.3x)), s = 0.5 + 0.1x (skewed regressor).

Rules h = C·min(sd, IQR/1.34)·n^(−1/7) with C ∈ {1.06, 1.5, 2, 2.5, 3},
Silverman's rule, and (n = 250 only) `cv_fast` and `cv_slope`, each with
`degree = 1` and `degree = 2`. Full table:
`tools/validation/bandwidth_calibration.md`. Mean cell RMSE relative to the
best of the 12 rule–degree combinations available at both sample sizes:

| Design, n | C = 2, quadratic (default) | C = 1.5, quadratic | C = 2.5, quadratic | C = 2, linear | Silverman, linear |
|---|---|---|---|---|---|
| D1 250 | 1.29 | 1.58 | 1.15 | 1.01 | 1.86 |
| D1 1000 | 1.28 | 1.61 | 1.10 | 1.28 | 2.30 |
| D2 250 | 1.12 | 1.00 | 1.27 | 1.81 | 1.21 |
| D2 1000 | 1.28 | 1.07 | 1.54 | 2.51 | 1.25 |
| D3 250 | 1.17 | 1.40 | 1.07 | 1.18 | 1.67 |
| D3 1000 | 1.18 | 1.46 | 1.05 | 1.61 | 2.12 |
| worst case | 1.29 | 1.61 | 1.54 | 2.51 | 2.30 |

Findings:

1. C = 2 with the local quadratic fit has the smallest worst-case loss
   (29 % above the best combination); it is the default
   (`h = "slope"`, `degree = 2`). No combination is best everywhere:
   smoother designs (D1, D3) favour wider bandwidths, the rapidly varying
   design D2 favours narrower ones.
2. At C = 2 the local quadratic fit reduces the mean absolute bias from
   0.57 to 0.07 (D1, n = 250), 0.58 to 0.29 (D2) and 0.08 to 0.02 (D3)
   relative to the local linear fit, at the cost of a larger standard
   deviation (0.97 vs 0.49, 0.23 vs 0.15, 0.10 vs 0.06). In D2 a bias of
   0.24–0.29 remains.
3. Silverman's rule with a local linear fit (the 0.3.0 default) is 21–130 %
   above the best combination.
4. `cv_slope` improved on `cv_fast` in all three designs (relative RMSE
   1.76 vs 2.56, 1.34 vs 1.62, 1.76 vs 2.51 with `degree = 1`) but was
   less accurate than the rule of thumb. Both are limited by the search
   range of `cv_fast` ([0.5, 2] × Silverman).

**Standard errors.** Script `tools/validation/size_default.R`; null designs
(Y = 1 + ε independent of X; N: X ~ N(0, 1), S: X ~ χ²(3)); 300
replications; default bandwidth and degree; nominal level 5 %; bootstrap
with R = 100. The last two columns are the share of fitted surfaces with
at least one significant cell before and after
`qqr_padjust(method = "BH")`:

| design | n | se_type | pooled_rejection_rate | max_cell_rate | any_cell_significant | any_cell_significant_BH |
|---|---|---|---|---|---|---|
| N | 250 | nid | 0.093 | 0.170 | 0.820 | 0.407 |
| N | 250 | ker | 0.035 | 0.090 | 0.483 | 0.060 |
| N | 250 | boot | 0.050 | 0.077 | 0.577 | 0.087 |
| N | 1000 | nid | 0.079 | 0.153 | 0.747 | 0.280 |
| N | 1000 | ker | 0.042 | 0.080 | 0.513 | 0.080 |
| N | 1000 | boot | 0.050 | 0.080 | 0.567 | 0.087 |
| S | 250 | nid | 0.076 | 0.157 | 0.740 | 0.263 |
| S | 250 | ker | 0.042 | 0.093 | 0.517 | 0.083 |
| S | 250 | boot | 0.050 | 0.067 | 0.500 | 0.097 |
| S | 1000 | nid | 0.074 | 0.150 | 0.767 | 0.237 |
| S | 1000 | ker | 0.044 | 0.100 | 0.563 | 0.097 |
| S | 1000 | boot | 0.052 | 0.073 | 0.613 | 0.103 |

Findings:

1. `"nid"` is oversized (7.4–9.3 % pooled, up to 17 % in single cells),
   `"ker"` is slightly conservative (3.5–4.4 %), `"boot"` is at the nominal
   level (5.0–5.2 %). The default is `se_type = "ker"`; `"boot"` is
   recommended for reported results.
2. Without adjustment, 48–82 % of null surfaces contain at least one
   significant cell (25 cells). After the Benjamini–Hochberg adjustment
   the share is 6–10 % for `"ker"` and `"boot"`, somewhat above 5 %, and
   24–41 % for `"nid"`: the adjustment cannot repair p-values that are too
   small.

Limitations: three calibration designs and two null designs, independent
observations, no controls, a 5 × 5 grid. The results do not establish the
behaviour for dependent (time-series) data, which is the typical QQR
application.

## 4. Inference: what is and is not supported

* `se_type = "ker"` (default), `"nid"`, `"boot"` are passed to
  `quantreg::summary.rq()` on the weighted local problem. Bootstrap
  results are reproducible with `set.seed()` (tested).
* No claim is made that these p-values are valid for the nonparametric
  surface: they ignore smoothing bias, data-driven bandwidths and
  multiplicity. The documentation (help page, README, both vignettes)
  says so.
* `qqr_padjust()` adjusts each p-value surface for multiplicity with
  `stats::p.adjust()` (Benjamini–Hochberg by default). It cannot repair
  p-values that are themselves too small; see section 3a.
* No uniform confidence bands or explicit bias correction are provided.

## 5. Dependence diagnostics (`qqr_cor_heatmap()`)

Spearman, Kendall and Pearson use `stats::cor()`; Somers' D(Y|X),
Blomqvist's β and co-exceedance use internal code checked against
brute-force definitions. Co-exceedance is now P̂(X ≥ Q̂<sub>X</sub>(τ<sub>i</sub>),
Y ≥ Q̂<sub>Y</sub>(τ<sub>j</sub>)) with benchmark (1 − τ<sub>i</sub>)(1 − τ<sub>j</sub>).
Within-cell correlations are computed on the observations falling into
each joint quantile cell. Because the cells truncate both variables, these
correlations are local, noisy summaries (cells have about n/K² observations).
They are not estimates of a population copula quantity.

## 6. Conclusion

Within the scope tested, the evidence supports the following:

* the implementation computes exactly the weighted local linear or local
  quadratic quantile regression it documents (deterministic equality with
  `quantreg`);
* with `grid = "quantile"`, `kernel_space = "cdf"` and a manual *h* it
  and `degree = 1` it implements the Sim and Zhou (2015) estimator; the
  default configuration is a local quadratic estimator evaluated at the
  regressor's quantiles with a slope-oriented bandwidth;
* `compare.qqr()` implements the Sim–Zhou robustness check;
* the estimator behaves as theory predicts (unbiased under an exact local
  model, boundary bias and bias–variance trade-off otherwise);
* p-values are approximate, pointwise and should be treated as descriptive.

The default bandwidth and degree were chosen on the simulation evidence of
section 3a. They are a robust compromise, not an optimum for every data
set: the three calibration designs cannot cover all applications, and in
the design with a rapidly varying slope the default still has a visible
bias. Results should always be accompanied by `qqr_sensitivity()`.

## References

Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B*, 57(1), 289–300. doi:10.1111/j.2517-6161.1995.tb02031.x

Fan, J. and Gijbels, I. (1996) *Local Polynomial Modelling and Its Applications*. London: Chapman & Hall.

Koenker, R. (2005) *Quantile Regression*. Cambridge University Press. doi:10.1017/CBO9780511754098

Sim, N. and Zhou, H. (2015) Oil prices, US stock return, and the dependence between their quantiles. *Journal of Banking & Finance*, 55, 1–8. doi:10.1016/j.jbankfin.2015.01.013

Somers, R. H. (1962) A new asymmetric measure of association for ordinal variables. *American Sociological Review*, 27(6), 799–811. doi:10.2307/2090408

Yu, K. and Jones, M. C. (1998) Local linear quantile regression. *Journal of the American Statistical Association*, 93(441), 228–237. doi:10.1080/01621459.1998.10474104
