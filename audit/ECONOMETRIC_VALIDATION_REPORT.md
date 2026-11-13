# Econometric validation report – qqrtools 0.3.0

This report states what `qqr()` estimates, checks it against the cited
method and against independent computations, and reports a Monte Carlo
study. All numbers come from code in this repository
(`tests/testthat/test-math-reference.R`, `tools/validation/simulation.R`);
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
spaced points over [min X, max X] (`grid = "range"`, the 0.2.0 behaviour).
The weighted problem is solved exactly by `quantreg::rq(method = "br")`.
P-values are those of `quantreg::summary.rq(se = se_type)` for the weighted
problem.

### Correspondence with Sim and Zhou (2015)

Sim and Zhou estimate β<sub>1</sub>(θ, τ) from a local linear quantile
regression of *Y* on X<sub>t</sub> − X<sup>θ</sup> with weights
K((F<sub>n</sub>(X<sub>t</sub>) − θ)/h), at θ on a grid, with h = 0.05,
and compare the θ-average of β̂<sub>1</sub>(θ, τ) with the QR slope.

| Element | Sim and Zhou (2015) | qqrtools 0.2.0 | qqrtools 0.3.0 |
|---|---|---|---|
| Location of local fits | sample θ-quantiles of X | equally spaced over range (labelled as quantiles) | θ-quantiles (default); range optional |
| Kernel distance | F<sub>n</sub>(X<sub>t</sub>) − θ | X<sub>t</sub> − x | X<sub>t</sub> − x (default) or F<sub>n</sub>(X<sub>t</sub>) − θ |
| Local regressor | X<sub>t</sub> − X<sup>θ</sup> | same | same |
| Controls | not in the original model (bivariate) | linear, local coefficients | linear, local coefficients |
| QR comparison | average over θ for each τ | average over τ, paired by index (wrong) | average over θ for each τ |

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

## 3. Monte Carlo study

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
3. **The default bandwidth gives noisy surfaces.** The density-oriented
   rules (Silverman ≈ 0.08–0.10 here, robust even smaller) give cell SDs
   around 1.0–1.4 at n = 250–1000 for slopes between about 0.7 and 3.3. The
   variance of a local slope scales like 1/(n h³), so rules of order
   n<sup>−1/5</sup> designed for density estimation under-smooth slopes. The
   θ-averaged slope used by `compare.qqr()` is much more precise.
   **Practical consequence:** single-cell estimates at the default bandwidth
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

## 4. Inference: what is and is not supported

* `se_type = "nid"`, `"ker"`, `"boot"` are passed to
  `quantreg::summary.rq()` on the weighted local problem. Bootstrap
  results are reproducible with `set.seed()` (tested).
* No claim is made that these p-values are valid for the nonparametric
  surface: they ignore smoothing bias, data-driven bandwidths and
  multiplicity. The documentation (help page, README, both vignettes) now
  says so.
* No uniform confidence bands, bias correction, or bandwidth selection
  tailored to QQR are provided.

## 5. Dependence diagnostics (`qqr_cor_heatmap()`)

Spearman, Kendall and Pearson use `stats::cor()`; Somers' D(Y|X),
Blomqvist's β and co-exceedance use internal code checked against
brute-force definitions. Somers' D was corrected (it previously computed
D(X|Y)). Co-exceedance is now P̂(X ≥ Q̂<sub>X</sub>(τ<sub>i</sub>),
Y ≥ Q̂<sub>Y</sub>(τ<sub>j</sub>)) with benchmark (1 − τ<sub>i</sub>)(1 − τ<sub>j</sub>).
Within-cell correlations are computed on the observations falling into
each joint quantile cell. Because the cells truncate both variables, these
correlations are local, noisy summaries (cells have about n/K² observations).
They are not estimates of a population copula quantity.

## 6. Conclusion

Within the scope tested, the evidence supports the following:

* the implementation computes exactly the weighted local linear quantile
  regression it documents (deterministic equality with `quantreg`);
* with `grid = "quantile"`, `kernel_space = "cdf"` and a manual *h* it
  implements the Sim and Zhou (2015) estimator; the default configuration
  is the closely related Yu–Jones local linear estimator evaluated at the
  regressor's quantiles;
* `compare.qqr()` now implements the Sim–Zhou robustness check;
* the estimator behaves as theory predicts (unbiased under an exact local
  model, boundary bias and bias–variance trade-off otherwise);
* p-values are approximate, pointwise and should be treated as descriptive.

The defaults favour transparency over accuracy. The main methodological
limitation is that the automatic bandwidth rules are not tuned to QQR
slopes. Results should always be accompanied by a sensitivity analysis
over *h*.

## References

Koenker, R. (2005) *Quantile Regression*. Cambridge University Press. doi:10.1017/CBO9780511754098

Sim, N. and Zhou, H. (2015) Oil prices, US stock return, and the dependence between their quantiles. *Journal of Banking & Finance*, 55, 1–8. doi:10.1016/j.jbankfin.2015.01.013

Somers, R. H. (1962) A new asymmetric measure of association for ordinal variables. *American Sociological Review*, 27(6), 799–811. doi:10.2307/2090408

Yu, K. and Jones, M. C. (1998) Local linear quantile regression. *Journal of the American Statistical Association*, 93(441), 228–237. doi:10.1080/01621459.1998.10474104
