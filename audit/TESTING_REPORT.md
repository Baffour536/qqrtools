# Testing report – qqrtools 0.4.0

Environment for all local results: Windows 11 Pro (build 26300), R 4.6.2
(ucrt), quantreg 6.1, testthat 3e, roxygen2 8.1.1, devtools 2.5.2, covr
(CRAN, November 2026). Date: 13 November 2026.

## 1. Unit and integration tests (`devtools::test()`, `NOT_CRAN=true`)

| File | Test blocks | Expectations | Failed | Skipped |
|---|---|---|---|---|
| `test-qqr-basic.R` (structure and workflow tests) | 15 | 137 | 0 | 0 |
| `test-math-reference.R` (independent reference calculations) | 10 | 96 | 0 | 0 |
| `test-validation.R` (inputs, edge cases, formula handling, selectors, SE types, legacy objects) | 13 | 75 | 0 | 0 |
| `test-plots-dependence.R` (graphics, colour scales, dependence measures) | 8 | 56 | 0 | 0 |
| `test-methods-export.R` (print/summary, tables, Excel/Word, bandwidth edge cases) | 9 | 63 | 0 | 0 |
| `test-bandwidth-degree.R` (local quadratic fit, slope bandwidth rules, sensitivity analysis, p-value adjustment) | 7 | 113 | 0 | 0 |
| **Total** | **62** | **540** | **0** | **0** |

Without `NOT_CRAN` (as on CRAN), the *np*- and *rgl*-based tests are
skipped.

Tests that compare against independent results use tolerances of 1e-8 (exact
algorithms) or explicitly justified statistical tolerances
(location-scale recovery).

Integration workflows covered: simulate → `qqr()` → `summary()` →
`qqr_surface()`/`rowave.qqr()` → `compare.qqr()` → plots → Excel/Word
export (files read back and contents checked).

## 2. Coverage (`covr::package_coverage(type = "tests")`, `NOT_CRAN=true`)

| File | Coverage |
|---|---|
| R/qqr-utils.R | 98.6 % |
| R/qqr-export.R | 98.0 % |
| R/qqr-robustness.R | 100.0 % |
| R/qqr-core.R | 94.6 % |
| R/qqr-methods.R | 93.5 % |
| R/qqr-bandwidth.R | 91.1 % |
| R/qqr-qr.R | 89.3 % |
| **Total** | **94.1 %** |

Uncovered lines are defensive branches that cannot be reached with valid
inputs or need an unavailable configuration: missing optional packages
(`requireNamespace()` failure messages), failure fallbacks of *np*
selectors, "all local fits non-unique" message, a legacy `"z"` coefficient alias, and `stop()` calls for degenerate design
matrices that `model.frame()` already prevents. All numerical paths of the
estimator, every exported function and all bandwidth rules except the
*np* failure branches are covered.

## 3. R CMD check

| Run | Command | Result |
|---|---|---|
| Version 0.4.0 | `rcmdcheck(args = "--as-cran")` incl. PDF manual, examples, tests, vignettes | **0 errors, 0 warnings, 0 notes** |

The hard-dependencies-only run (`_R_CHECK_DEPENDS_ONLY_=true`) was carried
out for version 0.3.0 and not repeated for 0.4.0; version 0.4.0 adds no
dependency (the new functions use only *stats* and *graphics*).

## 4. Other checks

| Check | Result |
|---|---|
| `devtools::document()` | NAMESPACE and 25 Rd files regenerated without warnings |
| README code blocks | all 3 executable blocks run in a fresh `Rscript` session |
| Examples | all run (`R CMD check`), file output to `tempdir()` |
| Vignettes | `qqrtools.Rmd` (fully evaluated) and `qqr-methodology.Rmd` build during check |
| Non-ASCII / portability | `R CMD check` OK; R sources ASCII |
| Undeclared dependencies | none (`checking dependencies in R code … OK`) |
| Developer paths | none in package files (`grep` for `C:/Users`, `OneDrive`, `/home/`) |
| Reproducibility | all tests and examples set seeds; bootstrap p-values reproducible with `set.seed()` (tested) |
| Monte Carlo study | `tools/validation/simulation.R`, 52 settings × 200 replications, results in `tools/validation/simulation_summary.md` |
| Calibration of defaults | `tools/validation/bandwidth_calibration.R` (3 designs × 2 sample sizes × up to 16 rule–degree combinations × 200 replications) and `tools/validation/size_default.R` (2 null designs × 2 sample sizes × 3 standard-error types × 300 replications); results in the `.md` files beside the scripts |

## 5. Performance (`tools/validation/benchmark.R`, single core)

m = K = 19 (361 local fits), one control, local linear fit (`degree = 1`);
timings measured for version 0.3.0 and not repeated. The default local
quadratic fit estimates one more coefficient per cell:

| n | h | se_type | seconds |
|---|---|---|---|
| 200 | silverman | none | 1.48* |
| 200 | silverman | nid | 0.82 |
| 1000 | silverman | none | 0.90 |
| 1000 | silverman | nid | 1.69 |
| 5000 | silverman | none | 4.84 |
| 5000 | silverman | nid | 10.05 |
| 1000 | cv_fast | none | 2.70 |
| 3000 | cv_fast | none | 18.98 |
| 1000 | sj | none | 0.80 |
| 1000 | yj | none | 0.83 |

\* first call includes package and namespace loading.

Cost is linear in m × K and roughly linear in n for the local fits;
`cv_fast` is quadratic in n. No optimisation was required.

## 6. Platforms

| Platform | Status |
|---|---|
| Windows 11, R 4.6.2 | all checks above executed locally |
| GitHub Actions `R-CMD-check` | runs on macOS (R release), Windows (R release) and Ubuntu (R devel, release, oldrel-1) for every push, together with a coverage job; see the Actions tab of the repository for the status of the current commit |
| R 4.1 on Ubuntu | not tested: current CRAN versions of the dependencies could not be installed on R 4.1 on the CI runner, so that job was removed |
| Clean environment | the CI runners start from an empty library; `remotes::install_github("Baffour536/qqrtools")` into an empty library followed by the README quick start succeeded on Windows for version 0.3.0 |

## 7. Not executed

* `lintr` was installed but not used as a release gate: the code base uses a
  consistent but non-tidyverse style, and style lints are not defects.
* Tests on R 4.1 (see section 6).
* Visual inspection of every generated figure.
