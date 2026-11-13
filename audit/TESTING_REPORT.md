# Testing report – qqrtools 0.3.0

Environment for all local results: Windows 11 Pro (build 26300), R 4.6.2
(ucrt), quantreg 6.1, testthat 3e, roxygen2 8.1.1, devtools 2.5.2, covr
(CRAN, November 2026). Date: 13 November 2026.

## 1. Unit and integration tests (`devtools::test()`, `NOT_CRAN=true`)

| File | Test blocks | Expectations | Failed | Skipped |
|---|---|---|---|---|
| `test-qqr-basic.R` (original suite, kept as regression tests) | 15 | 137 | 0 | 0 |
| `test-math-reference.R` (independent reference calculations) | 10 | 96 | 0 | 0 |
| `test-validation.R` (inputs, edge cases, formula handling, selectors, SE types, 0.2.0 objects) | 13 | 75 | 0 | 0 |
| `test-plots-dependence.R` (graphics, colour scales, dependence measures) | 8 | 56 | 0 | 0 |
| `test-methods-export.R` (print/summary, tables, Excel/Word, bandwidth edge cases) | 9 | 63 | 0 | 0 |
| **Total** | **55** | **427** | **0** | **0** |

Without `NOT_CRAN` (as on CRAN), the *np*- and *rgl*-based tests are
skipped.

The original suite (137 expectations) passed unchanged against the new
code, except that one error-message pattern was kept by adapting the new
message (`"must be numeric"`). No test was weakened to make it pass; tests
that compare against independent results use tolerances of 1e-8 (exact
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
| R/qqr-core.R | 94.1 % |
| R/qqr-methods.R | 93.4 % |
| R/qqr-bandwidth.R | 89.8 % |
| R/qqr-qr.R | 89.3 % |
| **Total** | **93.7 %** |

Uncovered lines are defensive branches that cannot be reached with valid
inputs or need an unavailable configuration: missing optional packages
(`requireNamespace()` failure messages), failure fallbacks of *np*
selectors, "all local fits non-unique" message, the `"z"` alias of objects
created by versions before 0.2.0, and `stop()` calls for degenerate design
matrices that `model.frame()` already prevents. All numerical paths of the
estimator, every exported function and all bandwidth rules except the
*np* failure branches are covered.

## 3. R CMD check

| Run | Command | Result |
|---|---|---|
| Baseline (GitHub/local 0.1.0 code before changes) | `rcmdcheck(args = c("--as-cran", "--no-manual"))` | 0 errors, 0 warnings, 1 note (invalid URLs, invalid DOI) |
| Release candidate | `rcmdcheck(args = "--as-cran")` incl. PDF manual, examples, tests, vignettes | **0 errors, 0 warnings, 1 note** |
| Release candidate, hard dependencies only | `_R_CHECK_DEPENDS_ONLY_=true`, `--no-tests --ignore-vignettes` | 0 errors, 0 warnings, 0 notes (examples run without any Suggests package) |

The remaining note (CRAN incoming feasibility) is:

* "New submission" – informational;
* "possibly invalid URLs" `https://github.com/Baffour536/qqrtools`,
  `…/issues` and the README badge – HTTP 404 because the repository was
  still named `qqr` and private when the check ran. **This note resolves
  when the repository is renamed and made public** (re-check listed in the
  release report).

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

## 5. Performance (`tools/validation/benchmark.R`, single core)

m = K = 19 (361 local fits), one control:

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
| macOS, Linux (R release, devel, oldrel-1, 4.1) | GitHub Actions workflow `R-CMD-check.yaml` added; **not yet executed** – it runs on the first push. Result to be recorded in `RELEASE_REPORT.md`. |
| Clean environment | the CI runners start from an empty library; the local hard-dependency-only check above is the local evidence. Installation from GitHub with `remotes::install_github()` is verified after publication (see release report). |

## 7. Not executed

* `lintr` was installed but not used as a release gate: the code base uses a
  consistent but non-tidyverse style, and style lints are not defects.
* Tests on macOS and Linux and on R 4.1 – pending the CI run.
* Visual inspection of every generated figure.
