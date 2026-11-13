# Release report – qqrtools 0.3.0

## Summary

`qqrtools` 0.3.0 is the first release of the package formerly called
`qqr`. It consolidates the unpublished local development (renaming and a
plotly visualisation suite) with a full review that corrected three
critical methodological defects. Those defects affected the estimation
grid, `compare.qqr()` and silently empty surfaces for factor controls and
transformed regressors. It also adds the Sim–Zhou empirical-CDF weighting,
input validation, 290 new test expectations and validation reports.

## Version choice

The GitHub repository contained `qqr` 0.2.0; the local development copy
had been set to `qqrtools` 0.1.0 but was never released. Because 0.3.0
changes default results (grid) and corrects `compare.qqr()`, a minor
version bump over the last published 0.2.0 is the appropriate semantic
version for a 0.x package; 0.3.0 also avoids a version that sorts below
the previously published one.

## Changes from qqr to qqrtools

* Package name `qqr` → `qqrtools` (DESCRIPTION, NAMESPACE, help pages,
  tests, vignettes, README, CITATION, RStudio project). Function names and
  the `"qqr"` class are unchanged, because "QQR" is the name of the method.
* Name availability: no CRAN package called `qqrtools` or `qqr` (checked
  against the CRAN package list on 13 Nov 2026). The related CRAN package
  `QuantileOnQuantile` is a different project.
* Migration: replace `library(qqr)` with `library(qqrtools)`; use
  `grid = "range"` to reproduce 0.2.0 grids; recompute `compare.qqr()`
  results.

## New features

`kernel_space = "cdf"` (Sim–Zhou weighting), `grid`/`thetas` arguments,
factor controls, transformed regressors/responses, models without
intercept, failure accounting (`n_failed`), Kish effective sample sizes,
co-exceedance independence benchmark, runnable examples, an evaluated
getting-started vignette, CI workflows.

## Fixed defects

Listed with severity in `TECHNICAL_AUDIT_REPORT.md` §3 (18 items, 3
critical, 7 major).

## Testing evidence

427 test expectations pass (0 failures); line coverage 93.7 %; local
`R CMD check --as-cran`: 0 errors, 0 warnings, 1 note (URLs not yet
reachable because the repository was private and not yet renamed). Details:
`TESTING_REPORT.md`.

## Validation findings

Deterministic equality with direct `quantreg` fits; correct limits;
Monte Carlo evidence of correct behaviour and of the main practical
limitation (noisy slopes at density-based default bandwidths; boundary
bias for wide bandwidths). Details: `ECONOMETRIC_VALIDATION_REPORT.md`.

## Documentation

New README (method, installation, quick start, interpretation, migration
table), corrected and extended methodology vignette, new getting-started
vignette, all help pages regenerated, NEWS with migration notes and errata,
verified references with DOIs.

## Security and release hygiene

* Git history (2 commits, 49 objects) scanned: no credentials, tokens,
  private keys, large binaries or local paths; author e-mail is a GitHub
  no-reply address.
* Committed tree contains only package source, documentation, tests,
  validation scripts and reports. Research data (`BRVM`, `Marshallian`,
  `BAFFOUR TEST DATA`), result workbooks, presentations, photos, PDFs,
  `.Rhistory`/`.RData` and library-repair scripts from the local folder
  were **not** committed.
* The maintainer e-mail in `DESCRIPTION` (also present in the already
  published 0.2.0 files) is the authors' chosen public contact address.
* Licence: MIT; all dependencies are GPL-compatible CRAN packages used
  through their public APIs.

## Publication status

| Step | Status |
|---|---|
| Commits on `master` (local clone) | done (see `git log`) |
| Push to `Baffour536/qqr` | **pending user confirmation** |
| Rename repository to `qqrtools` (GitHub settings) | **pending – web UI** (`gh` CLI is not signed in) |
| Make repository public | **pending – after the rename and a green CI run** |
| GitHub Actions R-CMD-check (Windows/macOS/Linux) | pending first push |
| Tag `v0.3.0` and GitHub release | pending |
| `remotes::install_github("Baffour536/qqrtools")` in a clean library | pending publication |
| Final commit SHA / release URL | to be recorded |

## Known limitations

See `TECHNICAL_AUDIT_REPORT.md` §6 and `ECONOMETRIC_VALIDATION_REPORT.md`
§4–6. In particular, p-values are pointwise and descriptive, the automatic
bandwidths are not tuned to QQR slopes, and only the Gaussian kernel is
available.

## Recommended future work

1. A QQR-specific bandwidth selector (e.g. cross-validated check loss of
   the local linear fit across τ) and default bandwidth guidance.
2. Bootstrap confidence bands for the surface and multiplicity-aware
   inference.
3. Local quadratic option to reduce boundary bias of the slopes.
4. Split `R/qqr-methods.R` into plotting, dependence and print modules.
5. CRAN submission once the URLs resolve (the check is otherwise clean).
