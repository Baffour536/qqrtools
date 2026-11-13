#' qqrtools: Tools for Quantile-on-Quantile Regression
#'
#' The **qqrtools** package provides tools for estimating and analysing
#' *quantile–on–quantile regression* (QQR) surfaces. QQR extends standard
#' quantile regression by modelling the response at multiple conditional
#' quantile levels of a *quantile regressor*, evaluated over a grid of
#' values of that regressor. This produces a surface of coefficients
#' indexed by quantile of the regressor (grid) and conditional quantile
#' level (tau).
#'
#' This package implements:
#'   * Local, kernel-weighted quantile regressions along a grid of the
#'     quantile regressor.
#'   * Storage of estimated coefficient and p-value surfaces for the
#'     intercept, the quantile regressor, and control variables.
#'   * Consistent user-facing naming of coefficients (no internal
#'     aliasing).
#'   * Tools to compare QQR surfaces with standard quantile regression
#'     (QR), including row means of QQR surfaces.
#'   * Flexible extraction and export utilities (wide/long surfaces,
#'     Excel, Word).
#'   * Static and interactive 3D visualisations of surfaces.
#'
#' ## Key Features
#'
#' ### Model Estimation
#' - **`qqr()`**: Fit a QQR model given a formula, data frame, and
#'   user-specified quantile regressor (`qvar`). The function returns
#'   an object of class `"qqr"` containing coefficient and p-value
#'   surfaces, grid values, quantile indices (`taus`), and the data
#'   used in the fit.
#' - Data used in the original fit is stored in the object (`fit$data`)
#'   to enable robust post-estimation comparisons and exports, avoiding
#'   scoping issues when re-evaluating calls.
#'
#' ### Naming Conventions
#' - Coefficients are named using **user-facing variable names**:
#'   - `"(Intercept)"` for the model intercept (if present)
#'   - The original `qvar` name for the quantile regressor
#'   - Control variables as given in the formula
#' - The internal local regressor (“z = x − z0”) is **never exposed** in
#'   public surfaces or comparisons.
#'
#' ### Comparison with Standard Quantile Regression
#' - **`compare.qqr()`**: Compares the row means of a QQR coefficient
#'   surface with coefficients from standard quantile regression
#'   (`quantreg::rq`). For each coefficient of interest, it returns
#'   a data frame with columns:
#'     * `tau`: quantile level
#'     * `qr`: coefficient from standard QR
#'     * `qqr_mean`: row-mean of the QQR surface
#' - Invalid coefficient names (e.g., internal aliases like `"z"`) are
#'   rejected with informative errors.
#' - Multi-coefficient comparisons (e.g., `coef = c("x", "z1")`) return
#'   a named list of comparison data frames.
#'
#' ### Surface Extraction
#' - **`qqr_surface()`**: Extract wide or long format tables for
#'   coefficient and p-value surfaces.
#' - Allows subsetting by grid points and quantile levels (`taus`).
#' - Supports optional rounding for display and export.
#'
#' - **`rowave.qqr()`**: Computes row means across quantile levels for
#'   a specified coefficient or all coefficients, with tidy labelling.
#'
#' ### Export Utilities
#' - **`excel.qqr()`**: Writes coefficient or p-value surfaces to an
#'   Excel workbook (one sheet per coefficient), with grid values and
#'   optional rounding.
#' - **`word.qqr()`**: Writes the same to a Word document with formatted
#'   tables (requires `officer` and `flextable`).
#'
#' ### Static and Interactive Plotting
#' - **`plot.qqr()`**: Static 3D visualisation of coefficient or p-value
#'   surfaces using multiple engines:
#'     * `"base"` — built-in R graphics
#'     * `"plot3D"` — enhanced 3D visualisation if installed
#' - **`Int.plot()`**: Interactive 3D plots (requires `rgl`).
#'
#' ## Typical Workflow
#' 1. Fit a QQR model:
#'    ```r
#'    fit <- qqr(
#'      y ~ qvar + controls,
#'      data = df,
#'      qvar = "qvar",
#'      taus = seq(0.05, 0.95, by = 0.05),
#'      m    = 19,
#'      h    = 0.05
#'    )
#'    ```
#' 2. View scalar summaries:
#'    ```r
#'    summary(fit)
#'    ```
#' 3. Compare with standard quantile regression:
#'    ```r
#'    cmp <- compare.qqr(fit, coef = "qvar")
#'    ```
#' 4. Extract surfaces for reporting or export:
#'    ```r
#'    surf <- qqr_surface(fit, coef = "qvar", long = TRUE)
#'    excel.qqr(fit, file = "qqr_surfaces.xlsx")
#'    word.qqr(fit, file = "qqr_surfaces.docx")
#'    ```
#' 5. Visualise:
#'    ```r
#'    plot(fit, coef = "qvar", what = "coef", engine = "plot3D")
#'    Int.plot(fit, coef = "qvar", what = "coef")
#'    ```
#'
#' ## Dependencies
#' - `quantreg` for quantile regression fitting
#' - `openxlsx` for Excel export
#' - `officer` and `flextable` for Word export
#' - `plot3D` and `rgl` (optional) for enhanced plotting
#'
#' @docType package
#' @name qqrtools
#' @importFrom stats setNames
#' @importFrom graphics par layout persp polygon text axis box image
#' @importFrom grDevices hcl.colors adjustcolor
"_PACKAGE"
