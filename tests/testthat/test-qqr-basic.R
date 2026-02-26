## tests/testthat/test-qqr-basic.R
##
## Updated tests to fix the two concrete failure modes you hit:
##
## (A) testthat::expect_* helpers in your installed testthat do NOT accept `info=`.
##     Replace `info=` with the supported `label=` argument (or omit entirely).
##
## (B) compare.qqr() reconstructs formula/data from x$call by evaluating
##     in parent.frame(). In testthat, this can resolve `data=` to an
##     unexpected object (for example stats::dnorm or stats::df) if the
##     symbol is not bound in the evaluation environment.
##
##     The robust fix in tests is:
##       - always bind the dataset to a simple name like `dat`
##       - call qqr(..., data = dat) using that symbol (not a literal df object)
##       - call compare.qqr() from the same test scope where `dat` exists
##     That ensures compare.qqr() can eval(cl$data) safely.
##
## If you want compare.qqr() to work even when data was passed as an
## expression, you should also store the model frame or data inside the
## qqr object at fit time; but the test fix below matches your current core.

library(testthat)
library(qqr)

## ----------------------------------------------------------------------
## Helper: deterministic DGP
## ----------------------------------------------------------------------
.make_dgp <- function(n = 200, with_controls = TRUE, seed = 123) {
  set.seed(seed)
  x <- rnorm(n)
  if (with_controls) {
    z1 <- rnorm(n)
    z2 <- rnorm(n)
    y  <- 0.5 + 0.8 * x + 0.3 * z1 - 0.2 * z2 + rnorm(n)
    data.frame(y = y, x = x, z1 = z1, z2 = z2)
  } else {
    y  <- 0.5 + 0.8 * x + rnorm(n)
    data.frame(y = y, x = x)
  }
}

## ----------------------------------------------------------------------
## 1. Core structure and dimensions of qqr() fit (user-facing naming)
## ----------------------------------------------------------------------

test_that("qqr() basic fit has correct structure and dimensions (user-facing naming)", {
  dat <- .make_dgp(n = 200, with_controls = TRUE, seed = 123)

  taus <- seq(0.1, 0.9, by = 0.2)
  m    <- 7

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,     # IMPORTANT: symbol `dat` is bound in this scope
      qvar  = "x",
      taus  = taus,
      m     = m,
      h     = 0.4
    )
  )

  expect_s3_class(fit, "qqr")

  ## Metadata
  expect_identical(fit$dependent,   "y")
  expect_identical(fit$independent, "x")
  expect_setequal(fit$controls, c("z1", "z2"))
  expect_equal(fit$n, nrow(dat))

  ## Grid and taus
  expect_length(fit$grid, m)
  expect_equal(fit$taus, taus)

  ## Naming policy: user-facing names only
  expect_true("(Intercept)" %in% fit$coef_names)
  expect_true("x"           %in% fit$coef_names)
  expect_true("z1"          %in% fit$coef_names)
  expect_true("z2"          %in% fit$coef_names)
  expect_false("z"          %in% fit$coef_names)

  ## Surfaces exist, are named correctly, and have correct dimensions
  expect_true(is.list(fit$coef_surfaces))
  expect_true(is.list(fit$pval_surfaces))

  expect_setequal(names(fit$coef_surfaces), fit$coef_names)
  expect_setequal(names(fit$pval_surfaces), fit$coef_names)

  for (nm in fit$coef_names) {
    mat_coef <- fit$coef_surfaces[[nm]]
    mat_pval <- fit$pval_surfaces[[nm]]

    expect_true(is.matrix(mat_coef), label = paste("coef surface is matrix:", nm))
    expect_true(is.matrix(mat_pval), label = paste("pval surface is matrix:", nm))

    expect_equal(dim(mat_coef), c(m, length(taus)), label = paste("coef dim:", nm))
    expect_equal(dim(mat_pval), c(m, length(taus)), label = paste("pval dim:", nm))

    ## Not all NA
    expect_lt(mean(is.na(mat_coef)), 1, label = paste("coef not all NA:", nm))
    expect_lt(mean(is.na(mat_pval)), 1, label = paste("pval not all NA:", nm))
  }

  ## qcoef must equal the quantile regressor surface under user-facing key
  expect_true(is.matrix(fit$qcoef))
  expect_equal(dim(fit$qcoef), c(m, length(taus)))
  expect_equal(fit$qcoef, fit$coef_surfaces[["x"]])
})

## ----------------------------------------------------------------------
## 2. Coefficient name mapping (.qqr_match_coef)
## ----------------------------------------------------------------------

test_that(".qqr_match_coef maps user-facing names correctly and rejects internal aliases", {
  dat <- .make_dgp(n = 120, with_controls = TRUE, seed = 222)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 5,
      h     = 0.5
    )
  )

  matcher <- getFromNamespace(".qqr_match_coef", "qqr")

  ## Intercept aliases accepted
  expect_identical(matcher(fit, "(Intercept)"), "(Intercept)")
  expect_identical(matcher(fit, "intercept"),   "(Intercept)")
  expect_identical(matcher(fit, "(intercept)"), "(Intercept)")

  ## Quantile regressor by its original name
  expect_identical(matcher(fit, "x"), "x")

  ## Controls
  expect_identical(matcher(fit, "z1"), "z1")
  expect_identical(matcher(fit, "z2"), "z2")

  ## Internal alias must be rejected
  expect_error(matcher(fit, "z"), regexp = "not recognised")

  ## Unknown coefficient rejected
  expect_error(matcher(fit, "bogus_name"), regexp = "not recognised")
})

## ----------------------------------------------------------------------
## 3. summary.qqr behaviour
## ----------------------------------------------------------------------

test_that("summary.qqr produces coherent summaries aligned with qcoef", {
  dat <- .make_dgp(n = 150, with_controls = TRUE, seed = 321)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = seq(0.1, 0.9, by = 0.2),
      m     = 6,
      h     = 0.4
    )
  )

  s <- summary(fit)

  expect_s3_class(s, "summary.qqr")
  expect_identical(s$dependent,   fit$dependent)
  expect_identical(s$independent, fit$independent)
  expect_setequal(s$controls,     fit$controls)
  expect_true(is.logical(s$intercept))
  expect_equal(s$taus, fit$taus)

  expect_true(is.numeric(s$grid_range))
  expect_length(s$grid_range, 2)
  expect_lte(s$grid_range[1], s$grid_range[2])

  expect_true(is.numeric(s$beta_min))
  expect_true(is.numeric(s$beta_max))
  expect_true(is.numeric(s$beta_mean))
  expect_true(is.numeric(s$beta_sd))
  expect_lte(s$beta_min, s$beta_max)

  expect_equal(length(s$row_means), nrow(fit$qcoef))
  expect_equal(length(s$col_means), ncol(fit$qcoef))

  ## Bandwidth fields added in updated summary.qqr
  expect_true("bandwidth" %in% names(s))
  expect_true("h_method"  %in% names(s))
  expect_true("h_details" %in% names(s))
  expect_true(is.numeric(s$bandwidth))
  expect_gt(s$bandwidth, 0)
  expect_true(s$h_method %in% c("silverman", "robust", "iqr", "cv", "manual"))
})

## ----------------------------------------------------------------------
## 4. Basic error handling and argument validation
## ----------------------------------------------------------------------

test_that("qqr() errors on invalid inputs and qvar", {
  set.seed(111)

  n <- 50
  x <- rnorm(n)
  z <- rnorm(n)
  y <- 0.3 + 0.7 * x + rnorm(n)

  dat <- data.frame(y = y, x = x, z = z)

  ## Non-numeric dependent
  dat_bad_y <- dat
  dat_bad_y$y <- factor(dat_bad_y$y > 0)
  expect_error(qqr(y ~ x, data = dat_bad_y, qvar = "x"), regexp = "must be numeric")

  ## qvar not in formula
  expect_error(suppressWarnings(qqr(y ~ x + z, data = dat, qvar = "wrong")), regexp = "not a variable in the formula|not an independent variable")

  ## Multiple regressors but qvar omitted
  expect_error(suppressWarnings(qqr(y ~ x + z, data = dat)), regexp = "qvar")

  ## Only dependent variable
  expect_error(qqr(y ~ 1, data = dat), regexp = "at least one independent variable")
})

## ----------------------------------------------------------------------
## 5. qqr_surface(): wide vs long, subsets, rounding
## ----------------------------------------------------------------------

test_that("qqr_surface() returns correctly shaped wide and long tables", {
  dat <- .make_dgp(n = 120, with_controls = TRUE, seed = 456)

  taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  m    <- 5

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = taus,
      m     = m,
      h     = 0.5
    )
  )

  surf_x <- qqr_surface(fit, coef = "x", what = "coef", long = FALSE, digits = 4)

  expect_true(is.data.frame(surf_x))
  expect_equal(nrow(surf_x), m)
  expect_true(all(c("Grid_Point", "Grid_Value") %in% names(surf_x)))

  tau_cols <- paste0("tau_", taus)
  expect_true(all(tau_cols %in% names(surf_x)))

  skip_if_not_installed("tidyr")

  surf_all <- qqr_surface(fit, coef = "all", what = "pvalue", long = TRUE, digits = 4)

  expect_true(is.data.frame(surf_all))
  expect_true(all(c("Coefficient", "Grid_Point", "tau", "value") %in% names(surf_all)))

  expect_true(all(surf_all$Coefficient %in% c("Intercept", "x", "z1", "z2")))
  expect_false(any(surf_all$Coefficient %in% "z"))
  expect_true(all(surf_all$tau %in% taus))

  surf_sub <- qqr_surface(
    fit,
    coef        = "x",
    what        = "coef",
    long        = FALSE,
    grid_subset = c(2, 4),
    tau_subset  = c(0.3, 0.7),
    digits      = 4
  )

  expect_equal(nrow(surf_sub), 2)
  expect_true(all(c("Grid_Point", "Grid_Value") %in% names(surf_sub)))
  expect_true(all(paste0("tau_", c(0.3, 0.7)) %in% names(surf_sub)))
})

## ----------------------------------------------------------------------
## 6. rowave.qqr(): row means across tau
## ----------------------------------------------------------------------

test_that("rowave.qqr() returns sensible row averages for all and single coefficients", {
  dat <- .make_dgp(n = 100, with_controls = TRUE, seed = 789)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 5,
      h     = 0.4
    )
  )

  ra_all <- rowave.qqr(fit, coef = "all", digits = 4)

  expect_true(is.data.frame(ra_all))
  expect_equal(nrow(ra_all), length(fit$grid))
  expect_true(all(c("grid_index", "grid_value") %in% names(ra_all)))
  expect_true(all(c("Intercept", "x", "z1", "z2") %in% names(ra_all)))
  expect_false(any(names(ra_all) %in% "z"))

  ra_x <- rowave.qqr(fit, coef = "x", digits = 4)
  expect_true(is.data.frame(ra_x))
  expect_true(all(c("grid_index", "grid_value") %in% names(ra_x)))
  expect_true("mean_x" %in% names(ra_x))
})

## ----------------------------------------------------------------------
## 7. Static plotting: plot.qqr()
## ----------------------------------------------------------------------

test_that("plot.qqr() static surfaces run without error", {
  dat <- .make_dgp(n = 80, with_controls = FALSE, seed = 135)

  fit <- suppressWarnings(
    qqr(
      y ~ x,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 5,
      h     = 0.5
    )
  )

  expect_silent(plot(fit, coef = "x", what = "coef", engine = "base"))
  expect_silent(plot(fit, coef = "x", what = "pvalue", engine = "base"))

  skip_if_not_installed("plot3D")
  expect_silent(plot(fit, coef = "x", what = "coef", engine = "plot3D"))
})

## ----------------------------------------------------------------------
## 8. Interactive plotting: Int.plot()
## ----------------------------------------------------------------------

test_that("Int.plot() runs when rgl is available", {
  skip_on_cran()
  skip_if_not_installed("rgl")

  dat <- .make_dgp(n = 60, with_controls = FALSE, seed = 246)

  fit <- suppressWarnings(
    qqr(
      y ~ x,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 4,
      h     = 0.5
    )
  )

  expect_invisible(Int.plot(fit, coef = "x", what = "coef", add_plane = TRUE, plane_z = 0))
})

## ----------------------------------------------------------------------
## 9. Excel export: excel.qqr()
## ----------------------------------------------------------------------

test_that("excel.qqr() writes a valid Excel file with coefficients and p-values", {
  skip_if_not_installed("openxlsx")

  dat <- .make_dgp(n = 80, with_controls = TRUE, seed = 975)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 5,
      h     = 0.4
    )
  )

  tmp <- tempfile(fileext = ".xlsx")

  expect_invisible(
    excel.qqr(
      fit,
      file               = tmp,
      which              = "both",
      digits             = 4,
      include_grid_value = TRUE
    )
  )

  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)

  wb <- openxlsx::loadWorkbook(tmp)
  sheets <- names(wb)

  expect_true(any(grepl("Intercept_coef", sheets)))
  expect_true(any(grepl("^x_coef", sheets)))
  expect_true(any(grepl("^z1_coef", sheets)) || any(grepl("^z2_coef", sheets)))
  expect_false(any(grepl("^z_coef", sheets)))
})

## ----------------------------------------------------------------------
## 10. Word export: word.qqr()
## ----------------------------------------------------------------------

test_that("word.qqr() writes a valid Word file with coefficient and p-value tables", {
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")

  dat <- .make_dgp(n = 80, with_controls = TRUE, seed = 864)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 5,
      h     = 0.4
    )
  )

  tmp <- tempfile(fileext = ".docx")

  expect_invisible(
    word.qqr(
      fit,
      file        = tmp,
      coef        = "x",
      which       = "both",
      grid_subset = NULL,
      digits      = 4
    )
  )

  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})

## ----------------------------------------------------------------------
## 11. se_type argument behaviour
## ----------------------------------------------------------------------

test_that("qqr() accepts se_type argument and stores it consistently", {
  dat <- .make_dgp(n = 80, with_controls = FALSE, seed = 4321)

  taus <- c(0.25, 0.5, 0.75)
  m    <- 4

  fit_default <- suppressWarnings(
    qqr(
      y ~ x,
      data  = dat,
      qvar  = "x",
      taus  = taus,
      m     = m,
      h     = 0.5
    )
  )
  expect_s3_class(fit_default, "qqr")
  expect_true("se_type" %in% names(fit_default))

  if ("se_type" %in% names(formals(qqr))) {
    fit_nid <- suppressWarnings(
      qqr(
        y ~ x,
        data    = dat,
        qvar    = "x",
        taus    = taus,
        m       = m,
        h       = 0.5,
        se_type = "nid"
      )
    )
    expect_s3_class(fit_nid, "qqr")
    expect_identical(fit_nid$se_type, "nid")
  }
})

## ----------------------------------------------------------------------
## 12. compare.qqr() and related helpers
## ----------------------------------------------------------------------

test_that("compare.qqr() compares QR coefficients with QQR row means (robust data binding)", {
  skip_if_not_installed("quantreg")

  ## CRITICAL: bind data to a stable symbol `dat`, and pass data = dat
  dat <- .make_dgp(n = 150, with_controls = TRUE, seed = 2025)

  taus <- seq(0.1, 0.9, by = 0.2)
  m    <- length(taus)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,     # ensures compare.qqr can eval(cl$data) -> dat
      qvar  = "x",
      taus  = taus,
      m     = m,
      h     = 0.4
    )
  )

  ## Single coefficient: quantile regressor
  dat <- .make_dgp(n = 150, with_controls = TRUE, seed = 2025)
  cmp_x <- compare.qqr(fit, coef = "x", digits = 6)
  expect_true(is.data.frame(cmp_x))
  expect_true(all(c("tau", "qr", "qqr_mean") %in% names(cmp_x)))
  expect_equal(cmp_x$tau, taus)
  expect_equal(nrow(cmp_x), m)
  expect_true(all(is.finite(cmp_x$qr)))
  expect_lt(mean(is.na(cmp_x$qqr_mean)), 1)

  ## Control coefficient
  cmp_z1 <- compare.qqr(fit, coef = "z1", digits = 6)
  expect_true(is.data.frame(cmp_z1))
  expect_equal(cmp_z1$tau, taus)

  ## Multiple coefficients returns a named list
  dat <- .make_dgp(n = 150, with_controls = TRUE, seed = 2025)
  cmp_many <- compare.qqr(fit, coef = c("x", "z1"), digits = 6)
  expect_true(is.list(cmp_many))
  expect_true(all(c("x", "z1") %in% names(cmp_many)))
  expect_true(is.data.frame(cmp_many[["x"]]))
  expect_true(is.data.frame(cmp_many[["z1"]]))

  ## Internal alias must not work
  expect_error(compare.qqr(fit, coef = "k"), regexp = "not found|not recognised|must be")
})

test_that("excel.compare.qqr() writes a valid Excel file (robust data binding)", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("openxlsx")

  dat <- .make_dgp(n = 120, with_controls = TRUE, seed = 3030)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 3,
      h     = 0.5
    )
  )

  tmp <- tempfile(fileext = ".xlsx")

  expect_invisible(
    excel.compare.qqr(
      fit,
      file   = tmp,
      coef   = c("x", "z1"),
      digits = 4
    )
  )

  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)

  wb <- openxlsx::loadWorkbook(tmp)
  sheets <- names(wb)
  expect_true(length(sheets) >= 2)
})

test_that("word.compare.qqr() writes a valid Word file (robust data binding)", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")

  dat <- .make_dgp(n = 120, with_controls = TRUE, seed = 4040)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.25, 0.5, 0.75),
      m     = 3,
      h     = 0.5
    )
  )

  tmp <- tempfile(fileext = ".docx")

  expect_invisible(
    word.compare.qqr(
      fit,
      file   = tmp,
      coef   = c("x", "z1"),
      digits = 4
    )
  )

  expect_true(file.exists(tmp))
  expect_gt(file.info(tmp)$size, 0)
})

test_that("plot.compare.qqr() returns a ggplot object and rejects multiple coefs (robust data binding)", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("ggplot2")

  dat <- .make_dgp(n = 150, with_controls = TRUE, seed = 5050)

  fit <- suppressWarnings(
    qqr(
      y ~ x + z1 + z2,
      data  = dat,
      qvar  = "x",
      taus  = c(0.1, 0.5, 0.9),
      m     = 3,
      h     = 0.4
    )
  )

  p <- plot.compare.qqr(fit, coef = "x", digits = 6)
  expect_true(inherits(p, "ggplot"))

  expect_error(
    plot.compare.qqr(fit, coef = c("x", "z1")),
    regexp = "single coefficient"
  )
})
