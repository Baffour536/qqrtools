## Print/summary methods, extraction helpers, exports and bandwidth edge
## cases.

.fit_mx <- function() {
  dat <- .ls_data(n = 150, seed = 41)
  qqr(y ~ x + w, data = dat, qvar = "x", taus = c(0.25, 0.5, 0.75), m = 4)
}

test_that("print and summary report the estimation settings", {
  fit <- .fit_mx()
  expect_output(print(fit), "empirical quantiles")
  expect_output(print(fit), "Slope-oriented")
  s <- summary(fit)
  expect_s3_class(s, "summary.qqr")
  expect_equal(s$col_means, colMeans(fit$qcoef))
  expect_equal(s$row_means, rowMeans(fit$qcoef))
  expect_equal(s$thetas, fit$thetas)
  out <- capture.output(print(s))
  expect_true(any(grepl("Selection method", out)))
  expect_true(any(grepl("quantile / x", out)))
  expect_true(any(grepl("Slope range", out)))
  fm <- .qfit(y ~ x, .ls_data(80, 2), taus = 0.5, m = 3, h = 0.3,
              grid = "range", kernel_space = "cdf")
  expect_output(print(fm), "equally spaced")
  expect_output(print(fm), "F_n")
  expect_output(print(summary(fm)), "Manual h")
  expect_identical(qqrtools:::.qqr_h_label(NULL), "unknown")
  expect_identical(qqrtools:::.qqr_h_label("other"), "other")
})

test_that("qqr_surface matches tau values despite floating-point noise", {
  dat <- .ls_data(n = 120, seed = 42)
  taus <- seq(0.05, 0.95, by = 0.05)[c(2, 6, 10)]   # 0.1, 0.3, 0.5
  fit <- qqr(y ~ x, dat, taus = taus, m = 3, se_type = "none")
  s1 <- qqr_surface(fit, coef = "x", tau_subset = c(0.1, 0.5))
  expect_identical(names(s1), c("Grid_Point", "Grid_Value", "tau_0.1",
                                "tau_0.5"))
  expect_equal(unname(as.matrix(s1[, 3:4])), unname(fit$qcoef[, c(1, 3)]))
  s2 <- qqr_surface(fit, coef = "x", tau_subset = 2:3, grid_subset = c(1, NA),
                    include_grid_value = FALSE)
  expect_equal(s2$Grid_Point, c(1, 3))
  expect_error(qqr_surface(fit, tau_subset = 0.42), "not in the fit")
  expect_error(qqr_surface(fit, tau_subset = "a"), "numeric")
  lng <- qqr_surface(fit, coef = c("x", "intercept"), what = "pvalue",
                     digits = 2)
  expect_equal(nrow(lng), 2 * 3 * 3)
  expect_identical(unique(lng$Coefficient), c("x", "Intercept"))
  ## long layout: coefficient, then grid point, then tau
  x_rows <- lng[lng$Coefficient == "x", ]
  expect_equal(x_rows$value, round(as.vector(t(fit$pval_surfaces$x)), 2))
  expect_error(qqr_surface(list()), "class 'qqr'")
})

test_that("rowave.qqr averages each grid point over tau", {
  fit <- .fit_mx()
  r <- rowave.qqr(fit, coef = "w", digits = 4)
  expect_equal(r$mean_w, round(rowMeans(fit$coef_surfaces$w), 4))
  ra <- rowave.qqr(fit, include_grid_value = FALSE)
  expect_identical(names(ra), c("grid_index", "Intercept", "x", "w"))
  r1 <- rowave.qqr(fit, coef = "x", include_grid_value = FALSE)
  expect_identical(names(r1), c("grid_index", "mean_x"))
  expect_error(rowave.qqr(list()), "class 'qqr'")
})

test_that("excel.qqr and excel.compare.qqr write the expected sheets", {
  skip_if_not_installed("openxlsx")
  fit <- .fit_mx()
  f <- tempfile(fileext = ".xlsx")
  excel.qqr(fit, file = f, which = "pvalue", include_grid_value = FALSE,
            digits = 3)
  expect_identical(openxlsx::getSheetNames(f),
                   c("Intercept_pval", "x_pval", "w_pval"))
  sh <- openxlsx::read.xlsx(f, sheet = "x_pval")
  expect_identical(names(sh)[1], "Grid_Point")
  expect_equal(unname(as.matrix(sh[, -1])),
               unname(round(fit$pval_surfaces$x, 3)), tolerance = 1e-8)
  f2 <- tempfile(fileext = ".xlsx")
  excel.qqr(fit, file = f2, which = "coef")
  expect_identical(openxlsx::getSheetNames(f2),
                   c("Intercept_coef", "x_coef", "w_coef"))
  f3 <- tempfile(fileext = ".xlsx")
  excel.compare.qqr(fit, file = f3, coef = c("x", "w"))
  expect_identical(openxlsx::getSheetNames(f3), c("coef_x", "coef_w"))
  expect_error(excel.qqr(list()), "class 'qqr'")
  unlink(c(f, f2, f3))
})

test_that("word.qqr writes star tables, with [NA] where p-values are missing", {
  skip_if_not_installed("officer")
  skip_if_not_installed("flextable")
  fit <- .fit_mx()
  fit$pval_surfaces$x[1, 1] <- NA
  f <- tempfile(fileext = ".docx")
  word.qqr(fit, file = f, which = "coef_stars", digits = 2,
           grid_subset = c(1, 2, NA))
  txt <- officer::docx_summary(officer::read_docx(f))$text
  expect_true(any(grepl("[NA]", txt, fixed = TRUE)))
  expect_true(any(grepl("significance stars", txt)))
  expect_true(any(grepl(formatC(fit$qcoef[2, 2], format = "f", digits = 2),
                        txt, fixed = TRUE)))
  f2 <- tempfile(fileext = ".docx")
  word.qqr(fit, file = f2, coef = "w", which = "pvalue",
           include_grid_value = FALSE)
  txt2 <- officer::docx_summary(officer::read_docx(f2))$text
  expect_true(any(grepl("p-values for w", txt2)))
  f3 <- tempfile(fileext = ".docx")
  word.qqr(fit, file = f3, which = "coef_stars", include_grid_value = FALSE)
  expect_true(file.exists(f3))
  expect_error(word.qqr(list()), "class 'qqr'")
  unlink(c(f, f2, f3))
})

test_that("plot.compare.qqr returns a ggplot with both series", {
  skip_if_not_installed("ggplot2")
  fit <- .fit_mx()
  p <- plot.compare.qqr(fit, coef = "w", title = "t")
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$method), c("QR", "QQR mean"))
  expect_error(plot.compare.qqr(fit, coef = c("x", "w")), "single")
  expect_error(compare.qqr(fit, coef = "nope"), "not recognised")
  expect_error(compare.qqr(fit, qr_args = 1), "list")
  expect_error(compare.qqr(fit, coef = "z"), "not a model coefficient")
  expect_message(compare.qqr(.qfit(y ~ x, .ls_data(60, 3), taus = 0.5, m = 3,
                                   h = 0.3, grid = "range",
                                   se_type = "none")),
                 "equally spaced")
})

test_that("bandwidth selectors fail loudly or fall back with a warning", {
  sel <- qqrtools:::.qqr_select_h
  expect_error(sel(c(1, 2, 3), NULL, "silverman"), "Too few")
  expect_error(sel(rep(1, 10), NULL, "silverman"), "zero variance")
  x <- c(rep(0, 20), 1:5)
  expect_error(sel(x, NULL, "iqr"), "IQR")
  expect_gt(sel(x, NULL, "robust")$h, 0)    # IQR = 0 -> sd is used
  expect_error(sel(stats::rnorm(20), NULL, "yj"), "requires the outcome")
  expect_error(sel(stats::rnorm(8), stats::rnorm(8), "yj"), "Too few")
  expect_error(sel(stats::rnorm(20), NULL, "nope"), "Unknown bandwidth")
  ## dpill() failure on a degenerate design falls back to Silverman
  skip_if_not_installed("KernSmooth")
  xx <- c(rep(0, 30), rep(1, 30))
  yy <- stats::rnorm(60)
  res <- suppressWarnings(sel(xx, yy, "yj"))
  expect_gt(res$h, 0)
})

test_that("quantreg failures inside a local fit are caught", {
  X <- cbind(1, c(-1, 0, 1, 2, 3, 4))
  ## rq() errors on infinite responses; the error must become a failed cell
  res <- qqrtools:::.qqr_local_fit(c(1, 2, Inf, 4, 5, 6), X, rep(1, 6), 0.5,
                                   "nid")
  expect_null(res$coef)
  ## rank-deficient local design
  res2 <- qqrtools:::.qqr_local_fit(1:6, cbind(1, rep(2, 6)), rep(1, 6),
                                    0.5, "nid")
  expect_null(res2$coef)
  expect_error(qqrtools:::.qqr_kernel(0, "epanechnikov"), "Unsupported")
  expect_equal(qqrtools:::.qqr_kish(c(0, 0)), 0)
  expect_error(qqrtools:::.qqr_match_coef(.fit_mx(), c("x", "w")),
               "single character")
})

test_that("np-based bandwidths are recorded (not on CRAN)", {
  skip_on_cran()
  skip_if_not_installed("np")
  set.seed(43)
  x <- stats::runif(50); y <- x + stats::rnorm(50, sd = 0.2)
  for (m in c("npcv", "sc")) {
    r <- suppressWarnings(suppressMessages(
      utils::capture.output(res <- qqrtools:::.qqr_select_h(x, y, m))))
    expect_gt(res$h, 0)
    expect_identical(res$method, m)
  }
})
