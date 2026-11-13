## Local polynomial degree, slope-oriented bandwidth rules, bandwidth
## sensitivity and multiplicity-adjusted p-values.

test_that("degree = 2 equals a direct weighted local quadratic rq() fit", {
  dat  <- .ls_data(n = 250, seed = 31)
  taus <- c(0.2, 0.5, 0.8)
  h    <- 0.2
  fit  <- .qfit(y ~ x + w, data = dat, qvar = "x", taus = taus, m = 4,
                h = h, degree = 2, se_type = "nid")
  expect_identical(fit$degree, 2L)
  expect_named(fit$coef_surfaces, c("(Intercept)", "x", "w"))
  for (i in seq_along(fit$grid)) {
    g  <- fit$grid[i]
    d2 <- data.frame(y = dat$y, xc = dat$x - g, w = dat$w,
                     xc2 = (dat$x - g)^2,
                     wt = stats::dnorm((dat$x - g) / h))
    for (j in seq_along(taus)) {
      ref <- quantreg::rq(y ~ xc + w + xc2, data = d2, weights = wt,
                          tau = taus[j])
      cf  <- stats::coef(ref)
      expect_equal(fit$coef_surfaces[["(Intercept)"]][i, j], unname(cf[1]),
                   tolerance = 1e-8)
      expect_equal(fit$qcoef[i, j], unname(cf[2]), tolerance = 1e-8)
      expect_equal(fit$coef_surfaces[["w"]][i, j], unname(cf[3]),
                   tolerance = 1e-8)
      pv <- suppressWarnings(summary(ref, se = "nid"))$coefficients[, 4]
      expect_equal(fit$pval_surfaces[["x"]][i, j], unname(pv[2]),
                   tolerance = 1e-8)
    }
  }
})

test_that("degree = 2 recovers the slope of a noiseless quadratic exactly", {
  set.seed(32)
  x   <- stats::runif(150)
  dat <- data.frame(x = x, y = 1 + 2 * x^2)
  th  <- c(0.1, 0.5, 0.9)
  f2  <- .qfit(y ~ x, data = dat, taus = c(0.25, 0.75), thetas = th,
               h = 0.3, se_type = "none", degree = 2)
  truth <- 4 * f2$grid
  expect_equal(f2$qcoef[, 1], truth, tolerance = 1e-6)
  expect_equal(f2$qcoef[, 2], truth, tolerance = 1e-6)
  ## the local linear slope is biased at the edges of the support
  f1 <- .qfit(y ~ x, data = dat, taus = 0.5, thetas = th, h = 0.3,
              se_type = "none", degree = 1)
  expect_gt(abs(f1$qcoef[1, 1] - truth[1]), 0.1)
  expect_gt(abs(f1$qcoef[3, 1] - truth[3]), 0.1)
})

test_that("degree is validated and printed", {
  dat <- .ls_data(n = 80, seed = 33)
  expect_error(qqr(y ~ x, data = dat, degree = 3), "degree")
  expect_error(qqr(y ~ x, data = dat, degree = c(1, 2)), "degree")
  expect_error(qqr(y ~ x, data = dat, degree = NA), "degree")
  f1 <- qqr(y ~ x, data = dat, taus = 0.5, m = 3, se_type = "none",
            degree = 1)
  f2 <- qqr(y ~ x, data = dat, taus = 0.5, m = 3, se_type = "none",
            degree = 2)
  expect_output(print(f1), "linear \\(degree 1\\)")
  expect_output(print(f2), "quadratic \\(degree 2\\)")
  expect_output(print(summary(f2)), "degree 2")
})

test_that("slope-oriented bandwidth rules match their definitions", {
  set.seed(34)
  n <- 180
  x <- stats::rexp(n)
  y <- x + stats::rnorm(n)
  sel <- qqrtools:::.qqr_select_h
  C   <- qqrtools:::.qqr_slope_constant()
  expect_equal(sel(x, NULL, "slope")$h,
               C * min(stats::sd(x), stats::IQR(x) / 1.34) * n^(-1 / 7))
  expect_equal(sel(x, y, "cv_slope")$h,
               sel(x, y, "cv_fast")$h * n^(1 / 5 - 1 / 7))
  expect_error(sel(x, NULL, "cv_slope"), "requires the outcome")
  ## wider than the density rule at the same scale
  expect_gt(sel(x, NULL, "slope")$h, sel(x, NULL, "robust")$h)
  expect_match(sel(x, NULL, "slope")$details, "n\\^\\(-1/7\\)")
  expect_match(qqrtools:::.qqr_h_label("slope"), "Slope")
  expect_match(qqrtools:::.qqr_h_label("cv_slope"), "slope rate")

  dat <- data.frame(x = x, y = y)
  for (meth in c("slope", "cv_slope")) {
    fit <- qqr(y ~ x, data = dat, taus = 0.5, m = 3, h = meth,
               se_type = "none")
    expect_identical(fit$h_method, meth)
    expect_equal(fit$bandwidth, sel(x, y, meth)$h)
  }
  ## zero IQR falls back to the standard deviation
  xz <- c(rep(0, 35), 1:5)
  expect_equal(sel(xz, NULL, "slope")$h, C * stats::sd(xz) * 40^(-1 / 7))
})

test_that("qqr_sensitivity refits the surface over bandwidths", {
  dat <- .ls_data(n = 200, seed = 35)
  fit <- qqr(y ~ x + w, data = dat, qvar = "x", taus = c(0.25, 0.5, 0.75),
             m = 4, se_type = "none")
  s <- qqr_sensitivity(fit)
  expect_s3_class(s, "qqr_sensitivity")
  expect_equal(s$h, c(0.5, 0.75, 1, 1.5, 2) * fit$bandwidth)
  expect_equal(nrow(s$summary), 5L)
  ## factor 1 reproduces the fit
  expect_equal(s$surfaces[[3]], fit$qcoef, tolerance = 1e-10)
  expect_equal(s$summary$cor_baseline[3], 1)
  expect_equal(s$summary$sign_agreement[3], 1)
  expect_equal(s$summary$mean_abs_diff[3], 0)
  expect_equal(unname(s$tau_means[3, ]), unname(colMeans(fit$qcoef)))
  ## another bandwidth equals a manual refit
  ref <- .qfit(y ~ x + w, data = dat, qvar = "x", taus = c(0.25, 0.5, 0.75),
               m = 4, se_type = "none", h = 2 * fit$bandwidth,
               degree = fit$degree)
  expect_equal(s$surfaces[[5]], ref$qcoef, tolerance = 1e-10)
  expect_equal(s$summary$mean[5], mean(ref$qcoef))

  s2 <- qqr_sensitivity(fit, h = c(0.4, 0.2), coef = "w")
  expect_equal(s2$h, c(0.2, 0.4))
  expect_identical(s2$coef, "w")
  expect_output(print(s), "bandwidth sensitivity")
  .with_null_device(expect_invisible(plot(s)))

  expect_error(qqr_sensitivity(list()), "'qqr' object")
  expect_error(qqr_sensitivity(fit, factors = c(-1, 1)), "positive")
  expect_error(qqr_sensitivity(fit, h = 0), "positive")
  expect_error(qqr_sensitivity(fit, coef = "nope"), "not recognised")
  old <- fit
  old$design <- NULL
  expect_error(qqr_sensitivity(old), "re-fit")
})

test_that("qqr_sensitivity keeps the kernel space and degree of the fit", {
  dat <- .ls_data(n = 150, seed = 36)
  fit <- .qfit(y ~ x, data = dat, taus = c(0.3, 0.7), m = 3, h = 0.2,
               kernel_space = "cdf", degree = 1, se_type = "none")
  s   <- qqr_sensitivity(fit, factors = c(1, 1.5))
  ref <- .qfit(y ~ x, data = dat, taus = c(0.3, 0.7), m = 3, h = 0.3,
               kernel_space = "cdf", degree = 1, se_type = "none")
  expect_equal(s$surfaces[[1]], fit$qcoef, tolerance = 1e-10)
  expect_equal(s$surfaces[[2]], ref$qcoef, tolerance = 1e-10)
})

test_that("qqr_padjust adjusts each p-value surface with p.adjust", {
  dat <- .ls_data(n = 200, seed = 37)
  fit <- qqr(y ~ x + w, data = dat, qvar = "x", taus = c(0.25, 0.5, 0.75),
             m = 4)
  adj <- qqr_padjust(fit)
  expect_identical(adj$p_adjust, "BH")
  for (nm in names(fit$pval_surfaces)) {
    expect_equal(as.numeric(adj$pval_surfaces[[nm]]),
                 stats::p.adjust(as.numeric(fit$pval_surfaces[[nm]]), "BH"))
    expect_identical(dim(adj$pval_surfaces[[nm]]), c(4L, 3L))
    expect_true(all(adj$pval_surfaces[[nm]] >= fit$pval_surfaces[[nm]]))
  }
  expect_identical(adj$coef_surfaces, fit$coef_surfaces)
  by <- qqr_padjust(fit, "BY")
  expect_true(all(by$pval_surfaces$x >= adj$pval_surfaces$x))
  expect_output(print(adj), "P-value adjustment")
  expect_output(print(summary(adj)), "P-value adjustment")
  expect_identical(qqr_padjust(fit, "none"), fit)
  ## the adjusted surface is what the extraction helpers return
  expect_equal(qqr_surface(adj, coef = "x", what = "pvalue",
                           include_grid_value = FALSE)[, -1],
               as.data.frame(adj$pval_surfaces$x), ignore_attr = TRUE)

  expect_error(qqr_padjust(adj), "already been adjusted")
  expect_error(qqr_padjust(fit, "nope"), "must be one of")
  expect_error(qqr_padjust(list()), "'qqr' object")
  nop <- qqr(y ~ x, data = dat, taus = 0.5, m = 3, se_type = "none")
  expect_error(qqr_padjust(nop), "no p-values")
})
