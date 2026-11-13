## Independent reference calculations for the QQR estimator.
## Each test recomputes a quantity from first principles (or with an
## established implementation) and compares it with qqr().

test_that("each surface cell equals a directly computed weighted rq() fit", {
  dat <- .ls_data(n = 250, seed = 11)
  taus <- c(0.2, 0.5, 0.8)
  h <- 0.15
  fit <- .qfit(y ~ x + w, data = dat, qvar = "x", taus = taus, m = 4, h = h)

  expect_equal(fit$grid,
               unname(stats::quantile(dat$x, fit$thetas, type = 7)))

  for (i in seq_along(fit$grid)) {
    g  <- fit$grid[i]
    wt <- stats::dnorm((dat$x - g) / h)
    d2 <- data.frame(y = dat$y, xc = dat$x - g, w = dat$w, wt = wt)
    for (j in seq_along(taus)) {
      ref <- quantreg::rq(y ~ xc + w, data = d2, weights = wt, tau = taus[j])
      cf  <- stats::coef(ref)
      expect_equal(fit$coef_surfaces[["(Intercept)"]][i, j], unname(cf[1]),
                   tolerance = 1e-8)
      expect_equal(fit$coef_surfaces[["x"]][i, j], unname(cf[2]),
                   tolerance = 1e-8)
      expect_equal(fit$coef_surfaces[["w"]][i, j], unname(cf[3]),
                   tolerance = 1e-8)
      pv <- suppressWarnings(summary(ref, se = "nid"))$coefficients[, 4]
      expect_equal(fit$pval_surfaces[["x"]][i, j], unname(pv[2]),
                   tolerance = 1e-8)
    }
  }
})

test_that("kernel_space = 'cdf' uses Sim-Zhou weights K((F_n(X) - theta)/h)", {
  dat <- .ls_data(n = 200, seed = 12)
  h <- 0.1
  fit <- .qfit(y ~ x, data = dat, taus = c(0.3, 0.7), thetas = c(0.25, 0.75),
               h = h, kernel_space = "cdf", se_type = "none")
  Fn <- stats::ecdf(dat$x)(dat$x)
  for (i in 1:2) {
    th <- fit$thetas[i]
    g  <- unname(stats::quantile(dat$x, th))
    wt <- stats::dnorm((Fn - th) / h)
    d2 <- data.frame(y = dat$y, xc = dat$x - g, wt = wt)
    for (j in 1:2) {
      ref <- quantreg::rq(y ~ xc, data = d2, weights = wt, tau = fit$taus[j])
      expect_equal(fit$qcoef[i, j], unname(stats::coef(ref)[2]),
                   tolerance = 1e-8)
    }
  }
})

test_that("with a huge bandwidth, QQR collapses to standard QR", {
  dat  <- .ls_data(n = 300, seed = 13)
  taus <- c(0.25, 0.5, 0.75)
  fit  <- .qfit(y ~ x + w, data = dat, qvar = "x", taus = taus, m = 5,
                h = 1e6, se_type = "none")
  qr <- stats::coef(quantreg::rq(y ~ x + w, data = dat, tau = taus))
  for (i in seq_len(5)) {
    expect_equal(unname(fit$qcoef[i, ]), unname(qr["x", ]), tolerance = 1e-6)
    expect_equal(unname(fit$coef_surfaces$w[i, ]), unname(qr["w", ]),
                 tolerance = 1e-6)
    ## local intercept = QR fitted quantile at x = grid point, w = 0
    expect_equal(unname(fit$coef_surfaces[["(Intercept)"]][i, ]),
                 unname(qr["(Intercept)", ] + qr["x", ] * fit$grid[i]),
                 tolerance = 1e-6)
  }
  cmp <- suppressMessages(compare.qqr(fit))
  expect_equal(cmp$qqr_mean, cmp$qr, tolerance = 1e-6)
})

test_that("a noiseless linear model is recovered exactly at every cell", {
  set.seed(14)
  x <- stats::rnorm(120)
  w <- stats::rnorm(120)
  dat <- data.frame(y = 2 - 3 * x + 0.5 * w, x = x, w = w)
  fit <- qqr(y ~ x + w, data = dat, qvar = "x", taus = c(0.1, 0.5, 0.9),
             m = 6, se_type = "none")
  expect_true(all(abs(fit$qcoef + 3) < 1e-8))
  expect_true(all(abs(fit$coef_surfaces$w - 0.5) < 1e-8))
})

test_that("the location-scale slope 2 + 0.5 qnorm(tau) is recovered", {
  ## The local linear model is exact for this DGP, so the estimator is
  ## unbiased for any h; a moderately wide h keeps the variance small.
  ## (With the density-oriented Silverman rule, h is ~0.06 here and the
  ## local slopes are much noisier; see ECONOMETRIC_VALIDATION_REPORT.md.)
  dat  <- .ls_data(n = 4000, seed = 15)
  taus <- c(0.25, 0.5, 0.75)
  fit  <- .qfit(y ~ x + w, data = dat, qvar = "x", taus = taus,
                thetas = c(0.25, 0.5, 0.75), h = 0.3, se_type = "none")
  truth <- matrix(2 + 0.5 * stats::qnorm(taus), nrow = 3, ncol = 3,
                  byrow = TRUE)
  expect_lt(max(abs(fit$qcoef - truth)), 0.35)
  expect_lt(max(abs(colMeans(fit$qcoef) - truth[1, ])), 0.25)
  expect_lt(max(abs(fit$coef_surfaces$w - 0.5)), 0.15)
})

test_that("compare.qqr averages the surface over theta for each tau", {
  dat  <- .ls_data(n = 300, seed = 16)
  taus <- c(0.2, 0.5, 0.8)
  fit  <- qqr(y ~ x + w, data = dat, qvar = "x", taus = taus, m = 7,
              se_type = "none")
  cmp  <- compare.qqr(fit)
  qr   <- stats::coef(quantreg::rq(y ~ x + w, data = dat, tau = taus))
  expect_equal(cmp$tau, taus)
  expect_equal(cmp$qqr_mean, unname(colMeans(fit$qcoef)))
  expect_equal(cmp$qr, unname(qr["x", ]), tolerance = 1e-8)
  expect_equal(cmp$difference, cmp$qqr_mean - cmp$qr)

  cmp_w <- compare.qqr(fit, coef = "w")
  expect_equal(cmp_w$qr, unname(qr["w", ]), tolerance = 1e-8)
  expect_equal(cmp_w$qqr_mean, unname(colMeans(fit$coef_surfaces$w)))

  all_cmp <- compare.qqr(fit, coef = "all")
  expect_named(all_cmp, c("(Intercept)", "x", "w"))
})

test_that("bandwidth rules match their closed-form definitions", {
  set.seed(17)
  x <- stats::rnorm(150)
  y <- x + stats::rnorm(150)
  n <- length(x)
  sel <- function(m) qqrtools:::.qqr_select_h(x, y, m)$h
  expect_equal(sel("silverman"), 1.06 * stats::sd(x) * n^(-1 / 5))
  expect_equal(sel("robust"),
               0.9 * min(stats::sd(x), stats::IQR(x) / 1.34) * n^(-1 / 5))
  expect_equal(sel("iqr"), stats::IQR(x) * n^(-1 / 5))
  expect_equal(sel("sj"), stats::bw.SJ(x))
  skip_if_not_installed("KernSmooth")
  expect_equal(sel("yj"), KernSmooth::dpill(x, y))
})

test_that("cv_fast minimises the leave-one-out check loss of y", {
  set.seed(18)
  x <- stats::runif(80)
  y <- sin(2 * pi * x) + 0.2 * stats::rnorm(80)
  res <- qqrtools:::.qqr_select_h(x, y, "cv_fast")
  hs  <- seq(0.5, 2, length.out = 15) * 1.06 * stats::sd(x) * 80^(-1 / 5)
  ## brute-force reference: weighted median of y_{-i}
  wmed <- function(v, w) {
    o <- order(v); v <- v[o]; w <- w[o]
    v[which(cumsum(w) / sum(w) >= 0.5)[1]]
  }
  loss <- sapply(hs, function(h) mean(sapply(seq_along(x), function(i) {
    wi <- stats::dnorm((x[-i] - x[i]) / h)
    u  <- y[i] - wmed(y[-i], wi)
    u * (0.5 - (u < 0))
  })))
  expect_equal(res$h, hs[which.min(loss)])
  ## a wiggly signal should not pick the widest candidate
  expect_lt(res$h, max(hs))
})

test_that("Kish effective sample sizes are stored per grid point", {
  dat <- .ls_data(n = 150, seed = 19)
  fit <- .qfit(y ~ x, data = dat, taus = c(0.3, 0.6), m = 3, h = 0.2,
               se_type = "none")
  for (i in 1:3) {
    w <- stats::dnorm((dat$x - fit$grid[i]) / 0.2)
    expect_equal(fit$effective_n[i, 1], sum(w)^2 / sum(w^2))
    expect_equal(fit$effective_n[i, 1], fit$effective_n[i, 2])
  }
})

test_that("grid = 'range' gives an equally spaced grid with ECDF thetas", {
  dat <- .ls_data(n = 120, seed = 20)
  fit <- .qfit(y ~ x, data = dat, taus = 0.5, m = 4, h = 0.2,
               grid = "range", se_type = "none")
  expect_equal(fit$grid, seq(min(dat$x), max(dat$x), length.out = 4))
  expect_equal(fit$thetas, sapply(fit$grid, function(g) mean(dat$x <= g)))
})
