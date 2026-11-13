## Input validation, edge cases and formula handling of qqr().

test_that("invalid taus, thetas, m and h are rejected with clear errors", {
  dat <- .ls_data(n = 60, seed = 1)
  expect_error(qqr(y ~ x, dat, taus = c(0, 0.5)), "strictly between 0 and 1")
  expect_error(qqr(y ~ x, dat, taus = c(0.5, 1.2)), "strictly between")
  expect_error(qqr(y ~ x, dat, taus = c(0.5, 0.3)), "strictly increasing")
  expect_error(qqr(y ~ x, dat, taus = c(0.5, NA)), "missing")
  expect_error(qqr(y ~ x, dat, taus = "a"), "numeric")
  expect_error(qqr(y ~ x, dat, m = 1), "'m'")
  expect_error(qqr(y ~ x, dat, m = 2.5), "'m'")
  expect_error(qqr(y ~ x, dat, thetas = c(0.2, 0.2)), "strictly increasing")
  expect_error(qqr(y ~ x, dat, thetas = 0.5), "at least two")
  expect_error(qqr(y ~ x, dat, thetas = c(0.2, 0.8), grid = "range"),
               "thetas")
  expect_error(qqr(y ~ x, dat, h = -1), "positive")
  expect_error(qqr(y ~ x, dat, h = c(0.1, 0.2)), "positive")
  expect_error(qqr(y ~ x, dat, h = "bogus"), "positive number or one of")
  expect_error(qqr(y ~ x, dat, se_type = "iid"), "should be one of")
  expect_error(qqr(y ~ x, dat, kernel_space = "rank"), "should be one of")
})

test_that("data problems give informative errors, not silent NAs", {
  dat <- .ls_data(n = 60, seed = 2)
  expect_error(qqr("y ~ x", dat), "formula")
  expect_error(qqr(y ~ x, as.list(dat)), "data frame")
  expect_error(qqr(y ~ x, dat[0, ]), "no rows")
  bad <- dat; bad$x <- 1
  expect_error(qqr(y ~ x, bad), "three distinct values")
  bad <- dat; bad$x <- as.character(round(bad$x, 1))
  expect_error(qqr(y ~ x, bad), "numeric")
  bad <- dat; bad$y[1] <- Inf
  expect_error(qqr(y ~ x, bad), "finite")
  expect_error(qqr(y ~ x * w, dat, qvar = "x"), "interaction")
  expect_error(qqr(y ~ x, dat[1:3, ], m = 2), "Too few observations|three")
  allna <- dat; allna$y <- NA_real_
  expect_error(qqr(y ~ x, allna), "No complete observations")
})

test_that("missing values are dropped via na.action", {
  dat <- .ls_data(n = 80, seed = 3)
  dat$w[c(2, 5)] <- NA
  dat$y[7] <- NA
  fit <- qqr(y ~ x + w, dat, qvar = "x", taus = 0.5, m = 3, se_type = "none")
  expect_equal(fit$n, 77)
  expect_equal(nrow(fit$data), 77)
  expect_error(qqr(y ~ x + w, dat, qvar = "x", na.action = stats::na.fail),
               "missing values")
})

test_that("factor controls are expanded and estimated (not silently NA)", {
  set.seed(4)
  n <- 300
  dat <- data.frame(x = stats::rnorm(n),
                    g = factor(sample(c("a", "b", "c"), n, TRUE)))
  dat$y <- dat$x + ifelse(dat$g == "b", 2, ifelse(dat$g == "c", -1, 0)) +
    stats::rnorm(n, sd = 0.3)
  fit <- qqr(y ~ x + g, dat, qvar = "x", taus = c(0.25, 0.75), m = 4,
             se_type = "none")
  expect_identical(fit$coef_names, c("(Intercept)", "x", "gb", "gc"))
  expect_identical(fit$controls, c("gb", "gc"))
  expect_false(anyNA(fit$coef_surfaces$gb))
  expect_lt(max(abs(fit$coef_surfaces$gb - 2)), 0.4)
  expect_lt(max(abs(fit$coef_surfaces$gc + 1)), 0.4)
  cmp <- compare.qqr(fit, coef = "gb")
  expect_false(anyNA(cmp$qr))
})

test_that("transformed quantile regressors and responses work", {
  dat <- .ls_data(n = 150, seed = 5)
  dat$ex <- exp(dat$x)
  f1 <- qqr(y ~ log(ex) + w, dat, qvar = "log(ex)", taus = c(0.3, 0.7),
            m = 3, se_type = "none")
  f2 <- qqr(y ~ x + w, dat, qvar = "x", taus = c(0.3, 0.7), m = 3,
            se_type = "none")
  expect_equal(unname(f1$qcoef), unname(f2$qcoef), tolerance = 1e-8)
  expect_identical(f1$independent, "log(ex)")
  f3 <- qqr(I(2 * y) ~ x, dat, taus = 0.5, m = 3, se_type = "none")
  f4 <- qqr(y ~ x, dat, taus = 0.5, m = 3, se_type = "none")
  expect_equal(f3$qcoef, 2 * f4$qcoef, tolerance = 1e-8)
})

test_that("models without an intercept are supported", {
  dat <- .ls_data(n = 150, seed = 6)
  fit <- qqr(y ~ 0 + x + w, dat, qvar = "x", taus = 0.5, m = 3,
             se_type = "none")
  expect_false(fit$intercept)
  expect_identical(fit$coef_names, c("x", "w"))
  expect_error(qqr_surface(fit, coef = "intercept"), "not recognised")
})

test_that("tiny bandwidths produce NA cells with a warning, not wrong values", {
  dat <- .ls_data(n = 60, seed = 7)
  expect_warning(
    fit <- .muffle_manual(suppressMessages(
      qqr(y ~ x + w, dat, qvar = "x", taus = 0.5, m = 5, h = 1e-4,
          se_type = "none"))),
    "could not be estimated"
  )
  expect_true(anyNA(fit$qcoef))
  expect_gt(fit$n_failed, 0)
})

test_that("perfectly collinear controls are reported as failed fits", {
  dat <- .ls_data(n = 80, seed = 8)
  dat$w2 <- 2 * dat$w
  expect_warning(
    fit <- qqr(y ~ x + w + w2, dat, qvar = "x", taus = 0.5, m = 3,
               se_type = "none"),
    "could not be estimated"
  )
  expect_true(all(is.na(fit$qcoef)))
})

test_that("a manual bandwidth triggers the reproducibility warning", {
  dat <- .ls_data(n = 60, seed = 9)
  expect_warning(qqr(y ~ x, dat, taus = 0.5, m = 3, h = 0.3,
                     se_type = "none"), "set manually")
  expect_warning(.muffle_manual(qqr(y ~ x, dat, taus = 0.5, m = 3, h = 2,
                                    kernel_space = "cdf", se_type = "none")),
                 "probability scale")
})

test_that("automatic selectors run inside qqr() and are recorded", {
  dat <- .ls_data(n = 120, seed = 10)
  for (meth in c("silverman", "robust", "iqr", "sj", "cv_fast")) {
    fit <- qqr(y ~ x, dat, taus = 0.5, m = 3, h = meth, se_type = "none")
    expect_identical(fit$h_method, meth)
    expect_gt(fit$bandwidth, 0)
    expect_type(fit$h_details, "character")
  }
  fit <- qqr(y ~ x, dat, taus = 0.5, m = 3, kernel_space = "cdf",
             h = "silverman", se_type = "none")
  Fn <- stats::ecdf(dat$x)(dat$x)
  expect_equal(fit$bandwidth, 1.06 * stats::sd(Fn) * 120^(-1 / 5))
  expect_match(fit$h_details, "empirical-CDF")
  skip_if_not_installed("KernSmooth")
  fit <- qqr(y ~ x, dat, taus = 0.5, m = 3, h = "yj", se_type = "none")
  expect_identical(fit$h_method, "yj")
})

test_that("np-based selectors work when np is installed", {
  skip_on_cran()
  skip_if_not_installed("np")
  dat <- .ls_data(n = 60, seed = 21)
  for (meth in c("npcv", "sc")) {
    fit <- suppressWarnings(suppressMessages(
      qqr(y ~ x, dat, taus = 0.5, m = 3, h = meth, se_type = "none")))
    expect_identical(fit$h_method, meth)
    expect_gt(fit$bandwidth, 0)
  }
})

test_that("p-value types ker and boot run; boot is reproducible with a seed", {
  dat <- .ls_data(n = 100, seed = 22)
  fk <- qqr(y ~ x, dat, taus = 0.5, m = 3, se_type = "ker")
  expect_true(all(fk$pval_surfaces$x >= 0 & fk$pval_surfaces$x <= 1,
                  na.rm = TRUE))
  set.seed(5); b1 <- qqr(y ~ x, dat, taus = 0.5, m = 2, se_type = "boot",
                         R = 50)
  set.seed(5); b2 <- qqr(y ~ x, dat, taus = 0.5, m = 2, se_type = "boot",
                         R = 50)
  expect_equal(b1$pval_surfaces$x, b2$pval_surfaces$x)
  fn <- qqr(y ~ x, dat, taus = 0.5, m = 2, se_type = "none")
  expect_true(all(is.na(fn$pval_surfaces$x)))
})

test_that("fits without stored thetas/design still work", {
  dat <- .ls_data(n = 100, seed = 23)
  fit <- .qfit(y ~ x + w, dat, qvar = "x", taus = c(0.3, 0.7), m = 3,
               grid = "range", h = 0.2)
  old <- fit
  old[c("thetas", "design", "y", "grid_type", "kernel_space", "formula",
        "n_failed")] <- NULL
  expect_equal(qqrtools:::.qqr_thetas(old), fit$thetas)
  expect_equal(suppressMessages(compare.qqr(old))$qr,
               suppressMessages(compare.qqr(fit))$qr)
  expect_output(print(old), "Quantile-on-Quantile")
  .with_null_device(expect_invisible(surface3d(old)))
})
