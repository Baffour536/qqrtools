## Visualisation functions and the quantile dependence heatmap.

.fit_small <- function() {
  dat <- .ls_data(n = 200, seed = 31)
  qqr(y ~ x + w, data = dat, qvar = "x", taus = c(0.2, 0.4, 0.6, 0.8),
      m = 4)
}

test_that("base-graphics 3D plots draw and return the fit invisibly", {
  fit <- .fit_small()
  .with_null_device({
    expect_invisible(surface3d(fit))
    expect_invisible(surface3d(fit, coef = "w", what = "pvalue",
                               clim = c(0, 1), add_zero_plane = TRUE,
                               legend = FALSE))
    expect_invisible(surface3d(fit, coef = "intercept", palette = "RdBu",
                               add_zero_plane = TRUE))
    expect_invisible(panel_surface(fit))
    expect_invisible(panel_surface(fit, coef = c("x", "w"), what = "pvalue",
                                   add_zero_plane = TRUE, legend = FALSE,
                                   clim = c(0, 1)))
    expect_invisible(plot(fit, engine = "base"))
    expect_invisible(plot(fit, coef = "w", what = "pvalue", engine = "base"))
  })
  .with_null_device({
    expect_error(surface3d(fit, coef = "nope"), "not recognised")
    expect_error(surface3d(fit, clim = c(1, 0)), "clim")
    expect_error(panel_surface(fit, clim = 3), "clim")
    expect_error(surface3d.qqr(list()), "class 'qqr'")
    expect_error(panel_surface.qqr(list()), "class 'qqr'")
  })
})

test_that("plot axes use the theta levels of the grid", {
  fit <- .fit_small()
  expect_equal(qqrtools:::.qqr_thetas(fit), seq(0.05, 0.95, length.out = 4))
})

test_that("plotly heatmap, contour and hub return plotly objects", {
  skip_if_not_installed("plotly")
  fit <- .fit_small()
  expect_s3_class(qqr_heatmap(fit), "plotly")
  expect_s3_class(qqr_heatmap(fit, coef = "all", sig_overlay = "stars"),
                  "plotly")
  expect_s3_class(qqr_heatmap(fit, coef = "intercept", sig_overlay = "dots",
                              palette = "viridis", grid_values = FALSE),
                  "plotly")
  expect_s3_class(qqr_heatmap(fit, what = "pvalue", palette = "YlOrRd",
                              center_zero = FALSE), "plotly")
  expect_s3_class(qqr_contour(fit, sig_mask = TRUE), "plotly")
  expect_s3_class(qqr_contour(fit, what = "pvalue", palette = "Plasma",
                              zero_contour = FALSE), "plotly")
  for (tp in c("surface", "heatmap", "contour"))
    expect_s3_class(qqr_plotly(fit, type = tp, colorscale = "RdBu",
                               show_zero_plane = TRUE, sig_overlay = TRUE),
                    "plotly")
  expect_error(qqr_heatmap(fit, coef = "nope"), "not recognised")
  expect_error(qqr_plotly(fit, clim = c(2, 1)), "clim")
  fit0 <- qqr(y ~ x, .ls_data(80, 2), taus = c(0.3, 0.6), m = 3,
              se_type = "none")
  expect_error(qqr_heatmap(fit0, sig_overlay = "stars"), "p-value surface")
  expect_error(qqr_contour(fit0, what = "pvalue"), "p-value surface")
})

test_that("diverging colour scales map negative values to blue", {
  cs <- qqrtools:::.qqr_colorscale("RdBu")
  expect_equal(cs[[1]][[1]], 0)
  expect_identical(cs[[1]][[2]], "#053061")      # dark blue at the minimum
  expect_identical(cs[[length(cs)]][[2]], "#67001f")  # dark red at maximum
  vir <- qqrtools:::.qqr_colorscale("plasma")
  expect_type(vir, "list")
  expect_identical(vir[[1]][[2]], grDevices::hcl.colors(11, "Plasma")[1])
  expect_identical(qqrtools:::.qqr_colorscale("Jet"), "Jet")
  expect_equal(qqrtools:::.qqr_zlim(c(-1, 3), NULL, TRUE), c(-3, 3))
  expect_equal(qqrtools:::.qqr_zlim(c(2, 2), NULL, FALSE), c(1.5, 2.5))
})

test_that("Somers' D matches a brute-force definition with ties", {
  set.seed(32)
  x <- sample(1:5, 40, TRUE)
  y <- sample(1:4, 40, TRUE)
  ref <- function(y, x) {
    C <- D <- TY <- 0
    n <- length(x)
    for (i in 1:(n - 1)) for (j in (i + 1):n) {
      sx <- sign(x[j] - x[i]); sy <- sign(y[j] - y[i])
      if (sx != 0 && sy != 0) { if (sx == sy) C <- C + 1 else D <- D + 1 }
      if (sx != 0 && sy == 0) TY <- TY + 1
    }
    (C - D) / (C + D + TY)
  }
  expect_equal(qqrtools:::.qqr_somers_d(y, x), ref(y, x))
  ## D(Y|X) differs from D(X|Y) when ties are present
  expect_false(isTRUE(all.equal(qqrtools:::.qqr_somers_d(y, x),
                                qqrtools:::.qqr_somers_d(x, y))))
  ## without ties it equals Kendall's tau-a
  u <- stats::rnorm(30); v <- u + stats::rnorm(30)
  expect_equal(qqrtools:::.qqr_somers_d(v, u),
               stats::cor(u, v, method = "kendall"))
})

test_that("Blomqvist's beta and co-exceedance match their definitions", {
  set.seed(33)
  x <- stats::rnorm(101); y <- x + stats::rnorm(101)
  mx <- stats::median(x); my <- stats::median(y)
  s  <- sign(x - mx) * sign(y - my)
  expect_equal(qqrtools:::.qqr_blomqvist(y, x),
               (sum(s > 0) - sum(s < 0)) / sum(s != 0))
  expect_equal(qqrtools:::.qqr_coexceedance(y, x, 0, 0),
               mean(x >= 0 & y >= 0))
})

test_that("qqr_cor_heatmap computes the documented matrices", {
  skip_if_not_installed("plotly")
  set.seed(34)
  n <- 600
  x <- stats::rnorm(n); y <- 0.5 * x + stats::rnorm(n)
  fit <- qqr(y ~ x, data.frame(x, y), taus = c(0.25, 0.5, 0.75), m = 3,
             se_type = "none")

  r <- qqr_cor_heatmap(fit)
  expect_s3_class(r, "qqr_cor_heatmap")
  expect_equal(dim(r$matrix), c(3, 3))
  expect_equal(unname(r$x_breaks), unname(stats::quantile(
    x, c(0, 0.375, 0.625, 1))))
  ## recompute one subsample cell by hand
  xb <- findInterval(x, r$x_breaks, rightmost.closed = TRUE)
  yb <- findInterval(y, r$y_breaks, rightmost.closed = TRUE)
  inc <- xb == 2 & yb == 3
  expect_equal(r$matrix[2, 3], stats::cor(x[inc], y[inc], method = "spearman"))
  expect_equal(sum(r$cell_n), n)

  ce <- qqr_cor_heatmap(fit, method = "coexceedance")
  qx <- stats::quantile(x, 0.75); qy <- stats::quantile(y, 0.25)
  expect_equal(ce$matrix[3, 1], mean(x >= qx & y >= qy))
  expect_equal(ce$independence[3, 1], 0.25 * 0.75)
  ## positive dependence => joint upper exceedance above independence
  expect_gt(ce$matrix[2, 2], ce$independence[2, 2])

  bi <- qqr_cor_heatmap(fit, mode = "binary", method = "pearson")
  xi <- as.numeric(x <= stats::quantile(x, 0.5))
  yi <- as.numeric(y <= stats::quantile(y, 0.25))
  expect_equal(bi$matrix[2, 1], stats::cor(xi, yi))

  for (meth in c("kendall", "somers_d", "pearson", "blomqvist"))
    expect_s3_class(qqr_cor_heatmap(fit, method = meth, n_bins = 2),
                    "qqr_cor_heatmap")
  expect_warning(qqr_cor_heatmap(fit, method = "blomqvist", mode = "binary"),
                 "Switching")
  expect_output(print(r), "Quantile Dependence Heatmap")
  expect_error(qqr_cor_heatmap(fit, n_bins = 1), "at least 2")
})

test_that("rgl interactive plot runs when rgl is available", {
  skip_on_cran()
  skip_if_not_installed("rgl")
  old <- options(rgl.useNULL = TRUE)
  on.exit(options(old))
  fit <- .fit_small()
  expect_invisible(Int.plot(fit, add_plane = TRUE))
  rgl::close3d()
})
