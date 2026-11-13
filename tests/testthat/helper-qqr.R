## Shared helpers for the qqrtools test suite.

## Deterministic data with a known conditional quantile function:
##   Q_Y(tau | x, w) = 1 + qnorm(tau) + (2 + 0.5 qnorm(tau)) x + 0.5 w
## for X ~ U(0, 1), i.e. a location-scale model whose QQR slope is
## beta(theta, tau) = 2 + 0.5 qnorm(tau) for every theta.
.ls_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  x <- stats::runif(n)
  w <- stats::rnorm(n)
  e <- stats::rnorm(n)
  y <- 1 + 2 * x + 0.5 * w + (1 + 0.5 * x) * e
  data.frame(y = y, x = x, w = w)
}

## Open a null graphics device for the duration of the calling test.
.with_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

## Quiet fit with a manual bandwidth (suppresses the reproducibility
## warning that is issued by design for manual h).
.qfit <- function(...) suppressWarnings(qqr(...))

## Muffle only the manual-bandwidth reminder, letting other warnings
## through to expect_warning().
.muffle_manual <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("set manually", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
}
