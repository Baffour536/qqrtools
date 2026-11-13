## =========================================================================
## qqr-robustness.R
## Bandwidth sensitivity analysis and multiplicity-adjusted p-values.
## =========================================================================

#' Bandwidth sensitivity of a QQR surface
#'
#' @description
#' Re-estimates a fitted quantile-on-quantile regression over several
#' bandwidths and summarises how much the surface of one coefficient
#' changes. All other settings (grid, quantile levels, kernel space, local
#' polynomial degree, controls) are those of the original fit.
#'
#' @details
#' No automatic rule is optimal for every data set, and a bandwidth chosen
#' after looking at the results invalidates naive inference. A transparent
#' alternative is to report the surface for a pre-specified set of
#' bandwidths around the selected one and to interpret only the features
#' that are stable across them.
#'
#' For each bandwidth the function reports the mean, standard deviation and
#' range of the surface and three comparisons with the baseline surface
#' (the surface of \code{x}): the Pearson correlation between the cells,
#' the share of cells with the same sign, and the mean absolute difference.
#' P-values are not recomputed.
#'
#' @param x An object of class \code{"qqr"}.
#' @param factors Positive multipliers applied to the bandwidth of
#'   \code{x}. Ignored when \code{h} is supplied. Default
#'   \code{c(0.5, 0.75, 1, 1.5, 2)}.
#' @param h Optional numeric vector of bandwidths to use instead of
#'   \code{factors} (on the scale of the fit: units of the regressor for
#'   \code{kernel_space = "x"}, probability scale for \code{"cdf"}).
#' @param coef Name of the coefficient. Defaults to the quantile regressor.
#'
#' @return An object of class \code{"qqr_sensitivity"}: a list with
#' \describe{
#'   \item{\code{summary}}{Data frame with one row per bandwidth: \code{h},
#'     \code{factor} (\code{h} divided by the baseline bandwidth),
#'     \code{mean}, \code{sd}, \code{min}, \code{max},
#'     \code{cor_baseline}, \code{sign_agreement}, \code{mean_abs_diff}
#'     and \code{n_na} (cells that could not be estimated).}
#'   \item{\code{surfaces}}{List of \eqn{m \times K} coefficient matrices,
#'     one per bandwidth.}
#'   \item{\code{tau_means}}{Matrix of surface averages over the grid for
#'     each tau (rows = bandwidths, columns = taus).}
#'   \item{\code{h}, \code{baseline_h}, \code{coef}, \code{thetas},
#'     \code{taus}}{Settings.}
#' }
#'
#' @seealso \code{\link{qqr}}, \code{\link{qqr_padjust}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(200))
#' dat$y <- 0.5 * dat$x + (1 + 0.3 * abs(dat$x)) * rnorm(200)
#' fit <- qqr(y ~ x, data = dat, taus = c(0.25, 0.5, 0.75), m = 5,
#'            se_type = "none")
#' sens <- qqr_sensitivity(fit)
#' sens
#' plot(sens)
#'
#' @export
qqr_sensitivity <- function(x, factors = c(0.5, 0.75, 1, 1.5, 2), h = NULL,
                            coef = NULL) {
  if (!inherits(x, "qqr")) stop("'x' must be a 'qqr' object.", call. = FALSE)
  if (is.null(x$y) || is.null(x$design))
    stop("'x' does not store the response and design matrix; re-fit it ",
         "with the current version of qqr().", call. = FALSE)
  if (is.null(coef)) coef <- x$independent
  nm <- .qqr_match_coef(x, coef)

  if (is.null(h)) {
    if (!is.numeric(factors) || length(factors) == 0L ||
        any(!is.finite(factors)) || any(factors <= 0))
      stop("'factors' must be positive numbers.", call. = FALSE)
    h <- sort(unique(factors)) * x$bandwidth
  } else {
    if (!is.numeric(h) || length(h) == 0L || any(!is.finite(h)) ||
        any(h <= 0))
      stop("'h' must be positive numbers.", call. = FALSE)
    h <- sort(unique(h))
  }

  xv     <- as.numeric(x$design[, x$independent])
  q_col  <- match(x$independent, colnames(x$design))
  kspace <- if (is.null(x$kernel_space)) "x" else x$kernel_space
  degree <- if (is.null(x$degree)) 1L else x$degree
  base   <- x$coef_surfaces[[nm]]

  surfaces <- lapply(h, function(hv) {
    .qqr_surfaces(y = x$y, X = x$design, x = xv, Fx = .qqr_ecdf_ranks(xv),
                  grid_vals = x$grid, thetas = .qqr_thetas(x),
                  taus = x$taus, h = hv, kernel = x$kernel,
                  kernel_space = kspace, degree = degree,
                  se_type = "none", q_col = q_col)$coef_surfaces[[nm]]
  })
  names(surfaces) <- paste0("h=", format(signif(h, 4L)))

  summ <- do.call(rbind, lapply(seq_along(h), function(k) {
    s  <- surfaces[[k]]
    ok <- is.finite(s) & is.finite(base)
    data.frame(
      h              = h[k],
      factor         = h[k] / x$bandwidth,
      mean           = mean(s, na.rm = TRUE),
      sd             = stats::sd(as.numeric(s), na.rm = TRUE),
      min            = min(s, na.rm = TRUE),
      max            = max(s, na.rm = TRUE),
      cor_baseline   = if (sum(ok) > 2L && stats::sd(s[ok]) > 0 &&
                           stats::sd(base[ok]) > 0)
        stats::cor(s[ok], base[ok]) else NA_real_,
      sign_agreement = if (any(ok)) mean(sign(s[ok]) == sign(base[ok])) else
        NA_real_,
      mean_abs_diff  = if (any(ok)) mean(abs(s[ok] - base[ok])) else NA_real_,
      n_na           = sum(!is.finite(s))
    )
  }))
  rownames(summ) <- NULL

  tau_means <- t(vapply(surfaces, colMeans, numeric(length(x$taus)),
                        na.rm = TRUE))
  colnames(tau_means) <- paste0("tau_", format(x$taus))

  structure(
    list(summary = summ, surfaces = surfaces, tau_means = tau_means,
         h = h, baseline_h = x$bandwidth, coef = coef,
         thetas = .qqr_thetas(x), taus = x$taus),
    class = "qqr_sensitivity"
  )
}

#' @rdname qqr_sensitivity
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (ignored by \code{print}; passed to
#'   \code{\link[graphics]{matplot}} by \code{plot}).
#' @method print qqr_sensitivity
#' @export
print.qqr_sensitivity <- function(x, digits = 3, ...) {
  cat("QQR bandwidth sensitivity for '", x$coef, "'\n", sep = "")
  cat("Baseline bandwidth: ", format(signif(x$baseline_h, 6L)), "\n\n",
      sep = "")
  print(format(x$summary, digits = digits), row.names = FALSE)
  cat("\ncor_baseline, sign_agreement and mean_abs_diff compare each ",
      "surface\nwith the baseline fit.\n", sep = "")
  invisible(x)
}

#' @rdname qqr_sensitivity
#' @description
#' The \code{plot} method draws the grid-averaged coefficient against tau,
#' one line per bandwidth.
#' @method plot qqr_sensitivity
#' @export
plot.qqr_sensitivity <- function(x, ...) {
  L <- nrow(x$tau_means)
  graphics::matplot(x$taus, t(x$tau_means), type = "b", lty = 1,
                    pch = seq_len(L), col = seq_len(L),
                    xlab = expression(tau),
                    ylab = paste0("Grid-averaged coefficient of ", x$coef),
                    ...)
  graphics::legend("topleft", legend = rownames(x$tau_means), lty = 1,
                   pch = seq_len(L), col = seq_len(L), bty = "n",
                   cex = 0.8)
  invisible(x)
}


#' Adjust QQR p-value surfaces for multiple testing
#'
#' @description
#' Replaces the pointwise p-values of a fitted \code{"qqr"} object by
#' p-values adjusted for the number of grid cells, using
#' \code{\link[stats]{p.adjust}}. Each coefficient's
#' \eqn{m \times K} surface is treated as one family of tests.
#'
#' @details
#' A QQR surface with the default settings contains \eqn{19 \times 19 = 361}
#' cells per coefficient, so about 18 cells are expected to be significant
#' at the 5\% level when there is no effect at all. The default
#' \code{method = "BH"} controls the false discovery rate (Benjamini and
#' Hochberg, 1995) and remains valid under the positive dependence that is
#' typical of neighbouring cells; \code{"BY"} (Benjamini and Yekutieli,
#' 2001) is valid under arbitrary dependence and more conservative;
#' \code{"holm"} and \code{"bonferroni"} control the family-wise error
#' rate.
#'
#' The adjustment addresses multiplicity only. The underlying p-values
#' still ignore smoothing bias and bandwidth selection. All plotting and
#' export functions use the adjusted p-values of the returned object.
#'
#' @param x An object of class \code{"qqr"} fitted with \code{se_type}
#'   other than \code{"none"}.
#' @param method Adjustment method; one of
#'   \code{\link[stats]{p.adjust.methods}}. Default \code{"BH"}.
#'
#' @return The \code{"qqr"} object with adjusted \code{pval_surfaces} and a
#'   component \code{p_adjust} naming the method.
#'
#' @references
#' Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery
#' rate: a practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society: Series B}, 57(1),
#' 289--300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Benjamini, Y. and Yekutieli, D. (2001) The control of the false
#' discovery rate in multiple testing under dependency. \emph{The Annals of
#' Statistics}, 29(4), 1165--1188. \doi{10.1214/aos/1013699998}
#'
#' @seealso \code{\link{qqr}}, \code{\link{qqr_sensitivity}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(200))
#' dat$y <- 0.5 * dat$x + rnorm(200)
#' fit <- qqr(y ~ x, data = dat, taus = c(0.25, 0.5, 0.75), m = 5)
#' fit_adj <- qqr_padjust(fit)
#' round(fit_adj$pval_surfaces$x, 3)
#'
#' @export
qqr_padjust <- function(x, method = "BH") {
  if (!inherits(x, "qqr")) stop("'x' must be a 'qqr' object.", call. = FALSE)
  if (!is.character(method) || length(method) != 1L ||
      !(method %in% stats::p.adjust.methods))
    stop("'method' must be one of: ",
         paste(stats::p.adjust.methods, collapse = ", "), call. = FALSE)
  if (!is.null(x$p_adjust))
    stop("The p-values of 'x' have already been adjusted (method '",
         x$p_adjust, "').", call. = FALSE)
  if (identical(x$se_type, "none") ||
      all(vapply(x$pval_surfaces, function(p) all(is.na(p)), logical(1L))))
    stop("'x' has no p-values. Re-fit with se_type other than \"none\".",
         call. = FALSE)
  if (identical(method, "none")) return(x)

  x$pval_surfaces <- lapply(x$pval_surfaces, function(p) {
    p[] <- stats::p.adjust(as.numeric(p), method = method)
    p
  })
  x$p_adjust <- method
  x
}
