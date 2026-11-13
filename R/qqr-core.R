## =========================================================================
## qqr-core.R
## Main fitting function qqr() and the local estimation helpers.
## =========================================================================

#' Quantile-on-Quantile Regression
#'
#' @description
#' Fits a quantile-on-quantile regression (QQR) model (Sim and Zhou, 2015)
#' by running a sequence of kernel-weighted local polynomial quantile
#' regressions of a dependent variable \eqn{Y} on a designated quantile
#' regressor \eqn{X} and optional control variables \eqn{W}. The local
#' regressions are evaluated at a grid of locations of \eqn{X} (by default
#' its empirical \eqn{\theta}-quantiles) and at a set of conditional
#' quantile levels \eqn{\tau}, producing coefficient and p-value surfaces
#' indexed by \eqn{(\theta, \tau)}.
#'
#' @details
#' \strong{Estimator.}
#' For each grid location \eqn{x^\theta} and each \eqn{\tau} in
#' \code{taus}, \code{qqr()} solves the weighted check-loss problem
#' \deqn{\min_{\alpha, \beta, \gamma} \sum_{t=1}^{n}
#'   \rho_\tau\{Y_t - \alpha - \beta (X_t - x^\theta) - \gamma^\top W_t\}
#'   \, K(d_t / h),}
#' where \eqn{\rho_\tau(u) = u\{\tau - I(u < 0)\}} and \eqn{K} is the
#' Gaussian kernel. The local slope \eqn{\hat\beta(\theta, \tau)} is the
#' QQR estimate of the effect of the \eqn{\theta}-quantile of \eqn{X} on
#' the \eqn{\tau}-quantile of \eqn{Y}; \eqn{\hat\alpha(\theta, \tau)} is
#' the local intercept (the fitted \eqn{\tau}-quantile of \eqn{Y} at
#' \eqn{X = x^\theta} and \eqn{W = 0}); and \eqn{\hat\gamma(\theta, \tau)}
#' are local control coefficients. Each problem is solved with
#' \code{\link[quantreg]{rq}} (Barrodale-Roberts simplex).
#'
#' \strong{Grid (\code{grid}).}
#' With \code{grid = "quantile"} (the default) the grid points are the
#' empirical quantiles \eqn{x^\theta = \hat Q_X(\theta)} at
#' \code{thetas} (by default \code{m} equally spaced levels from 0.05 to
#' 0.95), as in Sim and Zhou (2015). With \code{grid = "range"} the grid is
#' \code{m} equally spaced values between \code{min(X)} and \code{max(X)}.
#' In both cases
#' the stored \code{thetas} are the quantile levels of the grid points
#' (for \code{"range"}, the empirical CDF evaluated at the grid).
#'
#' \strong{Kernel distance (\code{kernel_space}).}
#' With \code{kernel_space = "x"} (default) the kernel distance is
#' \eqn{d_t = X_t - x^\theta}, i.e. the bandwidth is measured in the units
#' of \eqn{X}, as in local linear quantile regression (Yu and Jones,
#' 1998). With \code{kernel_space = "cdf"} the distance is
#' \eqn{d_t = F_n(X_t) - \theta}, where \eqn{F_n} is the empirical CDF of
#' \eqn{X}, exactly as in Sim and Zhou (2015); the bandwidth is then on the
#' probability scale (Sim and Zhou use \eqn{h = 0.05}) and automatic
#' selectors are applied to the ranks \eqn{F_n(X_t)}.
#'
#' \strong{Local polynomial degree (\code{degree}).}
#' With \code{degree = 1} the problem above (local linear) is solved, as
#' in Sim and Zhou (2015). With \code{degree = 2} (default) the term
#' \eqn{\delta (X_t - x^\theta)^2} is added to the local model (local
#' quadratic). The slope is a first derivative; the local linear slope has
#' a bias of order \eqn{h} near the boundary of the support of \eqn{X},
#' whereas the local quadratic slope has a bias of order \eqn{h^2}
#' everywhere, which is why an odd difference between polynomial degree
#' and derivative order is recommended (Fan and Gijbels, 1996). The
#' coefficient \eqn{\hat\delta} is not stored.
#'
#' \strong{Bandwidth (\code{h}).}
#' Either a single positive number (a reproducibility warning is issued) or
#' one of the automatic rules.
#' \itemize{
#'   \item \code{"slope"} (default): \eqn{h = C\,\hat\sigma\,n^{-1/7}}
#'     with \eqn{\hat\sigma = \min(\mathrm{sd}, \mathrm{IQR}/1.34)}. The
#'     variance of a local slope is of order \eqn{1/(n h^3)}, so its
#'     optimal bandwidth is of order \eqn{n^{-1/7}} (Fan and Gijbels,
#'     1996), wider than the \eqn{n^{-1/5}} of density and level
#'     estimation. The constant \eqn{C} was calibrated by simulation; see
#'     the methodology vignette.
#'   \item \code{"cv_slope"}: the cross-validated bandwidth of
#'     \code{"cv_fast"} multiplied by \eqn{n^{1/5 - 1/7}}, a rate
#'     adjustment from level to slope estimation.
#'   \item \code{"silverman"}, \code{"robust"}, \code{"iqr"},
#'     \code{"sj"}: density bandwidths of order \eqn{n^{-1/5}}.
#'   \item \code{"cv_fast"}, \code{"yj"}, \code{"npcv"}, \code{"sc"}:
#'     bandwidths for estimating the level of a regression function.
#' }
#' The \eqn{n^{-1/5}} rules under-smooth slopes and give noisy surfaces.
#' \code{"yj"} needs \pkg{KernSmooth}; \code{"npcv"} and \code{"sc"} need
#' \pkg{np}. No rule is optimal for every data set: use
#' \code{\link{qqr_sensitivity}} to check that conclusions are stable
#' across bandwidths.
#'
#' \strong{Inference.}
#' P-values are those reported by \code{\link[quantreg]{summary.rq}} for
#' the weighted local problem with \code{se = se_type}. In simulations
#' without any effect, 5\% tests based on \code{"ker"} (Powell kernel
#' sandwich, the default) rejected in 3.5--4.4\% of the cells, those based
#' on \code{"boot"} in about 5\%, and those based on \code{"nid"} in
#' 7--9\%; \code{"boot"} is therefore recommended for final results (see
#' the methodology vignette). All p-values are pointwise
#' and do not account for bandwidth selection, smoothing bias or the
#' multiplicity of grid cells; treat them as descriptive
#' (\code{\link{qqr_padjust}} adjusts them for multiplicity). Use
#' \code{se_type = "none"} to skip them (faster). With
#' \code{se_type = "boot"} results depend on the random number generator;
#' call \code{set.seed()} first.
#'
#' \strong{Formula support.}
#' Controls may be numeric or factors (factors are expanded with
#' treatment contrasts, so their surfaces are named after the dummy
#' columns, e.g. \code{"regionB"}). The quantile regressor must be numeric
#' and must enter the formula as a main effect only (it may be
#' transformed, e.g. \code{qvar = "log(x)"}, but must not appear in an
#' interaction).
#'
#' @param formula A model formula such as \code{y ~ x} or
#'   \code{y ~ x + w1 + w2}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param qvar Character string naming the quantile regressor exactly as it
#'   appears in the formula (e.g. \code{"x"} or \code{"log(x)"}). May be
#'   omitted when the formula has a single regressor.
#' @param taus Strictly increasing numeric vector of conditional quantile
#'   levels in \eqn{(0, 1)}. Default \code{seq(0.05, 0.95, by = 0.05)}.
#' @param m Integer \eqn{\ge 2}; number of grid points. Ignored when
#'   \code{thetas} is supplied. Default \code{19}.
#' @param h Bandwidth: a single positive number or the name of an automatic
#'   rule (see Details). Default \code{"slope"}.
#' @param kernel Kernel function. Only \code{"gaussian"} is implemented.
#' @param se_type Standard-error type passed to
#'   \code{\link[quantreg]{summary.rq}}: \code{"ker"} (default),
#'   \code{"nid"}, \code{"boot"} or \code{"none"}.
#' @param ... Further arguments passed to \code{\link[quantreg]{summary.rq}}
#'   (for example \code{R} when \code{se_type = "boot"}).
#' @param na.action Function for handling missing values. Default
#'   \code{stats::na.omit}.
#' @param grid Either \code{"quantile"} (default; grid at empirical
#'   quantiles of the quantile regressor) or \code{"range"} (equally spaced
#'   over its range).
#' @param thetas Optional strictly increasing vector of quantile levels in
#'   \eqn{(0, 1)} for the grid when \code{grid = "quantile"}. Defaults to
#'   \code{seq(0.05, 0.95, length.out = m)}.
#' @param kernel_space Either \code{"x"} (default; kernel on
#'   \eqn{X_t - x^\theta}) or \code{"cdf"} (kernel on
#'   \eqn{F_n(X_t) - \theta}, as in Sim and Zhou, 2015).
#' @param degree Degree of the local polynomial in the quantile regressor:
#'   \code{2} (local quadratic, default) or \code{1} (local linear).
#'
#' @return An object of class \code{"qqr"}: a list with components
#' \describe{
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{formula}}{The model formula.}
#'   \item{\code{dependent}, \code{independent}}{Names of \eqn{Y} and of
#'     the quantile regressor.}
#'   \item{\code{controls}}{Names of the control coefficients (dummy
#'     columns for factors).}
#'   \item{\code{grid}, \code{thetas}}{Grid values of \eqn{X} and their
#'     quantile levels (length \eqn{m}).}
#'   \item{\code{taus}}{The conditional quantile levels (length
#'     \eqn{K}).}
#'   \item{\code{coef_surfaces}, \code{pval_surfaces}}{Named lists of
#'     \eqn{m \times K} matrices (rows = grid points, columns = taus), one
#'     per coefficient.}
#'   \item{\code{qcoef}}{The slope surface of the quantile regressor.}
#'   \item{\code{intercept}}{Logical; whether the model has an intercept.}
#'   \item{\code{coef_names}}{Coefficient names.}
#'   \item{\code{bandwidth}, \code{h_method}, \code{h_details}}{The
#'     bandwidth, its selection method and a reporting description.}
#'   \item{\code{grid_type}, \code{kernel_space}, \code{degree},
#'     \code{kernel}, \code{se_type}}{Estimation settings.}
#'   \item{\code{effective_n}}{\eqn{m \times K} matrix of Kish effective
#'     local sample sizes \eqn{(\sum w)^2 / \sum w^2}.}
#'   \item{\code{n_failed}}{Number of local fits that failed (their cells
#'     are \code{NA}).}
#'   \item{\code{n}}{Number of observations used.}
#'   \item{\code{data}}{The model frame.}
#'   \item{\code{y}, \code{design}}{Response vector and (uncentred) design
#'     matrix, used by \code{\link{compare.qqr}}.}
#' }
#'
#' @seealso \code{\link{summary.qqr}}, \code{\link{qqr_surface}},
#'   \code{\link{compare.qqr}}, \code{\link{surface3d}},
#'   \code{\link{qqr_heatmap}}, \code{\link{excel.qqr}}
#'
#' @references
#' Fan, J. and Gijbels, I. (1996) \emph{Local Polynomial Modelling and Its
#' Applications}. London: Chapman & Hall.
#'
#' Koenker, R. and Bassett, G. (1978) Regression quantiles.
#' \emph{Econometrica}, 46(1), 33--50. \doi{10.2307/1913643}
#'
#' Sim, N. and Zhou, H. (2015) Oil prices, US stock return, and the
#' dependence between their quantiles. \emph{Journal of Banking & Finance},
#' 55, 1--8. \doi{10.1016/j.jbankfin.2015.01.013}
#'
#' Silverman, B. W. (1986) \emph{Density Estimation for Statistics and Data
#' Analysis}. London: Chapman & Hall.
#'
#' Yu, K. and Jones, M. C. (1998) Local linear quantile regression.
#' \emph{Journal of the American Statistical Association}, 93(441),
#' 228--237. \doi{10.1080/01621459.1998.10474104}
#'
#' @examples
#' set.seed(42)
#' n   <- 200
#' x   <- rnorm(n)
#' w   <- rnorm(n)
#' y   <- 0.5 + x + 0.3 * w + (1 + 0.5 * abs(x)) * rnorm(n)
#' dat <- data.frame(y = y, x = x, w = w)
#'
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.1, 0.25, 0.5, 0.75, 0.9), m = 5)
#' fit
#' summary(fit)
#' round(fit$qcoef, 2)
#'
#' ## Sim and Zhou (2015) weighting on the empirical-CDF scale
#' fit_sz <- suppressWarnings(
#'   qqr(y ~ x + w, data = dat, qvar = "x", taus = c(0.25, 0.5, 0.75),
#'       m = 5, kernel_space = "cdf", h = 0.05, degree = 1)
#' )
#'
#' @export
qqr <- function(formula,
                data,
                qvar         = NULL,
                taus         = seq(0.05, 0.95, by = 0.05),
                m            = 19L,
                h            = "slope",
                kernel       = c("gaussian"),
                se_type      = c("ker", "nid", "boot", "none"),
                ...,
                na.action    = stats::na.omit,
                grid         = c("quantile", "range"),
                thetas       = NULL,
                kernel_space = c("x", "cdf"),
                degree       = 2L) {

  cl           <- match.call()
  kernel       <- match.arg(kernel)
  se_type      <- match.arg(se_type)
  grid_type    <- match.arg(grid)
  kernel_space <- match.arg(kernel_space)

  ## -----------------------------------------------------------------------
  ## 1. Validate scalar arguments
  ## -----------------------------------------------------------------------
  if (!inherits(formula, "formula"))
    stop("'formula' must be a formula, e.g. y ~ x + w.", call. = FALSE)
  if (!is.data.frame(data))
    stop("'data' must be a data frame.", call. = FALSE)
  if (nrow(data) == 0L)
    stop("'data' has no rows.", call. = FALSE)

  taus <- .qqr_check_probs(taus, "taus")

  if (!is.numeric(degree) || length(degree) != 1L || is.na(degree) ||
      !(degree %in% c(1, 2)))
    stop("'degree' must be 1 (local linear) or 2 (local quadratic).",
         call. = FALSE)
  degree <- as.integer(degree)

  if (is.null(thetas)) {
    if (!is.numeric(m) || length(m) != 1L || !is.finite(m) || m < 2 ||
        m != round(m))
      stop("'m' must be a single whole number >= 2.", call. = FALSE)
    m <- as.integer(m)
  } else {
    if (grid_type == "range")
      stop("'thetas' can only be used with grid = \"quantile\".",
           call. = FALSE)
    thetas <- .qqr_check_probs(thetas, "thetas")
    m      <- length(thetas)
    if (m < 2L) stop("'thetas' must contain at least two values.",
                     call. = FALSE)
  }

  ## -----------------------------------------------------------------------
  ## 2. Model frame, response and design matrix
  ## -----------------------------------------------------------------------
  mf <- stats::model.frame(formula, data = data, na.action = na.action,
                           drop.unused.levels = TRUE)
  if (nrow(mf) == 0L)
    stop("No complete observations remain after handling missing values.",
         call. = FALSE)

  tt <- attr(mf, "terms")
  if (attr(tt, "response") != 1L)
    stop("'formula' must have a response on the left-hand side.",
         call. = FALSE)

  term_labels <- attr(tt, "term.labels")
  if (length(term_labels) == 0L)
    stop("Formula must have at least one independent variable.",
         call. = FALSE)

  dep_name <- names(mf)[1L]
  y        <- stats::model.response(mf)
  if (!is.numeric(y) || is.matrix(y))
    stop("The dependent variable must be numeric (a single numeric ",
         "column).", call. = FALSE)

  ## Quantile regressor
  if (is.null(qvar)) {
    if (length(term_labels) != 1L)
      stop("Formula has multiple regressors. Set 'qvar' to one of: ",
           paste(term_labels, collapse = ", "), call. = FALSE)
    qvar <- term_labels[1L]
  }
  if (!is.character(qvar) || length(qvar) != 1L)
    stop("'qvar' must be a single character string.", call. = FALSE)
  if (!(qvar %in% term_labels))
    stop("qvar '", qvar, "' is not an independent variable in the ",
         "formula. Available: ", paste(term_labels, collapse = ", "),
         call. = FALSE)

  ## qvar must not appear inside an interaction term
  fac <- attr(tt, "factors")
  if (qvar %in% rownames(fac)) {
    in_terms <- colnames(fac)[fac[qvar, ] > 0]
    if (length(setdiff(in_terms, qvar)) > 0L)
      stop("The quantile regressor '", qvar, "' appears in an interaction (",
           paste(setdiff(in_terms, qvar), collapse = ", "),
           "). Interactions with the quantile regressor are not supported.",
           call. = FALSE)
  }

  x <- mf[[qvar]]
  if (!is.numeric(x) || is.matrix(x))
    stop("The quantile regressor '", qvar, "' must be numeric.",
         call. = FALSE)
  if (any(!is.finite(x)) || any(!is.finite(y)))
    stop("The dependent variable and quantile regressor must be finite ",
         "(no Inf/NaN).", call. = FALSE)
  if (length(unique(x)) < 3L)
    stop("The quantile regressor must take at least three distinct values.",
         call. = FALSE)

  X <- stats::model.matrix(tt, mf)
  if (!(qvar %in% colnames(X)))
    stop("Could not locate '", qvar, "' in the design matrix.",
         call. = FALSE)
  if (any(!is.finite(X)))
    stop("The design matrix contains non-finite values.", call. = FALSE)

  has_intercept <- attr(tt, "intercept") == 1L
  int_col       <- if (has_intercept) "(Intercept)" else character(0L)
  control_cols  <- setdiff(colnames(X), c(int_col, qvar))
  coef_names    <- c(int_col, qvar, control_cols)
  X             <- X[, coef_names, drop = FALSE]
  n             <- nrow(X)
  p             <- ncol(X)

  if (n <= p + degree)
    stop("Too few observations (", n, ") for ", p + degree - 1L,
         " coefficients.",
         call. = FALSE)

  ## -----------------------------------------------------------------------
  ## 3. Grid
  ## -----------------------------------------------------------------------
  Fx <- .qqr_ecdf_ranks(x)
  if (grid_type == "quantile") {
    if (is.null(thetas)) thetas <- seq(0.05, 0.95, length.out = m)
    grid_vals <- as.numeric(stats::quantile(x, probs = thetas, type = 7,
                                            names = FALSE))
  } else {
    grid_vals <- seq(min(x), max(x), length.out = m)
    thetas    <- vapply(grid_vals, function(g) mean(x <= g), numeric(1L))
  }

  ## -----------------------------------------------------------------------
  ## 4. Bandwidth
  ## -----------------------------------------------------------------------
  s_var <- if (kernel_space == "x") x else Fx
  if (is.character(h)) {
    if (length(h) != 1L || !(h %in% .qqr_h_methods()))
      stop("h must be a positive number or one of: ",
           paste(sprintf("'%s'", .qqr_h_methods()), collapse = ", "),
           ". Got: ", deparse(h), call. = FALSE)
    hr <- .qqr_select_h(
      x      = s_var,
      y      = if (h %in% .qqr_h_needs_y()) y else NULL,
      method = h,
      kernel = kernel
    )
    h_val     <- hr$h
    h_method  <- hr$method
    h_details <- hr$details
    if (kernel_space == "cdf")
      h_details <- paste0(h_details,
                          " [applied to the empirical-CDF ranks of ", qvar,
                          "]")
  } else if (is.numeric(h) && length(h) == 1L && is.finite(h) && h > 0) {
    h_val     <- h
    h_method  <- "manual"
    h_details <- paste0("Manually specified h = ", format(h),
                        if (kernel_space == "cdf")
                          " (probability scale)" else "")
    warning("Bandwidth h was set manually to ", format(h), ". ",
            "Manual bandwidth choices affect the results; report them and ",
            "check sensitivity. Automatic rules: ",
            paste(.qqr_h_methods(), collapse = ", "), ".", call. = FALSE)
  } else {
    stop("h must be a single positive number or one of: ",
         paste(sprintf("'%s'", .qqr_h_methods()), collapse = ", "),
         ". Got: ", deparse(h), call. = FALSE)
  }
  if (kernel_space == "cdf" && h_val >= 1)
    warning("With kernel_space = \"cdf\" the bandwidth is on the ",
            "probability scale; h = ", format(h_val), " >= 1 gives almost ",
            "uniform weights.", call. = FALSE)

  ## -----------------------------------------------------------------------
  ## 5. Local fits
  ## -----------------------------------------------------------------------
  K     <- length(taus)
  q_col <- match(qvar, coef_names)
  sf <- .qqr_surfaces(y = y, X = X, x = x, Fx = Fx, grid_vals = grid_vals,
                      thetas = thetas, taus = taus, h = h_val,
                      kernel = kernel, kernel_space = kernel_space,
                      degree = degree, se_type = se_type, q_col = q_col, ...)
  coef_surfaces <- sf$coef_surfaces
  pval_surfaces <- sf$pval_surfaces
  effective_n   <- sf$effective_n
  n_failed      <- sf$n_failed
  n_nonuniq     <- sf$n_nonuniq

  if (n_failed > 0L)
    warning(n_failed, " of ", m * K, " local quantile regressions could not ",
            "be estimated (too little local information, e.g. a small ",
            "bandwidth or collinear controls). The corresponding cells are ",
            "NA.", call. = FALSE)
  if (n_nonuniq > 0L && n_nonuniq == m * K)
    message("Note: quantreg reported non-unique solutions in all local ",
            "fits; this is common with discrete data and is usually harmless.")

  if (se_type != "none") {
    n_pna <- sum(is.na(pval_surfaces[[q_col]]) &
                   !is.na(coef_surfaces[[q_col]]))
    if (n_pna > 0L)
      warning("P-values could not be computed for ", n_pna, " of ", m * K,
              " cells of '", qvar, "' (they are NA).", call. = FALSE)
  }

  structure(
    list(
      call          = cl,
      formula       = formula,
      dependent     = dep_name,
      independent   = qvar,
      controls      = control_cols,
      grid          = grid_vals,
      thetas        = thetas,
      taus          = taus,
      coef_surfaces = coef_surfaces,
      pval_surfaces = pval_surfaces,
      qcoef         = coef_surfaces[[qvar]],
      intercept     = has_intercept,
      coef_names    = coef_names,
      bandwidth     = h_val,
      h_method      = h_method,
      h_details     = h_details,
      grid_type     = grid_type,
      kernel_space  = kernel_space,
      degree        = degree,
      kernel        = kernel,
      se_type       = se_type,
      effective_n   = effective_n,
      n_failed      = n_failed,
      n             = n,
      data          = mf,
      y             = as.numeric(y),
      design        = X
    ),
    class = "qqr"
  )
}


## =========================================================================
## Internal helpers
## =========================================================================

## Validate a vector of probabilities.
## @noRd
.qqr_check_probs <- function(v, name) {
  if (!is.numeric(v) || length(v) == 0L || any(!is.finite(v)))
    stop("'", name, "' must be a non-empty numeric vector without missing ",
         "values.", call. = FALSE)
  if (any(v <= 0 | v >= 1))
    stop("All '", name, "' must lie strictly between 0 and 1.",
         call. = FALSE)
  if (length(v) > 1L && any(diff(v) <= 0))
    stop("'", name, "' must be strictly increasing (sorted, no duplicates).",
         call. = FALSE)
  as.numeric(v)
}

## Empirical CDF ranks F_n(X_t) = #{s : X_s <= X_t} / n.
## @noRd
.qqr_ecdf_ranks <- function(x) {
  stats::ecdf(x)(x)
}

## Kernel function.
## @noRd
.qqr_kernel <- function(u, kernel = "gaussian") {
  switch(kernel,
         gaussian = stats::dnorm(u),
         stop("Unsupported kernel '", kernel, "'.", call. = FALSE))
}

## Kish effective sample size (sum w)^2 / sum w^2.
## @noRd
.qqr_kish <- function(w) {
  sw2 <- sum(w^2)
  if (!is.finite(sw2) || sw2 <= 0) return(0)
  sum(w)^2 / sw2
}

## All local fits for one bandwidth. X is the uncentred design matrix whose
## column q_col holds the quantile regressor x; Fx are its ECDF ranks.
## @noRd
.qqr_surfaces <- function(y, X, x, Fx, grid_vals, thetas, taus, h, kernel,
                          kernel_space, degree, se_type, q_col, ...) {
  m <- length(grid_vals)
  K <- length(taus)
  p <- ncol(X)
  coef_names <- colnames(X)
  empty <- function() matrix(NA_real_, nrow = m, ncol = K)
  coef_surfaces <- stats::setNames(replicate(p, empty(), simplify = FALSE),
                                   coef_names)
  pval_surfaces <- stats::setNames(replicate(p, empty(), simplify = FALSE),
                                   coef_names)
  effective_n <- empty()
  n_failed  <- 0L
  n_nonuniq <- 0L

  for (i in seq_len(m)) {
    d <- if (kernel_space == "x") x - grid_vals[i] else Fx - thetas[i]
    w <- .qqr_kernel(d / h, kernel)
    effective_n[i, ] <- .qqr_kish(w)

    Xc <- X
    Xc[, q_col] <- x - grid_vals[i]
    ## local quadratic: add the squared deviation (its coefficient, half
    ## the local curvature, is not stored)
    if (degree == 2L) Xc <- cbind(Xc, (x - grid_vals[i])^2)

    for (j in seq_len(K)) {
      res <- .qqr_local_fit(y, Xc, w, taus[j], se_type, ...)
      n_nonuniq <- n_nonuniq + res$nonunique
      if (is.null(res$coef)) {
        n_failed <- n_failed + 1L
        next
      }
      for (k in seq_len(p)) {
        coef_surfaces[[k]][i, j] <- res$coef[k]
        pval_surfaces[[k]][i, j] <- res$pval[k]
      }
    }
  }
  list(coef_surfaces = coef_surfaces, pval_surfaces = pval_surfaces,
       effective_n = effective_n, n_failed = n_failed,
       n_nonuniq = n_nonuniq)
}

## One weighted local quantile regression.
## y: response; X: local (centred) design including intercept column;
## w: kernel weights. Returns list(coef, pval, nonunique) with coef = NULL
## on failure.
## @noRd
.qqr_local_fit <- function(y, X, w, tau, se_type, ...) {
  p    <- ncol(X)
  keep <- w > 0
  fail <- list(coef = NULL, pval = NULL, nonunique = 0L)
  if (sum(keep) <= p || .qqr_kish(w) < p + 1)
    return(fail)

  Xk <- X[keep, , drop = FALSE]
  wk <- w[keep]
  if (qr(Xk * sqrt(wk))$rank < p)
    return(fail)

  ## Safe column names so that any variable name (e.g. "log(x)" or factor
  ## dummies) can be used inside the model formula.
  vn <- paste0("v", seq_len(p))
  df <- as.data.frame(Xk)
  names(df) <- vn
  df$.y <- y[keep]
  df$.w <- wk
  form <- stats::as.formula(paste(".y ~ 0 +", paste(vn, collapse = " + ")))

  nonunique <- 0L
  fit <- tryCatch(
    withCallingHandlers(
      quantreg::rq(form, data = df, weights = .w, tau = tau, method = "br"),
      warning = function(wn) {
        if (grepl("nonunique", conditionMessage(wn))) {
          nonunique <<- 1L
          invokeRestart("muffleWarning")
        }
      }
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || length(fit$coefficients) != p ||
      any(!is.finite(fit$coefficients)))
    return(fail)

  cf   <- as.numeric(fit$coefficients[vn])
  pval <- rep(NA_real_, p)

  if (!identical(se_type, "none")) {
    sm <- tryCatch(
      suppressWarnings(quantreg::summary.rq(fit, se = se_type, ...)),
      error = function(e) NULL
    )
    if (!is.null(sm) && !is.null(sm$coefficients)) {
      cm   <- sm$coefficients
      pcol <- grep("^Pr\\(", colnames(cm))
      if (length(pcol) == 1L && all(vn %in% rownames(cm)))
        pval <- as.numeric(cm[vn, pcol])
    }
  }

  list(coef = cf, pval = pval, nonunique = nonunique)
}

## Map a user-facing coefficient name to the stored surface name.
## Accepts "intercept"/"(Intercept)" (case-insensitive) for the intercept.
## @noRd
.qqr_match_coef <- function(x, coef) {
  if (!is.character(coef) || length(coef) != 1L || is.na(coef))
    stop("'coef' must be a single character string.", call. = FALSE)
  all_names <- names(x$coef_surfaces)
  if (coef %in% all_names) return(coef)
  if (tolower(coef) %in% c("intercept", "(intercept)") &&
      "(Intercept)" %in% all_names)
    return("(Intercept)")
  if (identical(coef, x$independent) && "z" %in% all_names)
    return("z")   ## legacy alias
  stop("Coefficient '", coef, "' not recognised. Available: ",
       paste(.qqr_user_coefs(x), collapse = ", "), call. = FALSE)
}

## User-facing coefficient names.
## @noRd
.qqr_user_coefs <- function(x) {
  vapply(names(x$coef_surfaces), function(nm) {
    if (nm == "z") x$independent else nm
  }, character(1L), USE.NAMES = FALSE)
}

## Quantile levels of the grid points, also for objects fitted with
## an earlier development version (which did not store them).
## @noRd
.qqr_thetas <- function(x) {
  if (!is.null(x$thetas)) return(x$thetas)
  xv <- x$data[[x$independent]]
  if (is.null(xv)) return(seq(0.05, 0.95, length.out = length(x$grid)))
  vapply(x$grid, function(g) mean(xv <= g), numeric(1L))
}

## Surface matrix for a coefficient.
## @noRd
.qqr_get_surface <- function(x, coef, what = c("coef", "pvalue")) {
  what <- match.arg(what)
  nm   <- .qqr_match_coef(x, coef)
  mat  <- if (what == "coef") x$coef_surfaces[[nm]] else x$pval_surfaces[[nm]]
  if (is.null(mat))
    stop("No ", what, " surface stored for '", coef, "'.", call. = FALSE)
  if (what == "pvalue" && all(is.na(mat)))
    stop("The p-value surface for '", coef, "' is empty. Re-fit with ",
         "se_type other than \"none\".", call. = FALSE)
  mat
}
