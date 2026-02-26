## =========================================================================
## qqr-core.R
## Core fitting function and internal helpers for the qqr package.
## =========================================================================


## -------------------------------------------------------------------------
## Internal: bandwidth selector
## -------------------------------------------------------------------------
## Computes the kernel bandwidth h from the quantile-regressor vector x
## (and, where needed, the outcome vector y).
##
## method = "silverman"  : h = 1.06 * sd(x) * n^{-1/5}
## method = "robust"     : h = 0.9 * min(sd(x), IQR(x)/1.34) * n^{-1/5}
## method = "iqr"        : h = IQR(x) * n^{-1/5}
## method = "cv_fast"    : LOO quantile-loss, 15-point grid (no extra pkg)
## method = "sj"         : Sheather-Jones solve-the-equation plug-in
##                         (base R bw.SJ; density-based)
## method = "yj"         : Ruppert-Sheather-Wand direct plug-in
##                         (KernSmooth::dpill; regression-based, uses x & y)
## method = "npcv"       : npregbw LS cross-validation (np package)
## method = "sc"         : npscoefbw smooth-coefficient CV (np package)
##
## Returns a named list:
##   $h        numeric bandwidth
##   $method   character label used in print/summary output
##   $details  character description suitable for reporting
## @keywords internal
.qqr_select_h <- function(x, y = NULL, method, taus, m, kernel, se_type) {

  n   <- length(x)
  s   <- stats::sd(x,  na.rm = TRUE)
  iqr <- stats::IQR(x, na.rm = TRUE)

  if (n < 4L)
    stop("Too few observations to select bandwidth automatically.")

  ## ------------------------------------------------------------------
  ## 1. Silverman's rule of thumb
  ##    h = 1.06 * sd(x) * n^{-1/5}
  ##    Source: Silverman (1986, p. 48)
  ## ------------------------------------------------------------------
  if (method == "silverman") {
    h <- 1.06 * s * n^(-0.2)
    return(list(
      h       = h,
      method  = "silverman",
      details = paste0(
        "Silverman's rule of thumb  [h = 1.06 * sd(x) * n^(-1/5)]",
        "\n               Reference : Silverman (1986) Density Estimation",
        " for Statistics and Data Analysis, Chapman & Hall",
        "\n               sd(x) = ", round(s, 6L),
        ",  n = ", n,
        ",  h = ", round(h, 6L)
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 2. Robust rule of thumb  (Silverman with IQR correction)
  ##    h = 0.9 * min(sd(x), IQR(x)/1.34) * n^{-1/5}
  ##    Sources: Silverman (1986, p. 47); Scott (1992, p. 131)
  ## ------------------------------------------------------------------
  if (method == "robust") {
    scale <- min(s, iqr / 1.34)
    if (!is.finite(scale) || scale <= 0) scale <- s
    h <- 0.9 * scale * n^(-0.2)
    return(list(
      h       = h,
      method  = "robust",
      details = paste0(
        "Robust rule of thumb  [h = 0.9 * min(sd, IQR/1.34) * n^(-1/5)]",
        "\n               Reference : Silverman (1986); Scott (1992)",
        " Multivariate Density Estimation, Wiley",
        "\n               sd(x) = ", round(s, 6L),
        ",  IQR(x) = ", round(iqr, 6L),
        ",  scale = ", round(scale, 6L),
        ",  n = ", n,
        ",  h = ", round(h, 6L)
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 3. IQR-based rule
  ##    h = IQR(x) * n^{-1/5}
  ## ------------------------------------------------------------------
  if (method == "iqr") {
    if (!is.finite(iqr) || iqr <= 0)
      stop("IQR of the quantile regressor is zero or non-finite. ",
           "Use method = 'silverman' or 'robust' instead.")
    h <- iqr * n^(-0.2)
    return(list(
      h       = h,
      method  = "iqr",
      details = paste0(
        "IQR-based rule  [h = IQR(x) * n^(-1/5)]",
        "\n               IQR(x) = ", round(iqr, 6L),
        ",  n = ", n,
        ",  h = ", round(h, 6L)
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 4. Fast LOO cross-validation (no additional packages required)
  ##    Renamed from "cv" to "cv_fast"; minimises quantile check-loss
  ##    at tau = 0.5 over a 15-point grid in [0.5*h_sj, 2*h_sj].
  ##    Source: Yu and Jones (1998, JASA 93, 228-237)
  ## ------------------------------------------------------------------
  if (method == "cv_fast") {
    h_silver  <- 1.06 * s * n^(-0.2)
    h_grid    <- seq(0.5 * h_silver, 2.0 * h_silver, length.out = 15L)
    tau_cv    <- 0.5
    x_clean   <- x[is.finite(x)]
    n_cv      <- length(x_clean)

    cv_scores <- vapply(h_grid, function(h_cand) {
      errors <- vapply(seq_len(n_cv), function(i) {
        x_loo <- x_clean[-i]
        wx    <- switch(kernel,
                        gaussian = stats::dnorm((x_loo - x_clean[i]) / h_cand))
        eff_n <- sum(wx[is.finite(wx)], na.rm = TRUE)
        if (!is.finite(eff_n) || eff_n < 1e-10) return(NA_real_)
        ord   <- order(x_loo)
        w_ord <- wx[ord]
        w_cum <- cumsum(w_ord) / sum(w_ord, na.rm = TRUE)
        pred  <- x_loo[ord][which(w_cum >= tau_cv)[1L]]
        if (is.na(pred)) return(NA_real_)
        u <- x_clean[i] - pred
        ifelse(u >= 0, tau_cv * u, (tau_cv - 1) * u)
      }, numeric(1L))
      mean(errors, na.rm = TRUE)
    }, numeric(1L))

    if (all(is.na(cv_scores)))
      stop("cv_fast: cross-validation failed for all candidate bandwidths. ",
           "Use method = 'silverman' or 'robust' instead.")

    h_opt <- h_grid[which.min(cv_scores)]

    return(list(
      h       = h_opt,
      method  = "cv_fast",
      details = paste0(
        "Fast LOO cross-validation  [15-point grid, quantile loss at tau=0.5]",
        "\n               Reference : Yu and Jones (1998) Local linear quantile",
        " regression, JASA 93(441), 228-237",
        "\n               Search range: [", round(min(h_grid), 6L),
        ", ", round(max(h_grid), 6L), "]",
        "\n               Optimal h = ", round(h_opt, 6L),
        "  (minimises LOO quantile-loss at tau = 0.5)"
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 5. Sheather-Jones solve-the-equation plug-in  (h = "sj")
  ##    Uses base R bw.SJ() -- no additional packages required.
  ##    Solves iteratively for h satisfying:
  ##      h = [ R(K) / (n * psi4_hat(g(h))) ]^{1/5}
  ##    where psi4_hat is a pilot estimate of int[f''(x)]^2 dx and
  ##    g(h) is a pilot bandwidth for that estimate.
  ##    Sources: Sheather and Jones (1991, JRSS-B 53, 683-690)
  ##             Jones, Marron and Sheather (1996, JASA 91, 401-407)
  ## ------------------------------------------------------------------
  if (method == "sj") {
    h <- tryCatch(
      stats::bw.SJ(x[is.finite(x)]),
      error = function(e) {
        warning("SJ bandwidth selection failed (", conditionMessage(e),
                "). Falling back to Silverman rule.")
        1.06 * s * n^(-0.2)
      }
    )
    return(list(
      h       = h,
      method  = "sj",
      details = paste0(
        "Sheather-Jones solve-the-equation plug-in  [bw.SJ(x)]",
        "\n               References:",
        "\n                 Sheather & Jones (1991) A reliable data-based",
        " bandwidth selection method for kernel density estimation.",
        "\n                 J. Royal Statistical Society B, 53(3), 683-690.",
        "\n                 Jones, Marron & Sheather (1996) A brief survey of",
        " bandwidth selection for density estimation.",
        "\n                 J. American Statistical Association, 91(433),",
        " 401-407.",
        "\n               n = ", n,
        ",  h = ", round(h, 6L)
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 6. Ruppert-Sheather-Wand direct plug-in  (h = "yj")
  ##    Uses KernSmooth::dpill(x, y) -- KernSmooth is a recommended
  ##    package that ships with standard R installations.
  ##    Minimises the MISE of a local linear regression of y on x:
  ##      h_opt = [ sigma^2 int K^2 / (n f_X(x) [m''(x)]^2 mu2^2) ]^{1/5}
  ##    m''(x) is estimated from pilot fits using Mallows Cp block size.
  ##    This is the only method here that uses BOTH x and y and is
  ##    tuned to regression (not density) loss -- making it the most
  ##    structurally appropriate selector for QQR.
  ##    Source: Ruppert, Sheather & Wand (1995, JASA 90, 1257-1270)
  ## ------------------------------------------------------------------
  if (method == "yj") {
    if (is.null(y))
      stop("method = 'yj' requires the outcome vector y. ",
           "Make sure y is passed to .qqr_select_h().")

    if (!requireNamespace("KernSmooth", quietly = TRUE))
      stop("method = 'yj' requires the 'KernSmooth' package. ",
           "It is a recommended R package; install it with: ",
           "install.packages('KernSmooth')")

    x_fin <- x[is.finite(x) & is.finite(y)]
    y_fin <- y[is.finite(x) & is.finite(y)]
    n_fin <- length(x_fin)

    if (n_fin < 10L)
      stop("Too few finite observations for the YJ plug-in selector.")

    h <- tryCatch({
      tn <- n_fin - 1L
      KernSmooth::dpill(x_fin, y_fin, gridsize = tn)
    }, error = function(e) {
      warning("YJ (dpill) bandwidth selection failed (", conditionMessage(e),
              "). Falling back to Silverman rule.")
      1.06 * s * n^(-0.2)
    })

    if (!is.finite(h) || h <= 0) {
      warning("dpill returned a non-positive bandwidth. ",
              "Falling back to Silverman rule.")
      h <- 1.06 * s * n^(-0.2)
    }

    return(list(
      h       = h,
      method  = "yj",
      details = paste0(
        "Ruppert-Sheather-Wand direct plug-in  [KernSmooth::dpill(x, y)]",
        "\n               Reference: Ruppert, Sheather & Wand (1995) An",
        " effective bandwidth selector for local least squares regression.",
        "\n                 Journal of the American Statistical Association,",
        " 90(432), 1257-1270.",
        "\n               Note: Uses both x (quantile regressor) and y",
        " (outcome). Minimises regression MISE -- most appropriate for QQR.",
        "\n               n = ", n_fin,
        ",  h = ", round(h, 6L)
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 7. np full cross-validation  (h = "npcv")
  ##    Uses np::npregbw(xdat=x, ydat=y, regtype="ll")
  ##    Minimises the least-squares cross-validation criterion:
  ##      h_CV = argmin_h (1/n) sum_t [y_t - m_hat_{-t}(x_t; h)]^2
  ##    where m_hat_{-t} is the leave-one-out local-linear fit.
  ##    Supports mixed continuous and categorical data via generalised
  ##    product kernels (Li & Racine 2004).
  ##    Sources: Hayfield & Racine (2008, JSS 27(5), 1-32)
  ##             Li & Racine (2004, Econometrica 72, 1723-1762)
  ## ------------------------------------------------------------------
  if (method == "npcv") {
    if (is.null(y))
      stop("method = 'npcv' requires the outcome vector y.")

    if (!requireNamespace("np", quietly = TRUE))
      stop("method = 'npcv' requires the 'np' package. ",
           "Install it with: install.packages('np')")

    x_fin <- x[is.finite(x) & is.finite(y)]
    y_fin <- y[is.finite(x) & is.finite(y)]

    h <- tryCatch({
      npregbw_fn <- get("npregbw", envir = asNamespace("np"))
      bw_obj <- npregbw_fn(xdat = x_fin, ydat = y_fin, regtype = "ll")
      bw_obj$bw
    }, error = function(e) {
      warning("npcv (npregbw) bandwidth selection failed (",
              conditionMessage(e),
              "). Falling back to Silverman rule.")
      1.06 * s * n^(-0.2)
    })

    if (!is.finite(h) || h <= 0) {
      warning("npregbw returned a non-positive bandwidth. ",
              "Falling back to Silverman rule.")
      h <- 1.06 * s * n^(-0.2)
    }

    return(list(
      h       = h,
      method  = "npcv",
      details = paste0(
        "np full LS cross-validation  [np::npregbw(xdat, ydat, regtype='ll')]",
        "\n               References:",
        "\n                 Hayfield & Racine (2008) Nonparametric econometrics:",
        " The np package. J. Statistical Software, 27(5), 1-32.",
        "\n                 Li & Racine (2004) Cross-validated local linear",
        " nonparametric regression. Econometrica, 72(6), 1723-1762.",
        "\n               Note: Minimises LOO least-squares loss using",
        " generalised product kernels for mixed-type data.",
        "\n               n = ", length(x_fin),
        ",  h = ", round(h, 6L)
      )
    ))
  }

  ## ------------------------------------------------------------------
  ## 8. np smooth-coefficient cross-validation  (h = "sc")
  ##    Uses np::npscoefbw(xdat=x, ydat=y)
  ##    Selects bandwidth for a smooth varying-coefficient model:
  ##      Y_i = W_i' gamma(Z_i) + u_i
  ##    where gamma(Z_i) varies smoothly with the conditioning variable.
  ##    This is the closest structural analogue to QQR in least-squares
  ##    form: QQR estimates beta(tau, z_i), i.e. a coefficient that
  ##    varies with z_i -- the same structure as a varying-coefficient
  ##    model, but within a quantile rather than mean regression.
  ##    Sources: Li, Racine & Wooldridge (2002, JASA 97, 1129-1139)
  ##             Li & Racine (2010, Econometric Theory 26, 1-31)
  ##             Hayfield & Racine (2008, JSS 27(5), 1-32)
  ## ------------------------------------------------------------------
  if (method == "sc") {
    if (is.null(y))
      stop("method = 'sc' requires the outcome vector y.")

    if (!requireNamespace("np", quietly = TRUE))
      stop("method = 'sc' requires the 'np' package. ",
           "Install it with: install.packages('np')")

    x_fin <- x[is.finite(x) & is.finite(y)]
    y_fin <- y[is.finite(x) & is.finite(y)]

    h <- tryCatch({
      npscoefbw_fn <- get("npscoefbw", envir = asNamespace("np"))
      bw_obj <- npscoefbw_fn(xdat = x_fin, ydat = y_fin)
      bw_obj$bw
    }, error = function(e) {
      warning("sc (npscoefbw) bandwidth selection failed (",
              conditionMessage(e),
              "). Falling back to Silverman rule.")
      1.06 * s * n^(-0.2)
    })

    ## npscoefbw may return a vector (one bw per variable); take the mean
    if (length(h) > 1L) h <- mean(h, na.rm = TRUE)

    if (!is.finite(h) || h <= 0) {
      warning("npscoefbw returned a non-positive bandwidth. ",
              "Falling back to Silverman rule.")
      h <- 1.06 * s * n^(-0.2)
    }

    return(list(
      h       = h,
      method  = "sc",
      details = paste0(
        "np smooth-coefficient CV  [np::npscoefbw(xdat, ydat)]",
        "\n               References:",
        "\n                 Li, Racine & Wooldridge (2002) Semiparametric",
        " estimation of partially linear models for dependent data with",
        " generated regressors. Econometrica 70(5), 2123-2163.",
        "\n                 Li & Racine (2010) Smooth varying-coefficient",
        " estimation and inference for qualitative and quantitative data.",
        "\n                 Econometric Theory, 26(6), 1-31.",
        "\n                 Hayfield & Racine (2008) J. Statistical Software",
        " 27(5), 1-32.",
        "\n               Note: Bandwidth for Y_i = W_i' gamma(Z_i) + u_i --",
        " closest structural analogue to QQR in least-squares form.",
        "\n               n = ", length(x_fin),
        ",  h = ", round(h, 6L)
      )
    ))
  }

  stop(
    "Unknown h method: '", method, "'. ",
    "Use one of: 'silverman', 'robust', 'iqr', 'cv_fast', ",
    "'sj', 'yj', 'npcv', 'sc', or a positive number."
  )
}


## =========================================================================
## Main fitting function
## =========================================================================

#' Quantile-on-Quantile Regression
#'
#' @description
#' Fits a quantile-on-quantile regression (QQR) model by running a sequence
#' of local linear quantile regressions of a dependent variable on a
#' selected quantile regressor and optional control variables. The local
#' regressions are evaluated at a grid of locations of the quantile
#' regressor and at a set of conditional quantile levels \eqn{\tau},
#' producing a full coefficient surface indexed by both dimensions.
#'
#' @details
#' For each grid point \eqn{z_0} in the support of the quantile regressor
#' \eqn{X}, the function constructs a locally weighted sample using a
#' kernel centred at \eqn{z_0} and fits a weighted quantile regression of
#' the form
#' \deqn{
#'   Y = \beta_0(z_0, \tau) + \beta_1(z_0, \tau)(X - z_0)
#'       + \gamma(z_0, \tau)^\top W + u,
#' }
#' for each \eqn{\tau} in \code{taus}, where \eqn{W} denotes the optional
#' control variables. This produces a surface of coefficients indexed by
#' grid location and quantile level. P-values are obtained from
#' \code{\link[quantreg]{summary.rq}} using the chosen standard error type.
#'
#' \strong{Bandwidth selection.}
#' The bandwidth \code{h} governs the width of the kernel around each grid
#' point and is the most important tuning parameter in QQR. To promote
#' reproducibility and discourage post-hoc specification searching, the
#' default is \code{h = "silverman"}, which computes the bandwidth
#' automatically from the data. Eight automatic methods are available.
#' The first three require no additional packages; the next three
#' (\code{"sj"} and \code{"yj"}) use only recommended packages that ship
#' with standard R; the last two (\code{"npcv"} and \code{"sc"}) require
#' the \pkg{np} package.
#'
#' \describe{
#'   \item{\code{"silverman"}}{
#'     \strong{Silverman's rule of thumb.}
#'     \eqn{h = 1.06 \cdot \hat\sigma_x \cdot n^{-1/5}}.
#'     Optimal under normality. Fast, requires only the sample standard
#'     deviation, and is widely reported in the QQR literature.
#'     \emph{Silverman (1986), Density Estimation for Statistics and
#'     Data Analysis, Chapman & Hall.}}
#'   \item{\code{"robust"}}{
#'     \strong{Robust rule of thumb.}
#'     \eqn{h = 0.9 \cdot \min(\hat\sigma_x,\, \widehat{IQR}_x/1.34)
#'     \cdot n^{-1/5}}.
#'     The IQR correction protects against heavy tails and outliers.
#'     Preferred for financial returns, exchange rates, commodity prices,
#'     and other macroeconomic series.
#'     \emph{Silverman (1986); Scott (1992), Multivariate Density
#'     Estimation, Wiley.}}
#'   \item{\code{"iqr"}}{
#'     \strong{IQR-based rule.}
#'     \eqn{h = \widehat{IQR}_x \cdot n^{-1/5}}.
#'     Appropriate when the spread of the quantile regressor is better
#'     characterised by its interquartile range.}
#'   \item{\code{"cv_fast"}}{
#'     \strong{Fast leave-one-out cross-validation.}
#'     Minimises the LOO quantile check-loss at \eqn{\tau = 0.5} over a
#'     15-point grid spanning
#'     \eqn{[0.5\,h_{\rm Silverman},\, 2\,h_{\rm Silverman}]}.
#'     No additional packages required. Computationally heavier than the
#'     rule-of-thumb methods but fully data-driven.
#'     \emph{Yu and Jones (1998), JASA 93(441), 228-237.}}
#'   \item{\code{"sj"}}{
#'     \strong{Sheather-Jones solve-the-equation plug-in.}
#'     Uses \code{base::bw.SJ(x)} to solve iteratively for
#'     \eqn{h = [R(K) / (n\,\hat\psi_4(g(h)))]^{1/5}},
#'     where \eqn{\hat\psi_4} is a pilot estimate of
#'     \eqn{\int [f''(x)]^2\,dx} and \eqn{g(h)} is a pilot bandwidth.
#'     Widely regarded as the best density-based selector available.
#'     No additional packages required.
#'     \emph{Sheather and Jones (1991), JRSS-B 53(3), 683-690.}}
#'   \item{\code{"yj"}}{
#'     \strong{Ruppert-Sheather-Wand direct plug-in for local regression.}
#'     Uses \code{KernSmooth::dpill(x, y)} to minimise the MISE of a
#'     local linear regression of \eqn{Y} on \eqn{X}:
#'     \eqn{h = [\hat\sigma^2_\varepsilon \int K^2 /
#'     (n\,\hat f_X\,[m''(x)]^2\,\mu_2^2)]^{1/5}}.
#'     The only selector here that uses both \eqn{X} and \eqn{Y} and is
#'     tuned to regression loss, making it the most structurally
#'     appropriate for QQR. Requires \pkg{KernSmooth} (a recommended
#'     package that ships with standard R).
#'     \emph{Ruppert, Sheather and Wand (1995), JASA 90(432), 1257-1270.}}
#'   \item{\code{"npcv"}}{
#'     \strong{np full least-squares cross-validation.}
#'     Uses \code{np::npregbw(xdat, ydat, regtype = "ll")} to minimise
#'     the LOO squared-error loss via generalised product kernels
#'     (supports mixed continuous and categorical data).
#'     The most assumption-free selector; computationally heaviest.
#'     Requires the \pkg{np} package (listed in \code{Suggests}).
#'     \emph{Hayfield and Racine (2008), J. Statistical Software
#'     27(5), 1-32.}}
#'   \item{\code{"sc"}}{
#'     \strong{np smooth-coefficient cross-validation.}
#'     Uses \code{np::npscoefbw(xdat, ydat)} to select bandwidth for
#'     the varying-coefficient model
#'     \eqn{Y_i = W_i^\top \gamma(Z_i) + u_i},
#'     which is the closest structural analogue to QQR in least-squares
#'     form. Requires the \pkg{np} package (listed in \code{Suggests}).
#'     \emph{Li and Racine (2010), Econometric Theory 26, 1-31;
#'     Hayfield and Racine (2008), JSS 27(5), 1-32.}}
#' }
#'
#' If a positive numeric value is supplied for \code{h}, it is used
#' directly and a warning is issued reminding the researcher that manual
#' bandwidth selection may introduce specification bias.
#'
#' The selected bandwidth, the method label, and a full description
#' (including the formula, numerical inputs, and final value) are stored
#' in the returned object as \code{bandwidth}, \code{h_method}, and
#' \code{h_details}, and are printed by \code{\link{print.qqr}} and
#' \code{\link{summary.qqr}}.
#'
#' @param formula A model formula of the form \code{y ~ x} or
#'   \code{y ~ x + controls}, where \code{y} is numeric.
#' @param data A data frame containing all variables in \code{formula}.
#' @param qvar Optional character string naming the quantile regressor.
#'   Must be one of the right-hand side variables. If \code{NULL} and
#'   the formula contains exactly one regressor, that variable is used.
#'   An error is raised if there are multiple regressors and \code{qvar}
#'   is not supplied.
#' @param taus Numeric vector of conditional quantile levels in \eqn{(0,1)}
#'   at which the local regressions are evaluated. Default is
#'   \code{seq(0.05, 0.95, by = 0.05)}.
#' @param m Positive integer giving the number of equally spaced grid
#'   points over the observed range of the quantile regressor. Default
#'   is \code{19}.
#' @param h Bandwidth for the kernel weighting. Either a character string
#'   naming an automatic selection method (\code{"silverman"},
#'   \code{"robust"}, \code{"iqr"}, \code{"cv_fast"}, \code{"sj"},
#'   \code{"yj"}, \code{"npcv"}, or \code{"sc"}), or a single positive
#'   numeric value. Default is \code{"silverman"}. See Details.
#' @param kernel Character string specifying the kernel function. Currently
#'   only \code{"gaussian"} is supported.
#' @param se_type Character string selecting the standard error type passed
#'   to \code{\link[quantreg]{summary.rq}} for p-value computation.
#'   Options are \code{"nid"} (normal, independent design; default),
#'   \code{"ker"} (kernel), \code{"boot"} (bootstrap), or \code{"none"}
#'   (skip p-value computation entirely). When p-values cannot be computed
#'   at a grid-tau combination, the corresponding entry is set to \code{NA}.
#' @param na.action Function for handling missing values in the model
#'   frame. Default is \code{stats::na.omit}.
#' @param ... Additional arguments forwarded to
#'   \code{\link[quantreg]{summary.rq}} when \code{se_type = "boot"}
#'   (e.g. \code{R} for the number of bootstrap replications).
#'
#' @return
#' An object of class \code{"qqr"} containing:
#' \describe{
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{dependent}}{Name of the dependent variable.}
#'   \item{\code{independent}}{Name of the quantile regressor.}
#'   \item{\code{controls}}{Character vector of control variable names
#'     (length zero if none).}
#'   \item{\code{grid}}{Numeric vector of length \code{m}: the grid of
#'     quantile-regressor values used for local fitting.}
#'   \item{\code{taus}}{Numeric vector of quantile levels.}
#'   \item{\code{coef_surfaces}}{Named list of \eqn{m \times K} matrices,
#'     one per coefficient (intercept, quantile regressor, controls).}
#'   \item{\code{pval_surfaces}}{Named list of \eqn{m \times K} p-value
#'     matrices with the same structure.}
#'   \item{\code{qcoef}}{The \eqn{m \times K} coefficient matrix for the
#'     quantile regressor.}
#'   \item{\code{intercept}}{Logical; \code{TRUE} if the formula includes
#'     an intercept.}
#'   \item{\code{coef_names}}{Character vector of coefficient names.}
#'   \item{\code{bandwidth}}{The numeric bandwidth \code{h} used.}
#'   \item{\code{h_method}}{Character label of the selection method:
#'     \code{"silverman"}, \code{"robust"}, \code{"iqr"},
#'     \code{"cv_fast"}, \code{"sj"}, \code{"yj"}, \code{"npcv"},
#'     \code{"sc"}, or \code{"manual"}.}
#'   \item{\code{h_details}}{Full description of how \code{h} was
#'     selected, including the formula, numerical inputs, and references,
#'     suitable for inclusion in a paper's methods section.}
#'   \item{\code{kernel}}{Kernel name.}
#'   \item{\code{se_type}}{Standard error type used.}
#'   \item{\code{effective_n}}{Matrix of effective local sample sizes.}
#'   \item{\code{n}}{Number of observations used in the fit.}
#'   \item{\code{data}}{The model frame used for fitting.}
#' }
#'
#' @seealso
#' \code{\link{summary.qqr}}, \code{\link{print.qqr}},
#' \code{\link{plot.qqr}}, \code{\link{Int.plot}},
#' \code{\link{surface3d}}, \code{\link{panel_surface}},
#' \code{\link{qqr_surface}}, \code{\link{rowave.qqr}},
#' \code{\link{compare.qqr}}, \code{\link{excel.qqr}},
#' \code{\link{word.qqr}}
#'
#' @references
#' Hayfield, T. and Racine, J.S. (2008) Nonparametric econometrics: The np
#' package. \emph{Journal of Statistical Software}, 27(5), 1--32.
#'
#' Li, Q. and Racine, J.S. (2010) Smooth varying-coefficient estimation and
#' inference for qualitative and quantitative data.
#' \emph{Econometric Theory}, 26(6), 1--31.
#'
#' Ruppert, D., Sheather, S.J. and Wand, M.P. (1995) An effective bandwidth
#' selector for local least squares regression.
#' \emph{Journal of the American Statistical Association}, 90(432),
#' 1257--1270.
#'
#' Scott, D.W. (1992) \emph{Multivariate Density Estimation: Theory,
#' Practice, and Visualization}. New York: Wiley.
#'
#' Sheather, S.J. and Jones, M.C. (1991) A reliable data-based bandwidth
#' selection method for kernel density estimation.
#' \emph{Journal of the Royal Statistical Society Series B}, 53(3),
#' 683--690.
#'
#' Silverman, B.W. (1986) \emph{Density Estimation for Statistics and Data
#' Analysis}. London: Chapman & Hall.
#'
#' Sim, N. and Zhou, H. (2015) Oil prices, US stock return and the
#' dependence between their quantiles.
#' \emph{Journal of Banking & Finance}, 55, 1--8.
#'
#' Yu, K. and Jones, M.C. (1998) Local linear quantile regression.
#' \emph{Journal of the American Statistical Association}, 93(441),
#' 228--237.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' n   <- 300
#' x   <- rnorm(n)
#' u   <- rnorm(n)
#' y   <- 0.5 + x * stats::qnorm(stats::runif(n)) + 0.3 * u
#' dat <- data.frame(y = y, x = x, u = u)
#'
#' ## Default: Silverman rule (no extra packages)
#' fit1 <- qqr(y ~ x + u, data = dat, qvar = "x",
#'             taus = seq(0.05, 0.95, by = 0.05), m = 19)
#' print(fit1)
#' summary(fit1)
#'
#' ## Robust rule (heavy-tailed regressors)
#' fit2 <- qqr(y ~ x + u, data = dat, qvar = "x", h = "robust")
#'
#' ## Sheather-Jones plug-in (base R, no extra packages)
#' fit3 <- qqr(y ~ x + u, data = dat, qvar = "x", h = "sj")
#'
#' ## Ruppert-Sheather-Wand plug-in (uses x and y; KernSmooth)
#' fit4 <- qqr(y ~ x + u, data = dat, qvar = "x", h = "yj")
#'
#' ## Fast LOO cross-validation (no extra packages)
#' fit5 <- qqr(y ~ x + u, data = dat, qvar = "x", h = "cv_fast")
#'
#' ## np full LS cross-validation (requires np package)
#' fit6 <- qqr(y ~ x + u, data = dat, qvar = "x", h = "npcv")
#'
#' ## np smooth-coefficient CV (requires np package)
#' fit7 <- qqr(y ~ x + u, data = dat, qvar = "x", h = "sc")
#'
#' ## Manual bandwidth (triggers a reproducibility warning)
#' fit8 <- qqr(y ~ x + u, data = dat, qvar = "x", h = 0.10)
#' }
#'
#' @export
qqr <- function(formula,
                data,
                qvar      = NULL,
                taus      = seq(0.05, 0.95, by = 0.05),
                m         = 19L,
                h         = "silverman",
                kernel    = c("gaussian"),
                se_type   = c("nid", "ker", "boot", "none"),
                ...,
                na.action = stats::na.omit) {

  kernel  <- match.arg(kernel)
  se_type <- match.arg(se_type)

  ## -----------------------------------------------------------------------
  ## 1. Parse formula and data
  ## -----------------------------------------------------------------------
  mf <- stats::model.frame(formula, data = data, na.action = na.action)
  if (ncol(mf) < 2L)
    stop("Formula must have at least one independent variable.")

  terms_obj     <- stats::terms(formula)
  has_intercept <- as.logical(attr(terms_obj, "intercept"))

  dep_name    <- colnames(mf)[1L]
  indep_names <- colnames(mf)[-1L]

  ## -----------------------------------------------------------------------
  ## 2. Identify quantile regressor
  ## -----------------------------------------------------------------------
  if (length(indep_names) == 1L) {
    if (is.null(qvar)) {
      qvar <- indep_names[1L]
    } else if (!(qvar %in% indep_names)) {
      stop("qvar '", qvar, "' is not an independent variable in the formula.")
    }
  } else {
    if (is.null(qvar))
      stop("Formula has multiple regressors. ",
           "Set qvar to one of: ",
           paste(indep_names, collapse = ", "))
    if (!(qvar %in% indep_names))
      stop("qvar '", qvar, "' is not an independent variable in the formula.")
  }

  control_names <- setdiff(indep_names, qvar)

  y        <- mf[[dep_name]]
  x        <- mf[[qvar]]
  controls <- if (length(control_names)) mf[control_names] else NULL

  if (!is.numeric(y) || !is.numeric(x))
    stop("Both dependent variable and quantile regressor must be numeric.")
  if (length(y) != length(x))
    stop("Dependent and quantile regressor have different lengths.")

  ## -----------------------------------------------------------------------
  ## 3. Bandwidth selection
  ## -----------------------------------------------------------------------
  h_spec <- h   ## preserve original user input for storage

  ## Methods that require y to be passed explicitly
  y_methods <- c("yj", "npcv", "sc")

  if (is.character(h)) {

    h_result <- .qqr_select_h(
      x       = x,
      y       = if (h %in% y_methods) y else NULL,
      method  = h,
      taus    = taus,
      m       = m,
      kernel  = kernel,
      se_type = se_type
    )
    h         <- h_result$h
    h_method  <- h_result$method
    h_details <- h_result$details

  } else if (is.numeric(h) && length(h) == 1L && is.finite(h) && h > 0) {

    h_method  <- "manual"
    h_details <- paste0(
      "Manually specified  [h = ", h, "]",
      "\n               WARNING: Manual bandwidth selection may ",
      "introduce specification bias.",
      "\n               For reproducible research, consider h = ",
      "'silverman', 'robust', 'iqr', 'cv_fast', 'sj', 'yj', 'npcv', or 'sc'."
    )

    warning(
      "Bandwidth h was set manually to ", h, ". ",
      "Manual bandwidth selection may affect results and reduce ",
      "reproducibility. ",
      "For automatic selection, use h = 'silverman' (default), 'robust', ",
      "'iqr', 'cv_fast', 'sj', 'yj', 'npcv', or 'sc'.",
      call. = FALSE
    )

  } else {
    stop(
      "h must be a positive number or one of: ",
      "'silverman', 'robust', 'iqr', 'cv_fast', 'sj', 'yj', 'npcv', 'sc'. ",
      "Got: ", deparse(h_spec)
    )
  }

  ## -----------------------------------------------------------------------
  ## 4. Build grid
  ## -----------------------------------------------------------------------
  grid <- seq(min(x, na.rm = TRUE),
              max(x, na.rm = TRUE),
              length.out = m)

  ## -----------------------------------------------------------------------
  ## 5. Coefficient names
  ## -----------------------------------------------------------------------
  coef_names <- c(
    if (has_intercept) "(Intercept)" else NULL,
    qvar,
    control_names
  )

  n_coef <- length(coef_names)
  n_grid <- length(grid)
  n_tau  <- length(taus)

  ## -----------------------------------------------------------------------
  ## 6. Initialise surface storage
  ## -----------------------------------------------------------------------
  coef_surfaces <- stats::setNames(
    lapply(seq_len(n_coef), function(i)
      matrix(NA_real_, nrow = n_grid, ncol = n_tau)),
    coef_names
  )

  pval_surfaces <- stats::setNames(
    lapply(seq_len(n_coef), function(i)
      matrix(NA_real_, nrow = n_grid, ncol = n_tau)),
    coef_names
  )

  effective_n <- matrix(NA_real_, nrow = n_grid, ncol = n_tau)

  ## -----------------------------------------------------------------------
  ## 7. Fit local quantile regressions
  ## -----------------------------------------------------------------------
  for (j in seq_along(taus)) {
    tau_j <- taus[j]

    for (i in seq_len(n_grid)) {
      res <- .qqr_lprq_point(
        y          = y,
        x          = x,
        controls   = controls,
        grid_value = grid[i],
        h          = h,
        tau        = tau_j,
        kernel     = kernel,
        intercept  = has_intercept,
        qvar       = qvar,
        coef_names = coef_names,
        se_type    = se_type,
        ...
      )

      if (is.null(res)) next

      for (k in seq_len(n_coef)) {
        coef_surfaces[[k]][i, j] <- res$coef[k]
        pval_surfaces[[k]][i, j] <- res$pval[k]
      }
      effective_n[i, j] <- res$effective_n
    }
  }

  ## -----------------------------------------------------------------------
  ## 8. Assemble and return
  ## -----------------------------------------------------------------------
  structure(
    list(
      call          = match.call(),
      dependent     = dep_name,
      independent   = qvar,
      controls      = control_names,
      grid          = grid,
      taus          = taus,
      coef_surfaces = coef_surfaces,
      pval_surfaces = pval_surfaces,
      qcoef         = coef_surfaces[[qvar]],
      intercept     = has_intercept,
      coef_names    = coef_names,
      bandwidth     = h,
      h_method      = h_method,
      h_details     = h_details,
      kernel        = kernel,
      se_type       = se_type,
      effective_n   = effective_n,
      n             = nrow(mf),
      data          = mf
    ),
    class = "qqr"
  )
}


## =========================================================================
## Internal helpers
## =========================================================================

## -------------------------------------------------------------------------
## Map user-facing coefficient name to internal storage key
## @keywords internal
## -------------------------------------------------------------------------
.qqr_match_coef <- function(x, coef) {
  coef_names <- names(x$coef_surfaces)

  if (coef %in% coef_names) return(coef)

  if (tolower(coef) %in% c("intercept", "(intercept)")) {
    if ("(Intercept)" %in% coef_names) return("(Intercept)")
  }

  if (coef == x$independent) {
    if (x$independent %in% coef_names) return(x$independent)
    if ("z" %in% coef_names)           return("z")   ## backwards compat
  }

  available <- c(
    if ("(Intercept)" %in% coef_names) "intercept" else NULL,
    x$independent,
    x$controls
  )
  stop("Coefficient '", coef, "' not recognised. Available: ",
       paste(available, collapse = ", "))
}


## -------------------------------------------------------------------------
## One local kernel-weighted quantile regression
## @keywords internal
## -------------------------------------------------------------------------
.qqr_lprq_point <- function(y, x, controls, grid_value, h, tau,
                             kernel     = c("gaussian"),
                             intercept  = TRUE,
                             qvar,
                             coef_names,
                             se_type    = c("nid", "ker", "boot", "none"),
                             ...) {

  kernel  <- match.arg(kernel)
  se_type <- match.arg(se_type)

  z  <- x - grid_value
  wx <- switch(kernel,
               gaussian = stats::dnorm(z / h))

  effective_n <- sum(wx[is.finite(wx)], na.rm = TRUE)

  na_result <- function() {
    v <- stats::setNames(rep(NA_real_, length(coef_names)), coef_names)
    list(effective_n = effective_n, coef = v, pval = v)
  }

  if (!any(is.finite(wx)) || effective_n < .Machine$double.eps)
    return(na_result())

  ## Build local design frame
  if (is.null(controls) || ncol(controls) == 0L) {
    df_fit         <- data.frame(y = y)
    df_fit[[qvar]] <- z
    rhs            <- qvar
  } else {
    df_fit         <- data.frame(y = y)
    df_fit[[qvar]] <- z
    df_fit         <- cbind(df_fit, controls)
    rhs            <- paste(c(qvar, names(controls)), collapse = " + ")
  }

  form <- stats::as.formula(
    paste("y ~", if (!intercept) "0 +" else "", rhs)
  )

  fit <- tryCatch(
    quantreg::rq(form, data = df_fit, weights = wx, tau = tau, ci = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) return(na_result())

  ## Coefficients
  cf      <- stats::coef(fit)
  cf_full <- stats::setNames(rep(NA_real_, length(coef_names)), coef_names)
  cf_full[names(cf)] <- cf

  ## P-values
  pv_full <- stats::setNames(rep(NA_real_, length(coef_names)), coef_names)

  if (!identical(se_type, "none")) {
    summ <- tryCatch(
      quantreg::summary.rq(fit, se = se_type, ...),
      error = function(e) NULL
    )
    if (!is.null(summ) && !is.null(summ$coefficients)) {
      cm   <- summ$coefficients
      pcol <- grep("^Pr\\(", colnames(cm))
      if (length(pcol) >= 1L)
        pv_full[rownames(cm)] <- cm[, pcol[1L]]
    }
  }

  list(effective_n = effective_n, coef = cf_full, pval = pv_full)
}
