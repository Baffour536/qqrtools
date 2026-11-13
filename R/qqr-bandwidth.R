## =========================================================================
## qqr-bandwidth.R
## Automatic bandwidth selectors used by qqr().
## =========================================================================

## Silverman rule of thumb, used as a fallback by several selectors.
.qqr_h_silverman <- function(x) {
  1.06 * stats::sd(x) * length(x)^(-0.2)
}

## -------------------------------------------------------------------------
## Internal: bandwidth selector
## -------------------------------------------------------------------------
## Computes the kernel bandwidth h from the smoothing variable x (either the
## quantile regressor itself or its empirical CDF ranks, depending on
## qqr(kernel_space = )) and, where needed, the outcome vector y.
##
## method = "silverman" : h = 1.06 * sd(x) * n^(-1/5)
## method = "robust"    : h = 0.9 * min(sd(x), IQR(x)/1.34) * n^(-1/5)
## method = "iqr"       : h = IQR(x) * n^(-1/5)
## method = "cv_fast"   : leave-one-out check-loss CV of a kernel-weighted
##                        conditional median of y given x (15-point grid)
## method = "sj"        : Sheather-Jones solve-the-equation plug-in (bw.SJ)
## method = "yj"        : Ruppert-Sheather-Wand direct plug-in
##                        (KernSmooth::dpill; uses x and y)
## method = "npcv"      : np::npregbw least-squares CV (uses x and y)
## method = "sc"        : np::npscoefbw smooth-coefficient CV (uses x and y)
##
## Returns list(h, method, details).
## @keywords internal
## @noRd
.qqr_select_h <- function(x, y = NULL, method, kernel = "gaussian") {

  ok <- is.finite(x) & (if (is.null(y)) TRUE else is.finite(y))
  x  <- x[ok]
  if (!is.null(y)) y <- y[ok]
  n   <- length(x)

  if (n < 4L)
    stop("Too few observations (", n, ") to select a bandwidth automatically.",
         call. = FALSE)

  s   <- stats::sd(x)
  iqr <- stats::IQR(x)
  if (!is.finite(s) || s <= 0)
    stop("The smoothing variable has zero variance; a bandwidth cannot be ",
         "selected. Check that the quantile regressor is not constant.",
         call. = FALSE)

  h_silver <- .qqr_h_silverman(x)
  fmt <- function(v) format(signif(v, 6L))
  needs_y <- function() {
    if (is.null(y))
      stop("h = '", method, "' requires the outcome vector y.", call. = FALSE)
  }
  fallback <- function(what, e) {
    warning(what, " bandwidth selection failed (", conditionMessage(e),
            "). Falling back to Silverman's rule.", call. = FALSE)
    h_silver
  }

  out <- switch(
    method,

    ## 1. Silverman (1986, eq. 3.31 with the normal-reference constant)
    silverman = list(
      h = h_silver,
      details = paste0(
        "Silverman's rule of thumb [h = 1.06 * sd * n^(-1/5)]; ",
        "Silverman (1986). sd = ", fmt(s), ", n = ", n,
        ", h = ", fmt(h_silver))
    ),

    ## 2. Robust rule of thumb (Silverman 1986, eq. 3.31)
    robust = {
      scale <- min(s, iqr / 1.34)
      if (!is.finite(scale) || scale <= 0) scale <- s
      h <- 0.9 * scale * n^(-0.2)
      list(h = h, details = paste0(
        "Robust rule of thumb [h = 0.9 * min(sd, IQR/1.34) * n^(-1/5)]; ",
        "Silverman (1986). sd = ", fmt(s), ", IQR = ", fmt(iqr),
        ", n = ", n, ", h = ", fmt(h)))
    },

    ## 3. IQR rule
    iqr = {
      if (!is.finite(iqr) || iqr <= 0)
        stop("The IQR of the smoothing variable is zero; use h = ",
             "'silverman' or 'robust' instead.", call. = FALSE)
      h <- iqr * n^(-0.2)
      list(h = h, details = paste0(
        "IQR rule [h = IQR * n^(-1/5)]. IQR = ", fmt(iqr),
        ", n = ", n, ", h = ", fmt(h)))
    },

    ## 4. Leave-one-out check-loss cross-validation (tau = 0.5).
    ##    For each candidate h, y_i is predicted by the kernel-weighted
    ##    median of y_{-i} with weights K((x_j - x_i)/h) and the
    ##    criterion is the mean check loss rho_0.5(y_i - yhat_i).
    ##    The idea of check-loss cross-validation for kernel quantile
    ##    regression follows Yu and Jones (1998).
    cv_fast = {
      needs_y()
      h_grid <- seq(0.5 * h_silver, 2 * h_silver, length.out = 15L)
      ord    <- order(y)
      ys     <- y[ord]
      xs     <- x[ord]
      score  <- vapply(h_grid, function(hc) {
        loss <- vapply(seq_len(n), function(i) {
          w <- stats::dnorm((xs - xs[i]) / hc)
          w[i] <- 0
          sw <- sum(w)
          if (!is.finite(sw) || sw <= 0) return(NA_real_)
          pred <- ys[which(cumsum(w) / sw >= 0.5)[1L]]
          u <- ys[i] - pred
          u * (0.5 - (u < 0))
        }, numeric(1L))
        mean(loss, na.rm = TRUE)
      }, numeric(1L))
      if (all(!is.finite(score)))
        stop("cv_fast: cross-validation failed for all candidate ",
             "bandwidths.", call. = FALSE)
      h <- h_grid[which.min(score)]
      list(h = h, details = paste0(
        "Leave-one-out check-loss cross-validation at tau = 0.5 over 15 ",
        "candidates in [0.5, 2] x Silverman; Yu and Jones (1998). ",
        "Range = [", fmt(min(h_grid)), ", ", fmt(max(h_grid)), "], h = ",
        fmt(h)))
    },

    ## 5. Sheather-Jones (1991) solve-the-equation plug-in
    sj = {
      h <- tryCatch(stats::bw.SJ(x), error = function(e) fallback("SJ", e))
      list(h = h, details = paste0(
        "Sheather-Jones solve-the-equation plug-in [stats::bw.SJ]; ",
        "Sheather and Jones (1991). n = ", n, ", h = ", fmt(h)))
    },

    ## 6. Ruppert, Sheather and Wand (1995) direct plug-in
    yj = {
      needs_y()
      if (!requireNamespace("KernSmooth", quietly = TRUE))
        stop("h = 'yj' requires the 'KernSmooth' package.", call. = FALSE)
      if (n < 10L)
        stop("Too few observations for the 'yj' plug-in selector.",
             call. = FALSE)
      h <- tryCatch(KernSmooth::dpill(x, y),
                    error = function(e) fallback("YJ (dpill)", e))
      if (!is.finite(h) || h <= 0) {
        warning("dpill() returned an invalid bandwidth. Falling back to ",
                "Silverman's rule.", call. = FALSE)
        h <- h_silver
      }
      list(h = h, details = paste0(
        "Ruppert-Sheather-Wand direct plug-in for local linear regression ",
        "[KernSmooth::dpill(x, y)]; Ruppert, Sheather and Wand (1995). ",
        "n = ", n, ", h = ", fmt(h)))
    },

    ## 7. np least-squares cross-validation (local linear)
    npcv = {
      needs_y()
      if (!requireNamespace("np", quietly = TRUE))
        stop("h = 'npcv' requires the 'np' package.", call. = FALSE)
      h <- tryCatch({
        bw <- np::npregbw(xdat = x, ydat = y, regtype = "ll",
                          bwmethod = "cv.ls")
        bw$bw
      }, error = function(e) fallback("npcv (npregbw)", e))
      if (!is.finite(h) || h <= 0) {
        warning("npregbw() returned an invalid bandwidth. Falling back to ",
                "Silverman's rule.", call. = FALSE)
        h <- h_silver
      }
      list(h = h, details = paste0(
        "Least-squares cross-validation for local linear regression ",
        "[np::npregbw(regtype = 'll')]; Li and Racine (2004), Hayfield and ",
        "Racine (2008). n = ", n, ", h = ", fmt(h)))
    },

    ## 8. np smooth-coefficient cross-validation
    sc = {
      needs_y()
      if (!requireNamespace("np", quietly = TRUE))
        stop("h = 'sc' requires the 'np' package.", call. = FALSE)
      h <- tryCatch({
        bw <- np::npscoefbw(xdat = x, ydat = y)
        mean(bw$bw)
      }, error = function(e) fallback("sc (npscoefbw)", e))
      if (!is.finite(h) || h <= 0) {
        warning("npscoefbw() returned an invalid bandwidth. Falling back ",
                "to Silverman's rule.", call. = FALSE)
        h <- h_silver
      }
      list(h = h, details = paste0(
        "Smooth-coefficient cross-validation [np::npscoefbw]; Li, Huang, ",
        "Li and Fu (2002), Li and Racine (2010). n = ", n, ", h = ",
        fmt(h)))
    },

    ## 9. Slope-oriented rule of thumb: n^(-1/7) rate of local polynomial
    ##    first-derivative estimation (Fan and Gijbels 1996), robust scale.
    slope = {
      scale <- min(s, iqr / 1.34)
      if (!is.finite(scale) || scale <= 0) scale <- s
      C <- .qqr_slope_constant()
      h <- C * scale * n^(-1 / 7)
      list(h = h, details = paste0(
        "Slope-oriented rule of thumb [h = ", C,
        " * min(sd, IQR/1.34) * n^(-1/7)]; rate from Fan and Gijbels ",
        "(1996), constant calibrated by simulation. sd = ", fmt(s),
        ", IQR = ", fmt(iqr), ", n = ", n, ", h = ", fmt(h)))
    },

    ## 10. cv_fast inflated from the level rate n^(-1/5) to the slope rate
    ##     n^(-1/7).
    cv_slope = {
      needs_y()
      h0 <- .qqr_select_h(x, y, "cv_fast", kernel)$h
      h  <- h0 * n^(1 / 5 - 1 / 7)
      list(h = h, details = paste0(
        "Leave-one-out check-loss cross-validation at tau = 0.5 (h_cv = ",
        fmt(h0), ") multiplied by n^(1/5 - 1/7) to move from the level ",
        "rate to the slope rate; Yu and Jones (1998), Fan and Gijbels ",
        "(1996). n = ", n, ", h = ", fmt(h)))
    },

    stop("Unknown bandwidth method '", method, "'. Use one of: ",
         paste(.qqr_h_methods(), collapse = ", "),
         ", or a single positive number.", call. = FALSE)
  )

  if (!is.finite(out$h) || out$h <= 0)
    stop("Bandwidth selection produced an invalid value (", out$h, ").",
         call. = FALSE)

  list(h = out$h, method = method, details = out$details)
}

## Names of the automatic bandwidth methods.
.qqr_h_methods <- function() {
  c("slope", "cv_slope", "silverman", "robust", "iqr", "cv_fast", "sj",
    "yj", "npcv", "sc")
}

## Methods that use the outcome y.
.qqr_h_needs_y <- function() c("cv_slope", "cv_fast", "yj", "npcv", "sc")

## Constant of the slope-oriented rule of thumb (calibrated in
## tools/validation/bandwidth_calibration.R).
.qqr_slope_constant <- function() 2

## Human-readable label for print methods.
.qqr_h_label <- function(method) {
  labs <- c(
    slope     = "Slope-oriented rule of thumb",
    cv_slope  = "Check-loss CV, slope rate",
    silverman = "Silverman's rule of thumb",
    robust    = "Robust rule of thumb",
    iqr       = "IQR-based rule",
    cv_fast   = "Leave-one-out check-loss CV",
    cv        = "Leave-one-out check-loss CV",
    sj        = "Sheather-Jones plug-in",
    yj        = "Ruppert-Sheather-Wand plug-in",
    npcv      = "np least-squares CV",
    sc        = "np smooth-coefficient CV",
    manual    = "Manual (user-supplied)"
  )
  if (is.null(method)) return("unknown")
  if (method %in% names(labs)) labs[[method]] else method
}
