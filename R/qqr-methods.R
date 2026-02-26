#' Print method for QQR objects
#'
#' @description
#' Provides a compact console summary of a fitted \code{"qqr"} object,
#' reporting the key modelling choices and stored surfaces. The output
#' lists the call, dependent variable, quantile regressor, control
#' variables, intercept flag, sample size, grid size, quantile levels and
#' the set of coefficients for which surfaces are stored.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @details
#' This method is intended for a quick inspection of the fitted QQR model
#' structure. It does not compute any new statistics; it only reports
#' information already stored in the object.
#'
#' @return
#' The input object \code{x} is returned invisibly, so the method can be
#' used at the end of a pipe if desired.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{summary.qqr}}
#'
#' @method print qqr
#'
#' @export
#' @rawNamespace export(print.qqr)
print.qqr <- function(x, ...) {
  cat("Quantile-on-Quantile Regression\n")
  cat(strrep("-", 50L), "\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Dependent variable  : ", x$dependent, "\n")
  cat("Quantile regressor  : ", x$independent, "\n")
  if (length(x$controls))
    cat("Control variables   : ", paste(x$controls, collapse = ", "), "\n")
  cat("Intercept           : ", if (isTRUE(x$intercept)) "yes" else "no", "\n")
  cat("Observations        : ", x$n, "\n")
  cat("Grid points (m)     : ", length(x$grid), "\n")
  cat("Quantile levels     : ", paste(x$taus, collapse = " "), "\n")
  cat("Kernel              : ", x$kernel, "\n")
  cat("SE type             : ", x$se_type, "\n")
  ## Bandwidth
  h_lab <- switch(
    x$h_method,
    silverman = "Silverman's rule",
    robust    = "Robust rule",
    iqr       = "IQR-based rule",
    cv        = "Cross-validation (LOO)",
    manual    = "Manual (user-supplied)",
    x$h_method
  )
  cat(sprintf("Bandwidth h         :  %.6f  [%s]\n", x$bandwidth, h_lab))
  if (!is.null(x$h_method) && x$h_method == "manual")
    cat("  ** Manual h may introduce specification bias.",
        "Consider h = 'silverman' or 'robust'. **\n")
  cat(strrep("-", 50L), "\n")
  invisible(x)
}


#' Summary method for QQR objects
#'
#' @description
#' Computes and returns a compact set of global summaries for the
#' coefficient surface of the quantile regressor in a fitted QQR model.
#' These include the range, mean and standard deviation of the slope
#' surface, together with row and column means across grid points and
#' quantile levels.
#'
#' @param object An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @details
#' Let \eqn{\beta(z_i, \tau_j)} denote the coefficient on the quantile
#' regressor at grid location \eqn{z_i} and quantile level \eqn{\tau_j}.
#' The summary method extracts the full matrix \eqn{\beta} and reports:
#' \itemize{
#'   \item the observed range of the grid values used in the local
#'         regressions;
#'   \item the minimum and maximum of \eqn{\beta(z_i, \tau_j)} over all
#'         \eqn{i, j};
#'   \item the overall mean and standard deviation of the slope values;
#'   \item row means (averages over \eqn{\tau_j} for each grid location);
#'   \item column means (averages over \eqn{z_i} for each quantile level).
#' }
#' These summaries are useful for characterising the overall strength and
#' variability of the quantile-on-quantile relationship.
#'
#' @return
#' An object of class \code{"summary.qqr"} containing:
#' \itemize{
#'   \item \code{call} The original call.
#'   \item \code{dependent}, \code{independent}, \code{controls},
#'         \code{intercept} Names and flags describing the model.
#'   \item \code{taus} The vector of quantile levels.
#'   \item \code{grid_range} The range of grid values for the quantile
#'         regressor.
#'   \item \code{beta_min}, \code{beta_max} The minimum and maximum slope
#'         values across the grid and quantiles.
#'   \item \code{beta_mean}, \code{beta_sd} The overall mean and standard
#'         deviation of the slope surface.
#'   \item \code{row_means} Vector of row means of the slope matrix.
#'   \item \code{col_means} Vector of column means of the slope matrix.
#' }
#'
#' The associated \code{\link{print.summary.qqr}} method formats and
#' prints these components to the console.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{print.summary.qqr}}
#'
#' @export
#' @rawNamespace export(summary.qqr)
summary.qqr <- function(object, ...) {
  beta <- object$qcoef
  out <- list(
    call        = object$call,
    dependent   = object$dependent,
    independent = object$independent,
    controls    = object$controls,
    intercept   = object$intercept,
    taus        = object$taus,
    grid_range  = range(object$grid, na.rm = TRUE),
    beta_min    = min(beta,            na.rm = TRUE),
    beta_max    = max(beta,            na.rm = TRUE),
    beta_mean   = mean(beta,           na.rm = TRUE),
    beta_sd     = stats::sd(as.numeric(beta), na.rm = TRUE),
    row_means   = rowMeans(beta,       na.rm = TRUE),
    col_means   = colMeans(beta,       na.rm = TRUE),
    bandwidth   = object$bandwidth,
    h_method    = object$h_method,
    h_details   = object$h_details,
    kernel      = object$kernel,
    se_type     = object$se_type,
    n           = object$n
  )
  class(out) <- "summary.qqr"
  out
}



#' Print method for QQR summaries
#'
#' @description
#' Prints a concise textual summary of an object of class
#' \code{"summary.qqr"}, as created by \code{\link{summary.qqr}}.
#' The output reports the main variables, grid range, and basic statistics
#' for the slope surface in a readable format.
#'
#' @param x An object of class \code{"summary.qqr"}.
#' @param digits Integer number of significant digits to use when printing
#'   numeric values.
#' @param ... Further arguments passed to or from other methods (currently
#'   unused).
#'
#' @details
#' This is a pure print method that does not recompute any statistics. It
#' formats the information contained in the \code{"summary.qqr"} object
#' and writes it to the console.
#'
#' @return
#' The input object \code{x} is returned invisibly.
#'
#' @seealso
#' \code{\link{summary.qqr}}, \code{\link{qqr}}
#'
#' @method print summary.qqr
#' @export
#' @rawNamespace export(print.summary.qqr)
print.summary.qqr <- function(x, digits = 4, ...) {
  fmt <- function(v) format(round(v, digits), nsmall = digits)
  cat("Quantile-on-Quantile Regression -- Summary\n")
  cat(strrep("=", 55L), "\n")
  cat("Dependent variable  : ", x$dependent, "\n")
  cat("Quantile regressor  : ", x$independent, "\n")
  if (length(x$controls))
    cat("Control variables   : ", paste(x$controls, collapse = ", "), "\n")
  cat("Intercept           : ",
      if (isTRUE(x$intercept)) "yes" else "no", "\n")
  cat("Observations (n)    : ", x$n, "\n")
  cat("Quantile levels (K) : ", length(x$taus), "\n")
  cat("Kernel              : ", x$kernel, "\n")
  cat("SE type             : ", x$se_type, "\n")
  cat(strrep("-", 55L), "\n")
  ## Bandwidth block
  h_lab <- switch(
    x$h_method,
    silverman = "Silverman rule of thumb",
    robust    = "Robust rule of thumb",
    iqr       = "IQR-based rule",
    cv        = "Leave-one-out cross-validation",
    manual    = "Manual (user-supplied)",
    x$h_method
  )
  cat(sprintf("Bandwidth (h)       :  %s\n", fmt(x$bandwidth)))
  cat(sprintf("Selection method    :  %s\n", h_lab))
  if (!is.null(x$h_details))
    cat("Details             : ", x$h_details, "\n")
  if (!is.null(x$h_method) && x$h_method == "manual")
    cat("  ** Manual h -- results may be sensitive to this choice **\n")
  cat(strrep("-", 55L), "\n")
  ## Surface statistics
  cat("Grid range          : [",
      fmt(x$grid_range[1]), ", ", fmt(x$grid_range[2]), "]\n", sep = "")
  cat("Slope range         : [",
      fmt(x$beta_min), ", ", fmt(x$beta_max), "]\n", sep = "")
  cat("Slope mean          :  ", fmt(x$beta_mean), "\n")
  cat("Slope SD            :  ", fmt(x$beta_sd),   "\n")
  cat(strrep("=", 55L), "\n")
  invisible(x)
}

#' Static 3D plots of QQR coefficient or p-value surfaces
#'
#' @description
#' Draws static three-dimensional perspective plots of coefficient or
#' p-value surfaces from a fitted \code{"qqr"} object, using either the
#' \pkg{plot3D} package or base graphics. These plots show how the
#' conditional effect of a regressor varies over quantiles of both the
#' regressor and the response.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Character string identifying which coefficient surface to
#'   plot. This can be the quantile regressor (using its original name),
#'   the word \code{"intercept"} (or \code{"(Intercept)"}) for the
#'   intercept, or the name of a control variable.
#' @param what Character string indicating which surface to display:
#'   \code{"coef"} for coefficient values or \code{"pvalue"} for
#'   p-values.
#' @param engine Character string selecting the plotting backend. If
#'   \code{"plot3D"}, the function attempts to use
#'   \code{\link[plot3D]{persp3D}} from the \pkg{plot3D} package. If
#'   \code{"base"}, or if \pkg{plot3D} is not available, the base
#'   \code{\link[graphics]{persp}} function is used.
#' @param ... Additional graphical parameters passed to the underlying
#'   plotting function, for example \code{theta}, \code{phi}, or axis
#'   labels.
#'
#' @details
#' The rows of each surface correspond to a sequence of quantiles of the
#' quantile regressor, and the columns correspond to the quantile levels
#' \code{taus} of the dependent variable supplied when fitting the model.
#' The vertical axis shows either coefficient values or p-values at each
#' combination of regressor quantile and response quantile.
#'
#' These plots are useful for visualising how the strength and sign of a
#' relationship change across the joint distribution of the variables.
#'
#' @return
#' The input object \code{x} is returned invisibly. The function is called
#' for its side effect of producing a 3D plot.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{Int.plot}},
#' \code{\link[graphics]{persp}}, \code{\link[plot3D]{persp3D}}
#'
#' @method plot qqr
#'
#' @export
#' @rawNamespace export(plot.qqr)
plot.qqr <- function(x,
                     coef   = x$independent,
                     what   = c("coef", "pvalue"),
                     engine = c("plot3D", "base"),
                     ...) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")
  what   <- match.arg(what)
  engine <- match.arg(engine)

  surf_name <- .qqr_match_coef(x, coef)

  zmat <- if (what == "coef") {
    x$coef_surfaces[[surf_name]]
  } else {
    x$pval_surfaces[[surf_name]]
  }

  x_axis <- seq(0.05, 0.95, length.out = nrow(zmat))
  y_axis <- x$taus

  zlab <- if (what == "coef") "Coefficient" else "p-value"
  main <- paste0("QQR ", zlab, " surface: ",
                 x$dependent, " on ", coef)

  xlab <- paste0("Independent quantile (", coef, ")")
  ylab <- "Dependent quantile (tau)"

  if (engine == "plot3D" && requireNamespace("plot3D", quietly = TRUE)) {
    persp3D_fn <- get("persp3D", envir = asNamespace("plot3D"))
    persp3D_fn(
      x = x_axis,
      y = y_axis,
      z = zmat,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      main = main,
      ...
    )
  } else {
    graphics::persp(
      x = x_axis,
      y = y_axis,
      z = zmat,
      xlab = xlab,
      ylab = ylab,
      zlab = zlab,
      main = main,
      ...
    )
  }

  invisible(x)
}



#' Interactive 3D plots of QQR coefficient or p-value surfaces
#'
#' @description
#' \code{Int.plot} provides interactive three-dimensional visualisations of
#' coefficient or p-value surfaces from a fitted \code{"qqr"} object using
#' the \pkg{rgl} package. It is an S3 generic with a method for objects of
#' class \code{"qqr"}.
#'
#' @param x For the \code{qqr} method, an object of class \code{"qqr"}
#'   produced by \code{\link{qqr}}.
#' @param coef Character string identifying which coefficient surface to
#'   plot. This can be the name of the quantile regressor (as supplied in
#'   \code{qvar}), the word \code{"intercept"} (or \code{"(Intercept)"})
#'   for the intercept, or the name of one of the control variables.
#' @param what Character string specifying which surface to display.
#'   Either \code{"coef"} for coefficient values or \code{"pvalue"} for
#'   the associated p-values.
#' @param add_plane Logical flag indicating whether to add a horizontal
#'   reference plane to the plot, for example at zero for coefficients or
#'   at a chosen significance level for p-values.
#' @param plane_z Numeric value giving the height of the optional
#'   reference plane on the vertical axis.
#' @param ... Additional graphical parameters passed to
#'   \code{\link[rgl]{plot3d}}, such as \code{col}, \code{size}, or
#'   viewing options.
#'
#' @details
#' The method constructs a grid over the independent and dependent
#' quantile axes and plots the corresponding values of the selected
#' surface in an interactive \pkg{rgl} window. Users can rotate, zoom and
#' inspect the surface using standard \pkg{rgl} controls.
#'
#' The independent quantile axis corresponds to quantiles of the quantile
#' regressor, while the dependent quantile axis corresponds to the
#' \code{taus} used in \code{\link{qqr}}.
#'
#' @return
#' The function is called for its side effects, namely opening an
#' interactive 3D graphics device. It returns the input object
#' \code{x} invisibly.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{plot.qqr}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n  <- 200
#' x  <- rnorm(n)
#' y  <- x * qnorm(runif(n))
#' dat <- data.frame(y = y, x = x)
#'
#' fit <- qqr(y ~ x, data = dat, m = 15)
#'
#' # Interactive coefficient surface
#' Int.plot(fit, coef = "x", what = "coef", add_plane = TRUE, plane_z = 0)
#'
#' # Interactive p-value surface with 5% reference plane
#' Int.plot(fit, coef = "x", what = "pvalue", add_plane = TRUE, plane_z = 0.05)
#' }
#'
#' @export
Int.plot <- function(x, ...) {
  UseMethod("Int.plot")
}

#' @rdname Int.plot
#' @method Int.plot qqr
#'
#' @export
#' @rawNamespace export(Int.plot.qqr)
Int.plot.qqr <- function(x,
                         coef      = x$independent,
                         what      = c("coef", "pvalue"),
                         add_plane = FALSE,
                         plane_z   = 0,
                         ...) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop("Package 'rgl' is required for Int.plot.qqr().")
  }

  what <- match.arg(what)

  surf_name <- .qqr_match_coef(x, coef)

  zmat <- if (what == "coef") {
    x$coef_surfaces[[surf_name]]
  } else {
    x$pval_surfaces[[surf_name]]
  }

  x_axis <- seq(0.05, 0.95, length.out = nrow(zmat))
  y_axis <- x$taus

  zlab <- if (what == "coef") "Coefficient" else "p-value"
  main <- paste0("QQR ", zlab, " surface: ",
                 x$dependent, " on ", coef)

  xlab <- paste0("Independent quantile (", coef, ")")
  ylab <- "Dependent quantile (tau)"

  d     <- expand.grid(x = x_axis, y = y_axis)
  z_vec <- as.vector(zmat)

  rgl::open3d()
  rgl::plot3d(
    x = d$x,
    y = d$y,
    z = z_vec,
    xlab = xlab,
    ylab = ylab,
    zlab = zlab,
    main = main,
    ...
  )

  if (add_plane) {
    rgl::planes3d(a = 0, b = 0, c = 1, d = -plane_z,
                  col = "red", alpha = 0.5)
  }

  invisible(x)
}

## Internal helper: project 3D coordinates onto 2D using a persp() matrix.
## Equivalent to the graphics package projection helper.
## @keywords internal
.qqr_trans3d <- function(x, y, z, pmat) {
  tr <- cbind(x, y, z, 1) %*% pmat
  list(x = tr[, 1L] / tr[, 4L],
       y = tr[, 2L] / tr[, 4L])
}



#' Ribbon Panel 3D Surface Plot for QQR Coefficient Surfaces
#'
#' @description
#' \code{panel_surface.qqr} produces a publication-quality three-dimensional
#' ribbon panel plot of QQR coefficient surfaces, with one coloured flat panel
#' per coefficient. Each panel spans the quantile axis on one horizontal
#' dimension and displays the coefficient value on the vertical axis, with
#' panels arranged side-by-side along the second horizontal dimension. The
#' colour of each panel encodes the coefficient magnitude using a
#' continuous palette. This style closely matches the visualisation
#' convention introduced in applied QQR research for displaying multiple
#' coefficient surfaces simultaneously.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Character vector of coefficient names to include as panels.
#'   Defaults to all user-facing coefficients in the model. Each element
#'   must be a valid coefficient name: \code{"(Intercept)"} or
#'   \code{"intercept"}, the quantile regressor name, or a control
#'   variable name.
#' @param what Character string indicating whether to plot coefficient
#'   values (\code{"coef"}) or p-values (\code{"pvalue"}).
#' @param palette Character string passed to
#'   \code{\link[grDevices]{hcl.colors}} specifying the colour palette.
#'   Good choices include \code{"viridis"}, \code{"plasma"},
#'   \code{"YlGnBu"}, and \code{"RdBu"}. Default is \code{"viridis"}.
#' @param theta Numeric azimuthal viewing angle in degrees passed to
#'   \code{\link[graphics]{persp}}. Default is \code{-30}.
#' @param phi Numeric colatitude viewing angle in degrees passed to
#'   \code{\link[graphics]{persp}}. Default is \code{25}.
#' @param r Numeric distance of the viewpoint from the centre of the
#'   plot, passed to \code{\link[graphics]{persp}}. Default is \code{3}.
#' @param expand Numeric expansion factor applied to the z-axis passed to
#'   \code{\link[graphics]{persp}}. Default is \code{0.6}.
#' @param shade Numeric shading value between 0 and 1 passed to
#'   \code{\link[graphics]{persp}}. Higher values give stronger shading.
#'   Default is \code{0.3}.
#' @param border Colour of panel borders. Use \code{NA} for no border.
#'   Default is \code{"grey30"}.
#' @param add_zero_plane Logical flag indicating whether to add a
#'   translucent horizontal reference plane at \code{z = 0}. Useful for
#'   identifying sign changes in coefficient surfaces. Default is
#'   \code{FALSE}.
#' @param zero_plane_col Colour of the zero-reference plane when
#'   \code{add_zero_plane = TRUE}. Default is \code{"red"}.
#' @param main Optional character string for the overall plot title.
#'   If \code{NULL} a default title is constructed.
#' @param zlab Character string for the vertical axis label.
#'   Default is \code{"Coefficient"}.
#' @param cex_names Numeric character expansion factor for panel labels
#'   on the variable axis. Default is \code{0.75}.
#' @param digits Integer number of decimal places for the colour legend
#'   tick labels. Default is \code{2}.
#' @param legend Logical flag indicating whether to draw a colour scale
#'   legend. Default is \code{TRUE}.
#' @param ... Additional arguments passed to
#'   \code{\link[graphics]{persp}}.
#'
#' @details
#' For each requested coefficient, the function extracts the row mean of
#' the corresponding surface (averaging over all conditional quantile
#' levels \code{taus}) and plots it as a flat coloured ribbon spanning
#' the quantile grid. The ribbons are arranged side-by-side, with the
#' variable names on the x-axis, quantile levels on the y-axis, and
#' coefficient values on the z-axis.
#'
#' Colour is mapped to the coefficient value using the chosen palette and
#' a shared scale across all panels so that magnitudes are visually
#' comparable. When \code{add_zero_plane = TRUE}, a semi-transparent
#' plane at \eqn{z = 0} is superimposed to help identify where
#' effects change sign.
#'
#' This function uses only base R graphics (\pkg{grDevices} and
#' \pkg{graphics}) and has no additional package dependencies beyond
#' those already required by \pkg{qqr}.
#'
#' @return
#' The function is called for its side effect of producing a 3D plot.
#' It returns the input object \code{x} invisibly.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{plot.qqr}}, \code{\link{Int.plot}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n   <- 300
#' x   <- rnorm(n)
#' z1  <- rnorm(n)
#' z2  <- rnorm(n)
#' y   <- 0.5 + x * qnorm(runif(n)) + 0.3 * z1 - 0.2 * z2
#' dat <- data.frame(y = y, x = x, z1 = z1, z2 = z2)
#'
#' fit <- qqr(y ~ x + z1 + z2, data = dat, qvar = "x",
#'            taus = seq(0.05, 0.95, by = 0.05), m = 19)
#'
#' # Default plot - all coefficients
#' panel_surface(fit)
#'
#' # Selected coefficients with custom palette and zero plane
#' panel_surface(fit,
#'               coef           = c("x", "z1", "z2"),
#'               palette        = "RdBu",
#'               add_zero_plane = TRUE,
#'               theta          = -40,
#'               phi            = 30)
#' }
#'
#' @export
panel_surface <- function(x, ...) {
  UseMethod("panel_surface")
}

#' @rdname panel_surface
#' @method panel_surface qqr
#' @export
#' @rawNamespace export(panel_surface.qqr)
panel_surface.qqr <- function(x,
                               coef           = NULL,
                               what           = c("coef", "pvalue"),
                               palette        = "viridis",
                               theta          = -30,
                               phi            = 25,
                               r              = 3,
                               expand         = 0.6,
                               shade          = 0.3,
                               border         = "grey30",
                               add_zero_plane = FALSE,
                               zero_plane_col = "red",
                               main           = NULL,
                               zlab           = NULL,
                               cex_names      = 0.75,
                               digits         = 2,
                               legend         = TRUE,
                               ...) {

  if (!inherits(x, "qqr"))
    stop("Object must be of class 'qqr'.")

  what <- match.arg(what)

  ## -----------------------------------------------------------------------
  ## Resolve coefficients to plot
  ## -----------------------------------------------------------------------
  user_coefs <- c(
    if (isTRUE(x$intercept)) "(Intercept)" else character(0L),
    x$independent,
    x$controls
  )

  if (is.null(coef)) {
    coef <- user_coefs
  } else {
    coef <- vapply(coef, function(cn) {
      if (tolower(cn) %in% c("intercept", "(intercept)")) {
        if ("(Intercept)" %in% user_coefs) return("(Intercept)")
      }
      if (!(cn %in% user_coefs))
        stop("Coefficient '", cn, "' not found. Available: ",
             paste(user_coefs, collapse = ", "))
      cn
    }, character(1L))
  }

  n_coef <- length(coef)
  taus   <- x$taus
  n_tau  <- length(taus)

  ## -----------------------------------------------------------------------
  ## Extract surfaces and compute row means per coefficient
  ## -----------------------------------------------------------------------
  surfaces <- x$coef_surfaces
  if (what == "pvalue") surfaces <- x$pval_surfaces

  # Row means: one value per (coef, tau) combination
  # Each coefficient gets a flat ribbon: constant across grid dimension
  ribbon_mat <- matrix(NA_real_, nrow = n_coef, ncol = n_tau)
  rownames(ribbon_mat) <- coef
  colnames(ribbon_mat) <- round(taus, 3)

  for (k in seq_len(n_coef)) {
    surf <- surfaces[[coef[k]]]
    ribbon_mat[k, ] <- rowMeans(surf, na.rm = TRUE)
  }

  ## -----------------------------------------------------------------------
  ## Build the z-matrix for persp()
  ## x-axis: coefficients (panels), with padding between each
  ## y-axis: quantile levels (taus)
  ## -----------------------------------------------------------------------
  panel_width <- 1L          # width of each ribbon in x units
  gap_width   <- 0L          # gap between ribbons (0 = touching)
  x_step      <- panel_width + gap_width

  n_x_pts <- n_coef * x_step + 1L
  x_axis  <- seq(0, n_coef * x_step, length.out = n_x_pts)
  y_axis  <- taus

  z_full  <- matrix(NA_real_, nrow = n_x_pts, ncol = n_tau)

  for (k in seq_len(n_coef)) {
    x_start <- (k - 1L) * x_step
    x_end   <- k * x_step
    idx     <- which(x_axis >= x_start & x_axis <= x_end)
    for (pt in idx) {
      z_full[pt, ] <- ribbon_mat[k, ]
    }
  }

  ## -----------------------------------------------------------------------
  ## Colour mapping
  ## -----------------------------------------------------------------------
  z_range   <- range(z_full, na.rm = TRUE)
  z_min     <- z_range[1L]
  z_max     <- z_range[2L]
  n_colours <- 100L
  pal       <- grDevices::hcl.colors(n_colours, palette = palette)

  colour_from_z <- function(z_val) {
    if (is.na(z_val)) return(grDevices::adjustcolor("grey80", alpha.f = 0.5))
    idx <- as.integer((z_val - z_min) / (z_max - z_min) * (n_colours - 1L)) + 1L
    idx <- max(1L, min(n_colours, idx))
    pal[idx]
  }

  # Face colours for each panel
  facecol <- matrix(NA_character_,
                    nrow = nrow(z_full) - 1L,
                    ncol = ncol(z_full) - 1L)

  for (i in seq_len(nrow(facecol))) {
    for (j in seq_len(ncol(facecol))) {
      z_face <- mean(c(z_full[i, j], z_full[i+1L, j],
                       z_full[i, j+1L], z_full[i+1L, j+1L]),
                     na.rm = TRUE)
      facecol[i, j] <- colour_from_z(z_face)
    }
  }

  ## -----------------------------------------------------------------------
  ## Labels
  ## -----------------------------------------------------------------------
  if (is.null(zlab)) {
    zlab <- if (what == "coef") "Coefficient" else "P-value"
  }

  if (is.null(main)) {
    main <- paste0("QQR ", if (what == "coef") "Coefficient" else "P-value",
                   " Surfaces -- ", x$dependent)
  }

  # Tick positions and labels for the variable axis
  x_tick_pos    <- seq(x_step / 2, by = x_step, length.out = n_coef)
  x_tick_labels <- coef

  ## -----------------------------------------------------------------------
  ## Base plot
  ## -----------------------------------------------------------------------
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (legend) {
    graphics::layout(matrix(c(1L, 2L), nrow = 1L), widths = c(5, 1))
  }

  graphics::par(mar = c(2, 2, 3, 1))

  pmat <- graphics::persp(
    x       = x_axis,
    y       = y_axis,
    z       = z_full,
    col     = facecol,
    theta   = theta,
    phi     = phi,
    r       = r,
    expand  = expand,
    shade   = shade,
    border  = border,
    xlab    = "",
    ylab    = "Quantiles",
    zlab    = zlab,
    main    = main,
    xaxt    = "n",
    ticktype = "detailed",
    nticks  = 5L,
    ...
  )

  ## Add variable name labels on the x-axis
  for (k in seq_len(n_coef)) {
    x_mid <- (k - 0.5) * x_step
    pt    <- .qqr_trans3d(x_mid, min(y_axis) - 0.05,
                               z_min - 0.05 * (z_max - z_min),
                               pmat)
    graphics::text(pt$x, pt$y,
                   labels = x_tick_labels[k],
                   cex    = cex_names,
                   adj    = 0.5)
  }

  ## Optional zero plane
  if (add_zero_plane && z_min <= 0 && z_max >= 0) {
    corners <- rbind(
      c(min(x_axis), min(y_axis), 0),
      c(max(x_axis), min(y_axis), 0),
      c(max(x_axis), max(y_axis), 0),
      c(min(x_axis), max(y_axis), 0)
    )
    proj <- .qqr_trans3d(corners[, 1L], corners[, 2L],
                              corners[, 3L], pmat)
    graphics::polygon(proj$x, proj$y,
                      col    = grDevices::adjustcolor(zero_plane_col, alpha.f = 0.25),
                      border = grDevices::adjustcolor(zero_plane_col, alpha.f = 0.6))
  }

  ## -----------------------------------------------------------------------
  ## Colour legend
  ## -----------------------------------------------------------------------
  if (legend) {
    graphics::par(mar = c(4, 1, 4, 3))
    legend_vals <- seq(z_min, z_max, length.out = n_colours)
    legend_mat  <- matrix(legend_vals, nrow = 1L)

    graphics::image(
      x    = 1L,
      y    = legend_vals,
      z    = legend_mat,
      col  = pal,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )

    n_ticks  <- 5L
    tick_at  <- seq(z_min, z_max, length.out = n_ticks)
    tick_lab <- round(tick_at, digits)
    graphics::axis(4L, at = tick_at, labels = tick_lab,
                   las = 2L, cex.axis = 0.8, tick = TRUE)
    graphics::box()
  }

  invisible(x)
}


#' Publication-Quality 3D QQR Coefficient Surface Plot
#'
#' @description
#' \code{surface3d} produces a publication-quality three-dimensional surface
#' plot of a QQR coefficient surface for a selected regressor. The x-axis
#' represents quantiles of the quantile regressor (grid points), the y-axis
#' represents the conditional quantile levels (taus), and the z-axis shows
#' the estimated QQR coefficient at each grid-tau combination. Colour is
#' mapped to the coefficient magnitude. This directly visualises the full
#' QQR surface without any averaging or smoothing.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Character string giving the name of the coefficient surface
#'   to plot. This must be one of: the quantile regressor name, the word
#'   \code{"intercept"} (or \code{"(Intercept)"}), or a control variable
#'   name. Defaults to the quantile regressor (\code{x$independent}).
#' @param what Character string indicating whether to plot coefficient
#'   values (\code{"coef"}) or p-values (\code{"pvalue"}).
#'   Default is \code{"coef"}.
#' @param palette Character string passed to
#'   \code{\link[grDevices]{hcl.colors}} specifying the colour palette.
#'   Good choices include \code{"viridis"} (default), \code{"plasma"},
#'   \code{"RdBu"}, \code{"YlGnBu"}, \code{"Temps"}, and \code{"turbo"}.
#' @param theta Numeric azimuthal viewing angle in degrees. Default
#'   \code{-30}.
#' @param phi Numeric colatitude viewing angle in degrees. Default
#'   \code{25}.
#' @param r Numeric distance of the viewpoint from the centre.
#'   Default \code{3}.
#' @param expand Numeric expansion factor for the z-axis. Default
#'   \code{0.6}.
#' @param shade Numeric shading between 0 and 1. Default \code{0.25}.
#' @param border Colour of grid lines on the surface. Use \code{NA}
#'   to suppress grid lines. Default \code{"grey40"}.
#' @param alpha Numeric transparency between 0 (transparent) and 1
#'   (opaque). Default \code{0.85}.
#' @param add_zero_plane Logical. If \code{TRUE}, adds a translucent
#'   reference plane at \eqn{z = 0} to highlight sign changes.
#'   Default \code{FALSE}.
#' @param zero_plane_col Colour of the zero reference plane when
#'   \code{add_zero_plane = TRUE}. Default \code{"red"}.
#' @param clim Numeric vector of length 2 giving fixed colour scale
#'   limits \code{c(min, max)}. If \code{NULL} (default), limits are
#'   set from the data range.
#' @param main Optional character string for the plot title. If
#'   \code{NULL} a default title is constructed from the variable names.
#' @param xlab Character string for the x-axis label. Default uses the
#'   coefficient name.
#' @param ylab Character string for the y-axis label. Default is
#'   \code{"Quantile (tau)"}.
#' @param zlab Character string for the z-axis label. Default is
#'   \code{"Coefficient"} or \code{"P-value"} depending on \code{what}.
#' @param legend Logical. If \code{TRUE} (default), draws a colour
#'   scale legend on the right side of the plot.
#' @param digits Integer number of decimal places for legend tick
#'   labels. Default \code{2}.
#' @param nticks Integer number of axis tick marks. Default \code{5}.
#' @param ... Additional arguments passed to
#'   \code{\link[graphics]{persp}}.
#'
#' @details
#' The full \eqn{m \times K} coefficient matrix is plotted directly,
#' where \eqn{m} is the number of grid points (quantiles of the quantile
#' regressor) and \eqn{K} is the number of conditional quantile levels.
#' No averaging or smoothing is applied. This is equivalent to the
#' MATLAB \code{surf()} function applied to the raw QQR output matrix.
#'
#' The colour of each surface facet is determined by the average
#' coefficient value of its four corner vertices, mapped linearly onto
#' the chosen palette. When \code{clim} is supplied, the colour scale
#' is fixed across multiple plots for direct comparability -- the same
#' convention used in the MATLAB workflow.
#'
#' @return
#' The function is called for its side effect of producing a plot.
#' It returns the input object \code{x} invisibly.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{plot.qqr}}, \code{\link{Int.plot}},
#' \code{\link{panel_surface}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n   <- 300
#' x   <- rnorm(n)
#' z1  <- rnorm(n)
#' y   <- 0.5 + x * qnorm(runif(n)) + 0.3 * z1
#' dat <- data.frame(y = y, x = x, z1 = z1)
#'
#' fit <- qqr(y ~ x + z1, data = dat, qvar = "x",
#'            taus = seq(0.05, 0.95, by = 0.05), m = 19)
#'
#' # Plot coefficient surface for the quantile regressor
#' surface3d(fit, coef = "x")
#'
#' # Plot with fixed colour limits and zero plane
#' surface3d(fit, coef = "x",
#'           clim           = c(-2, 2),
#'           add_zero_plane = TRUE,
#'           palette        = "RdBu",
#'           theta          = 45,
#'           phi            = 30)
#'
#' # Plot p-value surface
#' surface3d(fit, coef = "x", what = "pvalue",
#'           clim    = c(0, 1),
#'           palette = "YlGnBu")
#'
#' # Plot control variable surface
#' surface3d(fit, coef = "z1", palette = "plasma")
#' }
#'
#' @export
surface3d <- function(x, ...) {
  UseMethod("surface3d")
}

#' @rdname surface3d
#' @method surface3d qqr
#' @export
#' @rawNamespace export(surface3d.qqr)
surface3d.qqr <- function(x,
                           coef           = x$independent,
                           what           = c("coef", "pvalue"),
                           palette        = "viridis",
                           theta          = -30,
                           phi            = 25,
                           r              = 3,
                           expand         = 0.6,
                           shade          = 0.25,
                           border         = "grey40",
                           alpha          = 0.85,
                           add_zero_plane = FALSE,
                           zero_plane_col = "red",
                           clim           = NULL,
                           main           = NULL,
                           xlab           = NULL,
                           ylab           = "Quantile (tau)",
                           zlab           = NULL,
                           legend         = TRUE,
                           digits         = 2,
                           nticks         = 5L,
                           ...) {

  ## -----------------------------------------------------------------------
  ## Input validation
  ## -----------------------------------------------------------------------
  if (!inherits(x, "qqr"))
    stop("Object must be of class 'qqr'.")

  what <- match.arg(what)

  ## Resolve coefficient name
  user_coefs <- c(
    if (isTRUE(x$intercept)) "(Intercept)" else character(0L),
    x$independent,
    x$controls
  )

  if (tolower(coef) %in% c("intercept", "(intercept)")) {
    coef <- "(Intercept)"
  }

  if (!(coef %in% user_coefs)) {
    stop(
      "Coefficient '", coef, "' not found. ",
      "Available coefficients: ",
      paste(user_coefs, collapse = ", ")
    )
  }

  ## -----------------------------------------------------------------------
  ## Extract the raw coefficient surface  (m x K matrix)
  ## -----------------------------------------------------------------------
  surf_mat <- if (what == "coef") {
    x$coef_surfaces[[coef]]
  } else {
    x$pval_surfaces[[coef]]
  }

  if (is.null(surf_mat))
    stop("Surface for coefficient '", coef, "' not found in the qqr object.")

  ## x-axis: quantile grid of the quantile regressor (m points)
  x_axis <- seq(0.05, 0.95, length.out = nrow(surf_mat))

  ## y-axis: conditional quantile levels (taus)
  y_axis <- x$taus

  ## -----------------------------------------------------------------------
  ## Colour mapping
  ## -----------------------------------------------------------------------
  z_min <- if (!is.null(clim)) clim[1L] else min(surf_mat, na.rm = TRUE)
  z_max <- if (!is.null(clim)) clim[2L] else max(surf_mat, na.rm = TRUE)

  if (z_min == z_max) {
    z_min <- z_min - 1e-6
    z_max <- z_max + 1e-6
  }

  n_col <- 100L
  pal   <- grDevices::hcl.colors(n_col, palette = palette)

  ## Function to map a z value to a colour
  map_col <- function(z_val) {
    if (is.na(z_val))
      return(grDevices::adjustcolor("grey90", alpha.f = alpha))
    idx <- as.integer(
      (z_val - z_min) / (z_max - z_min) * (n_col - 1L)
    ) + 1L
    idx <- max(1L, min(n_col, idx))
    grDevices::adjustcolor(pal[idx], alpha.f = alpha)
  }

  ## Face colours: average the four corner z-values per facet
  n_row <- nrow(surf_mat)
  n_col_mat <- ncol(surf_mat)

  facecol <- matrix(NA_character_,
                    nrow = n_row - 1L,
                    ncol = n_col_mat - 1L)

  for (i in seq_len(n_row - 1L)) {
    for (j in seq_len(n_col_mat - 1L)) {
      z_face <- mean(c(
        surf_mat[i,     j    ],
        surf_mat[i + 1L, j    ],
        surf_mat[i,     j + 1L],
        surf_mat[i + 1L, j + 1L]
      ), na.rm = TRUE)
      facecol[i, j] <- map_col(z_face)
    }
  }

  ## -----------------------------------------------------------------------
  ## Axis labels and title
  ## -----------------------------------------------------------------------
  if (is.null(zlab))
    zlab <- if (what == "coef") "Coefficient" else "P-value"

  if (is.null(xlab))
    xlab <- paste0("Quantile of ", coef)

  if (is.null(main))
    main <- paste0(
      "QQR ", if (what == "coef") "Coefficient" else "P-value",
      " Surface\n",
      x$dependent, " ~ ", coef
    )

  ## -----------------------------------------------------------------------
  ## Layout: surface + optional legend
  ## -----------------------------------------------------------------------
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (legend) {
    graphics::layout(
      matrix(c(1L, 2L), nrow = 1L),
      widths = c(5.5, 1)
    )
    graphics::par(mar = c(3, 3, 4, 1))
  } else {
    graphics::par(mar = c(3, 3, 4, 2))
  }

  ## -----------------------------------------------------------------------
  ## Main 3D surface
  ## -----------------------------------------------------------------------
  pmat <- graphics::persp(
    x        = x_axis,
    y        = y_axis,
    z        = surf_mat,
    col      = facecol,
    theta    = theta,
    phi      = phi,
    r        = r,
    expand   = expand,
    shade    = shade,
    border   = border,
    xlab     = xlab,
    ylab     = ylab,
    zlab     = zlab,
    main     = main,
    ticktype = "detailed",
    nticks   = nticks,
    ...
  )

  ## -----------------------------------------------------------------------
  ## Zero plane
  ## -----------------------------------------------------------------------
  if (add_zero_plane && z_min <= 0 && z_max >= 0) {
    corners_x <- c(min(x_axis), max(x_axis), max(x_axis), min(x_axis))
    corners_y <- c(min(y_axis), min(y_axis), max(y_axis), max(y_axis))
    corners_z <- rep(0, 4L)
    proj <- .qqr_trans3d(corners_x, corners_y, corners_z, pmat)
    graphics::polygon(
      proj$x, proj$y,
      col    = grDevices::adjustcolor(zero_plane_col, alpha.f = 0.25),
      border = grDevices::adjustcolor(zero_plane_col, alpha.f = 0.70),
      lwd    = 1.5
    )
  }

  ## -----------------------------------------------------------------------
  ## Colour legend
  ## -----------------------------------------------------------------------
  if (legend) {
    graphics::par(mar = c(5, 0.5, 5, 3.5))

    legend_z   <- seq(z_min, z_max, length.out = n_col)
    legend_mat <- matrix(legend_z, nrow = 1L)

    graphics::image(
      x    = 1L,
      y    = legend_z,
      z    = legend_mat,
      col  = pal,
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )

    tick_at  <- seq(z_min, z_max, length.out = nticks)
    tick_lab <- formatC(round(tick_at, digits),
                        format = "f", digits = digits)

    graphics::axis(
      side     = 4L,
      at       = tick_at,
      labels   = tick_lab,
      las      = 2L,
      cex.axis = 0.80,
      tick     = TRUE
    )
    graphics::box()
  }

  invisible(x)
}
