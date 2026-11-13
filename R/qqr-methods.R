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
  grid_type <- if (is.null(x$grid_type)) "range" else x$grid_type
  kspace    <- if (is.null(x$kernel_space)) "x" else x$kernel_space
  cat("Grid points (m)     : ", length(x$grid),
      if (grid_type == "quantile") "(empirical quantiles of the regressor)"
      else "(equally spaced over the regressor range)", "\n")
  cat("Quantile levels     : ", paste(x$taus, collapse = " "), "\n")
  cat("Kernel              : ", x$kernel,
      if (kspace == "cdf") "on F_n(x) - theta" else "on x - grid value", "\n")
  cat("Local polynomial    : ",
      if (identical(x$degree, 2L)) "quadratic (degree 2)" else
        "linear (degree 1)", "\n")
  cat("SE type             : ", x$se_type, "\n")
  if (!is.null(x$p_adjust))
    cat("P-value adjustment  : ", x$p_adjust, "\n")
  cat(sprintf("Bandwidth h         :  %.6g  [%s]\n", x$bandwidth,
              .qqr_h_label(x$h_method)))
  if (!is.null(x$h_method) && x$h_method == "manual")
    cat("  ** Manual bandwidth: report it and check sensitivity. **\n")
  if (!is.null(x$n_failed) && x$n_failed > 0L)
    cat("Failed local fits   : ", x$n_failed, "(cells are NA)\n")
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
#'   \item column means (averages over the grid for each quantile level
#'         \eqn{\tau_j}; with \code{grid = "quantile"} these are the
#'         \eqn{\theta}-averaged QQR slopes of Sim and Zhou, 2015, that
#'         \code{\link{compare.qqr}} sets against standard QR).
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
    thetas      = .qqr_thetas(object),
    grid_type   = if (is.null(object$grid_type)) "range" else object$grid_type,
    kernel_space = if (is.null(object$kernel_space)) "x" else
      object$kernel_space,
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
    degree      = if (is.null(object$degree)) 1L else object$degree,
    se_type     = object$se_type,
    p_adjust    = object$p_adjust,
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
  cat("Local polynomial    :  degree", if (is.null(x$degree)) 1L else
    x$degree, "\n")
  cat("SE type             : ", x$se_type, "\n")
  if (!is.null(x$p_adjust))
    cat("P-value adjustment  : ", x$p_adjust, "\n")
  cat(strrep("-", 55L), "\n")
  if (!is.null(x$grid_type))
    cat("Grid / kernel space : ", x$grid_type, "/", x$kernel_space, "\n")
  ## Bandwidth block
  cat(sprintf("Bandwidth (h)       :  %s\n", format(signif(x$bandwidth, 6L))))
  cat(sprintf("Selection method    :  %s\n", .qqr_h_label(x$h_method)))
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
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' plot(fit, coef = "x", engine = "base", theta = 40)
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

  zmat <- .qqr_get_surface(x, coef, what)

  x_axis <- .qqr_thetas(x)
  y_axis <- x$taus
  if (length(unique(x_axis)) < length(x_axis) || length(y_axis) < 2L)
    stop("A 3D surface needs at least two distinct grid quantiles and two ",
         "taus.", call. = FALSE)

  zlab <- if (what == "coef") "Coefficient" else "p-value"
  main <- paste0("QQR ", zlab, " surface: ",
                 x$dependent, " on ", coef)

  xlab <- paste0("Quantile of ", x$independent, " (theta)")
  ylab <- paste0("Quantile of ", x$dependent, " (tau)")

  if (engine == "plot3D" && requireNamespace("plot3D", quietly = TRUE)) {
    plot3D::persp3D(
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
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
#'   Int.plot(fit, coef = "x", add_plane = TRUE, plane_z = 0)
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

  zmat <- .qqr_get_surface(x, coef, what)

  x_axis <- .qqr_thetas(x)
  y_axis <- x$taus

  zlab <- if (what == "coef") "Coefficient" else "p-value"
  main <- paste0("QQR ", zlab, " surface: ",
                 x$dependent, " on ", coef)

  xlab <- paste0("Quantile of ", x$independent, " (theta)")
  ylab <- paste0("Quantile of ", x$dependent, " (tau)")

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
#' Plots one 3D coefficient (or p-value) surface per coefficient side by side
#' in a single figure, sharing a common colour scale and z-axis so that all
#' panels are directly visually comparable.
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
#' @param clim Numeric vector of length 2, \code{c(min, max)}, giving the
#'    shared colour and z-axis limits applied to every panel.  If
#'    \code{NULL} (default) the limits are computed automatically
#'    as the range of all finite values across all displayed
#'    surfaces, ensuring a common scale.  Supply a fixed value
#'    (e.g. \code{clim = c(-0.5, 0.5)}) to make plots from
#'    different models visually comparable.
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
#' For each requested coefficient the full \eqn{m \times K} surface is
#' drawn with \code{\link[graphics]{persp}} in its own panel (x-axis: grid
#' quantiles \eqn{\theta} of the quantile regressor; y-axis: \eqn{\tau};
#' z-axis: coefficient or p-value). The panels are placed side by side and
#' share one colour scale and one z-range.
#'
#' Colour is mapped to the coefficient value using the chosen palette and
#' a shared scale across all panels so that magnitudes are visually
#' comparable. When \code{add_zero_plane = TRUE}, a semi-transparent
#' plane at \eqn{z = 0} is superimposed to help identify where
#' effects change sign.
#'
#' This function uses only base R graphics (\pkg{grDevices} and
#' \pkg{graphics}) and has no additional package dependencies beyond
#' those already required by \pkg{qqrtools}.
#'
#' @return
#' The function is called for its side effect of producing a 3D plot.
#' It returns the input object \code{x} invisibly.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{plot.qqr}}, \code{\link{Int.plot}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' panel_surface(fit)
#' panel_surface(fit, coef = c("x", "w"), palette = "RdBu",
#'               add_zero_plane = TRUE)
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
                              clim           = NULL,
                              main           = NULL,
                              zlab           = NULL,
                              cex_names      = 0.95,
                              digits         = 2,
                              legend         = TRUE,
                              ...) {

  # ---------------------------------------------------------------------------
  # 0.  Input validation
  # ---------------------------------------------------------------------------
  if (!inherits(x, "qqr"))
    stop("Object must be of class 'qqr'.")

  what <- match.arg(what)

  # ---------------------------------------------------------------------------
  # 1.  Resolve which coefficients to plot
  # ---------------------------------------------------------------------------
  if (is.null(coef)) {
    coef <- names(x$coef_surfaces)
  } else {
    coef <- vapply(coef, function(cn) .qqr_match_coef(x, cn), character(1L),
                   USE.NAMES = FALSE)
  }

  n_coef <- length(coef)
  taus   <- x$taus       # outcome quantile levels (y-axis in each panel)

  # ---------------------------------------------------------------------------
  # 2.  Extract the correct surface list
  # ---------------------------------------------------------------------------
  # Each element is a matrix of dimension [m x J]:
  #   rows    = grid points (quantiles of the focal regressor)
  #   columns = tau levels  (quantiles of the outcome)
  if (what == "coef") {
    surfaces <- x$coef_surfaces
  } else {
    surfaces <- x$pval_surfaces
    if (is.null(surfaces))
      stop("No p-value surfaces found. Re-fit with se_type != 'none'.")
  }

  # Verify every requested coefficient has a surface
  missing_surfs <- setdiff(coef, names(surfaces))
  if (length(missing_surfs) > 0L)
    stop("Surface(s) not found for: ", paste(missing_surfs, collapse = ", "))

  # ---------------------------------------------------------------------------
  # 3.  Compute SHARED colour limits across ALL requested surfaces
  # ---------------------------------------------------------------------------
  # This is the key feature of panel_surface: every panel uses the same
  # colour mapping, so a dark red in panel 1 means the same magnitude as
  # a dark red in panel 3.  Without shared limits the panels are misleading.

  if (is.null(clim)) {
    all_vals <- unlist(lapply(coef, function(cn) as.vector(surfaces[[cn]])))
    all_vals <- all_vals[is.finite(all_vals)]
    if (length(all_vals) == 0L)
      stop("All surface values are non-finite. Cannot determine colour range.")
    z_min <- min(all_vals)
    z_max <- max(all_vals)
  } else {
    if (!is.numeric(clim) || length(clim) != 2L || clim[1L] >= clim[2L])
      stop("'clim' must be a numeric vector c(min, max) with min < max.")
    z_min <- clim[1L]
    z_max <- clim[2L]
  }

  # Avoid degenerate range (all identical values)
  if (isTRUE(all.equal(z_min, z_max))) {
    z_min <- z_min - 0.5
    z_max <- z_max + 0.5
  }

  z_range <- z_max - z_min

  # ---------------------------------------------------------------------------
  # 4.  Build the shared colour palette
  # ---------------------------------------------------------------------------
  n_colours <- 200L   # more colours = smoother gradient
  pal       <- grDevices::hcl.colors(n_colours, palette = palette, rev = FALSE)

  # Helper: map a numeric z value to a palette index
  z_to_col_idx <- function(z_val) {
    if (!is.finite(z_val)) return(NA_integer_)
    idx <- as.integer((z_val - z_min) / z_range * (n_colours - 1L)) + 1L
    max(1L, min(n_colours, idx))
  }

  # ---------------------------------------------------------------------------
  # 5.  Resolve axis labels and outer title
  # ---------------------------------------------------------------------------
  if (is.null(zlab))
    zlab <- if (what == "coef") "Coefficient" else "P-value"

  if (is.null(main))
    main <- paste0("QQR ",
                   if (what == "coef") "Coefficient" else "P-value",
                   " Surfaces -- ", x$dependent)

  # ---------------------------------------------------------------------------
  # 6.  Set up graphics layout
  # ---------------------------------------------------------------------------
  # Layout: [panel_1 | panel_2 | ... | panel_n | legend]
  # The legend column is narrow (width = 1) relative to each panel (width = 4).

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (legend) {
    layout_mat <- matrix(c(seq_len(n_coef), rep(n_coef + 1L, 1L)), nrow = 1L)
    col_widths  <- c(rep(4L, n_coef), 1L)
    graphics::layout(layout_mat, widths = col_widths)
  } else {
    graphics::layout(matrix(seq_len(n_coef), nrow = 1L))
  }

  # Outer margin: top = 3 lines (shared title), right = 0 (legend handles it)
  graphics::par(oma = c(1L, 0L, 3L, if (legend) 0L else 1L))

  # ---------------------------------------------------------------------------
  # 7.  Draw one persp() surface per coefficient
  # ---------------------------------------------------------------------------
  for (k in seq_len(n_coef)) {

    cn    <- coef[k]
    z_mat <- surfaces[[cn]]    # [m x J] matrix

    # x-axis: grid quantiles theta; y-axis: tau levels (shared)
    m_pts  <- nrow(z_mat)
    x_axis <- .qqr_thetas(x)
    y_axis <- taus
    if (m_pts < 2L || length(taus) < 2L)
      stop("panel_surface() needs at least two grid points and two taus.",
           call. = FALSE)

    # -- Compute face colours using the SHARED palette -----------------------
    # persp() colours the (m-1) x (J-1) faces.
    # We use the average z of the four corner vertices of each face.
    n_fx <- m_pts - 1L
    n_fy <- length(taus) - 1L
    facecol <- matrix(NA_character_, nrow = n_fx, ncol = n_fy)

    for (i in seq_len(n_fx)) {
      for (j in seq_len(n_fy)) {
        z_face <- mean(c(z_mat[i,     j],
                         z_mat[i + 1L, j],
                         z_mat[i,     j + 1L],
                         z_mat[i + 1L, j + 1L]),
                       na.rm = TRUE)
        col_idx       <- z_to_col_idx(z_face)
        facecol[i, j] <- if (is.na(col_idx))
          grDevices::adjustcolor("grey80", alpha.f = 0.4)
        else
          pal[col_idx]
      }
    }

    # -- Open panel and draw surface -----------------------------------------
    graphics::par(mar = c(2L, 1L, 3L, 1L))

    pmat <- graphics::persp(
      x        = x_axis,
      y        = y_axis,
      z        = z_mat,
      col      = facecol,
      theta    = theta,
      phi      = phi,
      r        = r,
      expand   = expand,
      shade    = shade,
      border   = border,
      xlab     = "theta",
      ylab     = "tau",
      zlab     = zlab,
      main     = cn,          # panel title = coefficient name
      cex.main = cex_names,   # user-controlled size for the panel title
      zlim     = c(z_min, z_max),   # SHARED z limits: critical for comparison
      ticktype = "detailed",
      nticks   = 4L,
      ...
    )

    # -- Optional zero-reference plane ---------------------------------------
    # Drawn only when zero lies within [z_min, z_max] for this panel
    if (add_zero_plane && z_min <= 0 && z_max >= 0) {
      corners_x <- c(min(x_axis), max(x_axis),
                     max(x_axis), min(x_axis))
      corners_y <- c(min(y_axis), min(y_axis),
                     max(y_axis), max(y_axis))
      corners_z <- rep(0, 4L)

      proj <- .qqr_trans3d(corners_x, corners_y, corners_z, pmat)

      graphics::polygon(
        proj$x,
        proj$y,
        col    = grDevices::adjustcolor(zero_plane_col, alpha.f = 0.25),
        border = grDevices::adjustcolor(zero_plane_col, alpha.f = 0.70)
      )
    }

  }  # end loop over coefficients

  # ---------------------------------------------------------------------------
  # 8.  Shared outer title
  # ---------------------------------------------------------------------------
  graphics::mtext(main,
                  side  = 3L,
                  outer = TRUE,
                  cex   = 1.15,
                  font  = 2L,
                  line  = 1L)

  # ---------------------------------------------------------------------------
  # 9.  Colour legend (right-hand column)
  # ---------------------------------------------------------------------------
  if (legend) {

    graphics::par(mar = c(4L, 1L, 4L, 3L))

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

    # Tick marks and labels on the right axis of the legend
    n_ticks  <- 5L
    tick_at  <- seq(z_min, z_max, length.out = n_ticks)
    tick_lab <- round(tick_at, digits)

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
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' surface3d(fit)
#' surface3d(fit, palette = "RdBu", add_zero_plane = TRUE, theta = 40)
#' surface3d(fit, what = "pvalue", clim = c(0, 1), palette = "YlGnBu")
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

  ## Resolve coefficient name and extract the raw m x K surface
  coef     <- .qqr_match_coef(x, coef)
  surf_mat <- .qqr_get_surface(x, coef, what)
  if (nrow(surf_mat) < 2L || ncol(surf_mat) < 2L)
    stop("surface3d() needs at least two grid points and two taus.",
         call. = FALSE)
  if (all(!is.finite(surf_mat)))
    stop("The surface contains no finite values.", call. = FALSE)

  ## x-axis: quantile levels theta of the grid points; y-axis: taus
  x_axis <- .qqr_thetas(x)
  y_axis <- x$taus

  ## -----------------------------------------------------------------------
  ## Colour mapping
  ## -----------------------------------------------------------------------
  if (!is.null(clim) &&
      (!is.numeric(clim) || length(clim) != 2L || !(clim[1L] < clim[2L])))
    stop("'clim' must be a numeric vector c(min, max) with min < max.",
         call. = FALSE)
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
    xlab <- paste0("Quantile of ", x$independent, " (theta)")

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

## ============================================================================
## qqr plotting functions -- plotly engine, Option-C axis convention
##
## Axis convention (consistent across ALL four functions):
##   x-axis = quantile regressor grid  z_i  (the independent/predictor axis)
##   y-axis = quantile levels  tau_j         (the dependent variable axis)
##   colour = estimated coefficient / p-value / dependence measure
##
## This places the predictor on the horizontal axis and the outcome quantile
## on the vertical axis, matching the standard regression display convention.
## qqr_heatmap and qqr_cor_heatmap share the same orientation so they can
## be placed side-by-side for direct visual comparison.
##
## Saving plotly output:
##   HTML (interactive) : htmlwidgets::saveWidget(p, "plot.html", selfcontained = TRUE)
##   PNG  (static)      : plotly::save_image(p, "plot.png", width = 1200, height = 800)
##                        Requires Python kaleido once:
##                        reticulate::py_install("kaleido")
## ============================================================================


## ============================================================================
## INTERNAL HELPERS (not exported, not documented)
## ============================================================================

## Manual RdBu colorscale  --  blue = negative, white = zero, red = positive
.qqr_rdbu <- function() {
  list(
    list(0.00, "#053061"), list(0.10, "#2166ac"),
    list(0.20, "#4393c3"), list(0.30, "#92c5de"),
    list(0.40, "#d1e5f0"), list(0.50, "#f7f7f7"),
    list(0.60, "#fddbc7"), list(0.70, "#f4a582"),
    list(0.80, "#d6604d"), list(0.90, "#b2182b"),
    list(1.00, "#67001f")
  )
}

## Map a palette name to an explicit plotly colorscale, always ordered from
## the lowest value (position 0) to the highest (position 1), so that no
## "reversescale" is needed. For "RdBu" (and other blue-red diverging
## palettes) negative values are blue and positive values red.
##   - "RdBu": the hand-coded ColorBrewer scale above;
##   - any grDevices::hcl.pals() name (case-insensitive), e.g. "viridis",
##     "plasma", "YlOrRd", "BrBG": built from grDevices::hcl.colors();
##   - anything else (e.g. "Jet", "Hot", "Portland") is passed to plotly
##     as a built-in colorscale name.
## @noRd
.qqr_colorscale <- function(palette, diverging_palettes = .qqr_diverging()) {
  if (is.list(palette)) return(palette)
  if (!is.character(palette) || length(palette) != 1L)
    stop("'palette' must be a single character string.", call. = FALSE)
  if (palette %in% c("RdBu", "RdBu_r"))
    return(.qqr_rdbu())
  pals <- grDevices::hcl.pals()
  hit  <- pals[tolower(pals) == tolower(palette)]
  if (length(hit)) {
    cols <- grDevices::hcl.colors(11L, palette = hit[1L])
    ## hcl.colors() returns red-to-blue for "RdYlBu"-type diverging
    ## palettes; flip those so that low = blue, high = red.
    if (hit[1L] %in% c("RdYlBu", "RdBu", "RdYlGn", "Spectral", "RdGy"))
      cols <- rev(cols)
    pos <- seq(0, 1, length.out = length(cols))
    return(lapply(seq_along(cols), function(i) list(pos[i], cols[i])))
  }
  palette
}

## Canonical list of diverging palette names recognised by this package
## @noRd
.qqr_diverging <- function() {
  c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr",
    "RdGy", "RdYlBu", "RdYlGn", "Spectral",
    "BlueRed", "RdBu_r", "Blue-Red", "Blue-Red 2", "Blue-Red 3",
    "Red-Green", "Purple-Green", "Purple-Brown", "Green-Brown",
    "Tropic", "Cork", "Vik", "Berlin", "Lisbon", "Tofino", "Broc")
}

## Resolve a "coef" argument that may be NULL (focal regressor), "all", or
## one or more coefficient names. Returns stored surface names.
## @noRd
.qqr_resolve_coefs <- function(x, coef, allow_all = TRUE) {
  if (is.null(coef)) return(.qqr_match_coef(x, x$independent))
  if (allow_all && length(coef) == 1L && identical(tolower(coef), "all"))
    return(names(x$coef_surfaces))
  vapply(coef, function(cn) .qqr_match_coef(x, cn), character(1L),
         USE.NAMES = FALSE)
}

## Colour limits, optionally symmetric around zero.
## @noRd
.qqr_zlim <- function(vals, clim, center_zero) {
  if (!is.null(clim)) {
    if (!is.numeric(clim) || length(clim) != 2L || !(clim[1L] < clim[2L]))
      stop("'clim' must be a numeric vector c(min, max) with min < max.",
           call. = FALSE)
    lim <- clim
  } else {
    vals <- vals[is.finite(vals)]
    if (!length(vals))
      stop("The surface contains no finite values to plot.", call. = FALSE)
    lim <- range(vals)
  }
  if (isTRUE(center_zero)) lim <- c(-1, 1) * max(abs(lim))
  if (lim[1L] == lim[2L]) lim <- lim + c(-0.5, 0.5)
  lim
}


## ============================================================================
## 1.  qqr_heatmap()
## ============================================================================

#' Heatmap of QQR Coefficient or P-value Surface
#'
#' @description
#' Produces an interactive 2-D heatmap of a quantile-on-quantile regression
#' coefficient or p-value surface using \pkg{plotly}.  Each cell at position
#' (\eqn{z_i}, \eqn{\tau_j}) is coloured according to the estimated
#' coefficient (or p-value) at that joint quantile state.  Optional
#' significance overlays (stars or dots) mark cells where the null hypothesis
#' is rejected at the chosen threshold.
#'
#' @details
#' \strong{Axis convention.}
#' The x-axis shows the quantile regressor grid values \eqn{z_i} (the
#' independent variable) and the y-axis shows the quantile levels \eqn{\tau_j}
#' of the dependent variable.  This places the predictor on the horizontal
#' axis, matching the standard convention for regression displays, and is
#' consistent with \code{\link{qqr_cor_heatmap}} so the two plots can be
#' placed side-by-side for direct comparison.
#'
#' \strong{Colour scale.}
#' The heatmap is the most compact and print-friendly way to display a QQR
#' surface.  Unlike a 3-D perspective plot, every cell value is unambiguous.
#' Diverging palettes (e.g.\ \code{"RdBu"}) with \code{center_zero = TRUE}
#' place white at exactly \eqn{\beta = 0}, making sign changes immediately
#' visible.  Blue cells indicate negative coefficients; red cells positive.
#'
#' \strong{Significance overlay.}
#' When \code{sig_overlay = "stars"}, star labels (\code{***}, \code{**},
#' \code{*}) are printed inside cells that fall below the corresponding
#' p-value threshold.  \code{sig_overlay = "dots"} prints a filled bullet
#' (\code{o}) in any significant cell.  Both options require a model
#' fitted with \code{se_type != "none"} so that p-value surfaces are
#' available.
#'
#' \strong{Multiple coefficients.}
#' When \code{coef = "all"}, a subplot grid is produced -- one panel per
#' coefficient (intercept, focal regressor, and controls) -- all sharing the
#' same colour scale for direct visual comparison.
#'
#' \strong{Saving output.}
#' The function returns a \pkg{plotly} object.  To persist it:
#' \itemize{
#'   \item Interactive HTML:
#'     \code{htmlwidgets::saveWidget(p, "heatmap.html", selfcontained = TRUE)}
#'   \item Static PNG (requires Python \code{kaleido}):
#'     \code{plotly::save_image(p, "heatmap.png", width = 1200, height = 800)}
#'   \item Install kaleido once:
#'     \code{reticulate::py_install("kaleido")}
#' }
#'
#' The concept of heatmap visualisation for QQR surfaces was introduced by
#' Roudane (2026) using \pkg{plotly}.  This implementation extends that idea
#' with significance overlays, symmetric zero-centred scaling, multi-coefficient
#' subplots, and a corrected axis convention.
#'
#' @param x           A fitted object of class \code{"qqr"}, as returned by
#'                    \code{\link{qqr}}.
#' @param coef        Character; name of the coefficient surface to plot.
#'                    Pass the exact variable name (e.g.\ \code{"oil"}),
#'                    \code{"intercept"} or \code{"(Intercept)"} for the
#'                    model intercept, or \code{"all"} to produce one subplot
#'                    per coefficient.  Default is the focal quantile regressor
#'                    (\code{x$independent}).
#' @param what        Character; which surface to display.
#'                    \code{"coef"} (default) for estimated coefficients,
#'                    \code{"pvalue"} for the associated p-values.  P-value
#'                    surfaces require the model to have been fitted with
#'                    \code{se_type != "none"}.
#' @param palette     Character; a plotly or RColorBrewer palette name.
#'                    Recommended diverging palettes (for coefficients):
#'                    \code{"RdBu"} (default; blue = negative, red = positive),
#'                    \code{"PuOr"}, \code{"BrBG"}.  Sequential palettes
#'                    (for p-values or one-sided effects): \code{"Viridis"},
#'                    \code{"Plasma"}, \code{"YlOrRd"}.  Any valid RColorBrewer
#'                    name or plotly built-in string is accepted (e.g.\
#'                    \code{"Jet"}, \code{"Hot"}).
#' @param center_zero Logical; if \code{TRUE} the colour midpoint is fixed at
#'                    zero so that blue and red are symmetrically opposite
#'                    around white.  Default: \code{TRUE} when
#'                    \code{what = "coef"} and a diverging palette is used;
#'                    \code{FALSE} otherwise.
#' @param clim        Numeric vector \code{c(min, max)} fixing the colour scale
#'                    limits across all cells.  Useful for comparing heatmaps
#'                    from different models on the same scale.  \code{NULL}
#'                    (default) uses the full range of the surface values.
#' @param sig_overlay Character; add a significance overlay inside each cell.
#'   \describe{
#'     \item{\code{"none"}}{No overlay (default).}
#'     \item{\code{"stars"}}{Prints \code{***} / \code{**} / \code{*} inside
#'       cells whose p-value falls below the corresponding threshold in
#'       \code{sig_levels}.  Thresholds are tested in order; the first one
#'       satisfied is used.}
#'     \item{\code{"dots"}}{Prints a filled bullet (\code{o}) inside
#'       any significant cell, regardless of how far below the threshold the
#'       p-value falls.}
#'   }
#'   Both \code{"stars"} and \code{"dots"} require \code{se_type != "none"}
#'   in the fitted model.
#' @param sig_levels  Named numeric vector of p-value thresholds for star
#'                    labels.  Names are the labels printed; values are the
#'                    cut-offs.  Default:
#'                    \code{c("***" = 0.01, "**" = 0.05, "*" = 0.10)}.
#' @param sig_size    Numeric; font size for significance annotations in
#'                    plotly units.  Default \code{11L}.
#' @param sig_colour  Character; text colour for significance annotations.
#'                    Default \code{"black"}.  Use \code{"white"} on dark
#'                    colour scales such as \code{"Viridis"}.
#' @param grid_values Logical; if \code{TRUE} (default) the x-axis shows the
#'                    actual quantile regressor values \eqn{z_i}.  If
#'                    \code{FALSE} the grid point index (1, 2, \ldots, m) is
#'                    shown instead -- useful when regressor values are long,
#'                    irregular, or confidential.
#' @param xlab        Character; x-axis label.  \code{NULL} auto-generates a
#'                    label using the regressor name and \eqn{z_i} notation.
#' @param ylab        Character; y-axis label.  \code{NULL} auto-generates
#'                    \code{"Quantile level (tau)"}.
#' @param title       Character; plot title.  \code{NULL} auto-generates a
#'                    title of the form \code{"QQR Coefficient Surface -- y"}.
#' @param digits      Integer; decimal places for axis tick labels and
#'                    hover-tooltip values.  Default \code{2L}.
#' @param snapshot    Logical; if \code{TRUE} and Python \code{kaleido} is
#'                    installed, saves a static PNG to \code{snapshot_file}
#'                    alongside the interactive widget.  Default \code{FALSE}.
#' @param snapshot_file Character; file path for the static PNG.  Used only
#'                    when \code{snapshot = TRUE}.  Default
#'                    \code{"qqr_heatmap.png"}.
#' @param ...         Currently unused; reserved for future arguments.
#'
#' @return A \code{plotly} object.  The object displays automatically in the
#'   RStudio Viewer pane, can be embedded in R Markdown / Quarto documents,
#'   and can be saved with \code{htmlwidgets::saveWidget()} or
#'   \code{plotly::save_image()}.
#'
#' @seealso
#'   \code{\link{qqr_contour}} for contour plots of the same surface;
#'   \code{\link{qqr_plotly}} for a unified 3-D / heatmap / contour hub;
#'   \code{\link{qqr_cor_heatmap}} for a pre-model quantile dependence
#'   diagnostic; \code{\link{panel_surface}} for a base-graphics ribbon
#'   panel plot.
#'
#' @references
#' Roudane, M. (2026). \emph{QuantileOnQuantile: Quantile-on-Quantile
#'   Regression Analysis}. R package version 1.0.3. CRAN.
#'   \url{https://cran.r-project.org/package=QuantileOnQuantile}.
#'   The heatmap visualisation concept is adapted from this work; the
#'   implementation has been substantially extended.
#'
#' Sim, N. and Zhou, H. (2015). Oil prices, US stock return and the dependence
#'   between their quantiles. \emph{Journal of Banking & Finance}, 55, 1--8.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   qqr_heatmap(fit)
#'   qqr_heatmap(fit, sig_overlay = "stars")
#'   qqr_heatmap(fit, coef = "all")
#' }
#' \dontrun{
#' p <- qqr_heatmap(fit, sig_overlay = "stars")
#' htmlwidgets::saveWidget(p, "heatmap.html", selfcontained = TRUE)
#' ## Static PNG (requires the Python 'kaleido' package)
#' qqr_heatmap(fit, snapshot = TRUE, snapshot_file = "heatmap.png")
#' }
#'
#' @export
qqr_heatmap <- function(x,
                        coef          = NULL,
                        what          = c("coef", "pvalue"),
                        palette       = "RdBu",
                        center_zero   = NULL,
                        clim          = NULL,
                        sig_overlay   = c("none", "stars", "dots"),
                        sig_levels    = c("***" = 0.01, "**" = 0.05,
                                          "*"   = 0.10),
                        sig_size      = 11L,
                        sig_colour    = "black",
                        grid_values   = TRUE,
                        xlab          = NULL,
                        ylab          = NULL,
                        title         = NULL,
                        digits        = 2L,
                        snapshot      = FALSE,
                        snapshot_file = "qqr_heatmap.png",
                        ...) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  if (!inherits(x, "qqr"))
    stop("'x' must be a fitted qqr object.")

  what        <- match.arg(what)
  sig_overlay <- match.arg(sig_overlay)

  ## -- 1. Resolve coefficient(s) ---------------------------------------------
  plot_coefs <- .qqr_resolve_coefs(x, coef, allow_all = TRUE)

  ## -- 2. Surface list -------------------------------------------------------
  surfaces <- stats::setNames(
    lapply(plot_coefs, function(cn) .qqr_get_surface(x, cn, what)),
    plot_coefs)
  if (sig_overlay != "none")
    for (cn in plot_coefs) .qqr_get_surface(x, cn, "pvalue")

  ## -- 3. Shared colour limits -----------------------------------------------
  diverging_palettes <- .qqr_diverging()
  is_diverging       <- is.character(palette) &&
    palette %in% diverging_palettes
  if (is.null(center_zero))
    center_zero <- (what == "coef" && is_diverging)

  lim   <- .qqr_zlim(unlist(surfaces), clim, center_zero)
  z_min <- lim[1L]; z_max <- lim[2L]

  ## -- 4. Axis labels --------------------------------------------------------
  ## Option C: x = z_i (regressor grid), y = tau (outcome quantile)
  if (is.null(xlab)) xlab <- if (grid_values)
    paste0(x$independent, " grid value (z\u1d62)") else "Grid point index"
  if (is.null(ylab)) ylab <- "Quantile level (\u03c4)"
  if (is.null(title)) {
    type_label <- if (what == "coef") "Coefficient" else "P-value"
    title <- paste0("QQR ", type_label, " Surface \u2014 ", x$dependent)
  }
  colorbar_title <- if (what == "coef") "Coeff." else "P-val."
  cs             <- .qqr_colorscale(palette, diverging_palettes)

  ## -- 5. Per-coefficient trace builder --------------------------------------
  ## mat is stored as  m x J  (rows = grid/X,  cols = tau/Y).
  ## Option C: x-axis = grid (length m), y-axis = tau (length J).
  ## plotly heatmap: z[row][col] -> (x[col], y[row]).
  ## So pass z = t(mat)  (J x m) so that t(mat)[j,i] at (x_vals[i], y_vals[j]).

  .make_trace <- function(cn) {
    mat    <- surfaces[[cn]]                                          ## m x J
    m      <- nrow(mat)
    x_vals <- if (grid_values) round(x$grid, digits) else seq_len(m) ## m
    y_vals <- round(x$taus, digits)                                   ## J
    z_plot <- t(mat)                                                  ## J x m

    ## Significance annotations
    anns <- list()
    if (sig_overlay != "none") {
      pmat <- x$pval_surfaces[[cn]]
      for (i in seq_len(m)) {
        for (j in seq_along(x$taus)) {
          pv  <- pmat[i, j]
          lbl <- ""
          for (nm in names(sig_levels))
            if (!is.na(pv) && pv < sig_levels[nm] && lbl == "") lbl <- nm
          if (sig_overlay == "dots" && lbl != "") lbl <- "\u2022"
          if (lbl != "")
            anns[[length(anns) + 1L]] <- list(
              x         = x_vals[i],   ## grid on x-axis
              y         = y_vals[j],   ## tau on y-axis
              text      = lbl,
              showarrow = FALSE,
              font      = list(size = sig_size, color = sig_colour)
            )
        }
      }
    }

    hover <- paste0(
      xlab,          ": %{x:.", digits,     "f}<br>",
      ylab,          ": %{y:.", digits,     "f}<br>",
      colorbar_title, ": %{z:.", digits+2L, "f}<extra>", cn, "</extra>"
    )

    list(
      trace = plotly::plot_ly(
        x             = x_vals,
        y             = y_vals,
        z             = z_plot,
        type          = "heatmap",
        colorscale    = cs,
        showscale     = TRUE,
        zmin          = z_min,
        zmax          = z_max,
        zmid          = if (center_zero) 0 else NULL,
        colorbar      = list(title      = colorbar_title,
                             tickformat = paste0(".", digits, "f"),
                             len        = 0.85),
        hovertemplate = hover
      ),
      anns  = anns,
      label = cn
    )
  }

  traces <- lapply(plot_coefs, .make_trace)

  ## -- 6. Assemble layout ----------------------------------------------------
  if (length(traces) == 1L) {

    p    <- traces[[1L]]$trace
    anns <- traces[[1L]]$anns
    p <- p |>
      plotly::layout(
        title         = list(text = paste0("<b>", title, "</b>"),
                             x = 0, font = list(size = 14L)),
        paper_bgcolor = "white",
        plot_bgcolor  = "white",
        margin        = list(l = 70, r = 20, b = 80, t = 70),
        xaxis         = list(title      = xlab,
                             tickformat = paste0(".", digits, "f"),
                             tickangle  = 45L,
                             tickfont   = list(size = 9L)),
        yaxis         = list(title      = ylab,
                             tickformat = paste0(".", digits, "f"),
                             tickfont   = list(size = 9L)),
        annotations   = anns
      )

  } else {
    n_cols    <- min(3L, length(traces))
    p_list    <- lapply(traces, `[[`, "trace")
    all_anns  <- lapply(traces, `[[`, "anns")
    subtitles <- lapply(traces, `[[`, "label")

    p <- plotly::subplot(p_list,
                         nrows  = ceiling(length(p_list) / n_cols),
                         shareX = FALSE, shareY = FALSE,
                         titleX = TRUE,  titleY = TRUE,
                         margin = 0.06)

    sub_title_anns <- lapply(seq_along(subtitles), function(k) {
      list(
        text      = paste0("<b>", subtitles[[k]], "</b>"),
        x         = ((k - 1L) %% n_cols + 0.5) / n_cols,
        y         = 1 - floor((k - 1L) / n_cols) /
          ceiling(length(traces) / n_cols),
        xref      = "paper", yref = "paper",
        showarrow = FALSE, font = list(size = 11L)
      )
    })

    p <- p |>
      plotly::layout(
        title         = list(text = paste0("<b>", title, "</b>"),
                             x = 0, font = list(size = 14L)),
        paper_bgcolor = "white",
        margin        = list(l = 60, r = 20, b = 60, t = 80),
        annotations   = c(unlist(all_anns, recursive = FALSE),
                          sub_title_anns)
      )
  }

  ## -- 7. Optional static snapshot -------------------------------------------
  if (snapshot) {
    tryCatch(
      plotly::save_image(p, snapshot_file),
      error = function(e) message(
        "Snapshot failed: ", e$message,
        "\nInstall kaleido: reticulate::py_install('kaleido')"
      )
    )
  }

  p
}


## ============================================================================
## 2.  qqr_contour()
## ============================================================================

#' Contour Plot of QQR Coefficient or P-value Surface
#'
#' @description
#' Produces an interactive filled contour plot of a quantile-on-quantile
#' regression coefficient or p-value surface using \pkg{plotly}.  Contour
#' lines connect points of equal coefficient value across the joint
#' (\eqn{z_i}, \eqn{\tau_j}) grid.  The zero contour is automatically
#' highlighted in bold black to make sign-change boundaries immediately
#' visible.
#'
#' @details
#' \strong{Axis convention.}
#' The x-axis shows the quantile regressor grid values \eqn{z_i} (the
#' independent variable) and the y-axis shows the quantile levels
#' \eqn{\tau_j} of the dependent variable.  This is consistent with
#' \code{\link{qqr_heatmap}} and \code{\link{qqr_cor_heatmap}} so all three
#' plots are directly comparable when placed side-by-side.
#'
#' \strong{Contour interpretation.}
#' Contour plots complement heatmaps: while a heatmap shows cell-level
#' values, contour lines emphasise regions of constant effect and make
#' gradient directions visible.  The bold zero contour (\code{zero_contour =
#' TRUE}) marks the exact boundary where the coefficient changes sign, which
#' is often the most economically meaningful feature of the surface.
#'
#' \strong{Significance mask.}
#' When \code{sig_mask = TRUE}, a semi-transparent grey overlay is drawn over
#' all cells whose p-value exceeds \code{sig_threshold}.  This visually
#' de-emphasises statistically noisy regions without discarding them.
#' Requires a model fitted with \code{se_type != "none"}.
#'
#' \strong{Saving output.}
#' \itemize{
#'   \item Interactive HTML:
#'     \code{htmlwidgets::saveWidget(p, "contour.html", selfcontained = TRUE)}
#'   \item Static PNG (requires Python \code{kaleido}):
#'     \code{plotly::save_image(p, "contour.png", width = 1000, height = 700)}
#' }
#'
#' The concept of contour visualisation for QQR surfaces is adapted from
#' Roudane (2026); the implementation has been substantially extended.
#'
#' @param x           A fitted object of class \code{"qqr"}, as returned by
#'                    \code{\link{qqr}}.
#' @param coef        Character; name of the coefficient surface to plot.
#'                    Pass the exact variable name (e.g.\ \code{"oil"}),
#'                    \code{"intercept"}, or \code{"(Intercept)"}.  Only one
#'                    coefficient at a time is supported (unlike
#'                    \code{\link{qqr_heatmap}} which accepts \code{"all"}).
#'                    Default: focal quantile regressor (\code{x$independent}).
#' @param what        Character; \code{"coef"} (default) for coefficient
#'                    values, \code{"pvalue"} for p-values.
#' @param n_levels    Integer; approximate number of contour bands.  More
#'                    levels give finer resolution but a busier plot.
#'                    Default \code{10L}.
#' @param palette     Character; plotly / RColorBrewer palette name.
#'                    Diverging default \code{"RdBu"}.  Sequential
#'                    alternatives: \code{"Viridis"}, \code{"Plasma"},
#'                    \code{"YlOrRd"}.
#' @param center_zero Logical; if \code{TRUE} (default when using a diverging
#'                    palette with \code{what = "coef"}) the colour midpoint
#'                    is fixed at zero.
#' @param clim        Numeric \code{c(min, max)}; fixed colour limits.
#'                    \code{NULL} uses the data range.
#' @param zero_contour Logical; if \code{TRUE} (default) a bold black
#'                    contour line is drawn at exactly \eqn{\beta = 0},
#'                    marking where the effect changes sign.  Only drawn when
#'                    zero falls within the colour range.
#' @param zero_lwd    Numeric; line width for the zero contour.
#'                    Default \code{2.5}.
#' @param show_labels Logical; if \code{TRUE} (default) contour level values
#'                    are printed as inline labels along the contour lines.
#' @param sig_mask    Logical; if \code{TRUE} a semi-transparent grey overlay
#'                    is drawn over cells where the p-value exceeds
#'                    \code{sig_threshold}, visually de-emphasising
#'                    non-significant regions.  Requires \code{se_type != "none"}.
#'                    Default \code{FALSE}.
#' @param sig_threshold Numeric; p-value threshold for the significance mask.
#'                    Cells with p > \code{sig_threshold} are greyed out when
#'                    \code{sig_mask = TRUE}.  Default \code{0.10}.
#' @param grid_values Logical; if \code{TRUE} (default) the x-axis shows
#'                    actual quantile regressor values \eqn{z_i}.  If
#'                    \code{FALSE} the grid index is shown.
#' @param xlab        Character; x-axis label.  \code{NULL} auto-generates
#'                    a label using the regressor name and \eqn{z_i} notation.
#' @param ylab        Character; y-axis label.  \code{NULL} auto-generates
#'                    \code{"Quantile level (tau)"}.
#' @param title       Character; plot title.  \code{NULL} auto-generates.
#' @param digits      Integer; decimal places for axis labels and hover
#'                    tooltips.  Default \code{2L}.
#' @param snapshot    Logical; save a static PNG snapshot.  Default
#'                    \code{FALSE}.
#' @param snapshot_file Character; file path for the PNG.  Default
#'                    \code{"qqr_contour.png"}.
#' @param ...         Currently unused.
#'
#' @return A \code{plotly} object.  Displays automatically in the RStudio
#'   Viewer pane and can be saved with \code{htmlwidgets::saveWidget()} or
#'   \code{plotly::save_image()}.
#'
#' @seealso
#'   \code{\link{qqr_heatmap}} for a tile-based display of the same surface;
#'   \code{\link{qqr_plotly}} for a unified 3-D / heatmap / contour hub;
#'   \code{\link{qqr_cor_heatmap}} for a pre-model dependence diagnostic.
#'
#' @references
#' Roudane, M. (2026). \emph{QuantileOnQuantile: Quantile-on-Quantile
#'   Regression Analysis}. R package version 1.0.3. CRAN.
#'   The contour visualisation concept is adapted from this work; the
#'   implementation has been substantially extended.
#'
#' Sim, N. and Zhou, H. (2015). Oil prices, US stock return and the dependence
#'   between their quantiles. \emph{Journal of Banking & Finance}, 55, 1--8.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   qqr_contour(fit)
#'   qqr_contour(fit, sig_mask = TRUE, sig_threshold = 0.05)
#' }
#' \dontrun{
#' p <- qqr_contour(fit)
#' htmlwidgets::saveWidget(p, "contour.html", selfcontained = TRUE)
#' }
#'
#' @export
qqr_contour <- function(x,
                        coef          = NULL,
                        what          = c("coef", "pvalue"),
                        n_levels      = 10L,
                        palette       = "RdBu",
                        center_zero   = NULL,
                        clim          = NULL,
                        zero_contour  = TRUE,
                        zero_lwd      = 2.5,
                        show_labels   = TRUE,
                        sig_mask      = FALSE,
                        sig_threshold = 0.10,
                        grid_values   = TRUE,
                        xlab          = NULL,
                        ylab          = NULL,
                        title         = NULL,
                        digits        = 2L,
                        snapshot      = FALSE,
                        snapshot_file = "qqr_contour.png",
                        ...) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  if (!inherits(x, "qqr"))
    stop("'x' must be a fitted qqr object.")

  what <- match.arg(what)

  ## -- 1. Resolve coefficient ------------------------------------------------
  coef <- .qqr_resolve_coefs(x, if (is.null(coef)) NULL else coef[1L],
                             allow_all = FALSE)

  ## -- 2. Extract surface ----------------------------------------------------
  mat <- .qqr_get_surface(x, coef, what)   ## m x J
  m   <- nrow(mat)
  if (m < 2L || ncol(mat) < 2L)
    stop("qqr_contour() needs at least two grid points and two taus.",
         call. = FALSE)

  ## Option C: x = grid (z_i),  y = tau
  x_vals <- if (grid_values) round(x$grid, digits) else seq_len(m)
  y_vals <- round(x$taus, digits)
  z_plot <- t(mat)          ## J x m  so z_plot[j,i] at (x_vals[i], y_vals[j])

  ## -- 3. Colour limits ------------------------------------------------------
  diverging_palettes <- .qqr_diverging()
  is_diverging       <- is.character(palette) &&
    palette %in% diverging_palettes
  if (is.null(center_zero))
    center_zero <- (what == "coef" && is_diverging)

  lim   <- .qqr_zlim(mat, clim, center_zero)
  z_min <- lim[1L]; z_max <- lim[2L]

  contour_size <- (z_max - z_min) / max(n_levels, 2L)
  cs           <- .qqr_colorscale(palette, diverging_palettes)

  ## -- 4. Labels -------------------------------------------------------------
  if (is.null(xlab)) xlab <- if (grid_values)
    paste0(x$independent, " grid value (z\u1d62)") else "Grid point index"
  if (is.null(ylab)) ylab <- "Quantile level (\u03c4)"
  z_label <- if (what == "coef") "Coefficient" else "P-value"
  if (is.null(title))
    title <- paste0("QQR ", z_label, " Contour \u2014 ",
                    x$dependent, " ~ ", coef)

  ## -- 5. Main filled contour trace ------------------------------------------
  p <- plotly::plot_ly(
    x             = x_vals,
    y             = y_vals,
    z             = z_plot,
    type          = "contour",
    colorscale    = cs,
    showscale     = TRUE,
    zmin          = z_min,
    zmax          = z_max,
    zmid          = if (center_zero) 0 else NULL,
    colorbar      = list(title      = z_label,
                         tickformat = paste0(".", digits, "f"),
                         len        = 0.85),
    contours      = list(
      start      = z_min,
      end        = z_max,
      size       = contour_size,
      showlabels = show_labels,
      labelfont  = list(size = 9L, color = "black"),
      coloring   = "heatmap"
    ),
    line          = list(smoothing = 0.85,
                         color     = "rgba(80,80,80,0.6)",
                         width     = 0.8),
    hovertemplate = paste0(
      xlab,    ": %{x:.", digits,     "f}<br>",
      ylab,    ": %{y:.", digits,     "f}<br>",
      z_label, ": %{z:.", digits+2L, "f}<extra></extra>"
    )
  )

  ## -- 6. Zero contour overlay -----------------------------------------------
  if (zero_contour && z_min <= 0 && z_max >= 0) {
    p <- p |>
      plotly::add_contour(
        x          = x_vals,
        y          = y_vals,
        z          = z_plot,
        showscale  = FALSE,
        contours   = list(start = 0, end = 0, size = 1,
                          showlabels = FALSE, coloring = "lines"),
        line       = list(color = "black", width = zero_lwd),
        hoverinfo  = "skip",
        name       = paste0(z_label, " = 0"),
        inherit    = FALSE
      )
  }

  ## -- 7. Significance mask --------------------------------------------------
  if (sig_mask) {
    pmat     <- .qqr_get_surface(x, coef, "pvalue")
    mask_mat <- t(ifelse(!is.na(pmat) & pmat > sig_threshold, 0.5, NA_real_))

    p <- p |>
      plotly::add_heatmap(
        x          = x_vals,
        y          = y_vals,
        z          = mask_mat,
        showscale  = FALSE,
        colorscale = list(list(0, "rgba(160,160,160,0.45)"),
                          list(1, "rgba(160,160,160,0.45)")),
        hoverinfo  = "skip",
        name       = paste0("p > ", sig_threshold),
        inherit    = FALSE
      )
  }

  ## -- 8. Layout -------------------------------------------------------------
  p <- p |>
    plotly::layout(
      title         = list(text = paste0("<b>", title, "</b>"),
                           x = 0, font = list(size = 14L)),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      margin        = list(l = 70, r = 20, b = 80, t = 70),
      xaxis         = list(title      = xlab,
                           tickformat = paste0(".", digits, "f"),
                           tickangle  = 45L,
                           tickfont   = list(size = 9L)),
      yaxis         = list(title      = ylab,
                           tickformat = paste0(".", digits, "f"),
                           tickfont   = list(size = 9L))
    )

  ## -- 9. Optional snapshot --------------------------------------------------
  if (snapshot) {
    tryCatch(
      plotly::save_image(p, snapshot_file),
      error = function(e) message("Snapshot failed: ", e$message)
    )
  }

  p
}


## ============================================================================
## 3.  qqr_plotly()
## ============================================================================

#' Interactive Plotly Visualisation Hub for QQR Objects
#'
#' @description
#' A unified interactive visualisation hub for quantile-on-quantile regression
#' surfaces using \pkg{plotly}.  Combines three display types in a single
#' function call:
#' \describe{
#'   \item{\code{"surface"}}{Interactive 3-D perspective surface.  The user
#'     can rotate, zoom, and hover-inspect any point.}
#'   \item{\code{"heatmap"}}{2-D colour grid.  Identical layout to
#'     \code{\link{qqr_heatmap}}.}
#'   \item{\code{"contour"}}{2-D filled iso-line plot.  Identical layout to
#'     \code{\link{qqr_contour}}.}
#' }
#' All three types share the same data extraction pipeline and the same
#' axis convention: the x-axis shows the quantile regressor grid (\eqn{z_i})
#' and the y-axis shows the quantile levels (\eqn{\tau}).
#'
#' @details
#' \strong{Axis convention (Option C).}
#' The x-axis shows the quantile regressor grid values \eqn{z_i} and the
#' y-axis shows the quantile levels \eqn{\tau_j}.  For the 3-D surface,
#' the regressor grid runs along the x-axis and tau along the y-axis.
#' This is consistent across all three plot types and with
#' \code{\link{qqr_heatmap}} and \code{\link{qqr_cor_heatmap}}, so any
#' combination of plots can be compared directly.
#'
#' \strong{Colour scales.}
#' \itemize{
#'   \item \code{"Viridis"} (default): perceptually uniform, colourblind-safe,
#'     prints well in greyscale.
#'   \item \code{"RdBu"}: diverging blue-white-red.  Best for coefficients
#'     that change sign.  Blue = negative, red = positive.
#'   \item \code{"Plasma"}: warm perceptually uniform palette.
#'   \item \code{"Jet"}: classic MATLAB rainbow.  Not colourblind-safe but
#'     familiar to many readers and good for presentations.
#'   \item Any valid plotly built-in name or RColorBrewer name is accepted.
#' }
#'
#' \strong{Zero plane (surface type only).}
#' When \code{show_zero_plane = TRUE}, a semi-transparent red plane is drawn
#' at \eqn{z = 0}, making it immediately clear where the coefficient surface
#' crosses from positive to negative territory.  Adjust opacity with
#' \code{zero_plane_opacity}.
#'
#' \strong{Saving output.}
#' \itemize{
#'   \item \code{htmlwidgets::saveWidget(p, "plot.html", selfcontained = TRUE)}
#'   \item \code{plotly::save_image(p, "plot.png")}
#'     (requires \code{reticulate::py_install("kaleido")})
#' }
#'
#' @param x               A fitted object of class \code{"qqr"}, as returned
#'                        by \code{\link{qqr}}.
#' @param coef            Character; coefficient name.  Default: focal
#'                        regressor (\code{x$independent}).
#' @param what            Character; \code{"coef"} (default) for coefficient
#'                        values or \code{"pvalue"} for p-values.
#' @param type            Character; plot type to produce.
#'   \describe{
#'     \item{\code{"surface"}}{Interactive 3-D perspective plot (default).}
#'     \item{\code{"heatmap"}}{2-D colour tile grid.}
#'     \item{\code{"contour"}}{2-D filled iso-line plot.}
#'   }
#' @param colorscale      Character or named list; plotly colour scale.
#'                        \code{"Viridis"} (default), \code{"RdBu"},
#'                        \code{"Plasma"}, \code{"Jet"}, \code{"Hot"}, or
#'                        any valid plotly / RColorBrewer name.
#' @param clim            Numeric \code{c(min, max)}; fixed colour / z-axis
#'                        limits.  \code{NULL} uses the data range.
#' @param center_zero     Logical; fix colour midpoint at zero.  Default
#'                        \code{NULL} auto-detects from the palette and
#'                        \code{what}.
#' @param show_zero_plane Logical; for \code{type = "surface"} only.  Draw a
#'                        semi-transparent red plane at \eqn{z = 0}.
#'                        Default \code{FALSE}.
#' @param zero_plane_opacity Numeric in \eqn{[0, 1]}; opacity of the zero
#'                        plane.  Default \code{0.25}.
#' @param show_contour    Logical; for \code{type = "surface"} only.  Show
#'                        grid lines projected onto the surface faces.
#'                        Default \code{TRUE}.
#' @param n_contours      Integer; number of contour levels for
#'                        \code{type = "contour"}.  Default \code{10L}.
#' @param show_labels     Logical; for \code{type = "contour"} only.  Print
#'                        inline contour value labels on lines.  Default
#'                        \code{TRUE}.
#' @param sig_overlay     Logical; for \code{type = "heatmap"} only.  Add
#'                        significance text annotations inside cells.
#'                        Requires \code{se_type != "none"} in the fitted
#'                        model.  Default \code{FALSE}.
#' @param sig_threshold   Numeric; p-value threshold for significance
#'                        annotation.  Cells with p < \code{sig_threshold}
#'                        receive the \code{sig_symbol} label.
#'                        Default \code{0.05}.
#' @param sig_symbol      Character; symbol to print in significant cells
#'                        when \code{sig_overlay = TRUE}.  Default \code{"*"}.
#' @param grid_values     Logical; if \code{TRUE} (default) axes show actual
#'                        \eqn{z_i} values; if \code{FALSE} shows the grid
#'                        point index.
#' @param aspect_ratio    Numeric vector of length 3 giving the x, y, z
#'                        scene aspect ratios for \code{type = "surface"}.
#'                        Default \code{c(1, 1, 0.7)}.
#' @param camera          List; camera position for \code{type = "surface"}.
#'                        Default provides a standard elevated angled view.
#' @param xlab            Character; x-axis label.  \code{NULL}
#'                        auto-generates using the regressor name and
#'                        \eqn{z_i} notation.
#' @param ylab            Character; y-axis label.  \code{NULL}
#'                        auto-generates \code{"Quantile level (tau)"}.
#' @param title           Character; plot title.  \code{NULL} auto-generates.
#' @param digits          Integer; decimal places for axis labels and hover
#'                        tooltips.  Default \code{2L}.
#' @param snapshot        Logical; save a static PNG snapshot.  Default
#'                        \code{FALSE}.
#' @param snapshot_file   Character; PNG file path.  Default
#'                        \code{"qqr_plotly.png"}.
#' @param ...             Currently unused.
#'
#' @return A \code{plotly} object.  Displays automatically in the RStudio
#'   Viewer pane and can be saved with \code{htmlwidgets::saveWidget()} or
#'   \code{plotly::save_image()}.
#'
#' @seealso
#'   \code{\link{qqr_heatmap}}, \code{\link{qqr_contour}},
#'   \code{\link{qqr_cor_heatmap}}, \code{\link{panel_surface}},
#'   \code{\link{Int.plot}}
#'
#' @references
#' Roudane, M. (2026). \emph{QuantileOnQuantile: Quantile-on-Quantile
#'   Regression Analysis}. R package version 1.0.3. CRAN.
#'   The 3-D, heatmap, and contour visualisation concepts are adapted from
#'   this work and substantially extended.
#'
#' Sievert, C. (2020). \emph{Interactive Web-Based Data Visualization with
#'   R, plotly, and shiny}. Chapman and Hall/CRC.
#'   \url{https://plotly-r.com}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   qqr_plotly(fit, type = "surface")
#'   qqr_plotly(fit, type = "heatmap", colorscale = "RdBu", sig_overlay = TRUE)
#'   qqr_plotly(fit, type = "contour")
#' }
#' \dontrun{
#' p <- qqr_plotly(fit, type = "surface", show_zero_plane = TRUE)
#' htmlwidgets::saveWidget(p, "surface.html", selfcontained = TRUE)
#' }
#'
#' @export
qqr_plotly <- function(x,
                       coef               = NULL,
                       what               = c("coef", "pvalue"),
                       type               = c("surface", "heatmap", "contour"),
                       colorscale         = "Viridis",
                       clim               = NULL,
                       center_zero        = NULL,
                       show_zero_plane    = FALSE,
                       zero_plane_opacity = 0.25,
                       show_contour       = TRUE,
                       n_contours         = 10L,
                       show_labels        = TRUE,
                       sig_overlay        = FALSE,
                       sig_threshold      = 0.05,
                       sig_symbol         = "*",
                       grid_values        = TRUE,
                       aspect_ratio       = c(1, 1, 0.7),
                       camera             = list(
                         eye = list(x = 1.4, y = 1.7, z = 1.2)
                       ),
                       xlab               = NULL,
                       ylab               = NULL,
                       title              = NULL,
                       digits             = 2L,
                       snapshot           = FALSE,
                       snapshot_file      = "qqr_plotly.png",
                       ...) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  if (!inherits(x, "qqr"))
    stop("'x' must be a fitted qqr object.")

  what <- match.arg(what)
  type <- match.arg(type)

  ## -- 1. Resolve coefficient ------------------------------------------------
  coef <- .qqr_resolve_coefs(x, if (is.null(coef)) NULL else coef[1L],
                             allow_all = FALSE)

  ## -- 2. Extract surface ----------------------------------------------------
  mat <- .qqr_get_surface(x, coef, what)   ## m x J (rows = grid, cols = tau)
  m   <- nrow(mat)
  if (m < 2L || ncol(mat) < 2L)
    stop("qqr_plotly() needs at least two grid points and two taus.",
         call. = FALSE)

  ## Option C: x = grid (z_i),  y = tau (J)
  x_vals <- if (grid_values) round(x$grid, digits) else seq_len(m)
  y_vals <- round(x$taus, digits)

  ## -- 3. Colour limits ------------------------------------------------------
  diverging_cs <- .qqr_diverging()
  is_div <- is.character(colorscale) && colorscale %in% diverging_cs
  if (is.null(center_zero))
    center_zero <- (what == "coef" && is_div)

  lim   <- .qqr_zlim(mat, clim, center_zero)
  z_min <- lim[1L]; z_max <- lim[2L]

  ## -- 4. Labels -------------------------------------------------------------
  if (is.null(xlab)) xlab <- if (grid_values)
    paste0(x$independent, " grid value (z\u1d62)") else "Grid index"
  if (is.null(ylab)) ylab <- "Quantile level (\u03c4)"
  z_label <- if (what == "coef") "Coefficient" else "P-value"

  if (is.null(title)) {
    type_label <- switch(type,
                         surface = "3D Surface",
                         heatmap = "Heatmap",
                         contour = "Contour Plot")
    title <- paste0("QQR ", type_label, " \u2014 ",
                    x$dependent, " ~ ", coef)
  }

  cs <- .qqr_colorscale(colorscale)

  hover_tmpl <- paste0(
    xlab,    ": %{x:.", digits,     "f}<br>",
    ylab,    ": %{y:.", digits,     "f}<br>",
    z_label, ": %{z:.", digits+2L, "f}<extra></extra>"
  )

  ## --------------------------------------------------------------------------
  ## 5a.  3-D SURFACE
  ## mat = m x J.  Need z_surf[j, i] = mat[i, j]  ? t(mat) shape J x m.
  ## x = x_vals (m),  y = y_vals (J),  z = t(mat) (J x m).
  ## --------------------------------------------------------------------------
  if (type == "surface") {

    z_surf <- t(mat)
    ## Spacing of the wire-frame lines; the quantile grid is generally not
    ## equally spaced in x, so use the average spacing.
    x_step <- (max(x_vals) - min(x_vals)) / (length(x_vals) - 1L)
    y_step <- (max(y_vals) - min(y_vals)) / (length(y_vals) - 1L)

    p <- plotly::plot_ly(
      x             = x_vals,
      y             = y_vals,
      z             = z_surf,
      type          = "surface",
      colorscale    = cs,
      showscale     = TRUE,
      cmin          = z_min,
      cmax          = z_max,
      colorbar      = list(title      = z_label,
                           tickformat = paste0(".", digits, "f")),
      contours      = list(
        x = list(show = show_contour, color = "rgba(0,0,0,0.25)",
                 start = min(x_vals), end = max(x_vals), size = x_step),
        y = list(show = show_contour, color = "rgba(0,0,0,0.25)",
                 start = min(y_vals), end = max(y_vals), size = y_step),
        z = list(show = FALSE)
      ),
      lighting      = list(ambient   = 0.55, diffuse  = 0.80,
                           specular  = 0.15, roughness = 0.90),
      lightposition = list(x = 60, y = 120, z = 80),
      hovertemplate = hover_tmpl
    )

    if (show_zero_plane && z_min <= 0 && z_max >= 0) {
      zero_mat <- matrix(0, nrow = length(y_vals), ncol = length(x_vals))
      p <- p |>
        plotly::add_surface(
          x          = x_vals,
          y          = y_vals,
          z          = zero_mat,
          colorscale = list(list(0, "rgba(200,0,0,0.3)"),
                            list(1, "rgba(200,0,0,0.3)")),
          showscale  = FALSE,
          opacity    = zero_plane_opacity,
          hoverinfo  = "skip",
          name       = "z = 0"
        )
    }

    p <- p |>
      plotly::layout(
        title         = list(text = paste0("<b>", title, "</b>"),
                             x = 0, font = list(size = 14L)),
        paper_bgcolor = "white",
        margin        = list(l = 10, r = 10, b = 10, t = 60),
        scene         = list(
          xaxis       = list(title      = xlab,
                             tickformat = paste0(".", digits, "f"),
                             showgrid   = TRUE, gridcolor = "#E6E6E6"),
          yaxis       = list(title      = ylab,
                             tickformat = paste0(".", digits, "f"),
                             showgrid   = TRUE, gridcolor = "#E6E6E6"),
          zaxis       = list(title      = z_label,
                             showgrid   = TRUE, gridcolor = "#F0F0F0"),
          aspectratio = list(x = aspect_ratio[1L],
                             y = aspect_ratio[2L],
                             z = aspect_ratio[3L]),
          camera      = camera
        )
      )
  }

  ## --------------------------------------------------------------------------
  ## 5b.  HEATMAP  -- z = t(mat), x = grid, y = tau
  ## --------------------------------------------------------------------------
  if (type == "heatmap") {

    p <- plotly::plot_ly(
      x             = x_vals,
      y             = y_vals,
      z             = t(mat),
      type          = "heatmap",
      colorscale    = cs,
      showscale     = TRUE,
      zmin          = z_min,
      zmax          = z_max,
      zmid          = if (center_zero) 0 else NULL,
      colorbar      = list(title      = z_label,
                           tickformat = paste0(".", digits, "f"),
                           len        = 0.85),
      hovertemplate = hover_tmpl
    )

    if (sig_overlay) {
      pmat <- .qqr_get_surface(x, coef, "pvalue")
      anns <- list()
      for (i in seq_len(m)) {
        for (j in seq_along(x$taus)) {
          pv <- pmat[i, j]
          if (!is.na(pv) && pv < sig_threshold)
            anns[[length(anns) + 1L]] <- list(
              x = x_vals[i], y = y_vals[j],
              text = sig_symbol, xref = "x", yref = "y",
              showarrow = FALSE,
              font = list(size = 11L, color = "black")
            )
        }
      }
      if (length(anns) > 0L)
        p <- p |> plotly::layout(annotations = anns)
    }

    p <- p |>
      plotly::layout(
        title         = list(text = paste0("<b>", title, "</b>"),
                             x = 0, font = list(size = 14L)),
        paper_bgcolor = "white",
        plot_bgcolor  = "white",
        margin        = list(l = 70, r = 20, b = 80, t = 70),
        xaxis         = list(title      = xlab,
                             tickformat = paste0(".", digits, "f"),
                             tickangle  = 45L,
                             tickfont   = list(size = 9L)),
        yaxis         = list(title      = ylab,
                             tickformat = paste0(".", digits, "f"),
                             tickfont   = list(size = 9L))
      )
  }

  ## --------------------------------------------------------------------------
  ## 5c.  CONTOUR  -- z = t(mat), x = grid, y = tau
  ## --------------------------------------------------------------------------
  if (type == "contour") {

    contour_size <- (z_max - z_min) / max(n_contours, 2L)

    p <- plotly::plot_ly(
      x             = x_vals,
      y             = y_vals,
      z             = t(mat),
      type          = "contour",
      colorscale    = cs,
      showscale     = TRUE,
      zmin          = z_min,
      zmax          = z_max,
      zmid          = if (center_zero) 0 else NULL,
      colorbar      = list(title      = z_label,
                           tickformat = paste0(".", digits, "f"),
                           len        = 0.85),
      contours      = list(
        start      = z_min,
        end        = z_max,
        size       = contour_size,
        showlabels = show_labels,
        labelfont  = list(size = 9L, color = "black"),
        coloring   = "heatmap"
      ),
      line          = list(smoothing = 0.85),
      hovertemplate = hover_tmpl
    )

    if (z_min <= 0 && z_max >= 0) {
      p <- p |>
        plotly::add_contour(
          x          = x_vals,
          y          = y_vals,
          z          = t(mat),
          showscale  = FALSE,
          contours   = list(start = 0, end = 0, size = 1,
                            showlabels = FALSE, coloring = "lines"),
          line       = list(color = "black", width = 2.5),
          hoverinfo  = "skip",
          name       = paste0(z_label, " = 0"),
          inherit    = FALSE
        )
    }

    p <- p |>
      plotly::layout(
        title         = list(text = paste0("<b>", title, "</b>"),
                             x = 0, font = list(size = 14L)),
        paper_bgcolor = "white",
        plot_bgcolor  = "white",
        margin        = list(l = 70, r = 20, b = 80, t = 70),
        xaxis         = list(title      = xlab,
                             tickformat = paste0(".", digits, "f"),
                             tickangle  = 45L,
                             tickfont   = list(size = 9L)),
        yaxis         = list(title      = ylab,
                             tickformat = paste0(".", digits, "f"),
                             tickfont   = list(size = 9L))
      )
  }

  ## -- Snapshot --------------------------------------------------------------
  if (snapshot) {
    tryCatch(
      plotly::save_image(p, snapshot_file),
      error = function(e) message("Snapshot failed: ", e$message)
    )
  }

  p
}


## ============================================================================
## 4.  qqr_cor_heatmap()
## ============================================================================

## Internal computation helpers (not exported) --------------------------------

## Somers' D(Y|X) = (C - D) / (number of pairs not tied on X)
##                = (C - D) / (C + D + T_Y),
## where T_Y counts pairs tied on Y but not on X (Somers, 1962).
## @noRd
.qqr_somers_d <- function(y, x) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n  <- length(x)
  if (n < 2L) return(NA_real_)
  C <- 0; D <- 0; TY <- 0
  for (i in seq_len(n - 1L)) {
    dx <- x[(i+1L):n] - x[i]; dy <- y[(i+1L):n] - y[i]
    C  <- C  + sum(dx * dy > 0)
    D  <- D  + sum(dx * dy < 0)
    TY <- TY + sum(dy == 0 & dx != 0)
  }
  denom <- C + D + TY
  if (denom == 0) return(NA_real_)
  (C - D) / denom
}

.qqr_blomqvist <- function(y, x) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n  <- length(x)
  if (n < 2L) return(NA_real_)
  mx    <- stats::median(x); my <- stats::median(y)
  same  <- sum((x > mx & y > my) | (x < mx & y < my))
  diff  <- sum((x > mx & y < my) | (x < mx & y > my))
  denom <- same + diff
  if (denom == 0L) return(NA_real_)
  (same - diff) / denom
}

.qqr_coexceedance <- function(y, x, q_x_lo, q_y_lo) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n  <- length(x)
  if (n == 0L) return(NA_real_)
  sum(x >= q_x_lo & y >= q_y_lo) / n
}

.qqr_cell_stat <- function(y_cell, x_cell, method,
                           q_x_lo = NULL, q_y_lo = NULL,
                           full_y = NULL,  full_x = NULL) {
  n <- sum(is.finite(x_cell) & is.finite(y_cell))
  if (n < 3L) return(NA_real_)
  switch(method,
         spearman     = stats::cor(x_cell, y_cell, use = "complete.obs",
                                   method = "spearman"),
         kendall      = stats::cor(x_cell, y_cell, use = "complete.obs",
                                   method = "kendall"),
         somers_d     = .qqr_somers_d(y_cell, x_cell),
         pearson      = stats::cor(x_cell, y_cell, use = "complete.obs",
                                   method = "pearson"),
         blomqvist    = .qqr_blomqvist(y_cell, x_cell),
         coexceedance = {
           if (is.null(q_x_lo) || is.null(q_y_lo) || is.null(full_y) ||
               is.null(full_x))
             stop("coexceedance requires q_x_lo, q_y_lo, full_y, full_x.")
           .qqr_coexceedance(full_y, full_x, q_x_lo, q_y_lo)
         },
         stop("Unknown method: '", method, "'.")
  )
}


#' Quantile Dependence Heatmap for QQR Objects
#'
#' @description
#' Computes and displays a two-dimensional quantile dependence heatmap for the
#' focal quantile regressor (\eqn{X}) and the dependent variable (\eqn{Y})
#' of a fitted \code{qqr} model using \pkg{plotly}.  Each cell at position
#' (\eqn{X\text{-quantile bin}_{i}}, \eqn{Y\text{-quantile bin}_{j}}) is
#' coloured by the chosen dependence measure computed within that joint
#' quantile stratum.
#'
#' Six dependence measures are supported:
#' \describe{
#'   \item{\code{"spearman"}}{Spearman rank correlation within each cell.
#'     The default.  Robust, assumption-free, and the most widely understood.
#'     Range \eqn{[-1, 1]}.}
#'   \item{\code{"kendall"}}{Kendall's \eqn{\tau_b} within each cell.  More
#'     robust than Spearman when cell sizes are small or ties are present.
#'     Range \eqn{[-1, 1]}.}
#'   \item{\code{"somers_d"}}{Somers' \eqn{D(Y|X)}: an asymmetric rank
#'     association measure \eqn{(C - D) / (C + D + T_Y)}, where \eqn{C}
#'     is concordant pairs, \eqn{D} discordant, and \eqn{T_Y} pairs tied on
#'     \eqn{Y} but not on \eqn{X} (so the denominator counts all pairs not
#'     tied on \eqn{X}).  Directional: \eqn{D(Y|X) \neq D(X|Y)} when there
#'     are ties.  Range \eqn{[-1, 1]}.}
#'   \item{\code{"pearson"}}{Pearson linear correlation within each cell.
#'     Assumes approximate linearity; rank-based measures are generally
#'     preferable in small cells.  Range \eqn{[-1, 1]}.}
#'   \item{\code{"blomqvist"}}{Blomqvist's \eqn{\beta} (medial correlation).
#'     Uses only the signs of deviations from cell medians, making it
#'     extremely robust for very small cells.  Range \eqn{[-1, 1]}.}
#'   \item{\code{"coexceedance"}}{Co-exceedance probability.  For each pair
#'     of levels \eqn{(\tau_i, \tau_j)} it computes the empirical joint
#'     exceedance probability
#'     \eqn{P(X \geq \hat Q_X(\tau_i),\; Y \geq \hat Q_Y(\tau_j))} over the
#'     full sample (not within cells; \code{mode} is ignored).  Under
#'     independence it is approximately \eqn{(1-\tau_i)(1-\tau_j)}.
#'     Range \eqn{[0, 1]}.}
#' }
#'
#' @details
#' \strong{Axis convention (Option C).}
#' The x-axis shows the X (regressor) quantile bins and the y-axis shows the
#' Y (dependent variable) quantile bins.  This places the predictor on the
#' horizontal axis, matching the natural regression convention and aligning
#' with \code{\link{qqr_heatmap}} so both plots can be compared side-by-side.
#'
#' \strong{Computation modes.}
#'
#' In \code{mode = "subsample"} (recommended), the data are partitioned into
#' \code{n_bins} x \code{n_bins} non-overlapping joint quantile cells and the
#' chosen measure is computed on observations within each cell.  Cell size
#' information is stored in \code{result$cell_n}.
#'
#' In \code{mode = "binary"}, the function replicates the method of Roudane
#' (2026): it correlates the binary indicators
#' \eqn{I(X \leq Q_X(\tau_i))} and \eqn{I(Y \leq Q_Y(\tau_j))} across the
#' \emph{full} sample for each cell.  This mode is retained for
#' reproducibility with the original paper but is methodologically inferior
#' to \code{"subsample"} for most applications.
#'
#' \strong{Diagnostic use.}
#' \code{qqr_cor_heatmap} is a \emph{pre-model diagnostic} -- it uses only the
#' raw data in \code{x$data}, not the estimated QQR coefficient surfaces.
#' Its purpose is to reveal the distributional dependence structure before
#' interpreting the fitted model.  If the heatmap is approximately uniform,
#' standard quantile regression may suffice; substantial variation across
#' cells motivates the full QQR analysis.
#'
#' \strong{Independence benchmark (co-exceedance).}
#' When \code{method = "coexceedance"} and \code{independence_line = TRUE},
#' a dashed contour marks where the observed co-exceedance equals the
#' independence benchmark \eqn{(1 - \tau_i)(1 - \tau_j)}.  The returned
#' list also contains \code{independence}, the benchmark matrix, so that
#' \code{result$matrix - result$independence} gives the excess joint
#' exceedance (positive values indicate positive dependence).
#'
#' \strong{Saving output.}
#' \itemize{
#'   \item \code{htmlwidgets::saveWidget(result$plot, "dep_heatmap.html",
#'     selfcontained = TRUE)}
#'   \item \code{plotly::save_image(result$plot, "dep_heatmap.png")}
#'     (requires \code{reticulate::py_install("kaleido")})
#' }
#'
#' The binary-indicator heatmap concept is adapted from Roudane (2026);
#' the subsample mode and all six dependence measures are original
#' contributions of this package.
#'
#' @param x           A fitted object of class \code{"qqr"}, as returned by
#'                    \code{\link{qqr}}.  The function uses \code{x$data},
#'                    \code{x$independent}, \code{x$dependent}, and
#'                    \code{x$taus}.
#' @param method      Character; dependence measure to compute within each
#'                    cell.  One of \code{"spearman"} (default),
#'                    \code{"kendall"}, \code{"somers_d"}, \code{"pearson"},
#'                    \code{"blomqvist"}, or \code{"coexceedance"}.  See
#'                    Description for full definitions and appropriate use cases.
#' @param mode        Character; computation strategy.
#'   \describe{
#'     \item{\code{"subsample"}}{(Default) Partitions data into non-overlapping
#'       joint quantile cells and computes the measure on the subset of
#'       observations within each cell.}
#'     \item{\code{"binary"}}{Correlates binary quantile indicators
#'       \eqn{I(X \leq Q_X(\tau_i))} and \eqn{I(Y \leq Q_Y(\tau_j))} over
#'       the full sample.  Only valid for \code{"spearman"}, \code{"kendall"},
#'       and \code{"pearson"}; other methods will switch automatically to
#'       \code{"subsample"} with a warning.}
#'   }
#' @param n_bins      Integer; number of quantile bins per axis.
#'                    \code{NULL} (default): uses the exact tau levels from
#'                    the fitted model (\code{x$taus}).  If the fit was run
#'                    with \code{taus = seq(0.05, 0.95, by = 0.05)}, both
#'                    axes will show \code{0.05, 0.10, ..., 0.95} and each
#'                    cell corresponds to one tau level of the fit.
#'                    Supply an integer (minimum 2) to override with a
#'                    custom uniform grid independent of the fit taus.
#' @param palette     Character; plotly / RColorBrewer palette name.  For
#'                    measures with range \eqn{[-1, 1]} the default is
#'                    \code{"RdBu"} (blue = negative, red = positive).  For
#'                    \code{"coexceedance"} the default is \code{"YlOrRd"}.
#'                    Any valid RColorBrewer or plotly name is accepted.
#' @param center_zero Logical; fix the colour midpoint at zero.  Default
#'                    \code{TRUE} for measures with range \eqn{[-1, 1]};
#'                    \code{FALSE} for \code{"coexceedance"}.
#' @param clim        Numeric \code{c(min, max)}; fixed colour limits.
#'                    Useful for comparing heatmaps from different variables
#'                    or samples on the same scale.  \code{NULL} uses the
#'                    observed data range.
#' @param annotate    Logical; if \code{TRUE} (default) the numeric value of
#'                    the dependence measure is printed inside each cell.
#'                    Cells with too few observations (\code{NA}) show an
#'                    em-dash (\code{"--"}).
#' @param ann_size    Numeric; font size for cell annotations in plotly units.
#'                    Default \code{11L}.
#' @param ann_thresh  Numeric in \eqn{[0, 1]}; controls the automatic
#'                    text colour switch.  Cells whose absolute value exceeds
#'                    \code{ann_thresh * max(|value|)} use white text; all
#'                    others use black.  Default \code{0.55}.
#' @param independence_line Logical; for \code{method = "coexceedance"} only.
#'                    Draw a dashed contour at the theoretical independence
#'                    benchmark \eqn{(1 - \tau_i)(1 - \tau_j)}.  Cells above
#'                    the contour exhibit positive tail dependence; cells
#'                    below exhibit negative tail dependence.  Default
#'                    \code{TRUE}.
#' @param grid_values Logical; currently unused in \code{qqr_cor_heatmap}.
#'                    Both axes always show the tau probability levels
#'                    (0.05, 0.10, ..., 0.95) so that cells are uniformly
#'                    spaced and the heatmap fills the full grid.  Kept
#'                    for signature consistency with the other plot functions.
#' @param xlab        Character; x-axis label.  \code{NULL} auto-generates
#'                    \code{"<X variable> quantile"}.
#' @param ylab        Character; y-axis label.  \code{NULL} auto-generates
#'                    \code{"<Y variable> quantile"}.
#' @param title       Character; plot title.  \code{NULL} auto-generates.
#' @param subtitle    Character; subtitle describing the method and mode.
#'                    \code{NULL} auto-generates.
#' @param digits      Integer; decimal places for cell annotations and axis
#'                    tick labels.  Default \code{2L}.
#' @param min_cell_n  Integer; minimum observations required in a cell for
#'                    it to be computed rather than returned as \code{NA}.
#'                    Cells below this threshold are shown in grey.
#'                    Default \code{3L}.
#' @param snapshot    Logical; save a static PNG snapshot of the plot.
#'                    Default \code{FALSE}.
#' @param snapshot_file Character; file path for the PNG.  Default
#'                    \code{"qqr_cor_heatmap.png"}.
#' @param ...         Currently unused.
#'
#' @return A named list of class \code{"qqr_cor_heatmap"} with the following
#'   elements:
#' \describe{
#'   \item{\code{plot}}{The \code{plotly} object.  Print directly, embed in
#'     R Markdown or Quarto, or save with \code{htmlwidgets::saveWidget()} /
#'     \code{plotly::save_image()}.}
#'   \item{\code{matrix}}{Numeric matrix of dimension
#'     \code{n_bins} x \code{n_bins} containing the computed dependence values.
#'     Rows correspond to X quantile bins; columns to Y quantile bins.}
#'   \item{\code{x_breaks}}{Numeric vector of length \code{n_bins + 1} giving
#'     the quantile break points for the X variable.}
#'   \item{\code{y_breaks}}{Numeric vector of length \code{n_bins + 1} giving
#'     the quantile break points for the Y variable.}
#'   \item{\code{cell_n}}{Integer matrix of observations per cell (only
#'     populated in \code{subsample} mode; otherwise all \code{NA}).}
#'   \item{\code{method}}{Character; the dependence method used.}
#'   \item{\code{mode}}{Character; the computation mode used.}
#' }
#' Use \code{print(result)} for a concise console summary.
#'
#' @seealso
#'   \code{\link{qqr_heatmap}} for the QQR coefficient surface heatmap;
#'   \code{\link{qqr_contour}} for contour plots of fitted surfaces;
#'   \code{\link{qqr_plotly}} for a unified interactive hub.
#'
#' @references
#' Blomqvist, N. (1950). On a measure of dependence between two random
#'   variables. \emph{Annals of Mathematical Statistics}, 21(4), 593--600.
#'
#' Roudane, M. (2026). \emph{QuantileOnQuantile: Quantile-on-Quantile
#'   Regression Analysis}. R package version 1.0.3. CRAN.
#'   The binary-indicator heatmap concept (\code{mode = "binary"}) is
#'   adapted from this work.
#'
#' Sim, N. and Zhou, H. (2015). Oil prices, US stock return and the
#'   dependence between their quantiles. \emph{Journal of Banking &
#'   Finance}, 55, 1--8.
#'
#' Somers, R. H. (1962). A new asymmetric measure of association for ordinal
#'   variables. \emph{American Sociological Review}, 27(6), 799--811.
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   res <- qqr_cor_heatmap(fit)
#'   res
#'   round(res$matrix, 2)
#'   ce <- qqr_cor_heatmap(fit, method = "coexceedance")
#'   round(ce$matrix - ce$independence, 3)
#' }
#'
#' @export
qqr_cor_heatmap <- function(x,
                            method            = c("spearman", "kendall",
                                                  "somers_d", "pearson",
                                                  "blomqvist", "coexceedance"),
                            mode              = c("subsample", "binary"),
                            n_bins            = NULL,
                            palette           = NULL,
                            center_zero       = NULL,
                            clim              = NULL,
                            annotate          = TRUE,
                            ann_size          = 11L,
                            ann_thresh        = 0.55,
                            independence_line = TRUE,
                            grid_values       = TRUE,
                            xlab              = NULL,
                            ylab              = NULL,
                            title             = NULL,
                            subtitle          = NULL,
                            digits            = 2L,
                            min_cell_n        = 3L,
                            snapshot          = FALSE,
                            snapshot_file     = "qqr_cor_heatmap.png",
                            ...) {

  if (!requireNamespace("plotly", quietly = TRUE))
    stop("Package 'plotly' is required. Install with: install.packages('plotly')")
  if (!inherits(x, "qqr"))
    stop("'x' must be a fitted qqr object.")

  method <- match.arg(method)
  mode   <- match.arg(mode)

  binary_ok <- c("spearman", "kendall", "pearson")
  if (mode == "binary" && !(method %in% binary_ok)) {
    warning("mode = 'binary' is not meaningful for method = '", method,
            "'. Switching to mode = 'subsample'.", call. = FALSE)
    mode <- "subsample"
  }

  ## -- 1. Extract data -------------------------------------------------------
  mf     <- x$data
  y_name <- x$dependent
  x_name <- x$independent

  if (!(y_name %in% names(mf)))
    stop("Dependent variable '", y_name, "' not found in x$data.")
  if (!(x_name %in% names(mf)))
    stop("Quantile regressor '", x_name, "' not found in x$data.")

  y_raw <- mf[[y_name]]; x_raw <- mf[[x_name]]
  ok    <- is.finite(x_raw) & is.finite(y_raw)
  y_raw <- y_raw[ok]; x_raw <- x_raw[ok]; n <- length(x_raw)

  if (n < 10L)
    stop("Too few complete observations (", n, ") to construct the heatmap.")

  ## -- 2. Grid size and tau levels --------------------------------------------
  ## Two cases:
  ##   (a) n_bins NULL (default): use the exact taus from the fit.
  ##       If fit used seq(0.05, 0.95, by=0.05) the axes show 0.05..0.95.
  ##   (b) n_bins supplied: build a fresh uniform grid of that size.

  if (is.null(n_bins)) {
    ## Case (a): honour the user's fit taus exactly
    tau_mids <- x$taus                        ## e.g. 0.05, 0.10 ... 0.95
    n_bins   <- length(tau_mids)
    if (n_bins < 2L)
      stop("The fit has a single tau; supply n_bins (>= 2).", call. = FALSE)
    ## Bin boundaries: 0, midpoints between consecutive taus, 1
    probs <- c(0, (tau_mids[-1L] + tau_mids[-n_bins]) / 2, 1)
  } else {
    ## Case (b): custom grid size, uniform spacing
    n_bins   <- as.integer(n_bins)
    if (n_bins < 2L) stop("n_bins must be at least 2.")
    probs    <- seq(0, 1, length.out = n_bins + 1L)
    tau_mids <- (probs[-1L] + probs[-length(probs)]) / 2
  }

  ## -- 3. Quantile break points -----------------------------------------------
  x_breaks <- stats::quantile(x_raw, probs = probs, na.rm = TRUE)
  y_breaks <- stats::quantile(y_raw, probs = probs, na.rm = TRUE)

  ## Axis tick values: always the tau probability levels (uniform spacing).
  ## Raw data quantile values are NOT used as axis coordinates because they
  ## are non-uniform and cause plotly to leave large empty regions.
  x_tick_vals <- round(tau_mids, digits)
  y_tick_vals <- round(tau_mids, digits)

  ## -- 4. Dependence matrix --------------------------------------------------
  ## dep_mat[i, j]: X-bin i (rows), Y-bin j (cols)
  dep_mat <- matrix(NA_real_,    nrow = n_bins, ncol = n_bins)
  cell_n  <- matrix(NA_integer_, nrow = n_bins, ncol = n_bins)
  rownames(dep_mat) <- paste0(x_name, "_Q", round(tau_mids, 2))
  colnames(dep_mat) <- paste0(y_name, "_Q", round(tau_mids, 2))

  if (method == "coexceedance") {
    ## Full-sample joint exceedance of the tau-quantiles (mode ignored)
    x_q <- stats::quantile(x_raw, probs = tau_mids, names = FALSE)
    y_q <- stats::quantile(y_raw, probs = tau_mids, names = FALSE)
    for (i in seq_len(n_bins))
      for (j in seq_len(n_bins))
        dep_mat[i, j] <- .qqr_coexceedance(y_raw, x_raw, x_q[i], y_q[j])
  } else if (mode == "subsample") {
    x_bin <- pmax(1L, pmin(n_bins,
                           findInterval(x_raw, x_breaks,
                                        rightmost.closed = TRUE)))
    y_bin <- pmax(1L, pmin(n_bins,
                           findInterval(y_raw, y_breaks,
                                        rightmost.closed = TRUE)))
    for (i in seq_len(n_bins)) {
      for (j in seq_len(n_bins)) {
        in_cell       <- (x_bin == i & y_bin == j)
        nc            <- sum(in_cell)
        cell_n[i, j]  <- nc
        if (nc < min_cell_n) { dep_mat[i, j] <- NA_real_; next }
        x_cell <- x_raw[in_cell]; y_cell <- y_raw[in_cell]
        ## Constant cells (e.g. ties) have no defined correlation
        if (length(unique(x_cell)) < 2L || length(unique(y_cell)) < 2L) {
          dep_mat[i, j] <- NA_real_; next
        }
        dep_mat[i, j] <- .qqr_cell_stat(y_cell, x_cell, method)
      }
    }
  } else {
    ## Binary mode: correlate I(X <= Q_X(tau_i)) with I(Y <= Q_Y(tau_j))
    ## across the FULL sample.  Thresholds are the exact tau quantile levels,
    ## not bin boundaries.  Using bin boundaries caused the last bin (tau=0.95)
    ## to use the 100th percentile as threshold -> all 1s -> var=0 -> NA.
    x_qtiles <- stats::quantile(x_raw, probs = tau_mids, na.rm = TRUE)
    y_qtiles <- stats::quantile(y_raw, probs = tau_mids, na.rm = TRUE)
    for (i in seq_len(n_bins)) {
      for (j in seq_len(n_bins)) {
        xi <- as.numeric(x_raw <= x_qtiles[i])
        yi <- as.numeric(y_raw <= y_qtiles[j])
        if (stats::var(xi) < .Machine$double.eps ||
            stats::var(yi) < .Machine$double.eps) {
          dep_mat[i, j] <- NA_real_; next
        }
        dep_mat[i, j] <- stats::cor(xi, yi, use = "complete.obs",
                                    method = method)
      }
    }
    cell_n[] <- NA_integer_
  }

  ## -- 5. Colour parameters --------------------------------------------------
  is_sequential <- (method == "coexceedance")
  if (is.null(palette))
    palette <- if (is_sequential) "YlOrRd" else "RdBu"
  if (is.null(center_zero))
    center_zero <- !is_sequential

  lim    <- .qqr_zlim(dep_mat, clim, center_zero && !is_sequential)
  z_min  <- lim[1L]; z_max <- lim[2L]
  z_span <- z_max - z_min

  ## -- 6. Labels -------------------------------------------------------------
  method_labels <- c(
    spearman     = "Spearman rank correlation",
    kendall      = "Kendall's tau-b",
    somers_d     = "Somers' D (Y|X)",
    pearson      = "Pearson linear correlation",
    blomqvist    = "Blomqvist's beta",
    coexceedance = "Co-exceedance probability"
  )
  legend_labels <- c(
    spearman     = "\u03c1<sub>S</sub>",
    kendall      = "\u03c4<sub>b</sub>",
    somers_d     = "D(Y|X)",
    pearson      = "r",
    blomqvist    = "\u03b2<sub>M</sub>",
    coexceedance = "P(joint)"
  )
  mode_labels <- c(
    subsample = "subsample mode",
    binary    = "binary indicator mode (Roudane 2026)"
  )

  ## Option C: x-axis = X quantile tau level,  y-axis = Y quantile tau level
  if (is.null(xlab))
    xlab <- paste0(x_name, " quantile level (\u03c4)")
  if (is.null(ylab))
    ylab <- paste0(y_name, " quantile level (\u03c4)")
  if (is.null(title))
    title <- paste0("Quantile Dependence Heatmap \u2014 ",
                    x_name, " vs ", y_name)
  if (is.null(subtitle))
    subtitle <- paste0(method_labels[method], "  |  ",
                       if (method == "coexceedance") "full sample"
                       else mode_labels[mode], "  |  ",
                       n_bins, " \u00d7 ", n_bins, " grid")

  ## -- 7. Cell annotations ---------------------------------------------------
  ## dep_mat[i, j] = value at (X-bin i, Y-bin j).
  ## Option C: x = X quantile,  y = Y quantile.
  ## plotly heatmap: z[row][col] -> (x[col], y[row]).
  ## Pass z = t(dep_mat)  so  t(dep_mat)[j, i] = dep_mat[i, j]
  ## ? placed at  (x_tick_vals[i],  y_tick_vals[j])  ?

  anns <- list()
  if (annotate) {
    for (i in seq_len(n_bins)) {
      for (j in seq_len(n_bins)) {
        val      <- dep_mat[i, j]
        ann_text <- if (is.na(val)) "\u2013"
        else formatC(val, digits = digits, format = "f")

        txt_col <- if (is.na(val) || z_span < .Machine$double.eps) {
          "black"
        } else if (abs(val - (z_min + z_max) / 2) / (z_span / 2) >= ann_thresh) {
          "white"
        } else {
          "black"
        }

        anns[[length(anns) + 1L]] <- list(
          x         = x_tick_vals[i],
          y         = y_tick_vals[j],
          text      = ann_text,
          xref      = "x", yref = "y",
          showarrow = FALSE,
          font      = list(size = ann_size, color = txt_col)
        )
      }
    }
  }

  ## -- 8. Colour scale -------------------------------------------------------
  cs <- .qqr_colorscale(palette, .qqr_diverging())

  ## -- 9. Main heatmap -------------------------------------------------------
  p <- plotly::plot_ly(
    x             = x_tick_vals,
    y             = y_tick_vals,
    z             = t(dep_mat),
    type          = "heatmap",
    colorscale    = cs,
    showscale     = TRUE,
    zmin          = z_min,
    zmax          = z_max,
    zmid          = if (center_zero && !is_sequential) 0 else NULL,
    colorbar      = list(title      = legend_labels[method],
                         tickformat = paste0(".", digits, "f"),
                         len        = 0.85),
    hovertemplate = paste0(
      xlab, ": %{x:.", digits,     "f}<br>",
      ylab, ": %{y:.", digits,     "f}<br>",
      method_labels[method], ": %{z:.", digits+2L, "f}<extra></extra>"
    )
  )

  ## -- 10. Independence benchmark contour (coexceedance only) ----------------
  if (method == "coexceedance" && independence_line) {
    ## Dashed line where the observed co-exceedance equals the
    ## independence benchmark (1 - tau_i)(1 - tau_j): excess = 0.
    ind_mat <- outer(1 - tau_mids, 1 - tau_mids, FUN = "*")
    excess  <- dep_mat - ind_mat
    p <- p |>
      plotly::add_contour(
        x          = x_tick_vals,
        y          = y_tick_vals,
        z          = t(excess),
        showscale  = FALSE,
        contours   = list(
          start      = 0,
          end        = 0,
          size       = 1,
          showlabels = FALSE,
          coloring   = "lines"
        ),
        line       = list(color = "black", dash = "dash", width = 1.5),
        hoverinfo  = "skip",
        name       = "Independence benchmark",
        inherit    = FALSE
      )
  }

  ## -- 11. Layout ------------------------------------------------------------
  p <- p |>
    plotly::layout(
      title         = list(
        text = paste0("<b>", title, "</b><br><sup>", subtitle, "</sup>"),
        x    = 0,
        font = list(size = 14L)
      ),
      paper_bgcolor = "white",
      plot_bgcolor  = "white",
      margin        = list(l = 70, r = 20, b = 90, t = 90),
      xaxis         = list(
        title     = xlab,
        tickvals  = x_tick_vals,
        ticktext  = as.character(x_tick_vals),
        tickfont  = list(size = 9L),
        tickangle = 45L
      ),
      yaxis         = list(
        title    = ylab,
        tickvals = y_tick_vals,
        ticktext = as.character(y_tick_vals),
        tickfont = list(size = 9L)
      ),
      annotations   = anns
    )

  ## -- 12. Snapshot ----------------------------------------------------------
  if (snapshot) {
    tryCatch(
      plotly::save_image(p, snapshot_file),
      error = function(e) message(
        "Snapshot failed: ", e$message,
        "\nInstall kaleido: reticulate::py_install('kaleido')"
      )
    )
  }

  ## -- 13. Return ------------------------------------------------------------
  result <- list(
    plot     = p,
    matrix   = dep_mat,
    x_breaks = x_breaks,
    y_breaks = y_breaks,
    cell_n   = cell_n,
    method   = method,
    mode     = if (method == "coexceedance") "full sample" else mode,
    independence = if (method == "coexceedance")
      outer(1 - tau_mids, 1 - tau_mids, FUN = "*") else NULL
  )
  class(result) <- c("qqr_cor_heatmap", "list")
  result
}


## ============================================================================
## print method for qqr_cor_heatmap
## ============================================================================

#' Print method for qqr_cor_heatmap objects
#'
#' @description
#' Prints a concise console summary of an object of class
#' \code{"qqr_cor_heatmap"}, as returned by \code{\link{qqr_cor_heatmap}}.
#' Reports the dependence method, computation mode, grid dimensions (and
#' which axis corresponds to X and Y), the value range, mean, and the count
#' of NA cells.  Also prints quick instructions for saving the embedded
#' \pkg{plotly} figure to HTML or PNG.
#'
#' @param x   An object of class \code{"qqr_cor_heatmap"}, as returned by
#'            \code{\link{qqr_cor_heatmap}}.
#' @param ... Currently unused; reserved for future arguments.
#'
#' @return Invisibly returns \code{x}, allowing use in pipes.
#'
#' @seealso \code{\link{qqr_cor_heatmap}}
#'
#' @method print qqr_cor_heatmap
#' @export
print.qqr_cor_heatmap <- function(x, ...) {
  method_labels <- c(
    spearman     = "Spearman rank correlation",
    kendall      = "Kendall's tau-b",
    somers_d     = "Somers' D(Y|X)",
    pearson      = "Pearson linear correlation",
    blomqvist    = "Blomqvist's beta",
    coexceedance = "Co-exceedance probability"
  )
  cat("Quantile Dependence Heatmap\n")
  cat(strrep("-", 50L), "\n")
  cat("Method :", method_labels[x$method], "\n")
  cat("Mode   :", x$mode, "\n")
  cat("Grid   :", nrow(x$matrix), "x", ncol(x$matrix),
      "bins  (rows = X quantile, cols = Y quantile)\n")
  vals <- x$matrix[is.finite(x$matrix)]
  if (length(vals))
    cat(sprintf("Range  : [%.4f, %.4f]   mean = %.4f\n",
                min(vals), max(vals), mean(vals)))
  n_na <- sum(!is.finite(x$matrix))
  if (n_na > 0L)
    cat("NA cells:", n_na,
        "(too few or constant observations in the cell)\n")
  cat(strrep("-", 50L), "\n")
  cat("$plot    -> plotly figure\n")
  cat("$matrix  -> raw dependence matrix (numeric)\n")
  cat("$cell_n  -> observations per cell (subsample mode only)\n")
  cat("Save HTML : htmlwidgets::saveWidget($plot, 'dep_heatmap.html')\n")
  cat("Save PNG  : plotly::save_image($plot, 'dep_heatmap.png')\n")
  invisible(x)
}
