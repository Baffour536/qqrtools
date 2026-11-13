
#' @title The generic compare function
#' @param x An object
#' @param ... Further arguments passed to specific methods. Currently unused in the generic.
#' @export
compare <- function(x, ...) {
  UseMethod("compare")
}




#' @title Compare QQR with standard quantile regression
#'
#' @description
#' Compares standard quantile regression (QR) coefficients with QQR
#' estimates averaged over the quantiles of the quantile regressor, the
#' robustness check proposed by Sim and Zhou (2015). For each
#' \eqn{\tau} in the fit, the QQR surface of the chosen coefficient is
#' averaged over the grid (\eqn{\theta}) dimension,
#' \deqn{\bar\beta(\tau) = \frac{1}{m} \sum_{i=1}^{m}
#'       \hat\beta(\theta_i, \tau),}
#' and set against the QR coefficient \eqn{\hat\beta^{QR}(\tau)} from
#' \code{\link[quantreg]{rq}} on the same data and design.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Coefficient(s) to compare: the quantile regressor name
#'   (default), a control name, \code{"intercept"}, a character vector of
#'   names, or \code{"all"}.
#' @param na.rm Logical; drop \code{NA} cells when averaging over the grid.
#'   Default \code{TRUE}.
#' @param digits Optional integer; round the returned numbers.
#' @param qr_args Optional named list of extra arguments for
#'   \code{\link[quantreg]{rq}} (e.g. \code{list(method = "fn")}).
#' @param ... Currently unused.
#'
#' @details
#' The comparison is meaningful when the grid represents the distribution
#' of the quantile regressor, i.e. with \code{grid = "quantile"} (the
#' default): the equally weighted average over \eqn{\theta_i} then
#' approximates the average effect over the distribution of \eqn{X}. With
#' \code{grid = "range"} the average weights all parts of the range of
#' \eqn{X} equally, including sparse tails, and a message is shown.
#'
#' For the slope of the quantile regressor and for control coefficients,
#' \eqn{\bar\beta(\tau)} and \eqn{\hat\beta^{QR}(\tau)} estimate comparable
#' quantities if the QR model is well specified. For the intercept they are
#' not directly comparable: the local intercept is the fitted
#' \eqn{\tau}-quantile of \eqn{Y} at \eqn{X = x^\theta}, while the QR
#' intercept is evaluated at \eqn{X = 0}.
#'
#' Earlier versions (\pkg{qqr} 0.2.0) instead averaged each surface over
#' \eqn{\tau} and paired grid point \eqn{j} with \eqn{\tau_j}, which does
#' not correspond to the Sim and Zhou comparison; results from those
#' versions should be recomputed.
#'
#' @return For a single coefficient, a data frame with columns \code{tau},
#'   \code{qr} (QR coefficient), \code{qqr_mean} (QQR average over the
#'   grid) and \code{difference} (\code{qqr_mean - qr}). For several
#'   coefficients, a named list of such data frames.
#'
#' @seealso \code{\link{plot.compare.qqr}}, \code{\link{rowave.qqr}},
#'   \code{\link[quantreg]{rq}}
#'
#' @references
#' Sim, N. and Zhou, H. (2015) Oil prices, US stock return, and the
#' dependence between their quantiles. \emph{Journal of Banking & Finance},
#' 55, 1--8. \doi{10.1016/j.jbankfin.2015.01.013}
#'
#' @examples
#' set.seed(1)
#' n   <- 200
#' x   <- rnorm(n)
#' z   <- rnorm(n)
#' y   <- 1 + x + 0.5 * z + rnorm(n)
#' dat <- data.frame(y = y, x = x, z = z)
#' fit <- qqr(y ~ x + z, data = dat, qvar = "x", m = 5,
#'            taus = c(0.25, 0.5, 0.75), se_type = "none")
#' compare.qqr(fit, digits = 3)
#' compare.qqr(fit, coef = "z", digits = 3)
#'
#' @export
#' @method compare qqr
#' @rawNamespace export(compare.qqr)
compare.qqr <- function(x,
                        coef    = x$independent,
                        na.rm   = TRUE,
                        digits  = NULL,
                        qr_args = list(),
                        ...) {

  if (!inherits(x, "qqr"))
    stop("Object must be of class 'qqr'.", call. = FALSE)
  if (!is.list(qr_args))
    stop("'qr_args' must be a list.", call. = FALSE)

  ## ---------------------------------------------------------------------
  ## Several coefficients -> list
  ## ---------------------------------------------------------------------
  if (length(coef) == 1L && identical(tolower(coef), "all"))
    coef <- .qqr_user_coefs(x)
  if (length(coef) > 1L) {
    if (!is.character(coef))
      stop("'coef' must be a character vector.", call. = FALSE)
    nms <- vapply(coef, function(cn) .qqr_match_coef(x, cn), character(1L),
                  USE.NAMES = FALSE)
    out <- lapply(nms, function(nm)
      compare.qqr(x, coef = nm, na.rm = na.rm, digits = digits,
                  qr_args = qr_args))
    return(stats::setNames(out, nms))
  }

  if (identical(coef, "z") && !("z" %in% names(x$coef_surfaces)))
    stop("Coefficient 'z' is not a model coefficient. Use the original ",
         "variable name: '", x$independent, "'.", call. = FALSE)
  nm <- .qqr_match_coef(x, coef)

  if (identical(x$grid_type, "range") || is.null(x$grid_type))
    message("compare.qqr(): the grid is equally spaced over the range of ",
            x$independent, ", so the average over the grid is not an ",
            "average over its distribution. Consider grid = \"quantile\".")

  ## ---------------------------------------------------------------------
  ## Standard QR on the same design
  ## ---------------------------------------------------------------------
  des <- .qqr_design(x)
  vn  <- paste0("v", seq_len(ncol(des$X)))
  df  <- as.data.frame(des$X)
  names(df) <- vn
  df$.y <- des$y
  form <- stats::as.formula(paste(".y ~ 0 +", paste(vn, collapse = " + ")))

  qr_fit <- tryCatch(
    suppressWarnings(do.call(quantreg::rq,
                             c(list(formula = form, data = df,
                                    tau = x$taus), qr_args))),
    error = function(e)
      stop("Standard quantile regression failed: ", conditionMessage(e),
           call. = FALSE)
  )
  cm <- matrix(stats::coef(qr_fit), nrow = length(vn),
               dimnames = list(colnames(des$X), NULL))
  stored <- if (nm == "z") x$independent else nm
  if (!(stored %in% rownames(cm)))
    stop("Coefficient '", coef, "' not found in the QR fit.", call. = FALSE)
  qr_vec <- as.numeric(cm[stored, ])

  ## ---------------------------------------------------------------------
  ## QQR average over the grid (theta) for each tau
  ## ---------------------------------------------------------------------
  surf    <- x$coef_surfaces[[nm]]
  qqr_vec <- as.numeric(colMeans(surf, na.rm = na.rm))
  qqr_vec[is.nan(qqr_vec)] <- NA_real_

  out <- data.frame(
    tau        = x$taus,
    qr         = qr_vec,
    qqr_mean   = qqr_vec,
    difference = qqr_vec - qr_vec
  )

  if (!is.null(digits))
    out[] <- lapply(out, function(v) if (is.numeric(v)) round(v, digits) else v)

  out
}

## Response and uncentred design matrix of a fit (also for qqr 0.2.0
## objects, which did not store them).
## @noRd
.qqr_design <- function(x) {
  if (!is.null(x$design) && !is.null(x$y))
    return(list(X = x$design, y = x$y))
  mf <- x$data
  if (is.null(mf))
    stop("The qqr object does not contain its data; re-fit it with the ",
         "current version of qqrtools.", call. = FALSE)
  tt <- attr(mf, "terms")
  X  <- stats::model.matrix(tt, mf)
  list(X = X, y = as.numeric(stats::model.response(mf)))
}



#' Export QR–QQR comparisons to Excel
#'
#' @description
#' Writes the output of \code{\link{compare.qqr}} to an Excel workbook
#' using \pkg{openxlsx}. When multiple coefficients are requested, each
#' coefficient is placed on a separate worksheet.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param file Character string giving the path to the Excel file to be
#'   created. Default is \code{"QQR_QR_comparison.xlsx"}.
#' @param coef As in \code{\link{compare.qqr}}: a character vector of
#'   coefficient names or \code{"all"}. The default is the quantile
#'   regressor, \code{x$independent}.
#' @param na.rm Logical; passed to \code{\link{compare.qqr}}.
#' @param digits Optional integer; passed to \code{\link{compare.qqr}} to
#'   control rounding for presentation.
#' @param qr_args Optional named list of additional arguments passed to
#'   \code{\link[quantreg]{rq}} (for example, \code{list(method = "fn")}).
#' @param ... Currently unused. Included for potential extensibility.
#'
#' @details
#' For each requested coefficient, the function creates a separate
#' worksheet in the Excel workbook and writes the comparison data frame
#' returned by \code{\link{compare.qqr}}. The sheet is named
#' \code{coef_<coefficient_name>}, truncated to 31 characters to
#' comply with Excel's sheet name length limit.
#'
#' @return
#' The path to the created Excel workbook, returned invisibly.
#'
#' @seealso
#' \code{\link{compare.qqr}}, \code{\link{word.compare.qqr}},
#' \code{\link{excel.qqr}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("openxlsx", quietly = TRUE)) {
#'   excel.compare.qqr(fit, file = tempfile(fileext = ".xlsx"), coef = "all")
#' }
#'
#' @export
excel.compare.qqr <- function(x,
                              file   = "QQR_QR_comparison.xlsx",
                              coef   = x$independent,
                              na.rm  = TRUE,
                              digits = NULL,
                              qr_args = list(),
                              ...) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for excel.compare.qqr().")
  }

  cmp <- compare.qqr(
    x      = x,
    coef   = coef,
    na.rm  = na.rm,
    digits = digits,
    qr_args = qr_args
  )

  wb <- openxlsx::createWorkbook()

  if (is.data.frame(cmp)) {
    openxlsx::addWorksheet(wb, "Comparison")
    openxlsx::writeData(wb, "Comparison", cmp)
  } else {
    for (nm in names(cmp)) {
      # Excel sheet name max length is 31
      sheet <- substr(paste0("coef_", nm), 1L, 31L)
      openxlsx::addWorksheet(wb, sheet)
      openxlsx::writeData(wb, sheet, cmp[[nm]])
    }
  }

  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(file)
}



#' Export QR–QQR comparisons to a Word document
#'
#' @description
#' Creates a Word document containing tables that compare standard
#' quantile regression (QR) coefficients with averaged QQR coefficients
#' for selected regressors, using \pkg{officer} and \pkg{flextable}.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param file Character string giving the path to the Word document to
#'   be created. Default is \code{"QQR_QR_comparison.docx"}.
#' @param coef As in \code{\link{compare.qqr}}: a character vector of
#'   coefficient names or \code{"all"}. The default is the quantile
#'   regressor, \code{x$independent}.
#' @param na.rm Logical; passed to \code{\link{compare.qqr}}.
#' @param digits Optional integer; passed to \code{\link{compare.qqr}}.
#' @param qr_args Optional named list of additional arguments passed to
#'   \code{\link[quantreg]{rq}} (for example, \code{list(method = "fn")}).
#' @param ... Currently unused. Included for potential extensibility.
#'
#'
#' @details
#' For each requested coefficient, the function inserts a short heading
#' and a table with the columns returned by \code{\link{compare.qqr}}
#' (\code{tau}, \code{qr}, \code{qqr_mean}, \code{difference}) into the
#' Word document. This layout is convenient
#' for reproducing the typical empirical comparison of QR and average
#' QQR coefficients across quantiles.
#'
#' @return
#' The path to the created Word document, returned invisibly.
#'
#' @seealso
#' \code{\link{compare.qqr}}, \code{\link{word.qqr}},
#' \code{\link[officer]{read_docx}},
#' \code{\link[flextable]{flextable}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("officer", quietly = TRUE) &&
#'     requireNamespace("flextable", quietly = TRUE)) {
#'   word.compare.qqr(fit, file = tempfile(fileext = ".docx"), digits = 3)
#' }
#'
#' @export
word.compare.qqr <- function(x,
                             file   = "QQR_QR_comparison.docx",
                             coef   = x$independent,
                             na.rm  = TRUE,
                             digits = NULL,
                             qr_args = list(),
                             ...) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")
  if (!requireNamespace("officer", quietly = TRUE) ||
      !requireNamespace("flextable", quietly = TRUE)) {
    stop("Packages 'officer' and 'flextable' are required for word.compare.qqr().")
  }

  cmp <- compare.qqr(
    x       = x,
    coef    = coef,
    na.rm   = na.rm,
    digits  = digits,
    qr_args = qr_args
  )

  doc <- officer::read_docx()

  add_one <- function(df, coef_label) {
    heading <- paste0(
      "Comparison of QR and averaged QQR coefficients for ", coef_label,
      " (dependent variable: ", x$dependent, ")"
    )
    doc <<- officer::body_add_par(doc, heading, style = "heading 2")

    ft <- flextable::flextable(df)
    ft <- flextable::autofit(ft)
    doc <<- flextable::body_add_flextable(doc, ft)
  }

  if (is.data.frame(cmp)) {
    add_one(cmp, paste(as.character(coef), collapse = ", "))
  } else {
    for (nm in names(cmp)) {
      add_one(cmp[[nm]], nm)
    }
  }

  print(doc, target = file)
  invisible(file)
}




#' Plot QR and averaged QQR coefficients across quantiles
#'
#' @description
#' Produces a \pkg{ggplot2} line plot comparing standard quantile
#' regression (QR) coefficients with averaged QQR coefficients for a
#' selected regressor, with quantile levels on the x-axis.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Character string giving the coefficient to plot. The
#'   default is the quantile regressor, \code{x$independent}. Only a
#'   single coefficient is allowed.
#' @param title Optional plot title. If \code{NULL}, a default title is
#'   constructed using the coefficient name and the dependent variable.
#' @param xlab Character string for the x-axis label. Default is
#'   \code{"Quantile level (tau)"}.
#' @param ylab Character string for the y-axis label. Default is
#'   \code{"Coefficient"}.
#' @param legend_title Character string for the legend title. Default is
#'   \code{""}.
#' @param na.rm Logical; passed to \code{\link{compare.qqr}}.
#' @param digits Optional integer; passed to \code{\link{compare.qqr}}
#'   for rounding the tabulated coefficients before plotting.
#' @param qr_args Optional named list of additional arguments passed to
#'   \code{\link[quantreg]{rq}} (for example, \code{list(method = "fn")}).
#' @param ... Currently unused. Included for potential extensibility.
#'
#'
#' @return
#' A \pkg{ggplot2} object, which can be further customised or printed.
#'
#' @seealso
#' \code{\link{compare.qqr}}, \code{\link{plot.qqr}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot.compare.qqr(fit)
#' }
#'
#' @export
#' @method plot compare.qqr
#' @rawNamespace export(plot.compare.qqr)
plot.compare.qqr <- function(x,
                             coef         = x$independent,
                             title        = NULL,
                             xlab         = "Quantile level (tau)",
                             ylab         = "Coefficient",
                             legend_title = "",
                             na.rm        = TRUE,
                             digits       = NULL,
                             qr_args      = list(),
                             ...) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")
  if (length(coef) != 1L) {
    stop("plot.compare.qqr() currently supports a single coefficient at a time.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot.compare.qqr().")
  }

  cmp <- compare.qqr(
    x       = x,
    coef    = coef,
    na.rm   = na.rm,
    digits  = digits,
    qr_args = qr_args
  )

  if (!is.data.frame(cmp) ||
      !all(c("tau", "qr", "qqr_mean") %in% names(cmp))) {
    stop("Unexpected structure returned by compare.qqr(); expected a single data frame.")
  }

  df_long <- data.frame(
    tau    = rep(cmp$tau, times = 2L),
    value  = c(cmp$qr, cmp$qqr_mean),
    method = rep(c("QR", "QQR mean"), each = nrow(cmp))
  )

  coef_label <- as.character(coef)
  default_title <- paste0(
    "QR vs averaged QQR coefficients for ", coef_label,
    " (", x$dependent, ")"
  )
  if (is.null(title)) title <- default_title

  ggplot2::ggplot(df_long, ggplot2::aes(x = tau, y = value, colour = method)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title  = title,
      x      = xlab,
      y      = ylab,
      colour = legend_title
    ) +
    ggplot2::theme_minimal()
}
