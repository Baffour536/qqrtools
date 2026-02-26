
#' @title The generic compare function
#' @param x An object
#' @param ... Further arguments passed to specific methods. Currently unused in the generic.
#' @export
compare <- function(x, ...) {
  UseMethod("compare")
}




#' @title Compare QQ Regression with Standard QR Coefficients
#' @description Compare QR coefficients with row means of QQR surfaces.
#' Fits a standard quantile regression (QR) using the same formula, data
#' and quantile levels used in the QQR fit, and compares the QR
#' coefficients with the **row means** of the corresponding QQR
#' coefficient surface.
#'
#' For each coefficient of interest, the function quietly runs
#' \code{\link[quantreg]{rq}} and then constructs a data frame with one
#' row per comparison point. When a single coefficient is requested, the
#' result has three columns: the quantile level \code{tau}, the QR
#' coefficient, and the QQR row-mean coefficient. The QQR row means are
#' computed with \code{\link{rowave.qqr}}.
#'
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Character string giving the coefficient to compare.
#'   This can be \code{"intercept"} or \code{"(Intercept)"}, the name of
#'   the quantile regressor, or the name of a control variable.
#'   At present \code{coef = "all"} is not supported for comparison; call
#'   the function separately for each coefficient of interest.
#' @param na.rm Logical flag indicating whether \code{NA} values should be
#'   removed when computing row means
#' @param digits Optional integer number of decimal places to which the
#'   numeric results should be rounded. If \code{NULL}, no rounding is
#'   applied.
#' @param qr_args Optional named list of additional arguments passed to
#'   \code{\link[quantreg]{rq}} (for example, \code{list(method = "fn")}).
#' @param ... Currently unused. Included for potential extensibility.
#'
#' @details
#' Suppose the QQR coefficient surface for a given regressor is stored as
#' a matrix of dimension \eqn{m \times K}, where \eqn{m} is the number of
#' grid points (quantiles of the quantile regressor) and \eqn{K} is the
#' number of conditional quantiles \code{taus}. The function first
#' computes the row means
#' \deqn{\bar\beta_j = \frac{1}{K} \sum_{k=1}^K \beta_{j,k}, \quad
#'       j = 1, \dots, m,}
#' using \code{\link{rowave.qqr}}.
#'
#' It then fits a standard QR,
#' \code{rq(formula, tau = x$taus, data = original_data)}, and extracts
#' the QR coefficients for the same regressor at each \code{tau}.
#'
#' To construct a one-to-one comparison, the function aligns the first
#' \code{n_pair = min(length(x$taus), m)} elements of the QR coefficient
#' vector with the first \code{n_pair} row means. In typical empirical
#' applications, the grid size \code{m} is chosen equal to the number of
#' quantile levels \code{length(x$taus)}, so this alignment is exact.
#'
#' @return
#' A data frame with \code{n_pair} rows and three columns:
#' \itemize{
#'   \item \code{tau} – the conditional quantile level.
#'   \item \code{qr} – the QR coefficient at that \code{tau}.
#'   \item \code{qqr_mean} – the corresponding QQR row-mean coefficient,
#'         computed across all conditional quantiles and aligned by index.
#' }
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{rowave.qqr}},
#' \code{\link[quantreg]{rq}}
#'
#' @examples
#' \dontrun{
#'   library(quantreg)
#'   set.seed(1)
#'   n  <- 200
#'   x  <- rnorm(n)
#'   z  <- rnorm(n)
#'   y  <- x * qnorm(runif(n)) + 0.5 * z
#'   dat <- data.frame(y = y, x = x, z = z)
#'
#'   fit <- qqr(y ~ x + z, data = dat, qvar = "x", m = 15,
#'              taus = seq(0.05, 0.95, by = 0.05))
#'   cmp_z <- compare.qqr(fit, coef = "z", digits = 4)
#' }
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

  ## ---------------------------------------------------------------------
  ## Basic validation
  ## ---------------------------------------------------------------------
  if (!inherits(x, "qqr")) {
    stop("Object must be of class 'qqr'.")
  }
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for compare.qqr().")
  }

  ## ---------------------------------------------------------------------
  ## Helper: normalise coefficient input (user-facing only)
  ## ---------------------------------------------------------------------
  normalise_coef <- function(coef_in) {
    if (!is.character(coef_in) || length(coef_in) != 1L) return(coef_in)
    if (identical(coef_in, "intercept")) return("(Intercept)")
    coef_in
  }

  ## User-facing coefficient universe
  user_coefs <- c(
    if (isTRUE(x$intercept)) "(Intercept)" else character(0L),
    x$independent,
    x$controls
  )

  ## ---------------------------------------------------------------------
  ## Handle coef = "all" or vector input
  ## ---------------------------------------------------------------------
  if (length(coef) > 1L || identical(coef, "all")) {

    if (identical(coef, "all")) {
      coef_vec <- user_coefs
    } else {
      if (!is.character(coef)) {
        stop("When 'coef' has length > 1, it must be a character vector.")
      }
      coef_vec <- vapply(coef, normalise_coef, character(1L))
    }

    out <- setNames(vector("list", length(coef_vec)), coef_vec)
    for (nm in coef_vec) {
      out[[nm]] <- compare.qqr(
        x       = x,
        coef    = nm,
        na.rm   = na.rm,
        digits  = digits,
        qr_args = qr_args
      )
    }
    return(out)
  }

  ## ---------------------------------------------------------------------
  ## Single coefficient path
  ## ---------------------------------------------------------------------
  if (!is.character(coef) || length(coef) != 1L) {
    stop("Argument 'coef' must be a single character string.")
  }

  coef <- normalise_coef(coef)

  ## Explicit rejection of internal aliases
  if (identical(coef, "z")) {
    stop(
      "Coefficient 'z' is an internal QQR alias and cannot be used here. ",
      "Use the original variable name: '", x$independent, "'."
    )
  }

  if (!coef %in% user_coefs) {
    stop(
      "Coefficient '", coef, "' was not found in the model. ",
      "Available coefficients are: ",
      paste(user_coefs, collapse = ", "),
      "."
    )
  }

  ## ---------------------------------------------------------------------
  ## Reconstruct formula and data robustly
  ## ---------------------------------------------------------------------
  cl <- x$call
  if (is.null(cl$formula)) {
    stop("The original qqr call does not contain a 'formula'.")
  }

  ## Formula
  call_env <- environment(cl$formula)
  if (is.null(call_env)) call_env <- parent.frame()
  qr_formula <- eval(cl$formula, envir = call_env)

  ## Data: prefer stored data, fallback to call$data
  if (!is.null(x$data)) {
    qr_data <- x$data
  } else if (!is.null(cl$data)) {
    qr_data <- eval(cl$data, envir = call_env)
  } else {
    stop(
      "Cannot reconstruct the data used in the qqr fit. ",
      "The qqr object should store the original data."
    )
  }

  if (!is.data.frame(qr_data)) {
    stop("Reconstructed data is not a data.frame.")
  }

  taus <- x$taus

  ## ---------------------------------------------------------------------
  ## Fit standard QR
  ## ---------------------------------------------------------------------
  rq_args <- c(
    list(
      formula = qr_formula,
      data    = qr_data,
      tau     = taus
    ),
    qr_args
  )

  qr_fit   <- do.call(quantreg::rq, rq_args)
  coef_mat <- coef(qr_fit)

  if (is.null(rownames(coef_mat))) {
    stop("QR coefficient matrix has no row names.")
  }

  if (!coef %in% rownames(coef_mat)) {
    stop(
      "Coefficient '", coef, "' was not found in the QR fit. ",
      "Available QR coefficients are: ",
      paste(rownames(coef_mat), collapse = ", "),
      "."
    )
  }

  qr_vec <- as.numeric(coef_mat[coef, ])

  ## ---------------------------------------------------------------------
  ## QQR row means
  ## ---------------------------------------------------------------------
  rm_df <- rowave.qqr(
    x,
    coef               = coef,
    na.rm              = na.rm,
    digits             = NULL,
    include_grid_value = FALSE
  )

  mean_col <- setdiff(names(rm_df), "grid_index")
  if (length(mean_col) != 1L) {
    stop("rowave.qqr() returned an unexpected structure.")
  }

  qqr_vec <- rm_df[[mean_col]]

  ## ---------------------------------------------------------------------
  ## Align QR and QQR vectors
  ## ---------------------------------------------------------------------
  n_pair <- min(length(taus), length(qqr_vec))

  out <- data.frame(
    tau      = taus[seq_len(n_pair)],
    qr       = qr_vec[seq_len(n_pair)],
    qqr_mean = qqr_vec[seq_len(n_pair)]
  )

  if (!is.null(digits)) {
    out[] <- lapply(out, function(v) {
      if (is.numeric(v)) round(v, digits) else v
    })
  }

  out
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
#' and a table with three columns (\code{tau}, \code{qr},
#' \code{qqr_mean}) into the Word document. This layout is convenient
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
    add_one(cmp, as.character(coef))
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
