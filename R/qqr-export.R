#' Export QQR coefficient and p-value surfaces to Excel
#'
#' @description
#' Writes coefficient and/or p-value surfaces from a fitted
#' \code{"qqr"} object to an Excel workbook using the \pkg{openxlsx}
#' package. Each coefficient is written to a separate worksheet, and
#' coefficients and p-values can be exported together or separately.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param file Character string giving the path to the Excel file to be
#'   created. The default is \code{"QQR_coefficients.xlsx"}.
#' @param which Character string indicating what to export: \code{"both"}
#'   for both coefficients and p-values, \code{"coef"} for coefficients
#'   only, or \code{"pvalue"} for p-values only.
#' @param overwrite Logical flag indicating whether an existing file at
#'   \code{file} should be overwritten. Passed to
#'   \code{\link[openxlsx]{saveWorkbook}}.
#' @param include_grid_value Logical flag indicating whether to include a
#'   column with the actual grid values (locations of the quantile
#'   regressor) in addition to the grid index.
#' @param digits Optional integer number of decimal places to which the
#'   exported numeric entries should be rounded. If \code{NULL}, values
#'   are written without rounding. Rounding only affects the exported
#'   workbook and does not alter the QQR object.
#'
#' @details
#' For each coefficient in \code{x$coef_names}, the corresponding matrix
#' is converted to a data frame with one row per grid point and one column
#' per quantile level. Columns are named \code{tau_0.05},
#' \code{tau_0.10}, and so on. Two leading columns give the grid index and
#' optionally the grid value. A separate worksheet is created for each
#' coefficient, with suffix \code{"_coef"} or \code{"_pval"} appended to a
#' readable label.
#'
#' This function requires the \pkg{openxlsx} package to be installed.
#'
#' @return
#' The path to the created workbook is returned invisibly.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{word.qqr}},
#' \code{\link[openxlsx]{writeData}}
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
#' # Export both coefficients and p-values
#' excel.qqr(fit, file = "All_QQR_results.xlsx",
#'           which = "both", digits = 4)
#' }
#'
#' @export
excel.qqr <- function(x,
                      file              = "QQR_coefficients.xlsx",
                      which             = c("both", "coef", "pvalue"),
                      overwrite         = TRUE,
                      include_grid_value = TRUE,
                      digits            = NULL) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for excel.qqr().")
  }

  which <- match.arg(which)

  wb <- openxlsx::createWorkbook()

  make_sheet_label <- function(internal_name) {
    if (internal_name == "(Intercept)") {
      "Intercept"
    } else if (internal_name == "z") {
      x$independent
    } else {
      internal_name
    }
  }

  round_df <- function(df) {
    # Round only numeric columns, leave Grid_Point as integer
    num_cols <- vapply(df, is.numeric, logical(1L))
    if ("Grid_Point" %in% names(df)) {
      num_cols[names(df) == "Grid_Point"] <- FALSE
    }
    if (!is.null(digits) && digits >= 0) {
      df[num_cols] <- lapply(df[num_cols], function(z) round(z, digits))
    }
    df
  }

  # Coefficients
  if (which %in% c("both", "coef")) {
    for (nm in x$coef_names) {
      mat <- x$coef_surfaces[[nm]]
      df  <- as.data.frame(mat)
      colnames(df) <- paste0("tau_", x$taus)
      df$Grid_Point <- seq_along(x$grid)
      if (include_grid_value) {
        df$Grid_Value <- x$grid
        df <- df[, c("Grid_Point", "Grid_Value",
                     paste0("tau_", x$taus))]
      } else {
        df <- df[, c("Grid_Point", paste0("tau_", x$taus))]
      }

      df <- round_df(df)  # format only the export

      sheet_name <- paste0(make_sheet_label(nm), "_coef")
      sheet_name <- gsub("[^A-Za-z0-9]+", "_", sheet_name)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, df)
    }
  }

  # P-values
  if (which %in% c("both", "pvalue")) {
    for (nm in x$coef_names) {
      mat <- x$pval_surfaces[[nm]]
      df  <- as.data.frame(mat)
      colnames(df) <- paste0("tau_", x$taus)
      df$Grid_Point <- seq_along(x$grid)
      if (include_grid_value) {
        df$Grid_Value <- x$grid
        df <- df[, c("Grid_Point", "Grid_Value",
                     paste0("tau_", x$taus))]
      } else {
        df <- df[, c("Grid_Point", paste0("tau_", x$taus))]
      }

      df <- round_df(df)  # format only the export

      sheet_name <- paste0(make_sheet_label(nm), "_pval")
      sheet_name <- gsub("[^A-Za-z0-9]+", "_", sheet_name)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet_name, df)
    }
  }

  openxlsx::saveWorkbook(wb, file, overwrite = overwrite)
  invisible(file)
}



#' Export selected QQR surfaces to a Word document
#'
#' @description
#' Creates a Word document containing nicely formatted tables of selected
#' coefficient and/or p-value surfaces from a fitted \code{"qqr"} object,
#' using the \pkg{officer} and \pkg{flextable} packages.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param file Character string giving the path to the Word document to be
#'   created. The default is \code{"QQR_results.docx"}.
#' @param coef Character string identifying which coefficient surface to
#'   export. This can be the quantile regressor name, \code{"intercept"}
#'   or \code{"(Intercept)"}, or the name of a control variable.
#' @param which Character string indicating which tables to include:
#'   \code{"both"} for both coefficients and p-values, \code{"coef"} for
#'   coefficients only, or \code{"pvalue"} for p-values only.
#' @param grid_subset Optional integer vector of grid indices to include
#'   in the tables. If \code{NULL}, all grid points are used. Values
#'   outside the valid range are truncated to the nearest valid index.
#'   If any \code{NA} values are supplied, they are replaced by the last
#'   grid index.
#' @param digits Optional integer number of decimal places to which
#'   numeric entries in the tables should be rounded. If \code{NULL},
#'   values are exported without rounding. Rounding does not modify the
#'   underlying QQR object.
#' @param include_grid_value Logical flag indicating whether the exported
#'   tables should include the \code{Grid_Value} column. When
#'   \code{TRUE} (the default), tables contain both the integer grid
#'   index (\code{Grid_Point}) and the corresponding numeric grid value.
#'   When \code{FALSE}, only \code{Grid_Point} and the \code{tau_*}
#'   columns are exported, which may be useful for compact reporting or
#'   when the grid values are documented elsewhere.
#'
#' @details
#' The function selects the requested coefficient and p-value matrices,
#' restricts them to the desired subset of grid points, and converts them
#' into data frames with one row per grid location and one column per
#' quantile level. Leading columns always include the grid index
#' (\code{Grid_Point}), and optionally include the grid value
#' (\code{Grid_Value}) when \code{include_grid_value = TRUE}.
#'
#' These data frames are then converted to \pkg{flextable} objects and
#' inserted into a Word document via \pkg{officer}, with simple headings
#' indicating the coefficient, dependent variable, and type of surface.
#'
#' This function requires the \pkg{officer} and \pkg{flextable} packages
#' to be installed.
#'
#' @return
#' The path to the created Word document is returned invisibly.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{excel.qqr}},
#' \code{\link[officer]{read_docx}}, \code{\link[flextable]{flextable}}
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
#' # Export coefficients and p-values for the quantile regressor,
#' # including grid values
#' word.qqr(fit,
#'          file              = "QQR_x_coef_pval.docx",
#'          coef              = "x",
#'          which             = "both",
#'          digits            = 4,
#'          include_grid_value = TRUE)
#'
#' # Export a compact version without the Grid_Value column
#' word.qqr(fit,
#'          file              = "QQR_x_coef_pval_compact.docx",
#'          coef              = "x",
#'          which             = "both",
#'          digits            = 4,
#'          include_grid_value = FALSE)
#' }
#'
#' @export
word.qqr <- function(x,
                     file               = "QQR_results.docx",
                     coef               = x$independent,
                     which              = c("both", "coef", "pvalue", "coef_stars"),
                     grid_subset        = NULL,
                     digits             = NULL,
                     include_grid_value = TRUE) {
  if (!inherits(x, "qqr"))
    stop("Object must be of class 'qqr'.")
  if (!requireNamespace("officer", quietly = TRUE) ||
      !requireNamespace("flextable", quietly = TRUE)) {
    stop("Packages 'officer' and 'flextable' are required for word.qqr().")
  }
  which <- match.arg(which)
  surf_name <- .qqr_match_coef(x, coef)
  coef_mat <- x$coef_surfaces[[surf_name]]
  pval_mat <- x$pval_surfaces[[surf_name]]
  n_grid <- nrow(coef_mat)
  # Grid selection
  if (is.null(grid_subset)) {
    idx <- seq_len(n_grid)
  } else {
    if (any(is.na(grid_subset))) {
      grid_subset[is.na(grid_subset)] <- n_grid
    }
    idx <- unique(pmax(1L,
                       pmin(n_grid,
                            as.integer(grid_subset))))
  }
  # -----------------------------
  # Significance stars helper
  # -----------------------------
  make_stars <- function(pval) {
    sig <- ifelse(
      is.na(pval), "\u2020",
      ifelse(pval < 0.01, "***",
             ifelse(pval < 0.05,  "**",
                    ifelse(pval < 0.10,  "*",
                           "")))
    )
    sig
  }
  # -----------------------------
  # Standard numeric table
  # -----------------------------
  make_df <- function(mat) {
    df <- as.data.frame(mat[idx, , drop = FALSE])
    colnames(df) <- paste0("tau_", x$taus)
    df$Grid_Point <- idx
    if (include_grid_value) {
      df$Grid_Value <- x$grid[idx]
      front <- c("Grid_Point","Grid_Value")
    } else {
      front <- "Grid_Point"
    }
    df <- df[, c(front,
                 paste0("tau_", x$taus))]
    if (!is.null(digits) && digits >= 0) {
      num_cols <- vapply(df,
                         is.numeric,
                         logical(1L))
      num_cols[names(df) == "Grid_Point"] <- FALSE
      df[num_cols] <- lapply(df[num_cols],
                             function(z)
                               round(z,digits))
    }
    df
  }
  # -----------------------------
  # Coef + Stars table
  # -----------------------------
  make_star_df <- function(coef_mat,
                         pval_mat) {
    coef_sub <- coef_mat[idx, , drop = FALSE]
    pval_sub <- pval_mat[idx, , drop = FALSE]
    out <- matrix("",
                  nrow=nrow(coef_sub),
                  ncol=ncol(coef_sub))
    for(j in seq_len(ncol(coef_sub))){
      coef_vec <- coef_sub[,j]
      pval_vec <- pval_sub[,j]
      if(!is.null(digits)){
        coef_vec <- round(coef_vec,digits)
      }
      stars <- make_stars(pval_vec)
      out[,j] <- paste0(coef_vec,stars)
    }
    df <- as.data.frame(out)
    colnames(df) <- paste0("tau_",x$taus)
    df$Grid_Point <- idx
    if(include_grid_value){
      df$Grid_Value <- x$grid[idx]
      front <- c("Grid_Point","Grid_Value")
    }else{
      front <- "Grid_Point"
    }
    df <- df[,c(front,
                paste0("tau_",x$taus))]
    df
  }
  doc <- officer::read_docx()
  # -----------------------------
  # COEF ONLY
  # -----------------------------
  if(which %in% c("both","coef")){
    df_coef <- make_df(coef_mat)
    title_coef <- paste0(
      "QQR coefficients for ",
      coef,
      " (dependent variable: ",
      x$dependent,
      ")"
    )
    doc <- officer::body_add_par(doc,
                                 title_coef,
                                 style="heading 2")

    ft <- flextable::flextable(df_coef)

    ft <- flextable::autofit(ft)

    doc <- flextable::body_add_flextable(doc,ft)
  }
  # -----------------------------
  # PVALUES ONLY
  # -----------------------------
  if(which %in% c("both","pvalue")){
    df_p <- make_df(pval_mat)
    title_p <- paste0(
      "QQR p-values for ",
      coef,
      " (dependent variable: ",
      x$dependent,
      ")"
    )
    doc <- officer::body_add_par(doc,
                                 title_p,
                                 style="heading 2")
    ft <- flextable::flextable(df_p)
    ft <- flextable::autofit(ft)
    doc <- flextable::body_add_flextable(doc,ft)
  }
  # -----------------------------
  # COEF WITH STARS
  # -----------------------------
  if(which == "coef_stars"){
    df_star <- make_star_df(coef_mat,
                            pval_mat)
    title_star <- paste0(
      "QQR coefficients with significance stars for ",
      coef,
      " (dependent variable: ",
      x$dependent,
      ")"
    )
    doc <- officer::body_add_par(doc,
                                 title_star,
                                 style="heading 2")
    ft <- flextable::flextable(df_star)
    ft <- flextable::autofit(ft)
    doc <- flextable::body_add_flextable(doc,ft)
    # Significance note
    doc <- officer::body_add_par(
      doc,
      "Notes: Entries report quantile-on-quantile regression coefficients. Significance stars denote p-values obtained from local quantile regressions: *** p<0.01, ** p<0.05, * p<0.10. \u2020 indicates that p-values could not be computed at the corresponding grid–quantile combination.",
      style="Normal"
    )
  }
  print(doc,
        target=file)
  invisible(file)
}

