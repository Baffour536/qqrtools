#' Extract QQR coefficient or p-value surfaces as tables
#'
#' @description
#' Extracts coefficient or p-value surfaces from a fitted \code{"qqr"}
#' object and returns them as data frames in either wide (matrix-like) or
#' long (tidy) format, optionally restricting the grid points, quantile
#' levels, and coefficients.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Either \code{"all"} to extract surfaces for all
#'   coefficients (intercept, quantile regressor, and any controls), or a
#'   character vector naming one or more coefficients. Coefficients may be
#'   referred to by \code{"intercept"} or \code{"(Intercept)"} for the
#'   intercept, by the name of the quantile regressor, or by control
#'   variable names.
#' @param what Character string indicating which surfaces to extract:
#'   \code{"coef"} for coefficient values or \code{"pvalue"} for
#'   p-values.
#' @param grid_subset Optional integer vector of grid indices to keep.
#'   If \code{NULL} (the default), all grid points are used. Values
#'   outside the valid range are truncated to the nearest valid index.
#' @param tau_subset Optional subset of quantile levels: either integer
#'   indices into \code{x$taus} (values \eqn{\ge 1}) or quantile values in
#'   \eqn{(0, 1)} that are matched to \code{x$taus} (with a tolerance of
#'   \code{1e-8}; unmatched values give an error). If \code{NULL}, all
#'   quantile levels are used.
#' @param include_grid_value Logical flag indicating whether to include a
#'   column containing the actual grid values (the locations of the
#'   quantile regressor) in addition to the grid index.
#' @param digits Optional integer number of decimal places to which
#'   numeric entries should be rounded for display. If \code{NULL},
#'   values are left on the original scale.
#' @param long Logical flag controlling the output shape. If \code{FALSE}
#'   and a single coefficient is requested, a wide matrix-like data frame
#'   is returned with one row per grid point and one column per quantile
#'   level. If \code{TRUE} or if multiple coefficients are requested, a
#'   tidy long data frame is returned with one row per
#'   coefficient-grid-quantile combination.
#'
#' @details
#' When \code{long = FALSE} and a single coefficient is requested, the
#' output has columns:
#' \itemize{
#'   \item \code{Grid_Point} - integer index of the grid location.
#'   \item \code{Grid_Value} - value of the quantile-regressor grid (if
#'         \code{include_grid_value = TRUE}).
#'   \item one column per selected quantile level, named \code{tau_0.05},
#'         \code{tau_0.10}, and so on.
#' }
#'
#' When a long format is returned, the output has columns:
#' \itemize{
#'   \item \code{Coefficient} - coefficient label.
#'   \item \code{Grid_Point} - grid index.
#'   \item \code{Grid_Value} - grid value (if requested).
#'   \item \code{tau} - quantile level.
#'   \item \code{value} - coefficient or p-value at that combination.
#' }
#'
#' Rounding through \code{digits} only affects the returned data frame and
#' does not modify the stored QQR object.
#'
#' @return
#' A data frame in wide or long format, depending on \code{long} and the
#' number of coefficients requested.
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{rowave.qqr}},
#' \code{\link{excel.qqr}}, \code{\link{word.qqr}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' ## Wide table: one row per grid point, one column per tau
#' qqr_surface(fit, coef = "x", digits = 3)
#' ## Long (tidy) table for all coefficients
#' head(qqr_surface(fit, coef = "all", long = TRUE))
#'
#' @export
qqr_surface <- function(x,
                        coef              = "all",
                        what              = c("coef", "pvalue"),
                        grid_subset       = NULL,
                        tau_subset        = NULL,
                        include_grid_value = TRUE,
                        digits            = NULL,
                        long              = NULL) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")

  what <- match.arg(what)

  # Helper to map internal coefficient name to user-facing label
  make_label <- function(internal_name) {
    if (internal_name == "(Intercept)") {
      "Intercept"
    } else if (internal_name == "z") {
      x$independent
    } else {
      internal_name
    }
  }

  # Decide which coefficients to include
  if (length(coef) == 1L && identical(coef, "all")) {
    coef_internal <- x$coef_names
  } else {
    # allow vector of names, including "intercept"
    coef_internal <- vapply(coef, function(cf) .qqr_match_coef(x, cf), character(1L))
    coef_internal <- unique(coef_internal)
  }

  # Decide output shape: if not specified, use long format when >1 coefficient
  if (is.null(long)) {
    long <- length(coef_internal) > 1L
  }

  # Grid subset indices
  # number of rows is the same for all coefficient surfaces
  n_grid <- nrow(x$coef_surfaces[[coef_internal[1L]]])
  if (is.null(grid_subset)) {
    idx_grid <- seq_len(n_grid)
  } else {
    grid_subset[is.na(grid_subset)] <- n_grid
    idx_grid <- unique(pmax(1L, pmin(n_grid, as.integer(grid_subset))))
  }

  # Tau subset indices
  n_tau <- length(x$taus)
  if (is.null(tau_subset)) {
    idx_tau <- seq_len(n_tau)
  } else {
    if (!is.numeric(tau_subset) || anyNA(tau_subset))
      stop("'tau_subset' must be numeric indices or tau values.",
           call. = FALSE)
    if (all(tau_subset < 1)) {
      # tau values: match to x$taus with a floating-point tolerance
      idx_tau <- vapply(tau_subset, function(tv) {
        hit <- which(abs(x$taus - tv) < 1e-8)
        if (length(hit) == 0L) NA_integer_ else hit[1L]
      }, integer(1L))
      if (anyNA(idx_tau))
        stop("tau value(s) not in the fit: ",
             paste(tau_subset[is.na(idx_tau)], collapse = ", "),
             ". Available: ", paste(x$taus, collapse = ", "), call. = FALSE)
    } else {
      idx_tau <- as.integer(tau_subset)
    }
    idx_tau <- unique(pmax(1L, pmin(n_tau, idx_tau)))
  }
  taus_sub <- x$taus[idx_tau]
  tau_cols <- paste0("tau_", taus_sub)

  # Helper to build a matrix for one coefficient
  extract_one <- function(nm) {
    mat_full <- if (what == "coef") {
      x$coef_surfaces[[nm]]
    } else {
      x$pval_surfaces[[nm]]
    }
    mat <- mat_full[idx_grid, idx_tau, drop = FALSE]
    colnames(mat) <- tau_cols
    mat
  }

  # Case 1: single coefficient and user wants wide matrix (long = FALSE)
  if (length(coef_internal) == 1L && !long) {
    mat <- extract_one(coef_internal[1L])
    df  <- as.data.frame(mat)
    df$Grid_Point <- idx_grid
    if (include_grid_value) {
      df$Grid_Value <- x$grid[idx_grid]
      df <- df[, c("Grid_Point", "Grid_Value", tau_cols)]
    } else {
      df <- df[, c("Grid_Point", tau_cols)]
    }

    # Optional rounding for display only
    if (!is.null(digits) && digits >= 0) {
      num_cols <- vapply(df, is.numeric, logical(1L))
      if ("Grid_Point" %in% names(df)) {
        num_cols[names(df) == "Grid_Point"] <- FALSE
      }
      df[num_cols] <- lapply(df[num_cols], function(z) round(z, digits))
    }

    return(df)
  }

  # Case 2: multiple coefficients, or long == TRUE: build tidy long table
  # in base R; one row per coefficient x grid point x tau, ordered by
  # coefficient, then grid point, then tau.
  pieces <- lapply(coef_internal, function(nm) {
    mat <- extract_one(nm)
    df <- data.frame(
      Coefficient = make_label(nm),
      Grid_Point  = rep(idx_grid, each = length(idx_tau)),
      stringsAsFactors = FALSE
    )
    if (include_grid_value)
      df$Grid_Value <- rep(x$grid[idx_grid], each = length(idx_tau))
    df$tau   <- rep(taus_sub, times = length(idx_grid))
    df$value <- as.vector(t(mat))
    df
  })

  long_df <- do.call(rbind, pieces)
  rownames(long_df) <- NULL

  # Optional rounding for display only
  if (!is.null(digits) && digits >= 0) {
    num_cols <- vapply(long_df, is.numeric, logical(1L))
    if ("Grid_Point" %in% names(long_df)) {
      num_cols[names(long_df) == "Grid_Point"] <- FALSE
    }
    long_df[num_cols] <- lapply(long_df[num_cols], function(z) round(z, digits))
  }

  long_df
}


#' Row means of QQR coefficient surfaces across quantiles
#'
#' @description
#' Computes row-wise averages of QQR coefficient surfaces across all
#' conditional quantile levels for each grid point. This can be done for a
#' single coefficient or for all coefficients stored in a \code{"qqr"}
#' object.
#'
#' @param x An object of class \code{"qqr"}, as returned by
#'   \code{\link{qqr}}.
#' @param coef Either \code{"all"} to compute row means for every
#'   coefficient surface (intercept, quantile regressor, and controls), or
#'   a single coefficient name. Coefficients may be referred to by
#'   \code{"intercept"} or \code{"(Intercept)"}, by the name of the
#'   quantile regressor, or by control variable names.
#' @param na.rm Logical flag indicating whether \code{NA} values should be
#'   removed when computing the row means. The default is \code{TRUE}.
#' @param digits Optional integer number of decimal places to which the
#'   numeric averages should be rounded for display. If \code{NULL},
#'   values are returned without rounding. Rounding does not modify the
#'   underlying QQR object.
#' @param include_grid_value Logical flag indicating whether the output
#'   should include the \code{grid_value} column (the underlying value of
#'   the quantile-regressor grid) alongside the integer grid index. The
#'   default is \code{TRUE}, matching the behaviour of the Excel export
#'   helpers. Setting this to \code{FALSE} drops the \code{grid_value}
#'   column and returns only the index and row means.
#' @param ... Currently unused. Included for potential method
#'   extensibility.
#'
#' @details
#' When \code{coef = "all"}, the function returns a data frame with one
#' row per grid point and one column per coefficient (after the grid
#' information), where each entry is the average of the corresponding
#' coefficient over all quantile levels. When a single coefficient is
#' requested, only one column of averages is returned.
#'
#' The \code{include_grid_value} argument controls whether the numeric
#' grid values are attached as a separate column. When
#' \code{include_grid_value = TRUE}, both \code{grid_index} (the position
#' of the grid point) and \code{grid_value} (the actual value of the
#' quantile regressor at that grid point) are included. When
#' \code{include_grid_value = FALSE}, only \code{grid_index} is returned,
#' which may be convenient when the grid values are not needed or when a
#' more compact table is desired.
#'
#' This summary is useful when a compact measure of the average effect of
#' each coefficient across the conditional quantile spectrum is required.
#'
#' @return
#' A data frame with columns:
#' \itemize{
#'   \item \code{grid_index} - integer index of the grid location.
#'   \item \code{grid_value} - value of the quantile-regressor grid at
#'         each index (included only when
#'         \code{include_grid_value = TRUE}).
#'   \item one or more columns of row means, either one per coefficient
#'         (when \code{coef = "all"}) or a single column for the chosen
#'         coefficient.
#' }
#'
#' @seealso
#' \code{\link{qqr}}, \code{\link{qqr_surface}}
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(x = rnorm(150), w = rnorm(150))
#' dat$y <- 1 + dat$x + 0.5 * dat$w + (1 + 0.3 * abs(dat$x)) * rnorm(150)
#' fit <- qqr(y ~ x + w, data = dat, qvar = "x",
#'            taus = c(0.25, 0.5, 0.75), m = 5)
#' rowave.qqr(fit, coef = "x", digits = 3)
#' rowave.qqr(fit, coef = "all", include_grid_value = FALSE)
#'
#' @export
rowave.qqr <- function(x,
                       coef               = "all",
                       na.rm              = TRUE,
                       digits             = NULL,
                       include_grid_value = TRUE,
                       ...) {

  if (!inherits(x, "qqr")) stop("Object must be of class 'qqr'.")

  # Helper to map internal coefficient name to a user-facing label
  make_label <- function(internal_name) {
    if (internal_name == "(Intercept)") {
      "Intercept"
    } else if (internal_name == "z") {
      x$independent
    } else {
      internal_name
    }
  }

  coef_names_internal <- x$coef_names
  grid_index <- seq_along(x$grid)

  # If coef = "all", compute row means for every coefficient surface
  if (length(coef) == 1L && identical(coef, "all")) {

    mean_mat_list <- lapply(coef_names_internal, function(nm) {
      mat <- x$coef_surfaces[[nm]]
      rowMeans(mat, na.rm = na.rm)
    })

    mean_mat <- do.call(cbind, mean_mat_list)
    colnames(mean_mat) <- vapply(coef_names_internal, make_label, character(1L))

    if (include_grid_value) {
      df <- data.frame(
        grid_index = grid_index,
        grid_value = x$grid,
        mean_mat,
        check.names = FALSE
      )
    } else {
      df <- data.frame(
        grid_index = grid_index,
        mean_mat,
        check.names = FALSE
      )
    }

  } else {
    # User requested a single coefficient (possibly with user-facing name)
    surf_name <- .qqr_match_coef(x, coef)
    mat       <- x$coef_surfaces[[surf_name]]
    avg       <- rowMeans(mat, na.rm = na.rm)

    lab <- make_label(surf_name)

    if (include_grid_value) {
      df <- data.frame(
        grid_index = grid_index,
        grid_value = x$grid,
        mean_coef  = avg,
        check.names = FALSE
      )
    } else {
      df <- data.frame(
        grid_index = grid_index,
        mean_coef  = avg,
        check.names = FALSE
      )
    }
    names(df)[names(df) == "mean_coef"] <- paste0("mean_", lab)
  }

  # Optional formatting only, does not touch the qqr object
  if (!is.null(digits)) {
    num_cols <- vapply(df, is.numeric, logical(1L))
    if ("grid_index" %in% names(df)) {
      num_cols[names(df) == "grid_index"] <- FALSE
    }
    df[num_cols] <- lapply(df[num_cols], function(z) round(z, digits))
  }

  df
}


