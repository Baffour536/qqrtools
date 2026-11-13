#' qqrtools: Tools for Quantile-on-Quantile Regression
#'
#' @description
#' \pkg{qqrtools} estimates quantile-on-quantile regression (QQR) surfaces
#' (Sim and Zhou, 2015): the effect of the \eqn{\theta}-quantile of a
#' regressor \eqn{X} on the \eqn{\tau}-quantile of a response \eqn{Y},
#' obtained from kernel-weighted local polynomial quantile regressions.
#'
#' @section Main functions:
#' \describe{
#'   \item{Estimation}{\code{\link{qqr}} fits the model and returns an object
#'     of class \code{"qqr"} with coefficient and p-value surfaces for the
#'     intercept, the quantile regressor and any controls.}
#'   \item{Summaries}{\code{print}, \code{\link{summary.qqr}},
#'     \code{\link{qqr_surface}} (wide or long tables) and
#'     \code{\link{rowave.qqr}} (averages over \eqn{\tau}).}
#'   \item{Robustness}{\code{\link{qqr_sensitivity}} refits the surface over
#'     a range of bandwidths; \code{\link{qqr_padjust}} adjusts the p-value
#'     surfaces for the number of cells tested.}
#'   \item{Comparison with QR}{\code{\link{compare.qqr}} averages each
#'     surface over \eqn{\theta} and compares it with standard quantile
#'     regression; \code{\link{plot.compare.qqr}},
#'     \code{\link{excel.compare.qqr}} and \code{\link{word.compare.qqr}}
#'     plot and export the comparison.}
#'   \item{Graphics}{\code{\link{surface3d}}, \code{\link{panel_surface}} and
#'     \code{\link{plot.qqr}} (base graphics); \code{\link{Int.plot}}
#'     (\pkg{rgl}); \code{\link{qqr_heatmap}}, \code{\link{qqr_contour}} and
#'     \code{\link{qqr_plotly}} (\pkg{plotly}).}
#'   \item{Data diagnostics}{\code{\link{qqr_cor_heatmap}} summarises the
#'     bivariate quantile dependence of \eqn{X} and \eqn{Y}.}
#'   \item{Export}{\code{\link{excel.qqr}} and \code{\link{word.qqr}}.}
#' }
#'
#' @references
#' Sim, N. and Zhou, H. (2015) Oil prices, US stock return, and the
#' dependence between their quantiles. \emph{Journal of Banking & Finance},
#' 55, 1--8. \doi{10.1016/j.jbankfin.2015.01.013}
#'
#' Koenker, R. (2005) \emph{Quantile Regression}. Cambridge University
#' Press. \doi{10.1017/CBO9780511754098}
#'
#' @keywords internal
#' @importFrom stats setNames
#' @importFrom graphics par layout persp polygon text axis box image
#' @importFrom grDevices hcl.colors adjustcolor
"_PACKAGE"
