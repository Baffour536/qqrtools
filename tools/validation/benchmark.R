## Timing benchmark for qqrtools::qqr(). Run from the package root:
##   Rscript tools/validation/benchmark.R
suppressPackageStartupMessages(library(qqrtools))
set.seed(1)
bench <- function(n, se, h = "silverman", m = 19, K = 19) {
  d <- data.frame(x = rnorm(n), w = rnorm(n))
  d$y <- d$x + 0.5 * d$w + rnorm(n)
  taus <- seq(0.05, 0.95, length.out = K)
  t <- system.time(suppressWarnings(
    qqr(y ~ x + w, d, qvar = "x", taus = taus, m = m, h = h, se_type = se)))
  data.frame(n = n, m = m, K = K, h = h, se_type = se,
             seconds = round(unname(t["elapsed"]), 2))
}
res <- rbind(
  bench(200, "none"), bench(200, "nid"),
  bench(1000, "none"), bench(1000, "nid"),
  bench(5000, "none"), bench(5000, "nid"),
  bench(1000, "none", h = "cv_fast"), bench(3000, "none", h = "cv_fast"),
  bench(1000, "none", h = "sj"), bench(1000, "none", h = "yj")
)
print(res, row.names = FALSE)
cat("\n", R.version.string, "|", Sys.info()[["sysname"]], "| quantreg",
    as.character(packageVersion("quantreg")), "\n")
