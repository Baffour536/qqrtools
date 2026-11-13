## =========================================================================
## Monte Carlo validation of qqrtools::qqr()
##
## Run from the package root after installing the package:
##   Rscript tools/validation/simulation.R
## Results are written to tools/validation/simulation_results.rds and a
## Markdown summary to tools/validation/simulation_summary.md.
##
## Designs (X ~ U(0,1), eps ~ N(0,1), theta = tau grid {.1,.25,.5,.75,.9}):
##   A  location-scale  Y = 1 + 2X + (0.5 + X) eps
##        true beta(theta, tau) = 2 + qnorm(tau)           (flat in theta)
##   B  curved          Y = 1 + 2X^2 + (0.5 + X) eps
##        true beta(theta, tau) = 4 theta + qnorm(tau)     (varies in both)
##   C  null            Y = 1 + eps (independent of X)
##        true beta = 0; used for the size of the pointwise p-values
## =========================================================================

suppressPackageStartupMessages(library(qqrtools))
library(parallel)

R_REP  <- as.integer(Sys.getenv("QQR_SIM_REPS", "200"))
SEED   <- 20261113L
grid_p <- c(0.1, 0.25, 0.5, 0.75, 0.9)

dgp <- function(design, n) {
  x <- stats::runif(n)
  e <- stats::rnorm(n)
  y <- switch(design,
              A = 1 + 2 * x + (0.5 + x) * e,
              B = 1 + 2 * x^2 + (0.5 + x) * e,
              C = 1 + e)
  data.frame(x = x, y = y)
}

truth <- function(design) {
  q <- stats::qnorm(grid_p)
  switch(design,
         A = matrix(2 + q, 5, 5, byrow = TRUE),
         B = outer(4 * grid_p, q, "+"),
         C = matrix(0, 5, 5))
}

## one replication: returns the slope surface (and p-values for C)
one_rep <- function(rep, design, n, h, space, se) {
  set.seed(SEED + rep)
  d   <- dgp(design, n)
  fit <- suppressWarnings(suppressMessages(
    qqr(y ~ x, data = d, taus = grid_p, thetas = grid_p, h = h,
        kernel_space = space, se_type = se)))
  list(coef = fit$qcoef, pval = fit$pval_surfaces$x, h = fit$bandwidth)
}

configs <- rbind(
  expand.grid(design = c("A", "B"), n = c(250L, 500L, 1000L),
              h = c("silverman", "robust", "0.1", "0.2", "0.3"),
              space = "x", se = "none", stringsAsFactors = FALSE),
  expand.grid(design = c("A", "B"), n = c(250L, 500L, 1000L),
              h = c("0.05", "0.1", "0.2"),
              space = "cdf", se = "none", stringsAsFactors = FALSE),
  expand.grid(design = "C", n = c(250L, 1000L),
              h = c("silverman", "0.2"), space = "x", se = "nid",
              stringsAsFactors = FALSE)
)

cl <- makeCluster(max(1L, detectCores() - 2L))
on.exit(stopCluster(cl), add = TRUE)
clusterEvalQ(cl, suppressPackageStartupMessages(library(qqrtools)))
clusterExport(cl, c("dgp", "one_rep", "grid_p", "SEED"))

t0 <- Sys.time()
results <- lapply(seq_len(nrow(configs)), function(k) {
  cf <- configs[k, ]
  hv <- if (cf$h %in% c("silverman", "robust")) cf$h else as.numeric(cf$h)
  reps <- parLapply(cl, seq_len(R_REP), one_rep, design = cf$design,
                    n = cf$n, h = hv, space = cf$space, se = cf$se)
  est  <- simplify2array(lapply(reps, `[[`, "coef"))      # 5 x 5 x R
  tr   <- truth(cf$design)
  bias <- apply(est, c(1, 2), mean, na.rm = TRUE) - tr
  sdv  <- apply(est, c(1, 2), stats::sd, na.rm = TRUE)
  rmse <- sqrt(apply(sweep(est, c(1, 2), tr)^2, c(1, 2), mean,
                     na.rm = TRUE))
  ## theta-averaged slope (compare.qqr) vs truth averaged over theta
  avg  <- apply(est, c(2, 3), mean)                         # 5 (tau) x R
  avg_rmse <- sqrt(mean(sweep(avg, 1, colMeans(tr))^2, na.rm = TRUE))
  size <- NA_real_
  if (cf$se != "none") {
    pv   <- simplify2array(lapply(reps, `[[`, "pval"))
    size <- mean(pv < 0.05, na.rm = TRUE)
    size_cell <- apply(pv < 0.05, c(1, 2), mean, na.rm = TRUE)
  } else size_cell <- NULL
  message(sprintf("[%2d/%2d] %s n=%4d h=%-9s %s  done (%.0fs)", k,
                  nrow(configs), cf$design, cf$n, cf$h, cf$space,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  list(config = cf, h_mean = mean(sapply(reps, `[[`, "h")),
       bias = bias, sd = sdv, rmse = rmse,
       mean_abs_bias = mean(abs(bias)), mean_sd = mean(sdv),
       mean_rmse = mean(rmse),
       interior_rmse = mean(rmse[2:4, 2:4]),
       avg_rmse = avg_rmse, size = size, size_cell = size_cell,
       na_share = mean(is.na(est)))
})

saveRDS(list(results = results, reps = R_REP, seed = SEED,
             session = list(R = R.version.string,
                            platform = R.version$platform,
                            qqrtools = as.character(packageVersion("qqrtools")),
                            quantreg = as.character(packageVersion("quantreg")))),
        "tools/validation/simulation_results.rds")

## ---- Markdown summary ---------------------------------------------------
tab <- do.call(rbind, lapply(results, function(r) data.frame(
  design = r$config$design, n = r$config$n, h = r$config$h,
  kernel_space = r$config$space, h_used = round(r$h_mean, 3),
  mean_abs_bias = round(r$mean_abs_bias, 3), mean_sd = round(r$mean_sd, 3),
  mean_rmse = round(r$mean_rmse, 3), interior_rmse = round(r$interior_rmse, 3),
  theta_avg_rmse = round(r$avg_rmse, 3),
  size_5pct = if (is.na(r$size)) "" else sprintf("%.3f", r$size),
  na_share = round(r$na_share, 3))))

md <- c("# qqrtools Monte Carlo summary", "",
        sprintf("Replications per design: %d; seed base %d; generated %s.",
                R_REP, SEED, format(Sys.time(), "%Y-%m-%d %H:%M")), "",
        "Averages over the 5 x 5 (theta, tau) grid {0.1, 0.25, 0.5, 0.75, 0.9}^2.",
        "`interior_rmse` averages the 3 x 3 interior cells; `theta_avg_rmse`",
        "is the RMSE of the theta-averaged slope used by compare.qqr();",
        "`size_5pct` is the rejection rate of H0: beta = 0 at 5% (design C,",
        "true beta = 0) pooled over cells.", "",
        paste0("| ", paste(names(tab), collapse = " | "), " |"),
        paste0("|", paste(rep("---", ncol(tab)), collapse = "|"), "|"),
        apply(tab, 1, function(r) paste0("| ", paste(r, collapse = " | "),
                                         " |")))
for (r in results) if (!is.null(r$size_cell)) {
  md <- c(md, "", sprintf("Cell rejection rates, design C, n = %d, h = %s",
                          r$config$n, r$config$h), "",
          "| theta \\ tau | 0.1 | 0.25 | 0.5 | 0.75 | 0.9 |",
          "|---|---|---|---|---|---|",
          sapply(1:5, function(i) paste0("| ", grid_p[i], " | ",
                 paste(sprintf("%.3f", r$size_cell[i, ]), collapse = " | "),
                 " |")))
}
for (r in results) if (r$config$design == "B" && r$config$n == 1000 &&
                       r$config$space == "x" &&
                       r$config$h %in% c("silverman", "0.2")) {
  md <- c(md, "", sprintf("Bias by cell, design B, n = 1000, h = %s",
                          r$config$h), "",
          "| theta \\ tau | 0.1 | 0.25 | 0.5 | 0.75 | 0.9 |",
          "|---|---|---|---|---|---|",
          sapply(1:5, function(i) paste0("| ", grid_p[i], " | ",
                 paste(sprintf("%+.3f", r$bias[i, ]), collapse = " | "),
                 " |")))
}
writeLines(md, "tools/validation/simulation_summary.md")
cat(md, sep = "\n")
