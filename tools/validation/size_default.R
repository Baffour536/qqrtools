## =========================================================================
## Size of the pointwise 5% test of beta(theta, tau) = 0 with the default
## bandwidth and degree (h = "slope", degree = 2) for the three
## standard-error types, and the effect of the Benjamini-Hochberg
## adjustment of qqr_padjust().
## Null designs: Y = 1 + eps, eps ~ N(0, 1), independent of
##   N: X ~ N(0, 1);   S: X ~ chi-square(3) (skewed).
##
## Run from the package root:  Rscript tools/validation/size_default.R
## =========================================================================

suppressPackageStartupMessages(library(qqrtools))

reps <- 300L
qs   <- c(0.1, 0.25, 0.5, 0.75, 0.9)
rx   <- list(N = function(n) rnorm(n), S = function(n) rchisq(n, 3))
out  <- list()
for (dn in names(rx)) for (n in c(250L, 1000L)) {
  for (se in c("nid", "ker", "boot")) {
    rej <- array(NA, c(reps, 5, 5))
    any_raw <- any_bh <- rep(NA, reps)
    for (r in seq_len(reps)) {
      set.seed(20261113L + r)
      dat <- data.frame(x = rx[[dn]](n), y = 1 + rnorm(n))
      f <- if (se == "boot")
        qqr(y ~ x, data = dat, taus = qs, thetas = qs, se_type = "boot",
            R = 100)
      else qqr(y ~ x, data = dat, taus = qs, thetas = qs, se_type = se)
      p <- f$pval_surfaces$x
      rej[r, , ] <- p < 0.05
      any_raw[r] <- any(p < 0.05, na.rm = TRUE)
      any_bh[r]  <- any(qqr_padjust(f)$pval_surfaces$x < 0.05, na.rm = TRUE)
    }
    cellrate <- apply(rej, c(2, 3), mean, na.rm = TRUE)
    out[[length(out) + 1L]] <- data.frame(
      design = dn, n = n, se_type = se,
      pooled_rejection_rate = mean(rej, na.rm = TRUE),
      max_cell_rate = max(cellrate),
      any_cell_significant = mean(any_raw),
      any_cell_significant_BH = mean(any_bh))
    cat("done", dn, n, se, "\n")
  }
}
res <- do.call(rbind, out)
num <- vapply(res, is.double, logical(1L))
res[num] <- lapply(res[num], round, 3)
lines <- c(
  "# Size of pointwise tests with the default bandwidth and degree",
  "",
  paste0("Null designs (no effect; N: normal X, S: skewed X), ", reps,
         " replications, 5 x 5 cells, h = \"slope\", degree = 2, nominal ",
         "level 5%; bootstrap with R = 100. The last two columns give the ",
         "share of fitted surfaces with at least one significant cell ",
         "before and after `qqr_padjust(method = \"BH\")`."),
  "",
  paste0("| ", paste(names(res), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(res)), collapse = "|"), "|"),
  apply(res, 1L, function(r) paste0("| ", paste(trimws(r), collapse = " | "),
                                    " |"))
)
writeLines(lines, "tools/validation/size_default.md")
cat(lines, sep = "\n")
