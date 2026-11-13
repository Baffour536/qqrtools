## =========================================================================
## Calibration of the slope-oriented bandwidth rule and of the default
## local polynomial degree.
##
## Rule: h = C * min(sd, IQR/1.34) * n^(-1/7). For each design, sample size,
## constant C and degree (1 = local linear, 2 = local quadratic) the script
## reports bias, SD and RMSE of the slope over a 5 x 5 (theta, tau) grid.
## Silverman's rule and the two cross-validation rules are included for
## reference.
##
## Run from the package root:  Rscript tools/validation/bandwidth_calibration.R
## =========================================================================

suppressPackageStartupMessages(library(qqrtools))

reps   <- 200L
qs     <- c(0.1, 0.25, 0.5, 0.75, 0.9)
ns     <- c(250L, 1000L)
consts <- c(1.06, 1.5, 2, 2.5, 3)

## Designs: Y = m(X) + s(X) * eps, eps ~ N(0, 1); the true QQR slope is
## m'(x) + s'(x) * qnorm(tau) at x = Q_X(theta).
designs <- list(
  ## D1: X uniform, slope linear in x (as design B of simulation.R)
  D1 = list(rx = function(n) runif(n), qx = function(p) p,
            m1 = function(x) 4 * x, s = function(x) 0.5 + x,
            s1 = function(x) 1 + 0 * x,
            m  = function(x) 1 + 2 * x^2),
  ## D2: X normal, sinusoidal mean (non-zero third derivative)
  D2 = list(rx = function(n) rnorm(n), qx = qnorm,
            m1 = function(x) 1.5 * cos(1.5 * x),
            s = function(x) 0.5 + 0.25 * abs(x),
            s1 = function(x) 0.25 * sign(x),
            m  = function(x) sin(1.5 * x)),
  ## D3: X skewed (chi-square 3), exponential-type mean
  D3 = list(rx = function(n) rchisq(n, 3), qx = function(p) qchisq(p, 3),
            m1 = function(x) 0.6 * exp(-0.3 * x),
            s = function(x) 0.5 + 0.1 * x,
            s1 = function(x) 0.1 + 0 * x,
            m  = function(x) 2 * (1 - exp(-0.3 * x)))
)

rob_scale <- function(x) min(sd(x), IQR(x) / 1.34)

settings <- rbind(
  expand.grid(rule = paste0("C=", consts), degree = 1:2,
              stringsAsFactors = FALSE),
  expand.grid(rule = c("silverman", "cv_fast", "cv_slope"), degree = 1:2,
              stringsAsFactors = FALSE)
)

out <- list()
for (dn in names(designs)) {
  d <- designs[[dn]]
  for (n in ns) {
    truth <- outer(d$qx(qs), qs, function(x, tau)
      d$m1(x) + d$s1(x) * qnorm(tau))
    est <- array(NA_real_, c(nrow(settings), reps, 5, 5))
    hs  <- matrix(NA_real_, nrow(settings), reps)
    for (r in seq_len(reps)) {
      set.seed(20261113L + r)
      x <- d$rx(n)
      y <- d$m(x) + d$s(x) * rnorm(n)
      dat <- data.frame(y = y, x = x)
      hcv <- NULL
      for (k in seq_len(nrow(settings))) {
        rl <- settings$rule[k]
        h  <- if (startsWith(rl, "C=")) {
          as.numeric(sub("C=", "", rl)) * rob_scale(x) * n^(-1 / 7)
        } else if (rl == "silverman") {
          1.06 * sd(x) * n^(-1 / 5)
        } else {
          if (n > 250L) next            # CV rules only at n = 250 (slow)
          if (is.null(hcv))
            hcv <- qqrtools:::.qqr_select_h(x, y, "cv_fast")$h
          if (rl == "cv_fast") hcv else hcv * n^(1 / 5 - 1 / 7)
        }
        hs[k, r] <- h
        f <- suppressWarnings(
          qqr(y ~ x, data = dat, taus = qs, thetas = qs, h = h,
              se_type = "none", degree = settings$degree[k]))
        est[k, r, , ] <- f$qcoef
      }
    }
    for (k in seq_len(nrow(settings))) {
      if (all(is.na(hs[k, ]))) next
      e    <- est[k, , , , drop = TRUE]
      bias <- apply(e, c(2, 3), mean, na.rm = TRUE) - truth
      sdv  <- apply(e, c(2, 3), sd, na.rm = TRUE)
      rmse <- sqrt(bias^2 + sdv^2)
      out[[length(out) + 1L]] <- data.frame(
        design = dn, n = n, rule = settings$rule[k],
        degree = settings$degree[k],
        h_mean = mean(hs[k, ], na.rm = TRUE),
        mean_abs_bias = mean(abs(bias)), mean_sd = mean(sdv),
        mean_rmse = mean(rmse),
        edge_rmse = mean(rmse[c(1, 5), ]),
        interior_rmse = mean(rmse[2:4, 2:4]),
        na_share = mean(is.na(e)))
    }
    cat("done", dn, n, "\n")
  }
}
res <- do.call(rbind, out)
saveRDS(res, "tools/validation/bandwidth_calibration.rds")

num <- vapply(res, is.numeric, logical(1L))
res[num] <- lapply(res[num], function(v) round(v, 3))
lines <- c(
  "# Bandwidth and degree calibration",
  "",
  paste0("Replications per setting: ", reps, "; seed base 20261113. ",
         "Rule `C=c` is h = c * min(sd, IQR/1.34) * n^(-1/7). ",
         "Averages over the 5 x 5 (theta, tau) grid ",
         "{0.1, 0.25, 0.5, 0.75, 0.9}^2; `edge_rmse` averages the rows ",
         "theta = 0.1 and 0.9. The cross-validation rules were run at ",
         "n = 250 only."),
  "",
  paste0("| ", paste(names(res), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(res)), collapse = "|"), "|"),
  apply(res, 1L, function(r) paste0("| ", paste(trimws(r), collapse = " | "),
                                    " |"))
)
writeLines(lines, "tools/validation/bandwidth_calibration.md")
cat("finished\n")
