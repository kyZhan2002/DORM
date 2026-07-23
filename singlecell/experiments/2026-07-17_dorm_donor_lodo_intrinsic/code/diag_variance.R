#!/usr/bin/env Rscript
## Airtight variance check: on the FIXED s1d1/CD72_1 combo, run standard vs intrinsic DORM over
## several fresh density-ratio realizations (same data & splits; only the internal cv.glmnet RNG
## varies, shared between the two priors per realization). Expect: standard ||P||/R^2 stable;
## intrinsic swinging tame <-> blow-up, depending on whether a source dominates the estimated ratio.
suppressPackageStartupMessages({ library(data.table) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

## fixed combo/data
D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, 200)
target <- "s1d1"; y <- "CD72_1"; K <- 20
set.seed(1)
src_donors <- setdiff(sort(unique(D$meta$batch)), target)
src_rows <- lapply(src_donors, function(b) sample_at_most(which(D$meta$batch == b), 400)); names(src_rows) <- src_donors
tgt_rows <- sample_at_most(which(D$meta$batch == target), 1200)
A <- screen_A(D, unlist(src_rows, use.names = FALSE), y, K); q <- K
ssp <- lapply(src_rows, split_half); tsp <- split_half(tgt_rows)
Xl <- lapply(ssp, function(s) make_X(D, s$a, A)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A))
Yl <- lapply(ssp, function(s) make_Y(D, s$a, y)); Yt <- lapply(ssp, function(s) make_Y(D, s$b, y))
X0 <- make_X(D, tsp$a, A); Y0 <- make_Y(D, tsp$a, y); X0t <- make_X(D, tsp$b, A); Y0t <- make_Y(D, tsp$b, y)
nl <- as.integer(lengths(Yl)); yv <- var(Y0)

tuned <- function(fit) {
  nt <- min(100, length(Y0t))
  tn <- tuning_Y(X0t[1:nt, , drop = FALSE], Y0t[1:nt], q, fit$beta_by_smax, SMAX_GRID)
  b <- fit$beta_by_smax[which(SMAX_GRID == tn$best_smax)[1], ]
  1 - get_err(X0, Y0, as.vector(b), q, length(Y0)) / yv
}

res <- rbindlist(lapply(1:5, function(k) {
  set.seed(200 + k); fs <- dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nl, q, penalty = TRUE, signed_prior = FALSE)
  set.seed(200 + k); fi <- dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nl, q, penalty = TRUE, signed_prior = TRUE)
  data.table(realization = k,
             Pnorm_std = round(sqrt(sum(fs$P^2)), 2), r2_std = round(tuned(fs), 3),
             Pnorm_int = round(sqrt(sum(fi$P^2)), 2), r2_int = round(tuned(fi), 3),
             rho_int_absmax = round(max(abs(fi$rho)), 2))
}))
cat("\n=== standard vs intrinsic across 5 fresh density-ratio realizations (same data) ===\n")
print(res, row.names = FALSE)
cat(sprintf("\nStandard:  ||P|| in [%.2f, %.2f], R^2 in [%.3f, %.3f]\n", min(res$Pnorm_std), max(res$Pnorm_std), min(res$r2_std), max(res$r2_std)))
cat(sprintf("Intrinsic: ||P|| in [%.2f, %.2f], R^2 in [%.3f, %.3f]\n", min(res$Pnorm_int), max(res$Pnorm_int), min(res$r2_int), max(res$r2_int)))
