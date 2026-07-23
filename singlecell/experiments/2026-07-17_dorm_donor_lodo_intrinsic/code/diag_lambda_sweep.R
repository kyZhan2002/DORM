#!/usr/bin/env Rscript
## Pick a stabilizing lambda for the intrinsic (signed-affine) prior on the known blow-up combo
## s1d1/CD72_1: for each lambda, run 3 fresh density-ratio realizations and report ||P||/R^2 spread.
suppressPackageStartupMessages({ library(data.table) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

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
tuned <- function(fit) { nt <- min(100, length(Y0t)); tn <- tuning_Y(X0t[1:nt,,drop=FALSE], Y0t[1:nt], q, fit$beta_by_smax, SMAX_GRID)
  1 - get_err(X0, Y0, as.vector(fit$beta_by_smax[which(SMAX_GRID==tn$best_smax)[1],]), q, length(Y0)) / yv }

out <- rbindlist(lapply(c(0.1, 0.5, 1, 2, 5, 20), function(lam) {
  SIGNED_LAMBDA <<- lam
  r <- rbindlist(lapply(1:3, function(k) {
    set.seed(300 + k); fi <- dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nl, q, penalty = TRUE, signed_prior = TRUE)
    data.table(Pn = sqrt(sum(fi$P^2)), r2 = tuned(fi), rho_absmax = max(abs(fi$rho)))
  }))
  data.table(lambda = lam, Pnorm_min = round(min(r$Pn),1), Pnorm_max = round(max(r$Pn),1),
             r2_min = round(min(r$r2),3), r2_max = round(max(r$r2),3), rho_absmax = round(max(r$rho_absmax),2))
}))
cat("\n=== intrinsic stability vs lambda (s1d1/CD72_1, 3 realizations each) ===\n")
print(out, row.names = FALSE)
cat("\n(standard DORM here: ||P||~0.9, R^2~0.73 for reference)\n")
