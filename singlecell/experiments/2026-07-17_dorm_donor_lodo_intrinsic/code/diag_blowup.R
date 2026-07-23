#!/usr/bin/env Rscript
## Targeted blow-up diagnostic on s1d1/CD72_1/seed1: expose signed rho, induced-mixture min,
## posterior weight range, and beta-norm for standard vs intrinsic (signed) prior.
suppressPackageStartupMessages({ library(data.table) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

## ---- reconstruct the exact combo (match run: HVG=200, MAXN_SRC=400, MAXN_TGT=1200, seed=1) ----
D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, 200)
target <- "s1d1"; y <- "CD72_1"; K <- 20
set.seed(1)
src_donors <- setdiff(sort(unique(D$meta$batch)), target)
src_rows <- lapply(src_donors, function(b) sample_at_most(which(D$meta$batch == b), 400)); names(src_rows) <- src_donors
src_rows <- src_rows[lengths(src_rows) >= 60]; src_donors <- names(src_rows)
tgt_rows <- sample_at_most(which(D$meta$batch == target), 1200)
A <- screen_A(D, unlist(src_rows, use.names = FALSE), y, K); q <- K
ssp <- lapply(src_rows, split_half); tsp <- split_half(tgt_rows)
Xl <- lapply(ssp, function(s) make_X(D, s$a, A)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A))
Yl <- lapply(ssp, function(s) make_Y(D, s$a, y))
X0 <- make_X(D, tsp$a, A); Y0 <- make_Y(D, tsp$a, y); X0t <- make_X(D, tsp$b, A)
L <- length(Xl); n0 <- nrow(X0)
cat(sprintf("Combo s1d1/CD72_1/seed1: L=%d sources, n0(eval)=%d, var(Y0)=%.3f\n\n", L, n0, var(Y0)))

## ---- replicate ONE cross-fit half's density ratios (as maximin_s_beta does), then both priors ----
models  <- or_estimation_logit(Xt, X0t, condA = FALSE, q = q)
dr_cond <- dratio_logit(models, X0, L, Nlist = sapply(Xt, nrow), n0 = nrow(X0t), normalize = FALSE)   # L x n0
drlist  <- lapply(seq_len(L), function(l) dratio_logit(models, Xl[[l]], L, Nlist = sapply(Xt, nrow), n0 = nrow(X0t), normalize = FALSE))

rho_std <- mulcvxr(L, dr_cond)
rho_int <- signed_affine_prior(dr_cond, c(list(dr_cond), drlist))

report <- function(tag, rho) {
  induced <- as.vector(t(rho) %*% (1 / dr_cond))          # sum_l rho_l * rbar_l(x) at each target point
  post <- posterior(rho, dr_cond, n0)                      # n0 x L posterior weights
  cat(sprintf("[%s prior]\n", tag))
  cat(sprintf("  rho: sum=%.3f  min=%.3f  max=%.3f  max|rho|=%.3f  #neg=%d\n",
              sum(rho), min(rho), max(rho), max(abs(rho)), sum(rho < 0)))
  cat(sprintf("  induced mixture ratio at target: min=%.3g  1%%q=%.3g  median=%.3g  (#<1e-2: %d, #<1e-3: %d)\n",
              min(induced), quantile(induced, .01), median(induced), sum(induced < 1e-2), sum(induced < 1e-3)))
  cat(sprintf("  posterior weights: min=%.3g  max=%.3g  max|.|=%.3g  #neg=%d\n\n",
              min(post), max(post), max(abs(post)), sum(post < 0)))
}
cat("rho values (standard vs intrinsic), per source donor:\n")
print(data.table(donor = src_donors, rho_std = round(rho_std, 3), rho_int = round(rho_int, 3)), row.names = FALSE)
cat("\n")
report("standard (simplex)", rho_std)
report("intrinsic (signed)", rho_int)

## ---- SOURCE-SIDE: this is where P/secondPmat lives (po_source * w, w = density ratio up to 1e6) ----
cat("Source-side induced-mixture min & posterior range (drives P/secondPmat):\n")
src_tab <- rbindlist(lapply(seq_len(L), function(l) {
  ind_s <- as.vector(t(rho_std) %*% (1 / drlist[[l]])); ind_i <- as.vector(t(rho_int) %*% (1 / drlist[[l]]))
  pl_s <- posterior(rho_std, drlist[[l]], ncol(drlist[[l]])); pl_i <- posterior(rho_int, drlist[[l]], ncol(drlist[[l]]))
  data.table(src = src_donors[l], max_w = max(drlist[[l]][l, ]),
             indmin_std = min(ind_s), post_absmax_std = max(abs(pl_s)),
             indmin_int = min(ind_i), post_absmax_int = max(abs(pl_i)))
}))
print(src_tab[, .(src, max_w = signif(max_w,3), indmin_std = signif(indmin_std,3), post_absmax_std = signif(post_absmax_std,3),
                  indmin_int = signif(indmin_int,3), post_absmax_int = signif(post_absmax_int,3))], row.names = FALSE)
cat(sprintf("\nWorst source-side posterior |weight|: standard=%.3g  intrinsic=%.3g\n\n",
            max(src_tab$post_absmax_std), max(src_tab$post_absmax_int)))

## ---- compute P directly from THIS realization via PQCalculation (the clincher) ----
Yt <- lapply(ssp, function(s) make_Y(D, s$b, y)); nl <- as.integer(lengths(Yl)); Nl <- sapply(Xl, nrow)
betalist <- lasso_imputation(Xt, Yt, nl)
pq_of <- function(rho) {
  post <- posterior(rho, dr_cond, n0)
  postlist <- lapply(seq_len(L), function(l) posterior(rho, drlist[[l]], Nl[l]))
  PQCalculation(X0, Xl, Yl, X0t, Xt, Yt, nl, post, postlist, drlist, betalist, q, penalty = TRUE, alpha = 0.5)
}
pq_s <- pq_of(rho_std); pq_i <- pq_of(rho_int)
cat(sprintf("This realization: ||P|| std=%.2f  int=%.2f   |  max||Q col|| std=%.2f int=%.2f\n",
            sqrt(sum(pq_s$P^2)), sqrt(sum(pq_i$P^2)),
            max(sqrt(colSums(pq_s$Q^2))), max(sqrt(colSums(pq_i$Q^2)))))
## how much of P is the source-side (secondPmat, ~ po_source * density-ratio w) vs target-side
pre_s <- DRcoef_calculation(X0, Xl, Yl, nl, posterior(rho_std, dr_cond, n0),
                            lapply(seq_len(L), function(l) posterior(rho_std, drlist[[l]], Nl[l])), drlist, betalist, q)
pre_i <- DRcoef_calculation(X0, Xl, Yl, nl, posterior(rho_int, dr_cond, n0),
                            lapply(seq_len(L), function(l) posterior(rho_int, drlist[[l]], Nl[l])), drlist, betalist, q)
cat(sprintf("||preP|| (pre-solve DR pooled coef) std=%.2f  int=%.2f\n", sqrt(sum(pre_s$preP^2)), sqrt(sum(pre_i$preP^2))))
