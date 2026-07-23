#!/usr/bin/env Rscript
## Does the density-ratio FEATURE DIMENSION explain why demixing finds no reduction (kappa*=0
## for all pairs, affine dim ~ L-1) and why intrinsic DORM blows up?
## Hypothesis: estimating the ratio on the FULL 200-dim gene space makes any two donors trivially
## separable (batch effects) => everything mutually irreducible => no latent structure. Re-estimate
## the ratio on the 20-dim A (the shared component DORM actually predicts on) and see if reducibility
## structure appears: some pairs with kappa* > 0, and a SMALLER affine effective dimension.
suppressPackageStartupMessages({ library(data.table) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R")); repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")
set.seed(1)
D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, 200); A <- head(D$rna_features, 20)
q <- length(A)
alldon <- sort(unique(D$meta$batch))
rows <- lapply(alldon, function(b) sample_at_most(which(D$meta$batch == b), 400))
Xl_full <- lapply(rows, function(r) make_X(D, r, A))               # [1, A(20), W(180)] = 201 cols
Xl_A    <- lapply(Xl_full, function(X) X[, 1:(q + 1), drop = FALSE])  # [1, A(20)]        =  21 cols

analyze <- function(Xl, tag) {
  R <- fit_source_density_ratios(Xl, dr_type = "logit")$Rmat        # L x M, ratio on this feature set
  L <- nrow(R)
  K <- matrix(0, L, L)
  for (i in 1:L) for (j in 1:L) if (i != j) K[i, j] <- kappa_coef(diag(L)[i, ], diag(L)[j, ], R)
  off <- K[row(K) != col(K)]
  Rc <- sweep(R, 2, colMeans(R)); sv <- svd(Rc)$d; pr <- (sum(sv)^2) / sum(sv^2)
  cat(sprintf("\n== %s (feature dim = %d) ==\n", tag, ncol(Xl[[1]]) - 1))
  cat(sprintf("  pairwise kappa*: min=%.3f mean=%.3f max=%.3f | %% pairs reducible (kappa*>0.05): %.0f%%\n",
              min(off), mean(off), max(off), 100 * mean(off > 0.05)))
  cat(sprintf("  affine effective dim (participation ratio) = %.2f / %d\n", pr, L - 1))
  cat(sprintf("  singular values rel s1: %s\n", paste(sprintf("%.2f", (sv / sv[1])[1:L]), collapse = " ")))
  ## demix at m=3 and m=6: is alpha still just donor-unit-vectors (no reduction), or does it blend?
  for (m in c(3, 6)) {
    rec <- tryCatch(recover_intrinsic_sources(Xl, m, oracle_Rmat = R, verbose = FALSE), error = function(e) NULL)
    if (is.null(rec)) { cat(sprintf("  m=%d: demix error\n", m)); next }
    ## "blend": max, over recovered Q_j, of the 2nd-largest |alpha| in its column (0 => pure donor vertex)
    blend <- max(apply(abs(rec$alpha), 2, function(c) sort(c, decreasing = TRUE)[2]))
    cat(sprintf("  m=%d: max pairwise kappa*=%.3f | vertices={%s} | max 2nd-|alpha| (blend)=%.2f\n",
                m, rec$max_kappa, paste(alldon[rec$sel], collapse = ","), blend))
  }
  invisible(list(K = K, pr = pr))
}

cat("Real donor data, L =", length(alldon), "donors, 400 cells each. Comparing density-ratio feature space.\n")
analyze(Xl_full, "FULL gene space")
analyze(Xl_A,    "A-only (shared component)")
cat("\nInterpretation: if A-only shows MORE reducible pairs and a SMALLER affine dim than the full space,\n")
cat("then the full-space ratio was hiding the latent structure -- demixing/intrinsic DORM should be run in A.\n")
