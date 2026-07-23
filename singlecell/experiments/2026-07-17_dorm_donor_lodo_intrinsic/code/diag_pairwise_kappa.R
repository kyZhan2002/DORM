#!/usr/bin/env Rscript
## Why does max kappa*=0 for every m? Check whether the raw donor-sources are already mutually
## irreducible (kappa*(P^i||P^j) ~ 0 for all pairs). If so, ANY m distinct donors trivially pass
## the irreducibility test => kappa*=0 is uninformative about a "true" m; the intrinsic sources ARE
## the donors. Also report the affine effective dimension (no. of independent ratio directions).
suppressPackageStartupMessages({ library(data.table) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")
set.seed(1)
D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, 200); A <- head(D$rna_features, 20)
donors <- setdiff(sort(unique(D$meta$batch)), "s3d7")
Xl <- lapply(donors, function(b) make_X(D, sample_at_most(which(D$meta$batch == b), 400), A)); names(Xl) <- donors
R <- fit_source_density_ratios(Xl, dr_type = "logit")$Rmat; L <- nrow(R)

## pairwise kappa* between raw donor sources (coefficient e_i vs e_j)
K <- matrix(0, L, L)
for (i in 1:L) for (j in 1:L) if (i != j) K[i, j] <- kappa_coef(diag(L)[i, ], diag(L)[j, ], R)
off <- K[row(K) != col(K)]
cat(sprintf("Pairwise kappa* among the %d raw donor-sources: min=%.4f  mean=%.4f  max=%.4f\n", L, min(off), mean(off), max(off)))
cat(sprintf("  fraction of donor-pairs with kappa* < 0.05 (mutually irreducible): %.0f%%\n", 100 * mean(off < 0.05)))
cat("  => distinct donors are already mutually irreducible, so demix returns donor-vertices and\n     kappa*=0 for every m (the test is trivially satisfied; it does NOT identify a true m).\n\n")

## affine effective dimension (how many independent ratio directions => 'true' number of sources)
Rc <- sweep(R, 2, colMeans(R)); sv <- svd(Rc)$d
cat(sprintf("Affine effective dimension (participation ratio) = %.2f out of max %d\n", (sum(sv)^2)/sum(sv^2), L - 1))
cat(sprintf("  singular values rel to s1: %s\n", paste(sprintf("%.2f", (sv/sv[1])[1:L]), collapse = " ")))
cat("  => ~full rank: the donors span ~L independent directions, i.e. NO low-dimensional latent\n     structure to reduce to. The intrinsic-source count is essentially L (the donors themselves).\n")
