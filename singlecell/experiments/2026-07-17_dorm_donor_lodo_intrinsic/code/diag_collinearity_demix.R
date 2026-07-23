#!/usr/bin/env Rscript
## DIAGNOSTIC (no DORM wiring): are the 11 donor-source density ratios affinely collinear,
## and what does Katz-Samuels demixing recover on this real data?
suppressPackageStartupMessages({ library(data.table) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

set.seed(1)
D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, 200)
target <- "s3d7"                                   # representative fold; sources = other 11 donors
src_donors <- setdiff(sort(unique(D$meta$batch)), target)
A <- head(D$rna_features, 20)                       # A/W split irrelevant for the covariate DR
Xlist <- lapply(src_donors, function(b) make_X(D, sample_at_most(which(D$meta$batch == b), 400), A))
names(Xlist) <- src_donors
L <- length(Xlist)
cat(sprintf("Sources = %d donors (target held out = %s); ~400 cells each, 200 HVG.\n\n", L, target))

## ---- density ratios vs pooled reference: R is L x M ----
fit <- fit_source_density_ratios(Xlist, dr_type = "logit")
R <- fit$Rmat
cat(sprintf("Density-ratio matrix R: %d sources x %d eval points. range=[%.3g, %.3g]\n", nrow(R), ncol(R), min(R), max(R)))

## (1) COLLINEARITY -------------------------------------------------------------
## pairwise correlation between the L source ratio vectors
cc <- cor(t(R))
cat(sprintf("\n[Collinearity] pairwise corr of source ratio vectors: mean=%.3f  min=%.3f  max(off-diag)=%.3f\n",
            mean(cc[upper.tri(cc)]), min(cc[upper.tri(cc)]), max(cc[upper.tri(cc)])))
## affine structure: center each source by the mean source, then SVD.
## affinely independent L sources => rank L-1; fast singular-value decay => near-dependent.
Rc <- sweep(R, 2, colMeans(R))                      # subtract mean source (row of length M)
sv <- svd(Rc)$d
sv_rel <- sv / sv[1]
cumvar <- cumsum(sv^2) / sum(sv^2)
cat(sprintf("[Affine SVD] singular values (rel to s1), first 11:\n  %s\n", paste(sprintf("%.3g", sv_rel[1:L]), collapse = " ")))
cat(sprintf("[Affine SVD] cumulative variance by top-k dims: top1=%.3f top2=%.3f top3=%.3f top4=%.3f\n",
            cumvar[1], cumvar[2], cumvar[3], cumvar[4]))
eff_rank <- (sum(sv)^2) / sum(sv^2)
cat(sprintf("[Affine SVD] participation ratio (effective affine dim) = %.2f  (max possible = %d)\n", eff_rank, L - 1))
cat(sprintf("[Affine SVD] condition number s1/s(L-1) = %.3g   (huge => near affinely-DEPENDENT => signed set ~unbounded)\n",
            sv[1] / sv[L - 1]))

## (2) DEMIXING -----------------------------------------------------------------
cat("\n[Demixing] Katz-Samuels recover_intrinsic_sources on the estimated ratios:\n")
cat("  max_kappa ~ 0 => recovered components mutually irreducible (clean latent sources);\n")
cat("  max_kappa -> ~1 => NOT separable (no clean intrinsic structure to demix).\n\n")
for (m in c(2, 3, 4, 6)) {
  rec <- tryCatch(recover_intrinsic_sources(Xlist, m, oracle_Rmat = R, verbose = FALSE),
                  error = function(e) { cat("  m=", m, " ERROR: ", conditionMessage(e), "\n"); NULL })
  if (is.null(rec)) next
  al <- round(rec$alpha, 2); rownames(al) <- src_donors
  cat(sprintf("--- m = %d --- max pairwise kappa* = %.4f ; demixed on donors {%s}\n",
              m, rec$max_kappa, paste(src_donors[rec$sel], collapse = ", ")))
  cat("  alpha (L x m; Q_j = sum_l alpha_lj P^(l)) column ranges:\n")
  cat(sprintf("   per-component alpha range: %s\n",
              paste(sapply(1:m, function(j) sprintf("Q%d:[%.2f,%.2f]", j, min(al[,j]), max(al[,j]))), collapse = "  ")))
  cat(sprintf("   max |alpha| = %.2f  (large => extreme signed combos => unstable recovery)\n\n", max(abs(al))))
}
cat("Done.\n")
