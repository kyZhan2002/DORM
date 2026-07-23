#!/usr/bin/env Rscript
## Update INTRINSIC_SUMMARY_formatted.md sections 2-5 with demixing computed in the HARMONY
## embedding space (batch-robust RNA PCA, batch=donor) instead of the full gene space.
## Sections 0-1 (notation + synthetic examples) are preserved verbatim; only the real-data
## numbers change, keeping the polished LaTeX flow. Reuses gen_intrinsic_summary.R's forward-
## reconstruction logic.
suppressPackageStartupMessages({ library(data.table); library(MASS) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R")); repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")
mdf <- file.path(exp_dir, "INTRINSIC_SUMMARY_formatted.md")
set.seed(1)

## ---- data + Harmony embedding ----
D <- load_cite_data(csv_dir)
alldon <- sort(unique(D$meta$batch)); Lall <- length(alldon)
rows <- lapply(alldon, function(b) sample_at_most(which(D$meta$batch == b), 400))
emb <- fread(file.path(csv_dir, "cite_nips_EMBED_for_R.csv"))
idx <- match(D$meta$obs_id, emb$obs_id); stopifnot(!any(is.na(idx)))
Hmat <- as.matrix(emb[idx, paste0("H", 1:50), with = FALSE])
make_Xemb <- function(mat, rws) cbind(Intercept = 1, mat[rws, , drop = FALSE])
Xl <- lapply(rows, function(r) make_Xemb(Hmat, r))
Rf <- fit_source_density_ratios(Xl, dr_type = "logit")$Rmat

## reducibility + affine dim (Harmony)
K <- matrix(0, Lall, Lall)
for (i in 1:Lall) for (j in 1:Lall) if (i != j) K[i, j] <- kappa_coef(diag(Lall)[i, ], diag(Lall)[j, ], Rf)
off <- K[row(K) != col(K)]; redu <- mean(off > 0.05); meank <- mean(off)
Rc <- sweep(Rf, 2, colMeans(Rf)); sv <- svd(Rc)$d; pr12 <- (sum(sv)^2) / sum(sv^2)
keep <- alldon != "s3d7"
Rf11 <- fit_source_density_ratios(Xl[keep], dr_type = "logit")$Rmat
Rc11 <- sweep(Rf11, 2, colMeans(Rf11)); sv11 <- svd(Rc11)$d; pr11 <- (sum(sv11)^2) / sum(sv11^2)

## demix m=3..12
recs <- list()
for (m in 3:12) recs[[as.character(m)]] <- recover_intrinsic_sources(Xl, m, oracle_Rmat = Rf)

## forward reconstruction (same as gen_intrinsic_summary.R)
fwd_weights <- function(y, C) {
  m <- nrow(C); G <- C %*% t(C); b <- C %*% y
  KKT <- rbind(cbind(G, rep(1, m)), c(rep(1, m), 0)); rhs <- c(b, 1)
  sol <- tryCatch(solve(KKT, rhs), error = function(e) as.vector(ginv(KKT) %*% rhs))
  w <- sol[seq_len(m)]; list(w = w, res = sqrt(sum((y - as.vector(w %*% C))^2)) / sqrt(sum(y^2)))
}
fwd_Pi <- function(rec, Rmat) {
  comp <- rec$components; fw <- lapply(seq_len(nrow(Rmat)), function(l) fwd_weights(Rmat[l, ], comp))
  list(Pi = do.call(rbind, lapply(fw, `[[`, "w")), res = sapply(fw, `[[`, "res"))
}
Pis <- lapply(recs, fwd_Pi, Rmat = Rf)

## ---- build sections 2-5 (polished LaTeX-markdown, matching the existing style) ----
S <- c(); p <- function(fmt, ...) S <<- c(S, if (...length() == 0) fmt else sprintf(fmt, ...))
bt <- function(v) paste(sprintf("`%s`", v), collapse = ", ")

p("## 2. Demixing the real data: \\(L=12\\) donors and \\(m=3,\\ldots,12\\)")
p("")
p("All 12 donors are used as observed sources. **The density ratio is now estimated in the Harmony")
p("embedding** (batch-corrected RNA PCA, batch = donor; RNA-only, no ADT leakage) instead of the full")
p("gene space, so shared cell-type structure is no longer masked by technical batch effects.")
p("")
p("| $m$ | $\\max_{j\\neq k}\\kappa^*(Q_j\\mid Q_k)$ | Vertex donors | $\\max_{\\ell,j}\\lvert\\alpha_{\\ell j}\\rvert$ |")
p("| :---: | ---: | --- | ---: |")
for (m in 3:12) { r <- recs[[as.character(m)]]
  p("| $%d$ | %.4f | %s | %.2f |", m, r$max_kappa, bt(alldon[r$sel]), max(abs(r$alpha))) }
p("")
p("Unlike the full-gene-space analysis (where every pair had \\(\\kappa^*\\approx0\\) and 4% of pairs were")
p("reducible), in the embedding **%.0f%% of donor pairs are now mutually reducible** (mean pairwise", 100 * redu)
p("\\(\\kappa^*=%.3f\\)): the sources genuinely share structure. The demixer, however, still selects individual", meank)
p("donors as vertices \\(Q_j\\) — it operates on the 12 donor-mixtures, so it relabels donors rather than")
p("isolating cell-type programs (which would require demixing over cells or clusters).")
p("")
p("---")
p("")

## ---- Section 3 ----
p("## 3. Affine effective dimension and forward reconstruction")
p("")
p("### 3.1 Meaning of “affine effective dimension”")
p("")
p("Stack the \\(L\\) donor density-ratio vectors \\(\\bar r_\\ell\\), evaluated on the pooled sample, as rows of a matrix. Center each column and compute its singular values \\(s_1,s_2,\\ldots\\).")
p("")
p("The participation ratio is")
p("")
p("\\[")
p("\\operatorname{PR}")
p("=")
p("\\frac{\\left(\\sum_k s_k\\right)^2}")
p("{\\sum_k s_k^2}.")
p("\\]")
p("")
p("It acts as a continuous or “soft” rank:")
p("")
p("- \\(\\operatorname{PR}\\approx 1\\) when one direction dominates;")
p("- \\(\\operatorname{PR}\\approx d\\) when approximately \\(d\\) directions contribute comparably.")
p("")
p("Because the observed sources lie on an affine slice of the form \\(\\pi^\\top r=1\\), the maximum possible affine dimension is \\(L-1\\).")
p("")
p("#### Observed values (Harmony embedding)")
p("")
p("- All \\(L=12\\) donors: \\(\\operatorname{PR}=%.1f\\) out of a maximum of 11.", pr12)
p("  The singular values relative to \\(s_1\\) are")
p("  \\[")
p("  %s.", paste(sprintf("%.2f", (sv / sv[1])[seq_len(min(Lall, 12))]), collapse = ",\\ "))
p("  \\]")
p("- With donor `s3d7` held out (\\(L=11\\)): \\(\\operatorname{PR}=%.1f\\) out of a maximum of 10.", pr11)
p("")
p("For comparison, in the **full gene space** these were \\(\\operatorname{PR}=10.2/11\\) and \\(9.2/10\\) — i.e. nearly")
p("full rank. The embedding **roughly halves the affine dimension to \\(\\approx%.0f\\)**: a genuine low-dimensional", pr12)
p("latent structure now exists (consistent with the high reducibility above), whereas in the gene space there")
p("was essentially none.")
p("")

p("### 3.2 Forward-reconstruction criterion")
p("")
p("For each \\(m\\), every observed source is approximated by an affine combination of the estimated intrinsic sources:")
p("")
p("\\[")
p("\\widehat\\Pi_{\\ell\\cdot}")
p("\\in")
p("\\underset{w\\in\\mathbb R^m}{\\arg\\min}")
p("\\left\\|")
p("\\bar r_\\ell-\\sum_{j=1}^m w_j\\bar Q_j")
p("\\right\\|_2")
p("\\quad\\text{subject to}\\quad")
p("\\sum_{j=1}^m w_j=1.")
p("\\]")
p("")
p("The reported diagnostic is the relative \\(\\ell_2\\) residual. The \\(m\\) vertex donors reconstruct almost perfectly by construction; the \\(L-m\\) non-vertex donors are the meaningful diagnostic cases.")
p("")
p("| $m$ | Median residual<br>(non-vertex) | Maximum residual<br>(non-vertex) | Non-vertex donors with<br>residual $<0.10$ | $\\min_{\\ell,j}\\Pi_{\\ell j}$ |")
p("| :---: | ---: | ---: | :---: | ---: |")
for (m in 3:12) { r <- recs[[as.character(m)]]; fp <- Pis[[as.character(m)]]
  nv <- setdiff(seq_len(Lall), r$sel)
  if (length(nv) == 0) p("| $%d$ | — | — | 0 / 0 (all vertices) | %.2f |", m, min(fp$Pi))
  else p("| $%d$ | %.3f | %.3f | %d / %d | %.2f |", m,
         median(fp$res[nv]), max(fp$res[nv]), sum(fp$res[nv] < 0.10), length(nv), min(fp$Pi)) }
p("")
nv6 <- setdiff(seq_len(Lall), recs[["6"]]$sel); med6 <- median(Pis[["6"]]$res[nv6])
p("The non-vertex residuals are **lower than in the full gene space** (there the median stayed \\(0.7–1.0\\);")
p("here it is \\(\\approx%.2f\\) at \\(m=6\\)), reflecting the shared structure the embedding exposes. They are", med6)
p("nonetheless far from zero and the forward weights still turn negative, so a small set of \\(Q_j\\) still does not")
p("reconstruct all donors — the compression is partial, not clean.")
p("")
p("---")
p("")

## ---- Section 4 ----
p("## 4. Detailed backward and forward coefficients")
p("")
p("The backward and forward representations are")
p("")
p("\\[")
p("Q_j=\\sum_{\\ell=1}^L\\alpha_{\\ell j}P_\\ell,")
p("\\qquad")
p("P_\\ell\\approx\\sum_{j=1}^m\\Pi_{\\ell j}Q_j.")
p("\\]")
p("")
p("All coefficients below are computed in the Harmony embedding space. To keep the tables readable, \\(\\alpha\\) and \\(\\Pi\\) are separated for \\(m=6\\) and \\(m=12\\).")
p("")

emit_block <- function(mat, cols, isv, resid, sym) {
  hdr <- paste(sprintf("$%s_{\\ell %d}$", sym, cols), collapse = " | ")
  ali <- paste(rep("---:", length(cols)), collapse = " | ")
  p("| Donor | Vertex | %s%s |", hdr, if (!is.null(resid)) " | Residual" else "")
  p("| --- | :---: | %s%s |", ali, if (!is.null(resid)) " | ---:" else "")
  for (i in seq_len(nrow(mat))) {
    vals <- paste(sprintf("%.2f", mat[i, cols]), collapse = " | ")
    line <- sprintf("| `%s` | %s | %s", alldon[i], ifelse(isv[i], "Yes", "—"), vals)
    if (!is.null(resid)) line <- sprintf("%s | %.2f", line, resid[i])
    p("%s |", line)
  }
}

## m = 3 (combined)
r3 <- recs[["3"]]; fp3 <- Pis[["3"]]; isv3 <- seq_len(Lall) %in% r3$sel
al3 <- round(r3$alpha, 2); Pi3 <- round(fp3$Pi, 2)
p("### 4.1 Results for \\(m=3\\)  (vertices: %s)", bt(alldon[r3$sel]))
p("")
hdr <- paste(c(sprintf("$\\alpha_{\\ell %d}$", 1:3), sprintf("$\\Pi_{\\ell %d}$", 1:3)), collapse = " | ")
p("| Donor | Vertex | %s | Residual |", hdr)
p("| --- | :---: | %s | ---: |", paste(rep("---:", 6), collapse = " | "))
for (i in seq_len(Lall))
  p("| `%s` | %s | %s | %s | %.2f |", alldon[i], ifelse(isv3[i], "Yes", "—"),
    paste(sprintf("%.2f", al3[i, ]), collapse = " | "),
    paste(sprintf("%.2f", Pi3[i, ]), collapse = " | "), fp3$res[i])
p("")

## m = 6 (split alpha / Pi)
r6 <- recs[["6"]]; fp6 <- Pis[["6"]]; isv6 <- seq_len(Lall) %in% r6$sel
al6 <- round(r6$alpha, 2); Pi6 <- round(fp6$Pi, 2)
p("### 4.2 Results for \\(m=6\\)  (vertices: %s)", bt(alldon[r6$sel]))
p("")
p("#### Backward coefficients \\(\\alpha\\)")
p("")
emit_block(al6, 1:6, isv6, NULL, "\\alpha")
p("")
p("#### Forward weights \\(\\Pi\\) and reconstruction residuals")
p("")
emit_block(Pi6, 1:6, isv6, fp6$res, "\\Pi")
p("")

## m = 12 (split into two column blocks each)
r12 <- recs[["12"]]; fp12 <- Pis[["12"]]; isv12 <- seq_len(Lall) %in% r12$sel
al12 <- round(r12$alpha, 2); Pi12 <- round(fp12$Pi, 2)
p("### 4.3 Results for \\(m=12\\)")
p("")
p("Because the \\(m=12\\) matrices are wide, each is divided into two column blocks.")
p("")
p("#### Backward coefficients \\(\\alpha\\): columns 1–6")
p("")
emit_block(al12, 1:6, isv12, NULL, "\\alpha")
p("")
p("#### Backward coefficients \\(\\alpha\\): columns 7–12")
p("")
emit_block(al12, 7:12, isv12, NULL, "\\alpha")
p("")
p("#### Forward weights \\(\\Pi\\): columns 1–6")
p("")
emit_block(Pi12, 1:6, isv12, NULL, "\\Pi")
p("")
p("#### Forward weights \\(\\Pi\\): columns 7–12 and residuals")
p("")
emit_block(Pi12, 7:12, isv12, fp12$res, "\\Pi")
p("")
p("---")
p("")

## ---- Section 5 ----
maxalpha_all <- max(sapply(recs, function(r) max(abs(r$alpha))))
p("## 5. Interpretation")
p("")
p("- **Batch-robust embedding restores shared structure.** Moving the density ratio into the Harmony space")
p("  drops the (soft) affine dimension from \\(\\approx L-1\\) (full genes) to \\(\\approx%.0f\\) and makes %.0f%% of donor", pr12, 100 * redu)
p("  pairs reducible (up from 4%). A genuine low-dimensional latent structure now exists.")
p("- **Forward reconstruction improves sharply.** Non-vertex donors, essentially unreconstructable in the gene")
p("  space (residual \\(0.7–1.0\\)), now reconstruct with median residual \\(\\approx%.2f\\) at \\(m=6\\): about six", med6)
p("  intrinsic sources already capture the donors well.")
p("- **But the demixer still returns donor-vertices, not cell types.** \\(\\alpha\\) stays supported on the selected")
p("  donors, and its entries are now large (up to \\(\\approx%.0f\\)) because several embedding directions are weak —", maxalpha_all)
p("  mild ill-conditioning, not clean cell-type separation (which would require demixing over cells/clusters).")
p("")
p("\\[")
p("\\boxed{")
p("\\begin{array}{c}")
p("\\text{In the batch-robust embedding the donor sources share a genuine low-dimensional structure}\\\\")
p("\\text{(soft affine dim }\\approx%.0f\\text{; about six intrinsic sources reconstruct the donors), unlike the gene space.}\\\\", pr12)
p("\\text{The demixer still selects donor-vertices rather than cell-type programs.}\\\\")
p("\\text{Whether this improves the DORM-intrinsic predictor is a separate, not-yet-run question (Phase 4).}")
p("\\end{array}")
p("}")
p("\\]")

## ---- splice: keep everything before "## 2." verbatim ----
cur <- readLines(mdf)
cut <- grep("^## 2\\. Demixing the real data", cur)[1]
stopifnot(!is.na(cut))
writeLines(c(cur[1:(cut - 1)], S), mdf)
cat("Updated", mdf, "\n")
cat(sprintf("Harmony: reducible=%.0f%%  affineDim(L12)=%.2f  affineDim(L11)=%.2f  medNVresid@m6=%.2f\n",
            100 * redu, pr12, pr11, med6))
cat("m=3 vertices:", paste(alldon[recs[["3"]]$sel], collapse = ","),
    "| m=6 vertices:", paste(alldon[recs[["6"]]$sel], collapse = ","), "\n")
