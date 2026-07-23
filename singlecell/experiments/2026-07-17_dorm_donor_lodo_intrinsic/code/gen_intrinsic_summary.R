#!/usr/bin/env Rscript
## Handover summary for the intrinsic-source / signed-affine machinery (no DORM prediction).
## Part 1: TWO synthetic examples with KNOWN mixing Pi (P_l = sum_j Pi_lj Q_j): (A) L=m=3 (independent
##   sources -> bounded enlarged rho set), (B) L=3,m=2 (dependent -> UNBOUNDED rho set). For each: the
##   explicit Pi, the enlarged rho set (figure), and demixing (alpha, kappa) via ESTIMATED and ORACLE
##   density ratios, compared to the true mixing.
## Part 2: real-data demixing L=12, m=3..12, with example alpha matrices.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(MASS) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R")); repo <- find_repo(code_dir); load_dorm_src(repo)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")
fig_dir <- file.path(exp_dir, "figures"); mdf <- file.path(exp_dir, "INTRINSIC_SUMMARY.md")
set.seed(1)
LN <- c(); add <- function(...) LN <<- c(LN, sprintf(...))

## grid on the plane rho = (u, v, 1-u-v)
g <- seq(-3, 4, length.out = 240); G <- as.matrix(expand.grid(u = g, v = g)); rho_grid <- cbind(G[,1], G[,2], 1-G[,1]-G[,2])
cell <- (g[2]-g[1])^2; areaSimplex <- sum((G[,1]>=0)&(G[,2]>=0)&(G[,1]+G[,2]<=1))*cell

## oracle density ratios rbar_l = P_l/P_pool at points X (Gaussian mixtures, diag cov sig^2 I; const cancels)
oracle_rbar <- function(Xev, Pi, mu, sig) {
  ph <- sapply(seq_len(nrow(mu)), function(j) exp(-rowSums(sweep(Xev, 2, mu[j, ])^2) / (2 * sig^2)))  # M x m
  pl <- ph %*% t(Pi)                                    # M x L (unnormalized P_l)
  ppool <- rowMeans(pl)                                 # equal source weights
  t(pl / ppool)                                         # L x M
}
allperm <- function(m) { P <- as.matrix(expand.grid(rep(list(seq_len(m)), m))); P[apply(P, 1, function(r) length(unique(r)) == m), , drop = FALSE] }
perm_err <- function(rec, tru) { P <- allperm(ncol(tru)); min(apply(P, 1, function(pm) sqrt(sum((rec[, pm, drop=FALSE] - tru)^2)))) }

run_example <- function(Pi, mu, sig, N, label) {
  Ln <- nrow(Pi); m <- ncol(Pi)
  Xl <- lapply(seq_len(Ln), function(l) { comp <- sample(seq_len(m), N, TRUE, Pi[l, ]); cbind(Intercept = 1, mu[comp, ] + matrix(rnorm(N * ncol(mu), sd = sig), N, ncol(mu))) })
  Xpool <- do.call(rbind, lapply(Xl, function(X) X[, -1, drop = FALSE]))
  Rorc <- oracle_rbar(Xpool, Pi, mu, sig)               # oracle L x M
  Rsub <- Rorc[, sort(sample(ncol(Rorc), min(700, ncol(Rorc))))]
  inR <- apply(rho_grid %*% Rsub, 1, min) >= 1e-3
  alpha_true <- t(ginv(Pi))                             # L x m true recovery coefs (Q_j = sum_l alpha_lj P_l)
  rec_o <- recover_intrinsic_sources(Xl, m, oracle_Rmat = Rorc)         # oracle-ratio demix
  rec_e <- recover_intrinsic_sources(Xl, m, dr_type = "logit")          # estimated-ratio demix
  list(label = label, Pi = Pi, m = m, inR = inR, areaR = sum(inR) * cell,
       hit = any(inR & (abs(G[,1]) > 3.9 | abs(G[,2]) > 3.9)), rng = range(rho_grid[inR, ]),
       alpha_true = alpha_true, rec_o = rec_o, rec_e = rec_e,
       err_o = perm_err(rec_o$alpha, alpha_true), err_e = perm_err(rec_e$alpha, alpha_true))
}
mu3 <- 4 * rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0)); Pi3 <- rbind(c(.6,.3,.1), c(.2,.6,.2), c(.1,.3,.6))
mu2 <- 4 * rbind(c(1,0,0,0), c(0,1,0,0));            Pi2 <- rbind(c(.8,.2), c(.5,.5), c(.2,.8))
exA <- run_example(Pi3, mu3, 1.0, 800, "A: L=3, m=3 (independent)")
exB <- run_example(Pi2, mu2, 1.0, 800, "B: L=3, m=2 (dependent)")

## ---- figure: two enlarged rho sets side by side ----
tri <- data.frame(u = c(0,1,0), v = c(0,0,1))
dt <- rbind(data.table(u=G[,1], v=G[,2], inR=exA$inR, cfg=exA$label), data.table(u=G[,1], v=G[,2], inR=exB$inR, cfg=exB$label))
vt <- rbind(data.table(u=exA$alpha_true[1,], v=exA$alpha_true[2,], cfg=exA$label),
            data.table(u=exB$alpha_true[1,], v=exB$alpha_true[2,], cfg=exB$label))
png(file.path(fig_dir, "intrinsic_signedset_L3.png"), 1150, 620, res = 120)
print(ggplot() + geom_raster(data=dt[inR==TRUE], aes(u,v), fill="#8B0000", alpha=0.5) +
  geom_polygon(data=tri, aes(u,v), fill=NA, color="red", linewidth=1.1) +
  geom_point(data=vt, aes(u,v), shape=17, size=3.2, color="black") +
  facet_wrap(~cfg) + coord_equal(xlim=c(-2.5,3.5), ylim=c(-2.5,3.5)) +
  labs(x=expression(rho[1]), y=expression(rho[2]),
       title="Enlarged rho set (dark) vs observed-source simplex (red outline); triangles = true latent Q_j") +
  theme_bw() + theme(strip.text=element_text(size=11), plot.title=element_text(size=11)))
dev.off()

pistr <- function(Pi) paste(sapply(seq_len(nrow(Pi)), function(l) sprintf("P%d = %s", l, paste(sprintf("%.1f Q%d", Pi[l,], seq_len(ncol(Pi))), collapse=" + "))), collapse="; ")
amat <- function(al) paste(apply(round(al,2), 2, function(c) sprintf("(%s)", paste(c, collapse=","))), collapse="  ")

add("# Intrinsic-source (signed-affine) machinery — handover summary\n")
add("Illustrations only; DORM prediction is not involved. eps (validity floor) = 1e-3.\n")
add("## 1. Enlarged rho set + demixing recovery (two synthetic examples with KNOWN mixing)\n")
add("rho = weights over the L=3 observed sources; lives on plane sum(rho)=1. Signed-affine set")
add("R = {sum(rho)=1, induced mixture sum_l rho_l*rbar_l(x) >= eps for all x}. Figure: `figures/intrinsic_signedset_L3.png`.\n")
for (ex in list(exA, exB)) {
  add("### Example %s\n", ex$label)
  add("How observed P are mixed from intrinsic Q (TRUE Pi):  %s", pistr(ex$Pi))
  add("- Enlarged rho set: area(R) = %.2f, **R/simplex = %.1fx**, rho range [%.2f, %.2f], bounded = **%s**.",
      ex$areaR, ex$areaR/areaSimplex, ex$rng[1], ex$rng[2], ifelse(ex$hit, "NO (unbounded)", "yes"))
  add("- Demixing kappa* (max pairwise): oracle-ratios = %.4f, estimated-ratios = %.4f  (~0 => recovered, irreducible).",
      ex$rec_o$max_kappa, ex$rec_e$max_kappa)
  add("- alpha (Q_j = sum_l alpha_lj P_l), columns = Q1..Q%d:", ex$m)
  add("    TRUE  t(pinv(Pi)) : %s", amat(ex$alpha_true))
  add("    ORACLE-recovered : %s   (Frobenius err vs true, best-perm = %.2f)", amat(ex$rec_o$alpha), ex$err_o)
  add("    EST-recovered    : %s   (err = %.2f)", amat(ex$rec_e$alpha), ex$err_e)
  add("")
}
add("**Takeaways.** (A) m=3=L: sources affinely INDEPENDENT => R is a BOUNDED triangle ~%.1fx the simplex,",
    exA$areaR/areaSimplex)
add("with the 3 true latent Q_j at its corners (outside the observed simplex = 'purer than any source');")
add("demixing recovers Pi (oracle err %.2f). (B) m=2 < L=3: sources affinely DEPENDENT => R is UNBOUNDED",
    exA$err_o)
add("(an infinite strip; a null direction over the 3 weights leaves the 2-latent mixture unchanged) — this")
add("non-compactness is the geometric root of the signed-prior instability at large L. Demixing still")
add("recovers the 2 latent (oracle err %.2f); note alpha is only unique up to that null direction, so the", exB$err_o)
add("recovered coefs (which use the two extreme sources) may differ from t(pinv(Pi)) while giving the same Q.\n")

## ---- Part 2: real-data demixing, L=12 (all donors), m=3..12 ----
D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, 200); A <- head(D$rna_features, 20)
alldon <- sort(unique(D$meta$batch))
Xl <- lapply(alldon, function(b) make_X(D, sample_at_most(which(D$meta$batch == b), 400), A)); names(Xl) <- alldon
Rf <- fit_source_density_ratios(Xl, dr_type = "logit")$Rmat
add("## 2. Demixing the real data, L = 12 (all donors), m = 3..12\n")
add("Sources = all %d donors. kappa* ~ 0 => recovered components mutually irreducible.\n", length(alldon))
add("| m | max pairwise kappa* | donors chosen as vertices | max|alpha| |")
add("|---|---|---|---|")
recs <- list()
for (m in 3:12) {
  rec <- tryCatch(recover_intrinsic_sources(Xl, m, oracle_Rmat = Rf, verbose = FALSE), error = function(e) NULL)
  if (is.null(rec)) { add("| %d | ERROR | - | - |", m); next }
  recs[[as.character(m)]] <- rec
  add("| %d | %.4f | %s | %.2f |", m, rec$max_kappa, paste(alldon[rec$sel], collapse = ", "), max(abs(rec$alpha)))
}
add("\n**kappa*=0 is uninformative here** (donors are pairwise mutually irreducible; mean pairwise kappa* 0.016,")
add("affine dim 9.2/10 -- see diag_pairwise_kappa.R), so the recovered 'intrinsic sources' are individual donors,")
add("NOT cell types, and the count does not reduce below ~L.\n")
## ---- Part 3: forward reconstruction (P_l = sum_j Pi_lj Q_j) + affine dimension ----
## Forward weights = affine regression of each source's density ratio onto the m estimated Q_j ratios
## (sum to 1, signs free). NB the purely-algebraic t(ginv(alpha)) is 0 on non-vertex rows (alpha's non-
## vertex rows are 0), so it cannot reconstruct the left-out donors; the ratio-space fit is the honest one.
fwd_weights <- function(y, C) {                          # y: length-M source ratio; C: m x M components
  m <- nrow(C); G <- C %*% t(C); b <- C %*% y
  KKT <- rbind(cbind(G, rep(1, m)), c(rep(1, m), 0)); rhs <- c(b, 1)
  sol <- tryCatch(solve(KKT, rhs), error = function(e) as.vector(ginv(KKT) %*% rhs))
  w <- sol[seq_len(m)]; yhat <- as.vector(w %*% C)
  list(w = w, res = sqrt(sum((y - yhat)^2)) / sqrt(sum(y^2)))
}
fwd_Pi <- function(rec, Rmat) {
  comp <- rec$components                                  # m x M estimated Q_j density ratios
  fw <- lapply(seq_len(nrow(Rmat)), function(l) fwd_weights(Rmat[l, ], comp))
  list(Pi = do.call(rbind, lapply(fw, `[[`, "w")), res = sapply(fw, `[[`, "res"))
}
Lall <- nrow(Rf)
Rc <- sweep(Rf, 2, colMeans(Rf)); sv <- svd(Rc)$d; pr <- (sum(sv)^2) / sum(sv^2)

add("## 3. Forward reconstruction (P_l = sum_j Pi_lj Q_j) and affine dimension\n")
add("### What \"affine effective dimension = 9.2 / 10\" means\n")
add("Stack the L donor density-ratio vectors rbar_l (evaluated on the pooled sample) as rows, center each")
add("column, and take the SVD. The **participation ratio** PR = (sum_k s_k)^2 / sum_k s_k^2 of the singular")
add("values s_k counts how many directions carry appreciable variance: PR ~ 1 if a single direction dominates,")
add("PR ~ d if d directions contribute about equally. It is a soft, continuous rank. Because the observed")
add("sources lie on the affine slice {prior^T r = 1}, the maximum possible is **L-1** (the sum constraint")
add("removes one degree of freedom).")
add("- diag_pairwise_kappa.R (L=11 sources, target s3d7 held out): **PR = 9.2 out of a max of 10**.")
add("- All L=%d donors here: **PR = %.1f out of a max of %d**; singular values rel. to s1: %s.",
    Lall, pr, Lall - 1, paste(sprintf("%.2f", (sv / sv[1])[seq_len(min(Lall, 12))]), collapse = " "))
add("PR ~ L-1 means the donors span **almost the full affine space**: there is essentially no low-dimensional")
add("latent structure to compress into a handful of Q_j. Together with pairwise kappa* ~ 0 (every donor is")
add("mutually irreducible), the genuine intrinsic-source count is ~L -- the Q_j the demixer returns are the")
add("donors themselves, not a smaller set of shared cell-type programs.\n")

add("### Forward-reconstruction fidelity by m\n")
add("For each m we express every observed source P_l as an affine mixture of the m estimated Q_j,")
add("Pi_l = argmin_w || rbar_l - sum_j w_j Qbar_j ||  s.t. sum_j w_j = 1, and report the relative L2 residual.")
add("The m *vertex* donors (the ones the demixer picked as Q_j) reconstruct ~perfectly by construction; the")
add("L-m non-vertex donors are the diagnostic.\n")
add("| m | median resid (non-vertex) | max resid (non-vertex) | # non-vertex recon well (<0.10) | most-negative Pi |")
add("|---|---|---|---|---|")
Pis <- list()
for (m in 3:12) {
  rec <- recs[[as.character(m)]]; if (is.null(rec)) { add("| %d | - | - | - | - |", m); next }
  fp <- fwd_Pi(rec, Rf); Pis[[as.character(m)]] <- fp
  nv <- setdiff(seq_len(Lall), rec$sel)
  if (length(nv) == 0) add("| %d | - | - | 0 / 0 (all donors are vertices) | %.2f |", m, min(fp$Pi))
  else add("| %d | %.3f | %.3f | %d / %d | %.2f |", m,
           median(fp$res[nv]), max(fp$res[nv]), sum(fp$res[nv] < 0.10), length(nv), min(fp$Pi))
}
add("\nResiduals stay high and the Pi weights go **negative** for the non-vertex donors until m approaches L,")
add("confirming that no small set of Q_j reconstructs the left-out donors -- exactly what PR ~ L-1 predicts.\n")

for (m in c(3, 6, 12)) {
  rec <- recs[[as.character(m)]]; if (is.null(rec)) next
  fp <- Pis[[as.character(m)]]; if (is.null(fp)) { fp <- fwd_Pi(rec, Rf); Pis[[as.character(m)]] <- fp }
  al <- round(rec$alpha, 2); Pi <- round(fp$Pi, 2); isv <- seq_len(Lall) %in% rec$sel
  add("### m = %d  (backward alpha: Q_j = sum_l alpha_lj P_l  |  forward Pi: P_l = sum_j Pi_lj Q_j)\n", m)
  add("| donor | vertex? | %s | %s | resid |",
      paste(sprintf("a:Q%d", 1:m), collapse = " | "), paste(sprintf("Pi:Q%d", 1:m), collapse = " | "))
  add("|%s", paste(rep("---|", 2 * m + 3), collapse = ""))
  for (i in seq_len(Lall))
    add("| %s | %s | %s | %s | %.2f |", alldon[i], ifelse(isv[i], "yes", ""),
        paste(sprintf("%.2f", al[i, ]), collapse = " | "),
        paste(sprintf("%.2f", Pi[i, ]), collapse = " | "), fp$res[i])
  add("")
}

add("### Interpretation (m = 3, 6, 12)\n")
add("- **alpha (backward)** is supported only on the m vertex donors: each Q_j ~ one donor (its alpha column is")
add("~ a unit vector) and the other donors contribute 0. The demixer does not blend donors into shared latent")
add("programs; it simply *re-labels* m of them as the 'pure' sources Q_j.")
add("- **Pi (forward)** shows the consequence. A vertex donor reconstructs as Pi_l ~ e_j (all weight on its own")
add("Q_j, resid ~ 0). A non-vertex donor needs a spread, often **signed**, combination of the Q_j and *still*")
add("leaves a sizeable residual -- it carries structure that none of the m chosen vertices span.")
add("- At **m = 12** every donor is its own vertex, so Pi ~ identity and residuals ~ 0: perfect but with zero")
add("compression. There is no intermediate m at which a few Q_j reconstruct all 12 donors. That is the concrete")
add("meaning of affine dim ~ L-1 and pairwise kappa* ~ 0: **on this data the intrinsic sources are the donors,**")
add("**so the signed-affine extension enlarges rho only marginally and buys no real demixing here.**\n")

writeLines(LN, mdf)
cat("Wrote", mdf, "\n")
