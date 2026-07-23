####################################################################################
##  Simu_intrinsic_recovery.R
##
##  Validation of the intrinsic-source recovery code (src/IntrinsicSources.R) via
##  the residue-demixing scheme. Two-stage design:
##    Stage 1 (oracle):    analytic density ratios from the known Gaussians are fed to
##                         the demixer -> isolates the ALGORITHM (should be near-exact).
##    Stage 2 (estimated): the full finite-sample pipeline with ML-estimated ratios
##                         (multinomial elastic-net logit, and random forest).
##
##  Checks, per stage, after aligning recovered components to the truth by permutation:
##    (a) alpha  vs true alpha* = t(pinv(Pi))            [square case L == m only]
##    (b) recovered component means vs true means        [headline, end-to-end]
##    (c) max pairwise kappa at convergence (-> 0)
##    (d) recovered-vs-true intrinsic density-ratio correlation
##    (e) NEGATIVE CONTROL: recovered means are closer to truth than any raw source mean
##
##  Run from anywhere:  Rscript simu/Simu_intrinsic_recovery.R
####################################################################################

find_repo_root = function() {
  d = normalizePath(getwd())
  for (i in 1:8) {
    if (file.exists(file.path(d, "src", "Functions3.R"))) return(d)
    parent = dirname(d); if (parent == d) break; d = parent
  }
  stop("Could not locate repo root (expected src/Functions3.R above the working dir).")
}
ROOT = find_repo_root()
source(file.path(ROOT, "src", "IntrinsicSources.R"))
source(file.path(ROOT, "src", "DataGen_intrinsic.R"))

set.seed(2026)

## --------------------------------------------------------------------------------
## Helpers
## --------------------------------------------------------------------------------

## Best permutation matching recovered component rows (of Mhat) to true rows (MU).
## Exhaustive over m! permutations (m is small); minimizes total matched distance.
best_perm = function(Mhat, MU) {
  m = nrow(MU)
  perms = if (m == 1) matrix(1) else do.call(rbind, .permn(seq_len(m)))
  costs = apply(perms, 1, function(pp) sum(rowSums((Mhat[pp, , drop = FALSE] - MU)^2)))
  perms[which.min(costs), ]
}
.permn = function(x) {                          # all permutations of x, as a list of vectors
  if (length(x) == 1) return(list(x))
  out = list()
  for (i in seq_along(x))
    for (p in .permn(x[-i])) out[[length(out) + 1]] = c(x[i], p)
  out
}

## True intrinsic-component density ratios q_j/pool on eval points (for check (d)).
oracle_component_ratios = function(eval_X, MU, sigma, Pi, Nlist, has_intercept = TRUE) {
  X = as.matrix(eval_X); if (has_intercept) X = X[, -1, drop = FALSE]
  m = nrow(MU)
  logphi = sapply(seq_len(m), function(j) -rowSums(sweep(X, 2, MU[j, ], "-")^2) / (2 * sigma^2))
  phi = exp(logphi - apply(logphi, 1, max))
  poolw = as.vector(t(Pi) %*% (Nlist / sum(Nlist)))   # pool = sum_k poolw_k Q_k
  pool = as.vector(phi %*% poolw)
  t(phi / pool)                                        # m x M, row j = q_j / pool
}

## Run recovery once (given a ratio matrix) and score it against the truth.
evaluate = function(rec, d, src_means, MU, Pi, L, m, comp_true) {
  mu_hat = t(rec$alpha) %*% src_means                  # m x p recovered component means
  perm = best_perm(mu_hat, MU)
  mu_hat = mu_hat[perm, , drop = FALSE]
  alpha_aligned = rec$alpha[, perm, drop = FALSE]
  comp_hat = rec$components[perm, , drop = FALSE]

  mean_err = mean(sqrt(rowSums((mu_hat - MU)^2)))       # (b) avg L2 mean error
  alpha_err = if (L == m) sqrt(sum((alpha_aligned - d$alpha_star)^2)) else NA  # (a)
  ratio_cor = mean(sapply(seq_len(m), function(j)       # (d) per-component ratio corr
    suppressWarnings(cor(comp_hat[j, ], comp_true[j, ]))))
  raw_err = mean(sapply(seq_len(m), function(j)         # (e) nearest raw source mean
    min(sqrt(rowSums(sweep(src_means, 2, MU[j, ], "-")^2)))))

  list(mean_err = mean_err, alpha_err = alpha_err, max_kappa = rec$max_kappa,
       ratio_cor = ratio_cor, raw_err = raw_err,
       sources = paste(rec$sel, collapse = "/"),
       mu_hat = mu_hat, perm = perm)
}

## --------------------------------------------------------------------------------
## Run one full setting (both stages, all requested ratio backends).
## --------------------------------------------------------------------------------
run_setting = function(name, L, m, p, MU, Pi, GAMMA, Nlist, nlist, n0, sigma,
                       seed, dr_types = c("logit", "rf")) {
  cat(sprintf("\n=============== Setting: %s  (L=%d, m=%d, p=%d) ===============\n",
              name, L, m, p))
  d = generate_intrinsic_data(L, m, p, Nlist, nlist, n0, Pi, MU, GAMMA,
                              rhoQ_target = rep(1 / m, m), sigma = sigma, seed = seed)
  poolX = do.call(rbind, d$Xlist)
  src_means = t(sapply(d$Xlist, function(x) colMeans(x[, -1, drop = FALSE])))  # L x p
  comp_true = oracle_component_ratios(poolX, MU, sigma, Pi, Nlist)             # m x M

  rows = list()

  ## Stage 1: oracle ratios
  Rmat_oracle = oracle_source_ratios(poolX, Pi, MU, sigma, Nlist)
  rec = recover_intrinsic_sources(d$Xlist, m, oracle_Rmat = Rmat_oracle)
  rows[["oracle"]] = evaluate(rec, d, src_means, MU, Pi, L, m, comp_true)

  ## Stage 2: estimated ratios
  recs_est = list(oracle = rec)
  for (dt in dr_types) {
    rec_dt = recover_intrinsic_sources(d$Xlist, m, dr_type = dt, eval_X = poolX)
    rows[[dt]] = evaluate(rec_dt, d, src_means, MU, Pi, L, m, comp_true)
    recs_est[[dt]] = rec_dt
  }

  ## Report
  cat(sprintf("%-10s %10s %10s %11s %11s %10s\n",
              "ratios", "mean_err", "alpha_err", "max_kappa*", "ratio_cor", "sources"))
  for (nm in names(rows)) {
    r = rows[[nm]]
    cat(sprintf("%-10s %10.4f %10s %11.5f %11.4f %10s\n", nm, r$mean_err,
                ifelse(is.na(r$alpha_err), "  -  ", sprintf("%.4f", r$alpha_err)),
                r$max_kappa, r$ratio_cor, r$sources))
  }
  raw = rows[["oracle"]]$raw_err
  cat(sprintf("negative control: nearest RAW source mean err = %.4f  (recovery should be << this)\n", raw))
  ## Gate on (i) the ALGORITHM being exact under oracle ratios and (ii) recovery with
  ## a reasonable density-ratio model (logit) clearly beating the raw baseline. The rf
  ## backend is reported as a stress test: rf multiclass probabilities are poorly
  ## calibrated here (low ratio_cor), so the kappa* estimates that drive Demix degrade --
  ## a density-ratio ESTIMATION limitation, not a demixing-algorithm failure.
  oracle_exact = rows[["oracle"]]$mean_err < 0.25 * raw && rows[["oracle"]]$max_kappa < 1e-2
  logit_beats  = "logit" %in% names(rows) && rows[["logit"]]$mean_err < raw
  pass = oracle_exact && logit_beats
  cat(sprintf("PASS (oracle exact & logit beats raw): %s   [oracle_exact=%s, logit_beats_raw=%s]\n",
              pass, oracle_exact, logit_beats))

  invisible(list(data = d, rows = rows, recs = recs_est, src_means = src_means,
                 poolX = poolX, raw_err = raw, pass = pass))
}

## --------------------------------------------------------------------------------
## Settings
## --------------------------------------------------------------------------------
## Setting A: square, p=2 (for the plot). Three well-separated Gaussian components.
MU_A = rbind(c(0, 5), c(5, -3), c(-5, -3))
GAMMA_A = cbind(0, matrix(c(1, 0,  0, 1,  -1, 1), nrow = 3, byrow = TRUE))
resA = run_setting("A: square L=m=3, p=2",
                   L = 3, m = 3, p = 2, MU = MU_A,
                   Pi = make_diag_dominant_Pi(3, 3, diag = 0.7),
                   GAMMA = GAMMA_A, Nlist = rep(2000, 3), nlist = rep(500, 3),
                   n0 = 2000, sigma = 1, seed = 11)

## Setting B: overcomplete L=4 > m=3 (affine reps non-unique; check means/ratios, not alpha).
MU_B = rbind(c(0, 6), c(6, -3), c(-6, -3))
GAMMA_B = cbind(0, matrix(c(1, 0,  0, 1,  -1, 1), nrow = 3, byrow = TRUE))
Pi_B = rbind(c(0.70, 0.15, 0.15),
             c(0.15, 0.70, 0.15),
             c(0.15, 0.15, 0.70),
             c(0.40, 0.40, 0.20))                  # 4th source: a distinct extra mixture
resB = run_setting("B: overcomplete L=4, m=3, p=2", L = 4, m = 3, p = 2, MU = MU_B,
                   Pi = Pi_B, GAMMA = GAMMA_B, Nlist = rep(2000, 4),
                   nlist = rep(500, 4), n0 = 2000, sigma = 1, seed = 22)

## Setting C: higher dimension p=8, square L=m=4.
set.seed(7)
MU_C = matrix(0, nrow = 4, ncol = 8)
for (j in 1:4) MU_C[j, (2 * j - 1):(2 * j)] = c(5, 5)   # each component distinct in its own 2 dims
GAMMA_C = cbind(0, matrix(rnorm(4 * 8), nrow = 4))
resC = run_setting("C: square L=m=4, p=8", L = 4, m = 4, p = 8, MU = MU_C,
                   Pi = make_diag_dominant_Pi(4, 4, diag = 0.6),
                   GAMMA = GAMMA_C, Nlist = rep(3000, 4), nlist = rep(600, 4),
                   n0 = 3000, sigma = 1, seed = 33)

## --------------------------------------------------------------------------------
## Plot for Setting A (p = 2): pooled sources + true vs recovered component means.
## --------------------------------------------------------------------------------
plot_path = file.path(ROOT, "simu", "Simu_intrinsic_recovery.png")   # .png is gitignored
tryCatch({
  png(plot_path, width = 1100, height = 520, res = 130)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  d = resA$data; L = 3; m = 3
  cols = c("#1b9e77", "#d95f02", "#7570b3")
  for (stage in c("oracle", "logit")) {
    plot(NA, xlim = range(resA$poolX[, 2]), ylim = range(resA$poolX[, 3]),
         xlab = "x1", ylab = "x2", main = paste0("Setting A: recovery (", stage, ")"))
    for (l in 1:L) {
      X = d$Xlist[[l]]
      points(X[, 2], X[, 3], col = adjustcolor(cols[((l - 1) %% m) + 1], 0.15), pch = 16, cex = 0.4)
    }
    points(d$MU[, 1], d$MU[, 2], pch = 4, cex = 2.2, lwd = 3, col = "black")   # true means (X)
    muh = resA$rows[[stage]]$mu_hat
    points(muh[, 1], muh[, 2], pch = 1, cex = 2.6, lwd = 3, col = "red")       # recovered (O)
    legend("topright", legend = c("true mean", "recovered"), pch = c(4, 1),
           col = c("black", "red"), pt.lwd = 3, bty = "n", cex = 0.8)
  }
  dev.off()
  cat(sprintf("\nPlot written to: %s\n", plot_path))
}, error = function(e) cat("\n(plot skipped:", conditionMessage(e), ")\n"))

cat("\n================= SUMMARY =================\n")
cat(sprintf("Setting A pass: %s\nSetting B pass: %s\nSetting C pass: %s\n",
            resA$pass, resB$pass, resC$pass))
