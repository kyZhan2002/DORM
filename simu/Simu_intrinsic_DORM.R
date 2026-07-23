####################################################################################
##  Simu_intrinsic_DORM.R
##
##  End-to-end demonstration of the intrinsic-source (signed-affine valid-mixture)
##  prior in DORM. NOTE: this is the projection-primary method -- it needs NO
##  intrinsic-source recovery, no irreducibility/separability, and no choice of m;
##  it just enables signed_prior = TRUE, which projects the target onto the
##  identifiable set of valid signed-affine mixtures (see signed_affine_prior).
##
##  Setting: L = 4 OBSERVED sources; the target is a mixture whose implied weights over
##  the 4 observed sources are SIGNED and lie OUTSIDE the observed-source simplex
##  Delta_L -- so standard DORM (prior confined to Delta_L) cannot represent the
##  target's mixture, while the signed-affine prior can. (The target is generated as a
##  mixture of m = 2 latent sources purely to place it outside Delta_L; the method does
##  not use or assume that latent structure.)
##
##  We compare:
##    - target mixture-weight recovery: || rho_hat - rho_bar_target_true ||  (HEADLINE)
##    - target prediction: MSE and EXCESS MSE over the best linear predictor
##      (target OLS floor), which isolates the part the mixture weights control.
##
##  Honest summary of what to expect (see comments at the bottom): the signed-affine
##  prior recovers the true signed mixture weights that standard DORM structurally
##  cannot, and gives a consistent (if modest) reduction in excess prediction error. The
##  outcome-MSE gain is bounded because source conditional means do not combine exactly
##  linearly across covariate mixtures; the clean, robust win is the mixture-structure
##  (weight) recovery. Averaged over several seeds for stability.
##
##  Run:  Rscript simu/Simu_intrinsic_DORM.R          (a few minutes; fits ML nuisances)
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
suppressMessages({
  source(file.path(ROOT, "src", "Functions3.R"))
  source(file.path(ROOT, "src", "IntrinsicSources.R"))
  source(file.path(ROOT, "src", "DataGen_intrinsic.R"))
})

## ---------------------------- configuration -------------------------------------
L = 4; m = 2; p = 3; q = p                 # low-dim A = full X (no separate W); penalty = FALSE
smax = 0.1                                 # small: intrinsic joint-mixing ~holds on target
n_seeds = as.integer(Sys.getenv("NSEEDS", 5))
sep = 4; sig = 1.3
Pi = rbind(c(0.85, 0.15),                  # 4 x 2 mixing, full column rank, diverse rows
           c(0.65, 0.35),
           c(0.40, 0.60),
           c(0.15, 0.85))
MU = sep * rbind(c(1, 0, 0), c(-1, 0.3, 0.2))          # 2 intrinsic component means in R^3
GAMMA = cbind(0, 3 * rbind(c(1, 1, -1), c(-1, 1, 1)))  # distinct component outcome models
rhoQ_target = c(0.8, 0.2)                  # target: mildly extrapolating latent mixture
Nlist = rep(2000, L); nlist = rep(1000, L); n0 = 2000

target_ols = function(d, q) {              # best achievable linear predictor on target
  A = d$X0[, 1:(q + 1)]
  as.vector(solve(t(A) %*% A + 1e-6 * diag(q + 1), t(A) %*% d$Y0))
}

## ---------------------------- one replication -----------------------------------
one_rep = function(seed) {
  d = generate_intrinsic_data(L, m, p, Nlist, nlist, n0, Pi, MU, GAMMA,
                              rhoQ_target = rhoQ_target, sigma = sig, seed = seed)

  r_std = get_DORM_beta(d$Xlist, d$Xtrainlist, d$Ylist, d$Ytrainlist, d$X0, d$X0train,
                        nlist, q, smax, dr_type = "logit")
  r_int = get_DORM_beta(d$Xlist, d$Xtrainlist, d$Ylist, d$Ytrainlist, d$X0, d$X0train,
                        nlist, q, smax, dr_type = "logit", signed_prior = TRUE)

  floor = get_err(d$X0, d$Y0, target_ols(d, q), q, n0)
  list(
    mse_std = get_err(d$X0, d$Y0, r_std$beta_star, q, n0),
    mse_int = get_err(d$X0, d$Y0, r_int$beta_star, q, n0),
    floor   = floor,
    wdist_std = sqrt(sum((r_std$rho - d$rho_bar_target)^2)),   # simplex weights vs true signed weights
    wdist_int = sqrt(sum((r_int$rho - d$rho_bar_target)^2)),
    rho_std = r_std$rho, rho_int = r_int$rho, rho_true = d$rho_bar_target
  )
}

## ---------------------------- run replications ----------------------------------
cat(sprintf("Intrinsic-DORM demo | L=%d observed, m=%d latent, p=%d, smax=%.2f, %d seeds\n",
            L, m, p, smax, n_seeds))
cat(sprintf("True target latent mix rhoQ = (%s); implied SIGNED weights over the %d observed sources:\n",
            paste(rhoQ_target, collapse = ", "), L))
cat("  rho_bar_target =", paste(round(as.vector(t(MASS::ginv(Pi)) %*% rhoQ_target), 3), collapse = ", "),
    " (outside Delta_L: has negative entries)\n\n")

res = lapply(seq_len(n_seeds), function(s) {
  cat(sprintf("  seed %d/%d ...\n", s, n_seeds)); one_rep(s)
})

grab = function(f) sapply(res, function(r) r[[f]])
tab = data.frame(seed = seq_len(n_seeds),
                 MSE_std   = round(grab("mse_std"), 3),
                 MSE_int   = round(grab("mse_int"), 3),
                 floor     = round(grab("floor"), 3),
                 excess_std = round(grab("mse_std") - grab("floor"), 3),
                 excess_int = round(grab("mse_int") - grab("floor"), 3),
                 wdist_std = round(grab("wdist_std"), 3),
                 wdist_int = round(grab("wdist_int"), 3))

cat("\n================= PER-SEED RESULTS =================\n")
print(tab, row.names = FALSE)

ex_std = mean(grab("mse_std") - grab("floor"))
ex_int = mean(grab("mse_int") - grab("floor"))
cat("\n================= SUMMARY (means) =================\n")
cat(sprintf("Target MSE:            standard = %.3f   intrinsic = %.3f   (OLS floor = %.3f)\n",
            mean(grab("mse_std")), mean(grab("mse_int")), mean(grab("floor"))))
cat(sprintf("Excess MSE over floor: standard = %.3f   intrinsic = %.3f   -> %.0f%% reduction\n",
            ex_std, ex_int, 100 * (ex_std - ex_int) / ex_std))
cat(sprintf("Dist. to TRUE signed target weights: standard = %.3f   intrinsic = %.3f   (%.1fx closer)\n",
            mean(grab("wdist_std")), mean(grab("wdist_int")),
            mean(grab("wdist_std")) / mean(grab("wdist_int"))))

cat("\nExample (seed 1) target weights over the observed sources:\n")
cat("  standard DORM (in Delta_L):", paste(round(res[[1]]$rho_std, 3), collapse = ", "), "\n")
cat("  intrinsic DORM (signed)   :", paste(round(res[[1]]$rho_int, 3), collapse = ", "), "\n")
cat("  TRUE rho_bar_target       :", paste(round(res[[1]]$rho_true, 3), collapse = ", "), "\n")

pass = (ex_int <= ex_std) && (mean(grab("wdist_int")) < 0.5 * mean(grab("wdist_std")))
cat(sprintf("\nPASS (intrinsic recovers signed weights & does not worsen excess MSE): %s\n", pass))

## ---------------------------------------------------------------------------------
## Interpretation:
##  * Standard DORM's prior is confined to the simplex Delta_L, so it collapses the
##    target onto the nearest convex combination of observed sources -- far from the
##    true SIGNED weights. Intrinsic-DORM recovers those signed weights closely: this
##    is the paper's "more precisely capture the mixture structure" (robust, clean).
##  * The prediction-MSE gain is consistent but modest: source conditional means do
##    not combine exactly linearly across covariate mixtures, so signed extrapolation
##    cannot fully reach a latent component's outcome model. Strong extrapolation
##    (target ~ a pure latent vertex) is harder still and can erode the MSE gain; the
##    mild-mixture regime here is where the benefit is stable.
## ---------------------------------------------------------------------------------
