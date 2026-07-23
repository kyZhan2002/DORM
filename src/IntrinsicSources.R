####################################################################################
##  IntrinsicSources.R
##
##  Recovery of intrinsic (latent) sources for DORM (mutual-contamination / demixing
##  framework of Katz-Samuels et al. 2019), following the appendix "Recovery of
##  intrinsic sources" of the DORM paper.
##
##  Model: each observed source P^(l) (l = 1..L) is a mixture of m <= L latent,
##  irreducible intrinsic sources {Q_j}:
##        P^(l) = sum_j pi_{lj} Q_j,   pi_{l.} in simplex,  Pi full column rank.
##  Because Pi has full column rank, each Q_j is an AFFINE (possibly signed)
##  combination of the observed sources:
##        Q_j = sum_l alpha_{lj} P^(l),   sum_l alpha_{lj} = 1.
##
##  We work with per-source density ratios r_l(x) = dP^(l)/dP^(ref) against a common
##  reference (the pooled sources), never densities. On the pooled reference every
##  ratio vector lies on the affine slice {r : prior^T r = 1} and decomposes as
##        r_l(x) = sum_j (Pi_{lj} / poolw_j) beta_j(x),   beta(x) in the m-simplex,
##  i.e. the source ratio vectors are convex combinations of m VERTICES
##  v_j = Pi[,j]/poolw_j. Recovering {Q_j} is therefore vertex/anchor identification.
##
##  ALGORITHM. The appendix's literal pairwise kappa*-residue sweep has a spurious
##  fixed point for m >= 3 (once a candidate sheds one component, kappa* against any
##  candidate still containing it is 0, so max-kappa can hit 0 at a wrong multi-
##  component mixture). We instead use the robust, provably-correct-under-separability
##  realization of the same "remove the explainable part until irreducible cores
##  remain" idea: the Successive Projection Algorithm (SPA) on the source ratio
##  vectors. SPA's deflation IS a (Euclidean) residue step; it finds the m anchor
##  points -> vertices v_j -> Pi (row-stochastic scaling) -> alpha = t(pinv(Pi)).
##  A post-hoc kappa* check confirms the recovered components are mutually irreducible.
##
##  Because alpha is REFERENCE-INVARIANT, recovery (pooled-source reference) and the
##  downstream DORM prior/posterior (target reference) may estimate their density
##  ratios independently; only the L x m matrix alpha is passed between them.
##
##  DORM path (projection-primary; needs no recovery, no irreducibility, no m):
##    signed_affine_prior()       - the prior DORM uses: KL-projection of the target
##                                  onto the identifiable valid signed-affine mixture set
##
##  Optional intrinsic-source recovery (Katz-Samuels demixing; for interpretation):
##    fit_source_density_ratios() - estimate forward ratios r_l = dP^(l)/dP^(pool)
##    kappa_hat()                 - finite-sample plug-in of kappa*(F || H)
##    residue_coef() / kappa_coef() - two-sample residue and kappa* in coefficient space
##    face_test_coef() / find_face_coef() - FaceTestHat / FindFaceHat (Katz Alg. 12-13)
##    demix_hat()                 - DemixHat (Katz Alg. 11): recursive vertex recovery
##    recover_intrinsic_sources() - orchestrator (estimated OR oracle ratios)
##
##  `signed_affine_prior()` calls `solve_robust()` from src/Functions3.R at run time;
##  source Functions3.R alongside this file before using the DORM hook.
####################################################################################

library(CVXR)
library(glmnet)
library(randomForest)
library(MASS)

## --------------------------------------------------------------------------------
## Density-ratio estimation against the pooled-source reference.
##
## We use a single multiclass probabilistic classifier P(source = l | x). With pool
## prior P(l) = N_l / sum_k N_k, Bayes' rule gives the forward ratio directly:
##        r_l(x) = dP^(l)/dP^(pool)(x) = P(l | x) / P(l).
## This avoids the label-overlap issue of "source-l vs pool" binary classifiers
## (the pool already contains source l) and, by construction, sum_l P(l) r_l = 1.
## --------------------------------------------------------------------------------
fit_source_density_ratios = function(Xlist, dr_type = "logit", eval_X = NULL,
                                     bound = 1e6, has_intercept = TRUE,
                                     alpha_dr = 0.5) {
  L = length(Xlist)
  Nlist = sapply(Xlist, nrow)
  prior = Nlist / sum(Nlist)

  drop_ic = function(M) {
    M = as.matrix(M)
    if (has_intercept) M[, -1, drop = FALSE] else M
  }

  Xtr = do.call(rbind, lapply(Xlist, drop_ic))
  lab = factor(rep(seq_len(L), times = Nlist), levels = as.character(seq_len(L)))

  Xev = if (is.null(eval_X)) Xtr else drop_ic(eval_X)
  Xev = as.matrix(Xev)

  if (dr_type == "logit") {
    cvfit = cv.glmnet(Xtr, lab, family = "multinomial", alpha = alpha_dr)
    predict_prob = function(newX) {
      p = predict(cvfit, newx = newX, s = "lambda.min", type = "response")
      p[, , 1, drop = TRUE]                      # M x L
    }
  } else if (dr_type == "rf") {
    rf = randomForest(x = as.data.frame(Xtr), y = lab)
    predict_prob = function(newX) {
      predict(rf, newdata = as.data.frame(newX), type = "prob")
    }
  } else {
    stop("fit_source_density_ratios: unsupported dr_type '", dr_type, "'")
  }

  prob_to_ratio = function(prob) {
    prob = prob[, as.character(seq_len(L)), drop = FALSE]   # enforce column order 1..L
    R = t(prob) / prior                                     # L x M (divides row l by prior[l])
    pmin(pmax(R, 1 / bound), bound)
  }

  Rmat = prob_to_ratio(predict_prob(Xev))

  predict_ratios = function(newX) prob_to_ratio(predict_prob(drop_ic(newX)))

  list(Rmat = Rmat, prior = prior, predict_ratios = predict_ratios,
       eval_X = Xev)
}

## --------------------------------------------------------------------------------
## Weighted tau-quantile (values already un-sorted, weights need not be normalized).
## --------------------------------------------------------------------------------
weighted_quantile = function(x, w, tau) {
  keep = is.finite(x) & is.finite(w) & (w > 0)
  if (!any(keep)) return(0)
  x = x[keep]; w = w[keep]
  o = order(x)
  x = x[o]; w = w[o]
  cw = cumsum(w) / sum(w)
  idx = which(cw >= tau)[1]
  if (is.na(idx)) idx = length(x)
  x[idx]
}

## --------------------------------------------------------------------------------
## kappa*(F || H) = essinf_x r_F(x)/r_H(x), approximated by a small tau-quantile of
## r_F/r_H UNDER H. Candidates stay valid distributions, so "under H" = reweight the
## (pooled) eval points by r_H >= 0. Capped at kappa_max for numerical stability.
##
## rF, rH: length-M vectors of the two candidates' density ratios at the eval points.
## --------------------------------------------------------------------------------
kappa_hat = function(rF, rH, tau = 0.02, kappa_max = 0.98, eps = 1e-8) {
  w = pmax(rH, 0)
  if (sum(w) <= 0) return(0)
  ratio = rF / pmax(rH, eps)
  ratio = pmax(ratio, 0)
  k = weighted_quantile(ratio, w, tau)
  min(max(k, 0), kappa_max)
}

## --------------------------------------------------------------------------------
## Katz-Samuels demixing (Demix / DemixHat, JMLR 2019) in coefficient space.
##
## Each candidate distribution is a coefficient vector a in R^L over the observed
## sources, F = sum_l a_l P^(l) with sum_l a_l = 1; its density ratio to the reference
## is a^T %*% Rmat. The two-sample operator kappa*(F | H) = essinf dF/dH is realized here
## by our density-ratio quantile (Katz's inf over a VC class becomes the low quantile of
## the ratio under H), and Res(F | H) = (F - kappa H)/(1 - kappa) is an affine update of
## the coefficient vector. DemixHat (Katz Alg. 11) is recursive: for K = 2 the two
## residues are a permutation of the two base distributions; for K > 2 it uses
## FindFace/FaceTest to place K-1 candidates on a common (K-1)-face, recurses on them,
## and peels off the last base by sequential two-sample residues.
## --------------------------------------------------------------------------------

## kappa*(F | H) for coefficient vectors (two-sample operator).
kappa_coef = function(aF, aH, Rmat, tau = 0.02, kappa_max = 0.98)
  kappa_hat(as.vector(aF %*% Rmat), as.vector(aH %*% Rmat), tau = tau, kappa_max = kappa_max)

## Residue(F | H) in coefficient space (Katz Alg. 1 / ResidueHat Alg. 10).
residue_coef = function(aF, aH, Rmat, tau = 0.02, kappa_max = 0.98) {
  k = kappa_coef(aF, aH, Rmat, tau, kappa_max)
  as.vector((aF - k * aH) / (1 - k))
}

## FaceTestHat (Katz Alg. 13): the rows of A lie on a common face iff every ordered pair
## is mutually reducible, i.e. all pairwise kappa* exceed the threshold eps.
face_test_coef = function(A, Rmat, eps, tau = 0.02) {
  K = nrow(A)
  if (K < 2) return(TRUE)
  for (i in 1:K) for (j in 1:K) if (i != j)
    if (kappa_coef(A[i, ], A[j, ], Rmat, tau) <= eps) return(FALSE)
  TRUE
}

## FindFaceHat (Katz Alg. 12): from a random Q in conv(S_2,...,S_K), push the K-1
## candidates Res((1/n) S_i + (1-1/n) Q | S_1) toward a common face as n grows, until
## FaceTest passes; randomized with a few restarts, falling back to the plain residues.
find_face_coef = function(A, Rmat, eps, tau = 0.02, max_iter = 60, restarts = 8) {
  K = nrow(A); rest = A[2:K, , drop = FALSE]
  for (rs in seq_len(restarts)) {
    w = rexp(K - 1); w = w / sum(w); Q = as.vector(w %*% rest)
    for (n in 2:max_iter) {
      R = t(sapply(2:K, function(i)
        residue_coef((1 / n) * A[i, ] + (1 - 1 / n) * Q, A[1, ], Rmat, tau)))
      if (face_test_coef(R, Rmat, eps, tau)) return(R)
    }
  }
  t(sapply(2:K, function(i) residue_coef(A[i, ], A[1, ], Rmat, tau)))   # fallback
}

## DemixHat (Katz Alg. 11): return a permutation of the base distributions as rows of
## coefficient vectors. A: K x L, row i = coefficients of the i-th input source.
demix_hat = function(A, Rmat, eps = 0.05, tau = 0.02, max_iter = 60, restarts = 8) {
  K = nrow(A)
  if (K == 2)
    return(rbind(residue_coef(A[1, ], A[2, ], Rmat, tau),
                 residue_coef(A[2, ], A[1, ], Rmat, tau)))
  R  = find_face_coef(A, Rmat, eps, tau, max_iter, restarts)
  Q  = demix_hat(R, Rmat, eps, tau, max_iter, restarts)         # K-1 recovered bases
  qK = colMeans(A)
  for (i in 1:(K - 1)) qK = residue_coef(qK, Q[i, ], Rmat, tau) # peel off the last base
  rbind(Q, qK)
}

## Greedy max-min selection of m maximally-separated observed sources, used when L > m to
## reduce to the square (K = m) demixing problem on a well-conditioned sub-mixing.
select_spread_sources = function(Rmat, m) {
  L = nrow(Rmat); if (m >= L) return(seq_len(L))
  D = matrix(0, L, L)
  for (i in 1:(L - 1)) for (j in (i + 1):L) { d = mean(abs(Rmat[i, ] - Rmat[j, ])); D[i, j] = d; D[j, i] = d }
  st = which(D == max(D), arr.ind = TRUE)[1, ]; sel = c(st[[1]], st[[2]])
  while (length(sel) < m) {
    rem = setdiff(seq_len(L), sel)
    sel = c(sel, rem[which.max(sapply(rem, function(r) min(D[r, sel])))])
  }
  sort(sel)[seq_len(m)]
}

## --------------------------------------------------------------------------------
## Orchestrator. Supply source covariates (with estimated ratios) OR a precomputed
## oracle ratio matrix (L x M) to isolate the demixing algorithm from nuisance error.
## --------------------------------------------------------------------------------
recover_intrinsic_sources = function(Xlist, m, dr_type = "logit", tau = 0.02,
                                     eval_X = NULL, oracle_Rmat = NULL, bound = 1e6,
                                     has_intercept = TRUE, eps = 0.05, verbose = FALSE) {
  L = length(Xlist)
  if (!is.null(oracle_Rmat)) {
    Rmat = oracle_Rmat
    fit = NULL
  } else {
    fit = fit_source_density_ratios(Xlist, dr_type = dr_type, eval_X = eval_X,
                                    bound = bound, has_intercept = has_intercept)
    Rmat = fit$Rmat
  }
  ## Square case K = m (Katz Section 4.4). If L > m, demix a maximally-separated m-subset.
  sel = if (L > m) select_spread_sources(Rmat, m) else seq_len(L)
  A0 = matrix(0, nrow = m, ncol = L)
  for (i in seq_len(m)) A0[i, sel[i]] = 1                  # each row = an observed source e_l
  A = demix_hat(A0, Rmat, eps = eps, tau = tau)            # m x L; rows = recovered base coeffs
  alpha = t(A)                                             # L x m; Q_j = sum_l alpha_lj P^(l)
  components = pmax(t(alpha) %*% Rmat, 0)                  # [.]_+ floored component density ratios

  ## Post-hoc irreducibility diagnostic: kappa*(Q_i || Q_j) ~ 0 for all pairs.
  max_kappa = 0
  if (m >= 2) for (i in 1:(m - 1)) for (j in (i + 1):m)
    max_kappa = max(max_kappa, kappa_hat(components[i, ], components[j, ], tau),
                    kappa_hat(components[j, ], components[i, ], tau))
  if (verbose) cat(sprintf("  Demix on sources {%s} | max pairwise kappa* = %.5f\n",
                           paste(sel, collapse = ", "), max_kappa))
  list(alpha = alpha, components = components, sel = sel,
       max_kappa = max_kappa, Rmat = Rmat, fit = fit)
}

## --------------------------------------------------------------------------------
## DORM hook (projection-primary): signed-affine prior via projection onto the
## valid-mixture set. This is the object DORM actually needs, and it requires NO
## intrinsic-source recovery, no irreducibility/separability, and no choice of m.
##
## DORM's default prior is the KL-projection of the target covariate law onto the
## observed-source simplex Delta_L (mulcvxr). The intrinsic-source extension enlarges
## the feasible set to the SIGNED-AFFINE valid-mixture set
##        R = { rho : 1'rho = 1,  sum_l rho_l * rbar_l(x) >= 0  for all x },
## i.e. affine weights whose induced mixture sum_l rho_l P^(l) is still a valid
## distribution (rbar_l = dP^(l)/dP^(ref)). R is identifiable from the observed sources
## alone and is compact iff the sources are affinely independent; it contains Delta_L,
## equals conv{Q_j} under separability, and lets the target be represented as a mixture
## that is "purer" than any single observed source. The prior is the projection of the
## target onto R:
##        rho = argmax_{rho in R} mean_x log( sum_l rho_l rbar_l(x) ).
##
## We solve exactly DORM's mulcvxr log-likelihood projection but (i) drop rho >= 0
## (signed weights allowed) and (ii) add the validity constraints sum_l rho_l rbar_l >= eps
## on the supplied evaluation points. Because the induced mixture ratio is then bounded
## below by eps > 0 wherever it is constrained, the posterior denominators stay positive
## and the standard posterior() is used downstream with no special handling.
##
## dr_obj:    L x n0 target-referenced ratios (dP^ref/dP^(l)); objective uses 1/dr_obj.
## dr_constr: list of L x . ratio matrices giving the points on which the induced mixture
##            must be a valid density (typically the target sample and the source samples).
## Returns the signed weight vector rho (length L, sums to 1).
## --------------------------------------------------------------------------------
## The ridge `lambda` here is the extrapolation regularizer: because signed weights have
## far higher variance than simplex weights, it must be substantially larger than
## mulcvxr's (default 0.1 vs 1e-4). It shrinks rho toward uniform mixing, so the prior
## extrapolates beyond Delta_L only when the target sample gives strong evidence for it;
## with estimated ratios this is what prevents overfitting into extreme signed weights.
signed_affine_prior = function(dr_obj, dr_constr = NULL, lambda = 0.1, eps = 1e-3) {
  L = nrow(dr_obj)
  ele = 1 / dr_obj                                        # L x n0, rbar_l on the target sample
  sc = median(ele); if (!is.finite(sc) || sc <= 0) sc = 1
  ele = ele / sc                                          # condition the log argument (as in mulcvxr)

  if (is.null(dr_constr)) dr_constr = list(dr_obj)
  Cmat = do.call(cbind, lapply(dr_constr, function(M) (1 / M) / sc))   # L x K validity points

  rho = Variable(L)
  obj = mean(-log(t(rho) %*% ele)) + lambda * sum_squares(rho)
  constraints = list(sum(rho) == 1, t(rho) %*% Cmat >= eps)            # signed rho; induced density >= eps
  problem = Problem(Minimize(obj), constraints = constraints)

  res = if (exists("solve_robust")) solve_robust(problem) else
        tryCatch(solve(problem), error = function(e) NULL)

  if (is.null(res)) {
    warning("signed_affine_prior: solver failed; falling back to the simplex prior (mulcvxr).")
    return(if (exists("mulcvxr")) mulcvxr(L, dr_obj) else rep(1 / L, L))
  }
  rho_val = as.vector(res$getValue(rho))
  rho_val / sum(rho_val)                                  # enforce sum-to-one exactly (signed retained)
}
