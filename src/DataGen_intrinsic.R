####################################################################################
##  DataGen_intrinsic.R
##
##  Data generator for the intrinsic-source setting: each observed source is a
##  mixture of m latent Gaussian components {Q_j}, and the target is a (possibly
##  extrapolating) mixture of the same components. Mirrors the covariate contract of
##  src/DataGeneration3.R (intercept in column 1; Xlist / Xtrainlist equal-size
##  splits; target X0 / X0train), so the output plugs straight into get_DORM_beta().
##
##  Intrinsic component j:   Q_j = N(MU[j, ], sigma^2 I_p).
##  Observed source l:       P^(l) = sum_j Pi[l, j] Q_j    (Pi rows on the simplex).
##  Outcome (component-tied): a point drawn from component j has
##                           Y = (1, x)^T GAMMA[j, ] + eps_Y * noise.
##
##  Identifiability (Katz-Samuels): Pi must have full column rank and the components
##  must be irreducible (distinct Gaussians => kappa* between any two is 0). A
##  diagonally dominant Pi (each source dominated by a distinct component) is the
##  clean recoverable regime.
####################################################################################

library(MASS)

## Build a diagonally dominant L x m mixing matrix (rows on the simplex): each source
## puts mass `diag` on its "own" component (l mod m) and spreads the rest evenly.
make_diag_dominant_Pi = function(L, m, diag = 0.7) {
  Pi = matrix((1 - diag) / (m - 1), nrow = L, ncol = m)
  for (l in seq_len(L)) {
    j = ((l - 1) %% m) + 1
    Pi[l, ] = (1 - diag) / (m - 1)
    Pi[l, j] = diag
  }
  Pi
}

## Sample a mixture of Gaussians: returns list(X = N x (1+p) with intercept, Z = labels).
.sample_mixture = function(N, weights, MU, sigma) {
  p = ncol(MU)
  Z = sample(seq_along(weights), N, replace = TRUE, prob = weights)
  Xmat = matrix(0, nrow = N, ncol = p)
  for (j in seq_along(weights)) {
    idx = which(Z == j)
    if (length(idx) > 0) {
      Xmat[idx, ] = mvrnorm(length(idx), mu = MU[j, ], Sigma = sigma^2 * diag(p))
    }
  }
  list(X = cbind(1, Xmat), Z = Z)
}

generate_intrinsic_data = function(L, m, p, Nlist, nlist, n0,
                                   Pi, MU, GAMMA, rhoQ_target,
                                   sigma = 1, eps_Y = 0.25, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (nrow(Pi) != L || ncol(Pi) != m) stop("Pi must be L x m")
  if (nrow(MU) != m || ncol(MU) != p) stop("MU must be m x p")
  if (nrow(GAMMA) != m || ncol(GAMMA) != p + 1) stop("GAMMA must be m x (p+1)")
  if (length(rhoQ_target) != m) stop("rhoQ_target must have length m")

  Xlist = vector("list", L); Xtrainlist = vector("list", L)
  Ylist = vector("list", L); Ytrainlist = vector("list", L)

  make_source = function(N, n, wl) {
    s = .sample_mixture(N, wl, MU, sigma)
    ## Outcome uses each point's own component coefficient row (component-tied Y).
    Yvec = rowSums(s$X[1:n, ] * GAMMA[s$Z[1:n], ]) + eps_Y * rnorm(n)
    list(X = s$X, Y = Yvec)
  }

  for (l in seq_len(L)) {
    a = make_source(Nlist[l], nlist[l], Pi[l, ])
    b = make_source(Nlist[l], nlist[l], Pi[l, ])       # independent auxiliary split
    Xlist[[l]] = a$X;      Ylist[[l]] = a$Y
    Xtrainlist[[l]] = b$X; Ytrainlist[[l]] = b$Y
  }

  ## Target: mixture over intrinsic components with weights rhoQ_target.
  make_target = function(N) {
    s = .sample_mixture(N, rhoQ_target, MU, sigma)
    Y = rowSums(s$X * GAMMA[s$Z, ]) + eps_Y * rnorm(N)
    list(X = s$X, Y = Y, Z = s$Z)
  }
  t0 = make_target(n0)
  t0train = make_target(n0)

  ## True structural quantities.
  alpha_star = t(ginv(Pi))                             # L x m, columns sum to 1
  ## Induced (signed) prior weights over observed sources implied by the target.
  rho_bar_target = as.vector(alpha_star %*% rhoQ_target)
  ## True target regression coefficient on (1, A) block is the rhoQ-weighted mixture
  ## of component coefficients (used only as a reference; DORM predicts on (1, A)).
  beta_target = as.vector(t(GAMMA) %*% rhoQ_target)

  list(Xlist = Xlist, Xtrainlist = Xtrainlist,
       Ylist = Ylist, Ytrainlist = Ytrainlist,
       X0 = t0$X, X0train = t0train$X, Y0 = t0$Y,
       Pi = Pi, alpha_star = alpha_star, MU = MU, GAMMA = GAMMA,
       rhoQ_target = rhoQ_target, rho_bar_target = rho_bar_target,
       beta_target = beta_target, sigma = sigma)
}

## --------------------------------------------------------------------------------
## Oracle density ratios r_l(x) = dP^(l)/dP^(pool)(x) for the Gaussian generator,
## evaluated on eval_X (rows are observations; intercept column optional). Used by the
## recovery simulation's Stage 1 to isolate the demixing algorithm from ML error.
##   P^(pool) = sum_l (N_l / sum N) P^(l),  P^(l) = sum_j Pi[l,j] N(MU_j, sigma^2 I).
## --------------------------------------------------------------------------------
oracle_source_ratios = function(eval_X, Pi, MU, sigma, Nlist, has_intercept = TRUE) {
  X = as.matrix(eval_X)
  if (has_intercept) X = X[, -1, drop = FALSE]
  L = nrow(Pi); m = nrow(MU); M = nrow(X); p = ncol(MU)

  ## Component densities phi_j(x) up to the common (2*pi*sigma^2)^(-p/2) factor,
  ## which cancels in all ratios.
  logphi = matrix(0, nrow = M, ncol = m)
  for (j in seq_len(m)) {
    d2 = rowSums(sweep(X, 2, MU[j, ], "-")^2)
    logphi[, j] = -d2 / (2 * sigma^2)
  }
  phi = exp(logphi - apply(logphi, 1, max))            # stabilize; per-row constant cancels
  ## Source densities (unnormalized-consistently): dens_l(x) = sum_j Pi[l,j] phi_j(x).
  dens = phi %*% t(Pi)                                 # M x L
  prior = Nlist / sum(Nlist)
  pool = as.vector(dens %*% prior)                     # M vector
  R = t(dens / pool)                                   # L x M, r_l = dens_l / pool
  R
}
