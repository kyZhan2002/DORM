library(Matrix)
library(MASS)

Generate_Simulation_Data_CondA = function(L, p, q, Nlist, nlist, n0, mixture,
                                          delta, s, MU_A, BETA, 
                                          ptemp = 6, MU_Acoef = 1, eps_A = 1,
                                          eps_W = 0.1, eps_beta = 0.1, eps_Y = 0.25,
                                          sd_WA = 1, tarmixA = TRUE) {
  
  ## This function generates simulation data with mixture conditional on A.
  ## Here the BETA should be of size L*(q+1+ptemp)
  ## mixture: a vector of length L specifying the mixture weights for target
  ## tarmixA: if TRUE, target A uses source-specific means; if FALSE, uses average of MU_A
  
  Xlist = list()
  Xtrainlist = list()
  Ylist = list()
  Ytrainlist = list()
  
  eps_S = 0.25
  
  if (length(mixture) != L || abs(sum(mixture) - 1) > 1e-6) {
    stop("Error: mixture must be a vector of length L that sums to 1.")
  }
  
  if (dim(MU_A)[1] != L || dim(MU_A)[2] != q) {
    stop("Error: MU_A matrix size must be L*q.")
  }
  MU_A = MU_Acoef * MU_A
  
  # Create source-specific WtoA matrices (different W|A models for each source)
  WtoA_list = list()
  base_WtoA = 1 * matrix(c(1, 0, -1, 0,
                           0, 1, 0, -1,
                           0, 0, 1, 0,
                          -1, 0, 2, 1,
                           0,2,-1,1,
                           2,-1,0,-2), nrow = ptemp, ncol = q, byrow = TRUE)
  
  for (l in 1:L) {
    # Add source-specific variation to WtoA
    WtoA_l = base_WtoA + matrix(rnorm(ptemp * q, mean = 0, sd = sd_WA), nrow = ptemp, ncol = q)
    zeros = matrix(rnorm((p-ptemp) * q, mean = 0, sd = 0.1), nrow = p - ptemp, ncol = q)
    WtoA_list[[l]] = rbind(WtoA_l, zeros)
  }
  
  if (dim(BETA)[1] != L || dim(BETA)[2] != q + 1 + ptemp) {
    stop("Error: BETA matrix size must be L*(q+1+ptemp).")
  }
  
  zeros = matrix(rnorm(L * (p - ptemp), mean = 0, sd = eps_beta), nrow = L)
  BETA_full = cbind(BETA, zeros)
  
  # Generate source data
  for (l in 1:L) {
    n = nlist[l]
    N = Nlist[l]
    A = mvrnorm(N, mu = MU_A[l,], Sigma = eps_A * diag(q))
    # Use source-specific W|A model
    W = A %*% t(WtoA_list[[l]]) + mvrnorm(N, mu = rep(0, p), Sigma = eps_W * diag(p))
    X = cbind(1, A, W)
    
    beta = BETA_full[l,]
    Y = X[1:n,] %*% beta + eps_Y * rnorm(n)
    
    Xlist[[l]] = X
    Ylist[[l]] = Y
  }
  
  # Generate auxiliary data
  for (l in 1:L) {
    n = nlist[l]
    N = Nlist[l]
    A = mvrnorm(N, mu = MU_A[l,], Sigma = eps_A * diag(q))
    W = A %*% t(WtoA_list[[l]]) + mvrnorm(N, mu = rep(0, p), Sigma = eps_W * diag(p))
    X = cbind(1, A, W)
    
    beta = BETA_full[l,]
    Y = X[1:n,] %*% beta + eps_Y * rnorm(n)
    
    Xtrainlist[[l]] = X
    Ytrainlist[[l]] = Y
  }
  
  # Generate target data
  # Assign each observation to a source based on mixture weights
  # This assignment determines which W|A model to use
  S = sample(1:L, n0, replace = TRUE, prob = mixture)
  
  # Generate A using source-specific means or average based on tarmixA
  A = matrix(NA, nrow = n0, ncol = q)
  if (tarmixA) {
    # Use source-specific means based on assignment S
    for (l in 1:L) {
      idx = which(S == l)
      if (length(idx) > 0) {
        A[idx,] = mvrnorm(length(idx), mu = MU_A[l,], Sigma = eps_A * diag(q))
      }
    }
  } else {
    # Use average of all source means
    MU_A_avg = colMeans(MU_A)
    A = mvrnorm(n0, mu = MU_A_avg, Sigma = eps_A * diag(q))
  }
  
  # Generate W using source-specific W|A models
  W = matrix(0, nrow = n0, ncol = p)
  for (l in 1:L) {
    idx = which(S == l)
    if (length(idx) > 0) {
      W[idx,] = A[idx,] %*% t(WtoA_list[[l]]) + mvrnorm(length(idx), mu = rep(0, p), Sigma = eps_W * diag(p))
    }
  }
  
  X0 = cbind(1, A, W)
  
  # Generate Y: mostly aligned with W|A assignment S, with small perturbation
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  
  Y0 = rep(0, n0)
  Assigns_Y0 = S  # Start with same assignment as W
  
  # For indi == 0, Y follows the same source as W
  Y0[indi == 0] = rowSums(X0[indi == 0,] * BETA_full[S[indi == 0],]) + eps_Y * rnorm(sum(indi == 0))
  
  # For indi == 1, Y has perturbation: use different source assignment
  SS = sample(1:L, sum(indi), replace = TRUE, prob = delta)
  Y0[indi == 1] = rowSums(X0[indi == 1,] * BETA_full[SS,]) + eps_Y * rnorm(sum(indi))
  Assigns_Y0[indi == 1] = SS
  
  Surrogate = 0.6 * Y0 + 0.05 * rowSums(X0[, 1:(q + 1)] * BETA_full[S, 1:(q + 1)]) + eps_S * rnorm(n0)
  
  # Generate training target data
  S_train = sample(1:L, n0, replace = TRUE, prob = mixture)
  
  A_train = matrix(NA, nrow = n0, ncol = q)
  if (tarmixA) {
    # Use source-specific means based on assignment S_train
    for (l in 1:L) {
      idx = which(S_train == l)
      if (length(idx) > 0) {
        A_train[idx,] = mvrnorm(length(idx), mu = MU_A[l,], Sigma = eps_A * diag(q))
      }
    }
  } else {
    # Use average of all source means
    MU_A_avg = colMeans(MU_A)
    A_train = mvrnorm(n0, mu = MU_A_avg, Sigma = eps_A * diag(q))
  }
  
  W_train = matrix(0, nrow = n0, ncol = p)
  for (l in 1:L) {
    idx = which(S_train == l)
    if (length(idx) > 0) {
      W_train[idx,] = A_train[idx,] %*% t(WtoA_list[[l]]) + mvrnorm(length(idx), mu = rep(0, p), Sigma = eps_W * diag(p))
    }
  }
  
  X0train = cbind(1, A_train, W_train)
  
  out = list('Xlist' = Xlist, 'Xtrainlist' = Xtrainlist,  
             'Ylist' = Ylist, 'Ytrainlist' = Ytrainlist,
             'X0' = X0, 'Y0' = Y0, 'S0' = Surrogate, 'X0train' = X0train,
             'Assigns' = S, 'Assigns_Y0' = Assigns_Y0, 
             'Assigns_train' = S_train, 'BETA' = BETA,
             'WtoA_list' = WtoA_list)
  
  return(out)
}

Generate_test_data_CondA = function(L, p, q, n0, mixture,
                                    delta, s, MU_A, BETA, WtoA_list,
                                    ptemp = 6, MU_Acoef = 1, eps_A = 1,
                                    eps_W = 0.1, eps_beta = 0.1, eps_Y = 0.25,
                                    sd_WA = 1, tarmixA = TRUE) {
  
  ## This function generates test X and Y with mixture conditional on A.
  ## WtoA_list: list of source-specific WtoA matrices (from training data generation)
  ## tarmixA: if TRUE, target A uses source-specific means; if FALSE, uses average of MU_A
  
  eps_S = 0.25
  
  if (length(mixture) != L || abs(sum(mixture) - 1) > 1e-6) {
    stop("Error: mixture must be a vector of length L that sums to 1.")
  }
  
  if (dim(MU_A)[1] != L || dim(MU_A)[2] != q) {
    stop("Error: MU_A matrix size must be L*q.")
  }
  MU_A = MU_Acoef * MU_A
  
  if (is.null(WtoA_list) || length(WtoA_list) != L) {
    stop("Error: WtoA_list must be provided and have length L.")
  }
  
  if (dim(BETA)[1] != L || dim(BETA)[2] != q + 1 + ptemp) {
    stop("Error: BETA matrix size must be L*(q+1+ptemp).")
  }
  
  # Pad BETA to full dimension
  zeros = matrix(rnorm(L * (p - ptemp), mean = 0, sd = eps_beta), nrow = L)
  BETA_full = cbind(BETA, zeros)
  
  # Generate test data
  S = sample(1:L, n0, replace = TRUE, prob = mixture)
  
  # Generate A using source-specific means or average based on tarmixA
  A = matrix(NA, nrow = n0, ncol = q)
  if (tarmixA) {
    # Use source-specific means based on assignment S
    for (l in 1:L) {
      idx = which(S == l)
      if (length(idx) > 0) {
        A[idx,] = mvrnorm(length(idx), mu = MU_A[l,], Sigma = eps_A * diag(q))
      }
    }
  } else {
    # Use average of all source means
    MU_A_avg = colMeans(MU_A)
    A = mvrnorm(n0, mu = MU_A_avg, Sigma = eps_A * diag(q))
  }
  
  # Generate W using source-specific W|A models
  W = matrix(0, nrow = n0, ncol = p)
  for (l in 1:L) {
    idx = which(S == l)
    if (length(idx) > 0) {
      W[idx,] = A[idx,] %*% t(WtoA_list[[l]]) + mvrnorm(length(idx), mu = rep(0, p), Sigma = eps_W * diag(p))
    }
  }
  
  X0 = cbind(1, A, W)
  
  # Generate Y with perturbation
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  
  Y0 = rep(0, n0)
  Y0[indi == 0] = rowSums(X0[indi == 0,] * BETA_full[S[indi == 0],]) + eps_Y * rnorm(sum(indi == 0))
  
  SS = sample(1:L, sum(indi), replace = TRUE, prob = delta)
  Y0[indi == 1] = rowSums(X0[indi == 1,] * BETA_full[SS,]) + eps_Y * rnorm(sum(indi))
  
  Surrogate = 0.6 * Y0 + 0.05 * rowSums(X0[, 1:(q + 1)] * BETA_full[S, 1:(q + 1)]) + eps_S * rnorm(n0)
  
  out = list('X0' = X0, 'Y0' = Y0, 'S0' = Surrogate, 'Assigns' = S, 'BETA' = BETA)
  
  return(out)
}