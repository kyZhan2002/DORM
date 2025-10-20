library(Matrix)
library(MASS)

Generate_Simulation_Data_contamination = function(L,p,q,Nlist,nlist,n0,
                                                  mixture, s, MU_A, BETA,  
                                                  ptemp = 6,MU_Acoef = 1, eps_A=0.25,
                                                  eps_W=0.01,eps_beta = 0.1,eps_Y = 0.25,contam_coef = 0.6){
  
  ## This function generates simulation data in a contamination model
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability rho, w.p. s we choose delta as model prior.
  
  Xlist = list() # Each X is N*(1+q+p)
  Xtrainlist = list()
  Ylist = list() # Each Y is n*1
  Ytrainlist = list()

  # MU_A = 1.0 * matrix(c(-2,0,1,0,
  #                      0,-1,2,1,
  #                      -1,2,0,1,
  #                      1,1,-1,2,
  #                      2,1,1,2),nrow=L,ncol=q,byrow=TRUE) # Each row is the mean of A_l
  
  if (dim(MU_A)[1] != L || dim(MU_A)[2] != q) {
    stop("Error: MU_A matrix size must be L*q.")
  }
  MU_A = MU_Acoef * MU_A
  
  WtoA = 0.5* matrix(c(1,0,-1,0,
                       0,1,0,-1,
                       0,0,1,0,
                       0,0,0,1,
                       0,0,0,1,
                       0,0,0,0),nrow=ptemp,ncol=q,byrow=TRUE)
  
  zeros = matrix(0,nrow=p-ptemp,ncol=q)
  WtoA = rbind(WtoA,zeros) # W = WtoA %*% A + Weps
  
  # BETA = 0.42 * matrix(c(8,7,10,-3,2,  0.5,-0.3,0,0,1,0,
  #                       5,2,5,2,-5, 0,0.5,-0.3,0,0,1,
  #                       10,9,8,-4,5,  0,0,0.3,-0.5,0,1,
  #                       -6,6,7,3,-4,  1,0,0,0.3,-0.5,0,
  #                       -4,10,9,-6,5, 0,0,0,1,0.3,-0.5),nrow=L,ncol=ptemp+q+1,byrow=TRUE)
  
  if (dim(BETA)[1] != L || dim(BETA)[2] != q+1+ptemp) {
    stop("Error: BETA matrix size must be L*(q+1+ptemp).")
  }
  
  zeros = matrix(rnorm(L*(p-ptemp),mean = 0,sd = eps_beta),nrow=L)
  BETA = cbind(BETA,zeros)
  
  for (l in 1:L) {
    n = nlist[l]
    N = Nlist[l]
    A = mvrnorm(N,mu = MU_A[l,], Sigma = eps_A * diag(q))  # mu_l = l
    W = A %*% t(WtoA) + mvrnorm(N,mu = rep(0,p),Sigma = eps_W * diag(p))
    X = cbind(1,A,W)  # a column of 1
    
    beta = BETA[l,]
    
    Y = X[1:n,] %*% beta + eps_Y *rnorm(n)
    # Surrogate = 0.6 * Y + 0.1 * X[1:n,1:(q+1)] %*% beta[1:(q+1)] + eps_S * rnorm(n)
    
    Xlist[[l]] = X
    Ylist[[l]] = Y
    # Slist[[l]] = Surrogate
  }
  
  # Generate Auxiliary Data
  
  for (l in 1:L) {
    n = nlist[l]
    N = Nlist[l]
    A = mvrnorm(N,mu = MU_A[l,], Sigma = eps_A * diag(q))  # mu_l = l
    W = A %*% t(WtoA) + mvrnorm(N,mu = rep(0,p),Sigma = eps_W * diag(p))
    X = cbind(1,A,W)  # a column of 1
    
    beta = BETA[l,]
    
    Y = X[1:n,] %*% beta + eps_Y *rnorm(n)
    
    Xtrainlist[[l]] = X
    Ytrainlist[[l]] = Y
  }
  
  X0 = matrix(data = NA, nrow = n0, ncol = p+q+1, byrow = FALSE)
  Y0 = rep(0,n0)
  
  # mixture is a vector with length L.
  # Vectorized sample for indi
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  n_contam = sum(indi)
  
  S = rep(-1, n0) # is model assignment,-1 means contaminated
  S[indi == 0] = sample(1:L, sum(indi == 0), replace = TRUE, prob = mixture)
  
  A = matrix(NA, nrow = n0, ncol = q)
  # not contaminated part
  for (l in 1:L) {
    n_l = sum(S == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A[S == l, ] = A_l
    }
  }
  # contaminated part
  A_con = mvrnorm(n_contam, mu= contam_coef*colSums(MU_A), Sigma = eps_A*diag(q))
  A[indi == 1,] = A_con
  
  W = A %*% t(WtoA) + mvrnorm(n0,mu=rep(0,p),Sigma = eps_W * diag(p))
  X0 = cbind(1,A,W)
  BETA_contam = contam_coef*colSums(BETA)
  Y0[indi == 0] = rowSums(X0[indi == 0, ] * BETA[S[indi == 0], ]) + eps_Y * rnorm(sum(indi == 0))
  # For indi == 1, use BETA_contam for all observations
  Y0[indi == 1] = X0[indi == 1, ] %*% BETA_contam + eps_Y * rnorm(sum(indi == 1))
  
  # Assign values to S0 and Assigns
  Assigns = S
  
  
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  n_contam = sum(indi)
  S_train = rep(-1, n0) # is model assignment,-1 means contaminated
  S_train[indi == 0] = sample(1:L, sum(indi == 0), replace = TRUE, prob = mixture)
  
  A_train = matrix(NA, nrow = n0, ncol = q)
  # not contaminated part
  for (l in 1:L) {
    n_l = sum(S_train == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A_train[S_train == l, ] = A_l
    }
  }
  # contaminated part
  A_con = mvrnorm(n_contam, mu= 0.8*colSums(MU_A), Sigma = eps_A*diag(q))
  A_train[indi == 1,] = A_con
  W_train = A_train %*% t(WtoA) + mvrnorm(n0,mu=rep(0,p),Sigma = eps_W * diag(p))
  X0train = cbind(1,A_train,W_train)
  Assigns_train = S_train
  
  out = list('Xlist'=Xlist, 'Xtrainlist' =Xtrainlist,  
             'Ylist' = Ylist, 'Ytrainlist' = Ytrainlist,
             'X0' = X0, 'Y0'=Y0 , 'X0train' = X0train,
             'Assigns' = Assigns, 'Assigns_train' = Assigns_train,'BETA' = BETA)
  ################################
  #  Nlist[l] = N_l, nlist[l] = n_l
  #  Xlist: a list of source covariates. Xlist[[l]] is N_l * (1+q+p) matrix, source l covariates.
  #  Xtrainlist: equal size to Xlist. Used for sample splitting.
  #  Ylist: a list of labels/outcomes. Ylist[[l]] is n_l vector. Ylist[[l]] is calculated by
  #         the first n_l rows of Xlist[[l]].
  #  X0: target covariates. A n_0 * (1+q+p) matrix.
  #  Y0: target labels/outcomes. A n_0-vector. We can use part of it or all of it.
  #  S0: Surrogate vector.
  #  Assigns(_train): n0-vector, tells us which site the i-th target data (X0[i],Y0[i]) is assigned to.
  #  
  #  ATTENTION: splitted data (Xlist,Xtrainlist) should be of equal size!
  ################################
  return(out)
}

Generate_test_data_contamination = function(L,p,q,n0,
                                            mixture, s, MU_A, BETA,  
                                            ptemp = 6,MU_Acoef = 1, eps_A=0.25,
                                            eps_W=0.01,eps_beta = 0.1,eps_Y = 0.25,contam_coef = 0.6){
  
  ## This function generates simulation data in a contamination model
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability rho, w.p. s we choose delta as model prior.
  
  # MU_A = 1.0 * matrix(c(-2,0,1,0,
  #                      0,-1,2,1,
  #                      -1,2,0,1,
  #                      1,1,-1,2,
  #                      2,1,1,2),nrow=L,ncol=q,byrow=TRUE) # Each row is the mean of A_l
  
  if (dim(MU_A)[1] != L || dim(MU_A)[2] != q) {
    stop("Error: MU_A matrix size must be L*q.")
  }
  MU_A = MU_Acoef * MU_A
  
  WtoA = 0.5* matrix(c(1,0,-1,0,
                       0,1,0,-1,
                       0,0,1,0,
                       0,0,0,1,
                       0,0,0,1,
                       0,0,0,0),nrow=ptemp,ncol=q,byrow=TRUE)
  
  zeros = matrix(0,nrow=p-ptemp,ncol=q)
  WtoA = rbind(WtoA,zeros) # W = WtoA %*% A + Weps
  
  # BETA = 0.42 * matrix(c(8,7,10,-3,2,  0.5,-0.3,0,0,1,0,
  #                       5,2,5,2,-5, 0,0.5,-0.3,0,0,1,
  #                       10,9,8,-4,5,  0,0,0.3,-0.5,0,1,
  #                       -6,6,7,3,-4,  1,0,0,0.3,-0.5,0,
  #                       -4,10,9,-6,5, 0,0,0,1,0.3,-0.5),nrow=L,ncol=ptemp+q+1,byrow=TRUE)
  
  if (dim(BETA)[1] != L || dim(BETA)[2] != q+1+ptemp) {
    stop("Error: BETA matrix size must be L*(q+1+ptemp).")
  }
  
  zeros = matrix(rnorm(L*(p-ptemp),mean = 0,sd = eps_beta),nrow=L)
  BETA = cbind(BETA,zeros)
  
  X0 = matrix(data = NA, nrow = n0, ncol = p+q+1, byrow = FALSE)
  Y0 = rep(0,n0)
  
  # mixture is a vector with length L.
  # Vectorized sample for indi
  # 1 for contaminated
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  n_contam = sum(indi)
  
  S = rep(-1, n0) # is model assignment,-1 means contaminated
  S[indi == 0] = sample(1:L, sum(indi == 0), replace = TRUE, prob = mixture)
  
  A = matrix(NA, nrow = n0, ncol = q)
  # not contaminated part
  for (l in 1:L) {
    n_l = sum(S == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A[S == l, ] = A_l
    }
  }
  # contaminated part
  A_con = mvrnorm(n_contam, mu= contam_coef*colSums(MU_A), Sigma = eps_A*diag(q))
  A[indi == 1,] = A_con
  
  W = A %*% t(WtoA) + mvrnorm(n0,mu=rep(0,p),Sigma = eps_W * diag(p))
  X0 = cbind(1,A,W)
  BETA_contam = contam_coef*colSums(BETA)
  #for(i in 1:n0){
  #  if(indi[i] == 0){
  #    Y0[i] = X0[i,] %*% BETA[S[i],] + eps_Y * rnorm(1)
  #  }else{
  #    Y0[i] = X0[i,] %*% BETA_contam + eps_Y * rnorm(1)
  #  }
  #}

  BETA_contam_matrix <- matrix(rep(BETA_contam, each = n_contam), nrow = n_contam, byrow = FALSE)
  
  # Vectorized computation to replace the loop
  Y0[indi == 0] <- rowSums(X0[indi == 0, ] * BETA[S[indi == 0], ]) + eps_Y * rnorm(n0 - n_contam)
  Y0[indi == 1] <- rowSums(X0[indi == 1, ] * BETA_contam_matrix) + eps_Y * rnorm(n_contam)
  # Assign values to S0, Y0, and Assigns
  Assigns = S
  
  out = list('X0' = X0, 'Y0'=Y0 ,'Assigns' = Assigns,'BETA' = BETA)
  ################################
  #  Nlist[l] = N_l, nlist[l] = n_l
  #  Xlist: a list of source covariates. Xlist[[l]] is N_l * (1+q+p) matrix, source l covariates.
  #  Xtrainlist: equal size to Xlist. Used for sample splitting.
  #  Ylist: a list of labels/outcomes. Ylist[[l]] is n_l vector. Ylist[[l]] is calculated by
  #         the first n_l rows of Xlist[[l]].
  #  X0: target covariates. A n_0 * (1+q+p) matrix.
  #  Y0: target labels/outcomes. A n_0-vector. We can use part of it or all of it.
  #  S0: Surrogate vector.
  #  Assigns(_train): n0-vector, tells us which site the i-th target data (X0[i],Y0[i]) is assigned to.
  #  
  #  ATTENTION: splitted data (Xlist,Xtrainlist) should be of equal size!
  ################################
  return(out)
}

Generate_test_Y_contamination = function(L,p,q,X0,Assigns,s,BETA, ptemp = 6,
                           eps_beta = 0.05, eps_Y = 0.25,contam_coef = 0.6){
  
  ## This function generates test Y given X.
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability, w.p. s we choose delta as model prior.

  
  n0 = nrow(X0)
  if(length(Assigns) != n0){
    stop("Error: Assigns must be equal size of X0.")
  }
  
  # BETA = 0.42 * matrix(c(8,7,10,-3,2,  0.5,-0.3,0,0,1,0,
  #                       5,2,5,2,-5, 0,0.5,-0.3,0,0,1,
  #                       10,9,8,-4,5,  0,0,0.3,-0.5,0,1,
  #                       -6,6,7,3,-4,  1,0,0,0.3,-0.5,0,
  #                       -4,10,9,-6,5, 0,0,0,1,0.3,-0.5),nrow=L,ncol=ptemp+q+1,byrow=TRUE)
  
  if (dim(BETA)[1] != L || dim(BETA)[2] != q+1+ptemp) {
    stop("Error: BETA matrix size must be L*(q+1+ptemp).")
  }
  
  zeros = matrix(rnorm(L*(p-ptemp),mean = 0,sd = eps_beta),nrow=L)
  BETA = cbind(BETA,zeros)
  
  Y0 = rep(0,n0)
  S = Assigns
  
  #indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  n_contam = sum(S == -1)
  
  BETA_contam = contam_coef*colSums(BETA)
  #for(i in 1:n0){
  #  if(indi[i] == 0){
  #    Y0[i] = X0[i,] %*% BETA[S[i],] + eps_Y * rnorm(1)
  #  }else{
  #    Y0[i] = X0[i,] %*% BETA_contam + eps_Y * rnorm(1)
  #  }
  #}
  
  # Vectorized computation to replace the loop
  Y0[S != -1] <- rowSums(X0[S != -1, ] * BETA[S[S != -1], ]) + eps_Y * rnorm(n0-n_contam)
  Y0[S == -1] <- rowSums(X0[S == -1, ] * BETA_contam) + eps_Y * rnorm(n_contam)
  
  
  out = list('Y0'=Y0)
  ################################
  #  Y0: target labels/outcomes. A n_0-vector. We can use part of it or all of it.
  ################################
  return(out)
}
