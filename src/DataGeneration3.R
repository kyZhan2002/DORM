library(Matrix)
library(MASS)

Generate_Simulation_Data = function(L,p,q,Nlist,nlist,n0,
                                    mixture,delta,s, MU_A, BETA, 
                                    ptemp = 6,MU_Acoef = 1, eps_A=0.25,
                                    eps_W=0.1,eps_beta = 0.1,eps_Y = 0.25){
  
  ## This function generates simulation data.
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability, w.p. s we choose delta as model prior.
  
  Xlist = list() # Each X is N*(1+q+p)
  Xtrainlist = list()
  Ylist = list() # Each Y is n*1
  Ytrainlist = list()
  # Slist = list() # Each S is n*1

  eps_S = 0.25
  
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
  S0 = rep(0,n0)
  Assigns = rep(0,n0)
  # mixture is a vector with length L.
  # Vectorized sample for indi
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  
  S = sample(1:L, n0, replace = TRUE, prob = mixture)
  A = matrix(NA, nrow = n0, ncol = q)
  for (l in 1:L) {
    n_l = sum(S == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A[S == l, ] = A_l
    }
  }
  W = A %*% t(WtoA) + mvrnorm(n0,mu=rep(0,p),Sigma = eps_W * diag(p))
  X0 = cbind(1,A,W)
  
  Y0 = rowSums(X0 * BETA[S,]) + eps_Y * rnorm(n0)
  SS = sample(1:L, sum(indi), replace = TRUE, prob = delta)
  Y0[indi == 1] = rowSums(X0[indi == 1,] * BETA[SS,]) + eps_Y * rnorm(sum(indi))
  Surrogate = 0.6 * Y0 + 0.05 * rowSums(X0[,1:(q+1)] * BETA[S,1:(q+1)]) + eps_S * rnorm(n0)
  
  # Assign values to S0, Y0, and Assigns
  S0 = Surrogate
  Assigns = S
  Assigns_Y0 = S
  Assigns_Y0[indi == 1] = SS
  
  S_train = sample(1:L, n0, replace = TRUE, prob = mixture)
  A_train = matrix(NA, nrow = n0, ncol = q)
  for (l in 1:L) {
    n_l = sum(S_train == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A_train[S_train == l, ] = A_l
    }
  }
  W_train = A_train %*% t(WtoA) + mvrnorm(n0, mu = rep(0, p), Sigma = eps_W * diag(p))
  X0train = cbind(1, A_train, W_train)
  Assigns_train = S_train

  out = list('Xlist'=Xlist, 'Xtrainlist' =Xtrainlist,  
             'Ylist' = Ylist, 'Ytrainlist' = Ytrainlist,
             'X0' = X0, 'Y0'=Y0 ,'S0'=S0,  'X0train' = X0train,
             'Assigns' = Assigns, 'Assigns_train' = Assigns_train, 'BETA' = BETA)
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

Generate_Simulation_Data_ExtraA = function(L,p,q,exA,Nlist,nlist,n0,
                                    mixture,delta,s, MU_A, BETA, 
                                    ptemp = 6,MU_Acoef = 1, eps_A=0.25,
                                    eps_W=0.01,eps_beta = 0.1){
  
  ## This function generates simulation data.
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability, w.p. s we choose delta as model prior.
  
  Xlist = list() # Each X is N*(1+q+p)
  Xtrainlist = list()
  Ylist = list() # Each Y is n*1
  Ytrainlist = list()
  # Slist = list() # Each S is n*1
  
  eps_Y = 0.25
  eps_S = 0.25
  
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
  
  zeros = matrix(rnorm(L*(p-ptemp+exA),mean = 0,sd = eps_beta),nrow=L)
  BETA = cbind(BETA,zeros)
  
  for (l in 1:L) {
    n = nlist[l]
    N = Nlist[l]
    A = mvrnorm(N,mu = MU_A[l,], Sigma = eps_A * diag(q))  # mu_l = l
    
    #########
    # The only difference is here: we add extra noise to A.
    #########
    
    extra = matrix(rnorm(N*exA, sd=0.5),nrow = N,ncol = exA)
    
    W = A %*% t(WtoA) + mvrnorm(N,mu = rep(0,p),Sigma = eps_W * diag(p))
    X = cbind(1,A,extra,W)  # a column of 1
    
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
    
    #########
    # The only difference is here: we add extra noise to A.
    #########
    
    extra = matrix(rnorm(N*exA, sd=0.5),nrow = N,ncol = exA)
    
    W = A %*% t(WtoA) + mvrnorm(N,mu = rep(0,p),Sigma = eps_W * diag(p))
    X = cbind(1,A,extra,W)  # a column of 1
    
    beta = BETA[l,]
    
    Y = X[1:n,] %*% beta + eps_Y *rnorm(n)
    
    Xtrainlist[[l]] = X
    Ytrainlist[[l]] = Y
  }
  
  X0 = matrix(data = NA, nrow = n0, ncol = p+q+1, byrow = FALSE)
  Y0 = rep(0,n0)
  S0 = rep(0,n0)
  Assigns = rep(0,n0)
  # mixture is a vector with length L.
  # Vectorized sample for indi
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  
  S = sample(1:L, n0, replace = TRUE, prob = mixture)
  A = matrix(NA, nrow = n0, ncol = q)
  for (l in 1:L) {
    n_l = sum(S == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A[S == l, ] = A_l
    }
  }
  extra = matrix(rnorm(n0*exA, sd=0.5),nrow = n0,ncol = exA)
  W = A %*% t(WtoA) + mvrnorm(n0,mu=rep(0,p),Sigma = eps_W * diag(p))
  X0 = cbind(1,A,extra,W)
  
  Y0 = rowSums(X0 * BETA[S,]) + eps_Y * rnorm(n0)
  SS = sample(1:L, sum(indi), replace = TRUE, prob = delta)
  Y0[indi == 1] = rowSums(X0[indi == 1,] * BETA[SS,]) + eps_Y * rnorm(sum(indi))
  Surrogate = 0.6 * Y0 + 0.05 * rowSums(X0[,1:(q+1)] * BETA[S,1:(q+1)]) + eps_S * rnorm(n0)
  
  # Assign values to S0, Y0, and Assigns
  S0 = Surrogate
  Assigns = S
  
  S_train = sample(1:L, n0, replace = TRUE, prob = mixture)
  A_train = matrix(NA, nrow = n0, ncol = q)
  for (l in 1:L) {
    n_l = sum(S_train == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A_train[S_train == l, ] = A_l
    }
  }
  extra = matrix(rnorm(n0*exA, sd=0.5),nrow = n0,ncol = exA)
  W_train = A_train %*% t(WtoA) + mvrnorm(n0, mu = rep(0, p), Sigma = eps_W * diag(p))
  X_train = cbind(1, A_train,extra, W_train)
  X0train = X_train
  Assigns_train = S_train
  
  out = list('Xlist'=Xlist, 'Xtrainlist' =Xtrainlist,  
             'Ylist' = Ylist, 'Ytrainlist' = Ytrainlist,
             'X0' = X0, 'Y0'=Y0 ,'S0'=S0,  'X0train' = X0train,
             'Assigns' = Assigns, 'Assigns_train' = Assigns_train)
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

Generate_test_data = function(L,p,q,n0,
                              mixture,delta,s, MU_A, BETA, 
                              ptemp = 6,MU_Acoef = 1, eps_A=0.25,
                              eps_W=0.01,  eps_beta = 0.1,eps_Y = 0.25){
  
  ## This function generates test X and Y.
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability, w.p. s we choose delta as model prior.
  
  eps_S = 0.25
  
  
  # MU_A = 1.0 * matrix(c(-2,0,1,0,
  #                      0,-1,2,1,
  #                      -1,2,0,1,
  #                      1,1,-1,2,
  #                      2,1,1,2),nrow=L,ncol=q,byrow=TRUE) # Each row is the mean of A_l
  
  if (dim(MU_A)[1] != L || dim(MU_A)[2] != q) {
    stop("Error: MU_A matrix size must be L*q.")
  }
  MU_A = MU_Acoef * MU_A
  
  WtoA = 0.3* matrix(c(1,0,-1,0,
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
  S0 = rep(0,n0)
  Assigns = rep(0,n0)
  # mixture is a vector with length L.
  # Vectorized sample for indi
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  
  S = sample(1:L, n0, replace = TRUE, prob = mixture)
  A = matrix(NA, nrow = n0, ncol = q)
  for (l in 1:L) {
    n_l = sum(S == l)
    if (n_l > 0) {
      A_l = mvrnorm(n_l, mu = MU_A[l,], Sigma = eps_A * diag(q))
      A[S == l, ] = A_l
    }
  }
  W = A %*% t(WtoA) + mvrnorm(n0,mu=rep(0,p),Sigma = eps_W * diag(p))
  X0 = cbind(1,A,W)
  
  Y0 = rowSums(X0 * BETA[S,]) + eps_Y * rnorm(n0)
  SS = sample(1:L, sum(indi), replace = TRUE, prob = delta)
  Y0[indi == 1] = rowSums(X0[indi == 1,] * BETA[SS,]) + eps_Y * rnorm(sum(indi))
  Surrogate = 0.6 * Y0 + 0.05 * rowSums(X0[,1:(q+1)] * BETA[S,1:(q+1)]) + eps_S * rnorm(n0)
  
  # Assign values to S0, Y0, and Assigns
  S0 = Surrogate
  Assigns = S
  
  out = list('X0' = X0, 'Y0'=Y0 ,'S0'=S0,'Assigns' = Assigns,'BETA'=BETA)
  ################################
  #  Nlist[l] = N_l, nlist[l] = n_l
  #  X0: target covariates. A n_0 * (1+q+p) matrix.
  #  Y0: target labels/outcomes. A n_0-vector. We can use part of it or all of it.
  #  S0: Surrogate vector.
  #  Assigns: n0-vector, tells us which site the i-th target data (X0[i],Y0[i]) is assigned to.
  ################################
  return(out)
}


Generate_test_Y = function(L,p,q,X0,Assigns,delta,s,BETA, ptemp = 6,
                           eps_beta = 0.05, eps_Y = 0.25){
  
  ## This function generates test Y given X.
  ## L means # of source sites, q is the length of A (without the column of 1)
  ## p is length of W, Nlist is the amount of data X in source and nlist is # of those with label Y
  ## n0 is num of target data X. No label.
  ## mixture is the prior probability, w.p. s we choose delta as model prior.
  
  eps_S = 0.25

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
  S0 = rep(0,n0)
  indi = sample(c(0, 1), n0, replace = TRUE, prob = c(1 - s, s))
  S = Assigns
  
  Y0 = rowSums(X0 * BETA[S,]) + eps_Y * rnorm(n0)
  SS = sample(1:L, sum(indi), replace = TRUE, prob = delta)
  Y0[indi == 1] = rowSums(X0[indi == 1,] * BETA[SS,]) + eps_Y * rnorm(sum(indi))
  Surrogate = 0.6 * Y0 + 0.05 * rowSums(X0[,1:(q+1)] * BETA[S,1:(q+1)]) + eps_S * rnorm(n0)
  S0 = Surrogate
  
  out = list('Y0'=Y0 ,'S0'=S0)
  ################################
  #  Y0: target labels/outcomes. A n_0-vector. We can use part of it or all of it.
  #  S0: Surrogate vector.
  ################################
  return(out)
}


