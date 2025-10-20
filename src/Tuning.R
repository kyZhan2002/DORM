library(Matrix)
library(MASS)
source("/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/Functions3.R")

tuning_Y = function(Xtest,Ytest,q,beta_array,smax_array = seq(0.1,1,0.1)){
  
  # This function tunes smax using Ytest.
  # Each row of matrix beta_array corresponds to each element in smax array.
  
  if( nrow(beta_array) != length(smax_array) ){
    stop("Rows of beta_array must be identical to the number of smax to be tuned.")
  }
  nsmax = length(smax_array)
    
  if( nrow(Xtest) != length(Ytest) ){
    stop("X and Y sample size should be the same.")
  }
  nn = length(Ytest)
  
  losses = rep(0,nsmax)
  for(i in 1:nsmax){
    beta = beta_array[i,]
    loss = get_err(Xtest,Ytest,beta,q,nn)
    losses[i] = loss
  }
  j = which.min(losses)
  
  out = list(losses = losses, best_smax = smax_array[j], best_loss = min(losses))
  return(out)
}
  

tuning_surrogate = function(X0,S0,q,beta_array,smax_array = seq(0.1,1,0.1)){
  
  # This function tunes smax using the surrogate S
  # Each row of matrix beta_array corresponds to each element in smax array.
  
  if( nrow(beta_array) != length(smax_array) ){
    stop("Rows of beta_array must be identical to the number of smax to be tuned.")
  }
  nsmax = length(smax_array)
  
  if( nrow(X0) != length(S0) ){
    stop("X and S sample size should be the same.")
  }
  
  correlation = rep(0,nsmax)
  
  for(i in 1:nsmax){
    beta = beta_array[i,]
    c = cor(S0, X0[,1:(q+1)] %*% beta)
    correlation[i] = c
  }
  j = which.max(correlation)
  
  out = list(correlation = correlation, best_smax = smax_array[j],
             best_correlation = max(correlation))
  return(out)
}


benchmarks = function(X0,Y0,q,beta_star,beta_MI,rho,P,Q){
  
  # This function calculates the performance of all the benchmarks.
  L = ncol(Q)
  if( nrow(X0) != length(Y0) ){
    stop("X and Y sample size should be the same.")
  }
  ny = length(Y0)
  err_star = get_err(X0,Y0,beta_star,q,ny)
  
  errSSQ = 1e8
  SSQvec = rep(0,L)
  minsiteQ = 1
  for(i in 1:L){
    err = get_err(X0,Y0,Q[,i],q,ny)
    SSQvec[i] = err
    if(err < errSSQ){
      minsiteQ = i
      errSSQ = err
    }
  }
  
  betaSA = rowMeans(Q)
  errSA = get_err(X0,Y0,betaSA,q,ny)
  
  betaRA = rep(0,q+1)
  for(l in 1:L){
    betaRA = betaRA + rho[l] * Q[,l]
  }
  errRA = get_err(X0,Y0,betaRA,q,ny)
  
  errMI = get_err(X0,Y0,beta_MI,q,ny)
  
  errRAP = get_err(X0,Y0,P,q,ny)
  
  errvec = c(err_star,errSSQ,errSA,errRA,errMI,errRAP)
  names(errvec) = c('REMIX',"Best single","SA","RA","MI","PA")
  
  losses = list(errvec = errvec, single_site = SSQvec)
  return(losses)
}

