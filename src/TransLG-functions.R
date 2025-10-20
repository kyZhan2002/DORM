source('/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/TransLasso-functions.R')
library(glmtrans)
library(glmnet)

List2TransLASSO = function(X0,Xlist,Y0,Ylist,ntar,q){
  
  L = length(Xlist)
  X = X0[1:ntar,1:(q+1)]
  y = Y0[1:ntar]
  nlist = rep(0,L)
  for(k in 1:L){
    y = c(y,Ylist[[k]])
    nn = length(Ylist[[k]])
    nlist[k] = nn
    X = rbind(X,Xlist[[k]][1:nn,1:(q+1)])
  }
  n.vec = c(ntar,nlist)
  
  return(list(X=X,y=y,n.vec=n.vec))
  
}

List2TransGLM = function(X0,Xlist,Y0,Ylist,ntar,q){
  
  L = length(Xlist)
  target = list(x = X0[1:ntar,1:(q+1)],y = Y0[1:ntar])
  source = list()
  for(k in 1:L){
    nn = length(Ylist[[k]])
    source[[k]] = list(x = Xlist[[k]][1:nn,1:(q+1)], y = Ylist[[k]])
  }
  
  return(list(source = source,target = target))
  
}

mytranslasso = function(X0,Xlist,Y0,Ylist,ntar,q,Inumd=5){
  # (q+1) is the dimension of X(A).
  
  DataTL = List2TransLASSO(X0,Xlist,Y0,Ylist,ntar,q)
  
  Inum = floor(ntar/Inumd)
  l1 = TRUE
  prop.re1 = Trans.lasso(DataTL$X, DataTL$y, DataTL$n.vec, I.til = 1:Inum, l1 = l1)
  prop.re2 = Trans.lasso(DataTL$X, DataTL$y, DataTL$n.vec, I.til = (ntar-Inum):ntar, l1=l1)
  #if(size.A0 > 0 & size.A0< M){ #Rank.re characterizes the performance of the sparsity index Rk
  #  Rank.re= ( sum(prop.re1$rank.pi[1:size.A0]<=size.A0) +
  #                     sum(prop.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
  #}else{ Rank.re = 1 }
  beta.prop = (prop.re1$beta.hat + prop.re2$beta.hat) / 2
  
  mse_TL = get_err(X0,Y0,beta.prop,length(beta.prop)-1,length(Y0))
  mse_TL
  
  return(list(beta_TL = beta.prop, mse_TL = mse_TL))
}

mytransglm = function(X0,Xlist,Y0,Ylist,ntar,q){
  D.training = List2TransGLM(X0,Xlist,Y0,Ylist,ntar,q)
  fit.gaussian = glmtrans(D.training$target, D.training$source,family = 'gaussian')
  y.pred.glmtrans = predict(fit.gaussian, X0[,1:(q+1)])
  mse_TG = mean((Y0 - y.pred.glmtrans)^2)
  
  return(list(beta_TG = fit.gaussian$beta, mse_TG = mse_TG ))
}


# Main PTL Algorithm
PTL_algorithm = function(Xlist, Ylist, X0, Y0, ntar) {

  X0 = X0[1:ntar,]
  Y0 = Y0[1:ntar]
  K = length(Xlist)
  
  # Step 1: Compute beta^(k) for each source dataset using cross-validation
  beta_k_list = lapply(1:K, function(k) {
    
    Y_k = Ylist[[k]]
    nn = length(Y_k)
    X_k = Xlist[[k]][1:nn,]
    cv_model = cv.glmnet(X_k, Y_k, alpha = 1)
    best_lambda = cv_model$lambda.min
    model = glmnet(X_k, Y_k, alpha = 1, lambda = best_lambda)
    as.vector(coef(model))[-1]  # Remove intercept
  })
  
  # Compute B matrix
  B = do.call(cbind, beta_k_list)
  # Compute Z_i
  Z_i = X0 %*% B
  
  # Step 2: Compute w_ptl
  w_ptl = solve(t(Z_i) %*% Z_i) %*% t(Z_i) %*% Y0
  
  # Step 3: Compute profiled responses
  e_i = Y0 - Z_i %*% w_ptl
  
  # Step 4: Compute delta_hat using cross-validation
  cv_model = cv.glmnet(X0, e_i, alpha = 1)
  best_lambda = cv_model$lambda.min
  model = glmnet(X0, e_i, alpha = 1, lambda = best_lambda)
  delta_hat = as.vector(coef(model))[-1]  # Remove intercept
  
  # Step 5: Compute beta_ptl
  beta_ptl = B %*% w_ptl + delta_hat
  
  mse_PTL = mean((Y0 - X0 %*% beta_ptl)^2)
  
  return(list(beta_PTL = beta_ptl, mse_PTL = mse_PTL ))
}


# These two functions are easy to implement:
# Just input the X0 and Xlist as structured in our method, and choose ntar as the target number size
## you want to use. 

# The small ntar is adverse for these two methods.

# NOTICE: (q+1) is the dimension of X(A).
# mse here is (Y-Xb)^2.

