library(Matrix)
library(MASS)
library(CVXR)
library(glmnet)
library(randomForest)
library(caret)
library(xgboost)
library(nnet)

simplex_uniform = function(dim){
  ## This function samples a dim vector from dim-simplex uniform distribution.
  vec = rexp(dim,rate = 1)
  vec = vec / sum(vec)
  return(vec)
}

or_estimation_ML = function(Xlist,X0,dr_type = 'rf', report = TRUE){
  
  ## This function estimates the odds ratio using random forest.
  
  L = length(Xlist)
  Nlist = sapply(Xlist, nrow)
  n0 = nrow(X0)
  modellist = list()
  
  for (l in 1:L) {
    # 0 stands for target and 1 stands for source l
    Y = c(rep(0, n0), rep(1, Nlist[l]))  
    Y = as.factor(Y)
    X = rbind(X0, Xlist[[l]])
    X = X[,-1] # Remove interception
    X = as.data.frame(X)
    
    # set.seed(123) # for reproducibility
    # Define the training control
    train_control = trainControl(method="cv", number=10)
    
    ####################################################
    # Set up a grid of tuning parameters for parameters.
    ## You can change the tuning grid if you want.
    ####################################################
    
    if(dr_type == 'rf'){
      tune_grid = expand.grid(
        # mtry = seq(5, sqrt(ncol(X)), by=2)
        mtry = floor(sqrt(ncol(X)))
      )
      rfmodel = train(
        x = X, y = Y, 
        method = "rf", 
        trControl = train_control, 
        tuneGrid = tune_grid,
      )
      # Save the tuned model and the best parameters
      modellist[[l]] = rfmodel
      
      
    }else if(dr_type == 'nnet'){
      tune_grid = expand.grid(
        #.size = c(3,5,7,9),   
        .size = 9,
        #.decay = c(0.1, 0.03,0.01) 
        .decay = 0.1,
      )
      
      nnet_model = train(
        x = X, y = Y,
        method = "nnet",
        trControl = train_control,
        tuneGrid = tune_grid,
        trace = FALSE,                 
        MaxNWts = 8000,                
      )
      
      # Save the tuned model and the best parameters
      modellist[[l]] = nnet_model
      
    }else if(dr_type == 'XGB'){
      tune_grid = expand.grid(
        nrounds = c(50,100),
        max_depth = c(6,9),
        eta = 0.3,
        gamma = c(0,0.1),
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
      )
      xgb_model = train(
        x = X, y = Y,
        method = "xgbTree",
        trControl = train_control,
        tuneGrid = tune_grid,
        verbosity = 0
      )
      
      # Save the tuned model and the best parameters
      modellist[[l]] = xgb_model
    }else{
      stop("Wrong dr_type.")
    }
    
    if(report){
    # Print the best parameters
    print(paste("Best parameters for model", l, ":", sep=" "))
    print(modellist[[l]]$bestTune)
    }
  }
  
  return(modellist)
} 

or_estimation_logit = function(Xlist,Xtar,alpha_dr = 0.5){
  
  ## This function estimates logistic model for odds ratio = exp(x'b)
  
  L = length(Xlist)
  modellist = list()
  for(l in 1:L){
    Xsou = Xlist[[l]]
    y1 = rep(1, nrow(Xsou))  # source for 1
    y2 = rep(0, nrow(Xtar))  # target for 0
    y = c(y1, y2)
    X = rbind(Xsou,Xtar)
    X = X[,-1] ## remove the intercept which already exist in X.
    
    cv_fit = cv.glmnet(X, y, alpha = alpha_dr, family = "binomial")
  
    best_lambda = cv_fit$lambda.min
    
    final_model = glmnet(X, y, alpha = alpha_dr, lambda = best_lambda, family = "binomial")
    modellist[[l]] = final_model
  }
  return(modellist)
  
}


dratio_ML = function(modellist,X0,L,Nlist,n0,bound = 1e6,normalize = FALSE){
  
  ## The function output a L*n0 density ratio matrix, each row is density ratio:dP(target X0) to dP(source)
  ## We multiply nl/n0 here!!! This is density ratio
  
  X0 = X0[,-1] # Remove interception
  dr = matrix(0,nrow=L,ncol=nrow(X0))
  for(l in 1:L){
    X0 = as.data.frame(X0)
    pred = predict(modellist[[l]],X0,type = "prob")
    ratio = ifelse(pred[,2] == 0, bound, pred[,1] / pred[,2])
    ratio[ratio==0] = 1/bound
    
    # truncation
    ratio = pmax(1/bound, pmin(ratio, bound))
    
    dr[l,]=Nlist[l] / n0 * ratio
    
    # mean of (1/dr[l,]) should be 1.
    if(normalize){
      nor_const = mean(1/dr[l,])
      dr[l,] = nor_const * dr[l,]
    }
  }
  return(dr)
}


dratio_logit = function(modellist,X0,L,Nlist,n0,bound=1e6,normalize = FALSE){
  
  # use logistic regression to get density ratio!
  
  dr = matrix(0,nrow=L,ncol=nrow(X0))
  for(l in 1:L){
    model = modellist[[l]]
    new_X = X0[,-1] ## remove the intercept which already exist in X.
    log_density_ratio = predict(model, new_X, type = "link")
    density_ratio = exp(log_density_ratio)
    density_ratio = pmax(1/bound, pmin(density_ratio, bound))
    dr[l,] = Nlist[l]/n0 * 1/density_ratio
    
    # mean of (1/dr[l,]) should be 1.
    if(normalize){
      nor_const = mean(1/dr[l,])
      dr[l,] = nor_const * dr[l,]
    }
    
  }
  return(dr)
}

lasso_imputation = function(Xlist,Ylist,nlist){
  
  ## This function conducts a lasso regression to get imputation model.
  
  L = length(Xlist)
  betalist = list()
  for (l in 1:L){
    n = nlist[l]
    X = Xlist[[l]][1:n,-1]
    Y = Ylist[[l]]
    cvfit = cv.glmnet(x=X, y=Y, alpha = 1, nfolds = 10)
    best_lambda = cvfit$lambda.min
    lasso_fit = glmnet(x=X, y=Y, alpha = 1, lambda = best_lambda)
    betalist[[l]] = as.vector(coef(lasso_fit))
  }
  return(betalist)
}

mulcvxr = function(L, dr){
  
  # This function uses CVXR to solve the optimization for prior rho.
  
  #dr = matrix(0, nrow = nrow(or), ncol = ncol(or))
  #for(l in 1:L){
  #  dr[l,] = Nlist[l] / n0 * or[l,] # odds ratio * r = density ratio
  #}
  
  ele = 1 / dr
  r = Variable(L)
  g = function(r){
    p = t(r) %*% ele
    return(mean(-log(p)))
  }
  
  # Adding a small penalty term to the objective function
  penalty = 1e-4 * sum_squares(r)
  
  EPS = 1e-8
  constraints = list(sum(r) == 1, r >= EPS, r <= 1 - EPS)
  
  problem = Problem(Minimize(g(r) + penalty), constraints = constraints)
  result = solve(problem)
  opt.sol = result$getValue(r)
  
  # cat("Optimal solution:\n", opt.sol)
  opt.sol
}

posterior = function(rho,dr,n0){
  
  # Calculate posterior eta from rho
  
  L =length(rho)
  post = matrix(data = NA,nrow = n0, ncol = L)
  for(l in 1:L){
    post[,l]= rho[l] / dr[l,] 
  }
  for(i in 1:n0){
    post[i,]=post[i,]/sum(post[i,])
  }
  return(post)
}

DRcoef_calculation = function(X0,Xlist,Ylist,nlist,post,postlist,drlist,
                            betalist,q){
  
  # Calculate the doubly robust coefficients
  
  L = length(Xlist)
  Nlist = sapply(Xlist, nrow)
  n0 = nrow(X0)
  firstPmat = matrix(0,nrow=L,ncol=q+1)
  secondPmat = matrix(0,nrow=L,ncol=q+1)
  Q = matrix(data = 0,nrow = q+1,ncol = L)
  
  for(l in 1:L){
    sum = rep(0, q+1)
    Ql2 = rep(0, q+1)
    res = Ylist[[l]] - Xlist[[l]][1:nlist[l],] %*% betalist[[l]]
    for(j in 1:nlist[l]){
      w = drlist[[l]][l,j] 
      r = res[j]
      po = postlist[[l]][j,l]
      temp = w * r * Xlist[[l]][j,1:(q+1)]
      Ql2 = Ql2 + temp
      sum = sum + po * temp
    }
    sum = (1/nlist[l]) * sum
    Ql2 = (1/nlist[l]) * Ql2
    secondPmat[l,] = sum
    
    sum = rep(0, q+1)
    Ql1 = rep(0, q+1)
    for(i in 1:n0){
      po = post[i,l]
      m = (betalist[[l]] %*% X0[i,])[1]
      temp = m * X0[i,1:(q+1)]
      Ql1 = Ql1 + temp
      sum = sum + po * temp
    }
    sum = (1/n0) * sum
    Ql1 = (1/n0) * Ql1
    firstPmat[l,] = sum
    Q[,l] = Ql1 + Ql2 
  }
  
  P = colSums(firstPmat) + colSums(secondPmat)
  
  out = list(preP = P,preQ = Q)
  return(out)
  
}

pseudo_solver = function(A,b,trc = 1e-10){
  
  # Solve for pseudo X and Y, which are used in high-dimensional glmnet
  
  p = nrow(A)
  if(ncol(A) != p){
    stop("X'X should be square matrix.")
  }
  if(length(b) != p){
    stop("X'y has wrong length.")
  }
  
  M = matrix(0,nrow = p+1, ncol = p+1)
  M[1:p,1:p] = A
  M[p+1,1:p] = b
  M[1:p,p+1] = b
  M[p+1,p+1] = b %*% b
  
  Evalue = eigen(M)$values
  Evalue[Evalue < trc] = 0
  
  #if(length(which(Evalue != 0)) > p){
  #  stop("Equation system unsolvable.")
  #}
  
  Dhalf = diag(sqrt(Evalue))
  U = eigen(M)$vectors
  
  C = U[,1:p] %*% Dhalf[1:p,1:p]
  C = t(C)
  
  pX = C[,1:p]
  pY = C[,p+1]
  
  out = list(pX = pX, pY= pY)
  return(out)
  
}

PQCalculation = function(X0,Xlist,Ylist,X0train,Xtrainlist,Ytrainlist,nlist,post,postlist,drlist,
                         betalist,q,penalty = FALSE,alpha = 0.5){
  
  ## This function calculates P(beta_mix) and Q(beta_l) in the algorithm with penalty.
  
  n0 = nrow(X0)
  n0t = nrow(X0train)
  Nlist = sapply(Xlist, nrow)
  L = length(Xlist)
  pre = DRcoef_calculation(X0,Xlist,Ylist,nlist,post,postlist,drlist,
                           betalist,q)
  preP = pre$preP
  preQ = pre$preQ
  Sigma = (1/n0) * t(X0[,1:(q+1)]) %*% X0[,1:(q+1)]
  
  if(penalty == TRUE){
  
    Q = matrix(data = 0,nrow = q+1,ncol = L)
    
    lenl = 40
    lambda_seq = 10^seq(-3, 0, length.out = lenl)
    DRloss_seq = rep(0,lenl)
    Sigmat = (1/n0t) * t(X0train[,1:(q+1)]) %*% X0train[,1:(q+1)]
    tpre = DRcoef_calculation(X0train,Xtrainlist,Ytrainlist,nlist,post,postlist,drlist,
                              betalist,q)
    tpreP = tpre$preP
    tpreQ = tpre$preQ
    
    for(l in 1:L){
      pse = pseudo_solver(Sigma,preQ[,l])
      pX = pse$pX
      pY = pse$pY
      
      for(k in 1:lenl){
        test_fit = glmnet(x = pX,y = pY, alpha = alpha, lambda = lambda_seq[k],intercept = FALSE)
        b = as.vector(coef(test_fit)[-1])
        DRloss_seq[k] = b %*% Sigmat %*% b - 2 * b %*% tpreQ[,l]
      }
      best_index = which.min(DRloss_seq)
      best_lambda = lambda_seq[best_index]
      eln_fit = glmnet(x= pX, y= pY, alpha = alpha, lambda = best_lambda, intercept = FALSE)
      Q[,l] = as.vector(coef(eln_fit)[-1])
    }
    
    pse = pseudo_solver(Sigma,preP)
    pX = pse$pX
    pY = pse$pY
    
    for(k in 1:lenl){
      test_fit = glmnet(x = pX,y = pY, alpha = alpha, lambda = lambda_seq[k],intercept = FALSE)
      b = as.vector(coef(test_fit)[-1])
      DRloss_seq[k] = b %*% Sigmat %*% b - 2 * b %*% tpreP
    }
    best_index = which.min(DRloss_seq)
    best_lambda = lambda_seq[best_index]
    eln_fit = glmnet(x= pX, y= pY, alpha = alpha, lambda = best_lambda, intercept = FALSE)
    P = as.vector(coef(eln_fit)[-1])
    
    out = list(P=P,Q=Q)
    
  }else{
    Sigma_inv = solve(Sigma)  
    P = Sigma_inv %*% preP
    Q = Sigma_inv %*% preQ
    out = list(P=P,Q=Q)
  }
  
  return(out)
}


opt_delta = function(P,Q,Sigma0,s){
  
  ## This function optimizes quadratic form of delta. (NOT USED IN NEWEST VERSION)
  
  L = ncol(Q)
  opt.weight=rep(NA, L)
  v = Variable(L)
  
  Gamma = s*t(Q) %*% Sigma0 %*% Q
  
  linterm = 2*(1-s)*t(P) %*% Sigma0 %*% Q
  
  Diag.matrix=diag(eigen(Gamma)$values)
  for(ind in 1:L){
    Diag.matrix[ind,ind]=max(Diag.matrix[ind,ind],1e-6)
  }
  
  Gamma.positive=eigen(Gamma)$vectors%*%Diag.matrix%*%t(eigen(Gamma)$vectors)
  objective = Minimize(quad_form(v,Gamma.positive) + linterm %*% v)
  
  constraints = list(sum(v)== 1, v>=0 )
  prob.weight= Problem(objective, constraints)
  result= solve(prob.weight)
  opt.status=result$status
  opt.sol=result$getValue(v)
  for(l in 1:L){
    opt.weight[l]=opt.sol[l]*(abs(opt.sol[l])>1e-8)
  }
  return(opt.weight)
}

opt_s_delta = function(P,Q,Sigma0,smax){
  
  ## This function optimizes quadratic form of delta and smax.
  
  L = ncol(Q)
  opt.weight=rep(NA, L+1)
  sd = Variable(L+1)
  
  QP = cbind(Q,P)
  Gamma = t(QP) %*% Sigma0 %*% QP
  
  Diag.matrix=diag(eigen(Gamma)$values)
  for(ind in 1:L){
    Diag.matrix[ind,ind]=max(Diag.matrix[ind,ind],1e-6)
  }
  Gamma.positive=eigen(Gamma)$vectors%*%Diag.matrix%*%t(eigen(Gamma)$vectors)
  eps = 1e-9
  constraints = list(
    sum(sd) == 1, 
    sd >= eps,       
    sd[L+1] >= 1-smax-eps
  )
  
  objective = Minimize(quad_form(sd,Gamma.positive))
  
  prob.weight= Problem(objective, constraints)
  result= solve(prob.weight)
  opt.status=result$status
  opt.sol=result$getValue(sd)
  for(l in 1:(L+1)){
    opt.weight[l]=opt.sol[l] * (abs(opt.sol[l])>1e-8)
  }
  opt.weight[L+1] = 1-opt.weight[L+1]
  opt.weight[1:L] = opt.weight[1:L]/opt.weight[L+1]
  return(opt.weight)
  
}

get_err = function(X,Y,beta,q,nn){
  # q is length(beta)-1
  errvec = Y[1:nn] - X[1:nn,1:(q+1)] %*% beta
  return( mean(errvec^2) )
}

extract_half_rows1 = function(mat) {
  half_rows = floor(nrow(mat) / 2)
  return(mat[1:half_rows, ])
}

extract_half_rows2 = function(mat) {
  half_rows = floor(nrow(mat) / 2)
  return(mat[(half_rows+1):nrow(mat), ])
}

maximin_s_beta = function(Xlist,Xtrainlist,Ylist,Ytrainlist,X0,X0train,nlist,q,smax,penalty = FALSE, alpha = 0.5,
                          rho_pseudo = "NA", dr_type = "rf",normalize = FALSE){
  
  # This is the main algorithm of calculating REMIX beta.
  
  L = length(Xlist)
  n0 = nrow(X0)
  Nlist = sapply(Xlist, nrow)
  lam = n0 ** (-1/3)
  
  ## Cross-fitting: using auxiliary data to train nuisance models
  
  if(dr_type == 'logit'){
    modellist = or_estimation_logit(Xtrainlist,X0train) # density ratio model fit by aux data
    dr = dratio_logit(modellist,X0,L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
    drlist = list()
    for(l in 1:L){
      drlist[[l]] = dratio_logit(modellist,Xlist[[l]],L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
    }
  }else{
    modellist = or_estimation_ML(Xtrainlist,X0train,dr_type = dr_type) # density ratio model fit by aux data
    dr = dratio_ML(modellist,X0,L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
    drlist = list()
    for(l in 1:L){
      drlist[[l]] = dratio_ML(modellist,Xlist[[l]],L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
    }
  }
  
  if(rho_pseudo == 'max'){
    
    maxind = which.max(Nlist)
    cat('dr_pseudo choose site with max sample size',maxind,'\n')
    Xtrainmax = extract_half_rows2(Xtrainlist[[maxind]])
    Xtrainlist_pseudo = Xtrainlist
    Xtrainlist_pseudo[[maxind]] = extract_half_rows1(Xtrainlist[[maxind]])
    
    if(dr_type == 'logit'){
      modellist = or_estimation_logit(Xtrainlist_pseudo,Xtrainmax) # density ratio model fit by aux data
      dr_new = dratio_logit(modellist,X0,L,Nlist = sapply(Xtrainlist_pseudo,nrow),n0 = nrow(Xtrainmax),normalize = normalize)
    }else{
      modellist = or_estimation_ML(Xtrainlist_pseudo,Xtrainmax,dr_type = dr_type) # density ratio model fit by aux data
      dr_new = dratio_ML(modellist,X0,L,Nlist = sapply(Xtrainlist_pseudo,nrow),n0 = nrow(Xtrainmax),normalize = normalize)
    }
    
    rho = mulcvxr(L,dr_new)
    
  }else if(rho_pseudo == 'pool'){

    cat('dr_pseudo use pooled source training data','\n')
    half_Xtrainlist = lapply(Xtrainlist, extract_half_rows1)
    lasthalf_Xtrainlist = lapply(Xtrainlist,extract_half_rows2)
    Xpool = do.call(rbind,lasthalf_Xtrainlist)
    
    if(dr_type == 'logit'){
      modellist = or_estimation_logit(half_Xtrainlist,Xpool) # density ratio model fit by aux data
      dr_new = dratio_logit(modellist,X0,L,Nlist = sapply(half_Xtrainlist,nrow),n0 = nrow(Xpool),normalize = normalize)
    }else{
      modellist = or_estimation_ML(half_Xtrainlist,Xpool,dr_type = dr_type) # density ratio model fit by aux data
      dr_new = dratio_ML(modellist,X0,L,Nlist = sapply(half_Xtrainlist,nrow),n0 = nrow(Xpool),normalize = normalize)
    }
    
    rho = mulcvxr(L,dr_new)
    
  }else{
    # use true density ratio against target to get rho
    rho = mulcvxr(L,dr)
  }
  
  betalist = lasso_imputation(Xtrainlist,Ytrainlist,nlist) # imputation model fit by aux data
  post = posterior(rho,dr,n0)
  postlist = list()
  for(l in 1:L){
    postlist[[l]] = posterior(rho,drlist[[l]],Nlist[[l]])
  }
  
  Sigma0 = (1/n0) * t(X0[,1:(q+1)]) %*% X0[,1:(q+1)]
  
  PQ = PQCalculation(X0,Xlist,Ylist,X0train,Xtrainlist,Ytrainlist,nlist,post,postlist,
                     drlist,betalist,q,penalty,alpha)
  P = PQ$P
  Q = PQ$Q
  
  dts_star = opt_s_delta(P,Q,Sigma0,smax)
  s = dts_star[L+1]
  opt.beta = (1-s)*P + s* Q %*% dts_star[1:L]
  # cat("\n delta* is",dts_star[1:L]," and s* is",s)
  
  dt_MI = opt_delta(P,Q,Sigma0,1)
  beta_MI = Q %*% dt_MI
  # cat("\n delta_MI is",dt_MI,"\n")
  
  # convert odds ratio to density ratio for the outcome dr.
  #dr = matrix(0,nrow = nrow(or),ncol=ncol(or))
  #for(l in 1:L){
  #  dr[l,] = Nlist[l]/n0 * or[l,] # odds ratio * r = density ratio
  #}
  
  out = list(beta_star = opt.beta,
             beta_MI = beta_MI,
             beta_RAP = P,
             DoublyR = Q,
             deltas = dts_star,
             rho = rho,
             dr = dr)
  return(out)
}

get_REMIX_beta = function(Xlist,Xtrainlist,Ylist,Ytrainlist,X0,X0train,nlist,q,smax,penalty = FALSE, alpha = 0.5,
                          rho_pseudo = "NA", dr_type = "rf",normalize = FALSE){
  
  # Use sample splitting: beta = 0.5 (betaA + betaB)
  
  out1 = maximin_s_beta(Xlist,Xtrainlist,Ylist,Ytrainlist,X0,X0train,nlist,q,smax,penalty,alpha,rho_pseudo,dr_type,normalize = normalize)
  out2 = maximin_s_beta(Xtrainlist,Xlist,Ytrainlist,Ylist,X0train,X0,nlist,q,smax,penalty,alpha,rho_pseudo,dr_type,normalize = normalize)
  opt.beta = 0.5 * (out1$beta_star + out2$beta_star)
  beta_MI = 0.5 * (out1$beta_MI + out2$beta_MI)
  P = 0.5 * (out1$beta_RAP + out2$beta_RAP)
  Q = 0.5 * (out1$DoublyR + out2$DoublyR)
  dts_star = 0.5 * (out1$deltas + out2$deltas)
  rho = 0.5 * (out1$rho + out2$rho)
  # dr = 0.5 * (out1$dr + out2$dr)
  dr = out1$dr
  
  out = list(beta_star = opt.beta,
             beta_MI = beta_MI,
             beta_RAP = P,
             DoublyR = Q,
             rho = rho,
             deltas = dts_star,
             DR_dP0bydPl = dr)
  return(out)
}





