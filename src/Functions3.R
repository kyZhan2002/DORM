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

or_estimation_ML = function(Xlist,X0,dr_type = 'rf', report = TRUE, condA = FALSE, q = NULL){
  
  ## This function estimates the odds ratio using random forest.
  
  L = length(Xlist)
  Nlist = sapply(Xlist, nrow)
  n0 = nrow(X0)
  modellist = list()
  
  # If condA is TRUE, also fit models for A only
  if(condA){
    if(is.null(q)){
      stop("q must be provided when condA = TRUE")
    }
    modellist_A = list()
  }
  
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
      
      # Fit model for A only if condA = TRUE
      if(condA){
        X_A = X[, 1:q]
        tune_grid_A = expand.grid(
          mtry = floor(sqrt(ncol(X_A)))
        )
        rfmodel_A = train(
          x = X_A, y = Y, 
          method = "rf", 
          trControl = train_control, 
          tuneGrid = tune_grid_A,
        )
        modellist_A[[l]] = rfmodel_A
      }
      
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
      
      # Fit model for A only if condA = TRUE
      if(condA){
        X_A = X[, 1:q]
        nnet_model_A = train(
          x = X_A, y = Y,
          method = "nnet",
          trControl = train_control,
          tuneGrid = tune_grid,
          trace = FALSE,                 
          MaxNWts = 8000,                
        )
        modellist_A[[l]] = nnet_model_A
      }
      
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
      
      # Fit model for A only if condA = TRUE
      if(condA){
        X_A = X[, 1:q]
        xgb_model_A = train(
          x = X_A, y = Y,
          method = "xgbTree",
          trControl = train_control,
          tuneGrid = tune_grid,
          verbosity = 0
        )
        modellist_A[[l]] = xgb_model_A
      }
    }else{
      stop("Wrong dr_type.")
    }
    
    if(report){
      # Print the best parameters
      print(paste("Best parameters for model", l, ":", sep=" "))
      print(modellist[[l]]$bestTune)
    }
  }
  
  if(condA){
    return(list(modellist = modellist, modellist_A = modellist_A))
  }else{
    return(modellist)
  }
} 

or_estimation_logit = function(Xlist,Xtar,alpha_dr = 0.5, condA = FALSE, q = NULL){
  
  ## This function estimates logistic model for odds ratio = exp(x'b)
  
  L = length(Xlist)
  modellist = list()
  
  # If condA is TRUE, also fit models for A only
  if(condA){
    if(is.null(q)){
      stop("q must be provided when condA = TRUE")
    }
    modellist_A = list()
  }
  
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
    
    # Fit model for A only if condA = TRUE
    if(condA){
      X_A = X[, 1:q]
      cv_fit_A = cv.glmnet(X_A, y, alpha = alpha_dr, family = "binomial")
      best_lambda_A = cv_fit_A$lambda.min
      final_model_A = glmnet(X_A, y, alpha = alpha_dr, lambda = best_lambda_A, family = "binomial")
      modellist_A[[l]] = final_model_A
    }
  }
  
  if(condA){
    return(list(modellist = modellist, modellist_A = modellist_A))
  }else{
    return(modellist)
  }
}


dratio_ML = function(modellist,X0,L,Nlist,n0,bound = 1e6,normalize = FALSE, condA = FALSE, modellist_A = NULL, q = NULL){
  
  ## The function output a L*n0 density ratio matrix, each row is density ratio:dP(target X0) to dP(source)
  ## We multiply nl/n0 here!!! This is density ratio
  
  X0 = X0[,-1] # Remove interception
  X0 = as.data.frame(X0)  # Convert once outside the loop
  
  dr = matrix(0,nrow=L,ncol=nrow(X0))
  for(l in 1:L){
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
  
  # If condA = TRUE, compute conditional density ratio AND return both
  if(condA){
    if(is.null(modellist_A)){
      stop("modellist_A must be provided when condA = TRUE")
    }
    if(is.null(q)){
      stop("q must be provided when condA = TRUE for dratio_ML")
    }
    
    dr_A = matrix(0,nrow=L,ncol=nrow(X0))
    X0_A = X0[, 1:q]
    
    for(l in 1:L){
      pred_A = predict(modellist_A[[l]],X0_A,type = "prob")
      ratio_A = ifelse(pred_A[,2] == 0, bound, pred_A[,1] / pred_A[,2])
      ratio_A[ratio_A==0] = 1/bound
      ratio_A = pmax(1/bound, pmin(ratio_A, bound))
      dr_A[l,] = Nlist[l] / n0 * ratio_A
    }
    
    # Conditional density ratio: dr(X) / dr(A)
    dr_cond = dr / dr_A
    
    if(normalize){
      for(l in 1:L){
        nor_const = mean(1/dr_cond[l,])
        dr_cond[l,] = nor_const * dr_cond[l,]
      }
    }
    
    # Return both unconditional and conditional
    return(list(dr_uncond = dr, dr_cond = dr_cond))
  }
  
  return(dr)
}


dratio_logit = function(modellist,X0,L,Nlist,n0,bound=1e6,normalize = FALSE, condA = FALSE, modellist_A = NULL, q = NULL){
  
  # use logistic regression to get density ratio!
  
  dr = matrix(0,nrow=L,ncol=nrow(X0))
  for(l in 1:L){
    model = modellist[[l]]
    new_X = X0[,-1] ## remove the intercept which already exist in X.
    log_density_ratio = predict(model, new_X, type = "link")
    # exp(log_density_ratio) = P(Y=1|X)/P(Y=0|X) = dP_source/dP_target
    # We want dP_target/dP_source, so we take 1/exp(log_density_ratio)
    density_ratio = 1 / exp(log_density_ratio)
    density_ratio = pmax(1/bound, pmin(density_ratio, bound))
    dr[l,] = Nlist[l]/n0 * density_ratio
    
    # mean of (1/dr[l,]) should be 1.
    if(normalize){
      nor_const = mean(1/dr[l,])
      dr[l,] = nor_const * dr[l,]
    }
  }
  
  # If condA = TRUE, compute conditional density ratio AND return both
  if(condA){
    if(is.null(modellist_A)){
      stop("modellist_A must be provided when condA = TRUE")
    }
    if(is.null(q)){
      stop("q must be provided when condA = TRUE for dratio_logit")
    }
    
    dr_A = matrix(0,nrow=L,ncol=nrow(X0))
    new_X_A = new_X[, 1:q]
    
    for(l in 1:L){
      model_A = modellist_A[[l]]
      log_density_ratio_A = predict(model_A, new_X_A, type = "link")
      density_ratio_A = 1 / exp(log_density_ratio_A)
      density_ratio_A = pmax(1/bound, pmin(density_ratio_A, bound))
      dr_A[l,] = Nlist[l]/n0 * density_ratio_A
    }
    
    # Conditional density ratio: dr(X) / dr(A)
    dr_cond = dr / dr_A
    
    if(normalize){
      for(l in 1:L){
        nor_const = mean(1/dr_cond[l,])
        dr_cond[l,] = nor_const * dr_cond[l,]
      }
    }
    
    # Return both unconditional and conditional
    return(list(dr_uncond = dr, dr_cond = dr_cond))
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
  
  # Calculate posterior eta from rho (vectorized version)
  
  L = length(rho)
  # Vectorized computation: rho[l] / dr[l,] for all l at once
  # dr is L x n, we want to multiply row l by rho[l]
  post = sweep(1/dr, 1, rho, "*")  # MARGIN=1 to multiply rows by rho
  post = t(post)  # transpose to get n0 x L
  
  # Normalize rows to sum to 1
  post = post / rowSums(post)
  
  return(post)
}

DRcoef_calculation = function(X0,Xlist,Ylist,nlist,post,postlist,drlist,
                              betalist,q){
  
  # Calculate the doubly robust coefficients (vectorized version)
  
  L = length(Xlist)
  Nlist = sapply(Xlist, nrow)
  n0 = nrow(X0)
  firstPmat = matrix(0,nrow=L,ncol=q+1)
  secondPmat = matrix(0,nrow=L,ncol=q+1)
  Q = matrix(data = 0,nrow = q+1,ncol = L)
  
  for(l in 1:L){
    # Vectorized computation for source data
    n = nlist[l]
    X_source = Xlist[[l]][1:n, 1:(q+1)]
    res = Ylist[[l]] - Xlist[[l]][1:n,] %*% betalist[[l]]
    w = drlist[[l]][l, 1:n]
    po_source = postlist[[l]][1:n, l]
    
    # Vectorized outer product sum: sum_j w_j * r_j * po_j * X_j
    # Convert res to vector and use proper vectorization
    res_vec = as.vector(res)
    weighted_factor = w * res_vec
    weighted_res = sweep(X_source, 1, weighted_factor, "*")  # multiply each row of X_source
    Ql2 = colSums(weighted_res) / n
    
    weighted_factor_po = po_source * w * res_vec
    secondPmat[l,] = colSums(sweep(X_source, 1, weighted_factor_po, "*")) / n
    
    # Vectorized computation for target data
    X_target = X0[, 1:(q+1)]
    m = X0 %*% betalist[[l]]
    po_target = post[, l]
    
    m_vec = as.vector(m)
    weighted_m = sweep(X_target, 1, m_vec, "*")
    Ql1 = colSums(weighted_m) / n0
    firstPmat[l,] = colSums(sweep(X_target, 1, po_target * m_vec, "*")) / n0
    
    Q[,l] = Ql1 + Ql2 
  }
  
  P = colSums(firstPmat) + colSums(secondPmat)
  
  out = list(preP = P, preQ = Q)
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
  M[p+1,p+1] = sum(b^2)  # More efficient than b %*% b
  
  eig = eigen(M, symmetric = TRUE)
  Evalue = eig$values
  Evalue[Evalue < trc] = 0
  U = eig$vectors
  
  # Only compute for non-zero eigenvalues
  nonzero_idx = which(Evalue > 0)
  if(length(nonzero_idx) > p){
    nonzero_idx = nonzero_idx[1:p]
  }
  
  Dhalf = sqrt(Evalue[nonzero_idx])
  C = t(U[, nonzero_idx] %*% diag(Dhalf))
  
  pX = C[, 1:p]
  pY = C[, p+1]
  
  out = list(pX = pX, pY = pY)
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
    
    lenl = 30
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
  
  ## This function optimizes quadratic form of delta.
  
  L = ncol(Q)
  opt.weight=rep(NA, L)
  v = Variable(L)
  
  Gamma = s*t(Q) %*% Sigma0 %*% Q
  
  linterm = 2*(1-s)*t(P) %*% Sigma0 %*% Q
  
  # Make Gamma positive definite more efficiently
  eig = eigen(Gamma, symmetric = TRUE)
  eig$values = pmax(eig$values, 1e-6)
  Gamma.positive = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  objective = Minimize(quad_form(v,Gamma.positive) + linterm %*% v)
  
  constraints = list(sum(v)== 1, v>=0 )
  prob.weight= Problem(objective, constraints)
  result= solve(prob.weight)
  opt.status=result$status
  opt.sol=result$getValue(v)
  opt.weight = ifelse(abs(opt.sol) > 1e-8, opt.sol, 0)
  return(opt.weight)
}

opt_s_delta = function(P,Q,Sigma0,smax){
  
  ## This function optimizes quadratic form of delta and smax.
  
  L = ncol(Q)
  opt.weight=rep(NA, L+1)
  sd = Variable(L+1)
  
  QP = cbind(Q,P)
  Gamma = t(QP) %*% Sigma0 %*% QP
  
  # Make Gamma positive definite more efficiently
  eig = eigen(Gamma, symmetric = TRUE)
  eig$values = pmax(eig$values, 1e-6)
  Gamma.positive = eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
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
  opt.weight = ifelse(abs(opt.sol) > 1e-8, opt.sol, 0)
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
                          rho_pseudo = "NA", dr_type = "rf",normalize = FALSE, condA = FALSE){
  
  # This is the main algorithm of calculating DORM beta.
  
  L = length(Xlist)
  n0 = nrow(X0)
  Nlist = sapply(Xlist, nrow)
  lam = n0 ** (-1/3)
  
  ## Cross-fitting: using auxiliary data to train nuisance models
  
  if(dr_type == 'logit'){
    models = or_estimation_logit(Xtrainlist,X0train, condA = condA, q = q)
    if(condA){
      modellist = models$modellist
      modellist_A = models$modellist_A
      # Compute BOTH unconditional and conditional density ratios in one call
      dr_result = dratio_logit(modellist,X0,L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
      dr_uncond = dr_result$dr_uncond
      dr_cond = dr_result$dr_cond
      
      drlist_uncond = list()
      drlist_cond = list()
      for(l in 1:L){
        dr_result_l = dratio_logit(modellist,Xlist[[l]],L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
        drlist_uncond[[l]] = dr_result_l$dr_uncond
        drlist_cond[[l]] = dr_result_l$dr_cond
      }
    }else{
      modellist = models
      dr_uncond = dratio_logit(modellist,X0,L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
      dr_cond = dr_uncond  # Same as unconditional when condA = FALSE
      
      drlist_uncond = list()
      drlist_cond = list()
      for(l in 1:L){
        drlist_uncond[[l]] = dratio_logit(modellist,Xlist[[l]],L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
        drlist_cond[[l]] = drlist_uncond[[l]]
      }
    }
  }else{
    models = or_estimation_ML(Xtrainlist,X0train,dr_type = dr_type, condA = condA, q = q)
    if(condA){
      modellist = models$modellist
      modellist_A = models$modellist_A
      # Compute BOTH unconditional and conditional density ratios in one call
      dr_result = dratio_ML(modellist,X0,L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
      dr_uncond = dr_result$dr_uncond
      dr_cond = dr_result$dr_cond
      
      drlist_uncond = list()
      drlist_cond = list()
      for(l in 1:L){
        dr_result_l = dratio_ML(modellist,Xlist[[l]],L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
        drlist_uncond[[l]] = dr_result_l$dr_uncond
        drlist_cond[[l]] = dr_result_l$dr_cond
      }
    }else{
      modellist = models
      dr_uncond = dratio_ML(modellist,X0,L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
      dr_cond = dr_uncond  # Same as unconditional when condA = FALSE
      
      drlist_uncond = list()
      drlist_cond = list()
      for(l in 1:L){
        drlist_uncond[[l]] = dratio_ML(modellist,Xlist[[l]],L,Nlist = sapply(Xtrainlist,nrow),n0 = nrow(X0train),normalize = normalize)
        drlist_cond[[l]] = drlist_uncond[[l]]
      }
    }
  }
  
  if(rho_pseudo == 'max'){
    
    maxind = which.max(Nlist)
    cat('dr_pseudo choose site with max sample size',maxind,'\n')
    Xtrainmax = extract_half_rows2(Xtrainlist[[maxind]])
    Xtrainlist_pseudo = Xtrainlist
    Xtrainlist_pseudo[[maxind]] = extract_half_rows1(Xtrainlist[[maxind]])
    
    if(dr_type == 'logit'){
      modellist_new = or_estimation_logit(Xtrainlist_pseudo,Xtrainmax, condA = condA, q = q)
      if(condA){
        modellist = modellist_new$modellist
        modellist_A = modellist_new$modellist_A
        dr_result_new = dratio_logit(modellist,X0,L,Nlist = sapply(Xtrainlist_pseudo,nrow),n0 = nrow(Xtrainmax),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
        dr_new = dr_result_new$dr_cond  # Use conditional for rho estimation
      }else{
        modellist = modellist_new
        dr_new = dratio_logit(modellist,X0,L,Nlist = sapply(Xtrainlist_pseudo,nrow),n0 = nrow(Xtrainmax),normalize = normalize)
      }
    }else{
      modellist_new = or_estimation_ML(Xtrainlist_pseudo,Xtrainmax,dr_type = dr_type, condA = condA, q = q)
      if(condA){
        modellist = modellist_new$modellist
        modellist_A = modellist_new$modellist_A
        dr_result_new = dratio_ML(modellist,X0,L,Nlist = sapply(Xtrainlist_pseudo,nrow),n0 = nrow(Xtrainmax),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
        dr_new = dr_result_new$dr_cond  # Use conditional for rho estimation
      }else{
        modellist = modellist_new
        dr_new = dratio_ML(modellist,X0,L,Nlist = sapply(Xtrainlist_pseudo,nrow),n0 = nrow(Xtrainmax),normalize = normalize)
      }
    }
    
    rho = mulcvxr(L,dr_new)
    
  }else if(rho_pseudo == 'pool'){
    
    cat('dr_pseudo use pooled source training data','\n')
    half_Xtrainlist = lapply(Xtrainlist, extract_half_rows1)
    lasthalf_Xtrainlist = lapply(Xtrainlist,extract_half_rows2)
    Xpool = do.call(rbind,lasthalf_Xtrainlist)
    
    if(dr_type == 'logit'){
      modellist_new = or_estimation_logit(half_Xtrainlist,Xpool, condA = condA, q = q)
      if(condA){
        modellist = modellist_new$modellist
        modellist_A = modellist_new$modellist_A
        dr_result_new = dratio_logit(modellist,X0,L,Nlist = sapply(half_Xtrainlist,nrow),n0 = nrow(Xpool),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
        dr_new = dr_result_new$dr_cond  # Use conditional for rho estimation
      }else{
        modellist = modellist_new
        dr_new = dratio_logit(modellist,X0,L,Nlist = sapply(half_Xtrainlist,nrow),n0 = nrow(Xpool),normalize = normalize)
      }
    }else{
      modellist_new = or_estimation_ML(half_Xtrainlist,Xpool,dr_type = dr_type, condA = condA, q = q)
      if(condA){
        modellist = modellist_new$modellist
        modellist_A = modellist_new$modellist_A
        dr_result_new = dratio_ML(modellist,X0,L,Nlist = sapply(half_Xtrainlist,nrow),n0 = nrow(Xpool),normalize = normalize, condA = TRUE, modellist_A = modellist_A, q = q)
        dr_new = dr_result_new$dr_cond  # Use conditional for rho estimation
      }else{
        modellist = modellist_new
        dr_new = dratio_ML(modellist,X0,L,Nlist = sapply(half_Xtrainlist,nrow),n0 = nrow(Xpool),normalize = normalize)
      }
    }
    
    rho = mulcvxr(L,dr_new)
    
  }else{
    # use density ratio against target to get rho
    # Use conditional dr if condA = TRUE, otherwise unconditional
    rho = mulcvxr(L,dr_cond)
  }
  
  betalist = lasso_imputation(Xtrainlist,Ytrainlist,nlist) # imputation model fit by aux data
  
  # Use conditional dr for posterior if condA = TRUE
  post = posterior(rho,dr_cond,n0)
  postlist = list()
  for(l in 1:L){
    postlist[[l]] = posterior(rho,drlist_cond[[l]],Nlist[[l]])
  }
  
  Sigma0 = (1/n0) * t(X0[,1:(q+1)]) %*% X0[,1:(q+1)]
  
  # Always use unconditional dr for DRcoef and PQ calculation
  PQ = PQCalculation(X0,Xlist,Ylist,X0train,Xtrainlist,Ytrainlist,nlist,post,postlist,
                     drlist_uncond,betalist,q,penalty,alpha)
  P = PQ$P
  Q = PQ$Q
  
  dts_star = opt_s_delta(P,Q,Sigma0,smax)
  s = dts_star[L+1]
  opt.beta = (1-s)*P + s* Q %*% dts_star[1:L]
  # cat("\n delta* is",dts_star[1:L]," and s* is",s)
  
  dt_MI = opt_delta(P,Q,Sigma0,1)
  beta_MI = Q %*% dt_MI
  # cat("\n delta_MI is",dt_MI,"\n")
  
  out = list(beta_star = opt.beta,
             beta_MI = beta_MI,
             beta_RAP = P,
             DoublyR = Q,
             deltas = dts_star,
             rho = rho,
             dr = dr_uncond)  # Return unconditional dr
  return(out)
}

get_DORM_beta = function(Xlist,Xtrainlist,Ylist,Ytrainlist,X0,X0train,nlist,q,smax,penalty = FALSE, alpha = 0.5,
                         rho_pseudo = "NA", dr_type = "rf",normalize = FALSE, condA = FALSE){
  
  # Use sample splitting: beta = 0.5 (betaA + betaB)
  
  cat("Fitting density ratio model using",dr_type,"with conditional on A=",condA)
  if(penalty){
    cat("Fitting high dimensional DORM using elastic net, alpha=",alpha)
  }
  
  out1 = maximin_s_beta(Xlist,Xtrainlist,Ylist,Ytrainlist,X0,X0train,nlist,q,smax,penalty,alpha,rho_pseudo,dr_type,normalize = normalize, condA = condA)
  out2 = maximin_s_beta(Xtrainlist,Xlist,Ytrainlist,Ylist,X0train,X0,nlist,q,smax,penalty,alpha,rho_pseudo,dr_type,normalize = normalize, condA = condA)
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
