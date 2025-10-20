source('/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/TransLG-functions.R')
source('/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/DataGeneration_contam.R')
source('/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/Functions3.R')
source('/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/Tuning.R')
source('/n/data1/hsph/biostat/celehs/lab/kez641/REMIX/SourceCode/Evaluation.R')

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# The first argument is the random_seed, the second argument is r
seed <- as.numeric(args[1])
r <- as.numeric(args[2])

# Now you can use `random_seed` and `r` in your R code
print(paste("Random seed:", seed))
print(paste("r value:", r))

set.seed(seed)

output_dir <- "results/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


#########################
# r = 0.5
mixture = c(r,0,1-r,0,0)    # How target is mixed by sources. Here 1/2 from site 1 and 1/2 from site 3.
# mixture from r = 0 to r = 1 to give a complete profile.
#########################

L = 5     # number of source sites
q = 4     # dim of low-d A = q+1
p = 195   # rest dim of covariate X. Total dim of X is (1+q+p).
ptemp = 6

Nlist = rep(1000,5)      # Source sample size of X. 
nlist = rep(500,5)      # Source sample size of Y.
n0 = 500                # Target sample size of X0.
# You can set these sizes larger in real simulation.

trues = 0.1                   # Probability that target use another mixture delta.
delta = c(0,0,1,0,0)

# Each covariate X = (1,A,W). First few elements of W is correlated with A, while the rest is noise.
MU_A = 0.5 * matrix(c(-2,0,1,0,
                      0,-1,2,1,
                      -1,2,0,1,
                      1,1,-1,2,
                      2,1,1,2),nrow=L,ncol=q,byrow=TRUE) # Each row is the mean of A_l

# True setting of BETA. Y = X %*% BETA + eps.
BETA =     1  *  matrix(c(0,4,4,3,-3,  0.5,-0.2,0,0,1,0,
                          -3,3,-1,3,-3, 0,0.5,-0.2,0,0,1,
                          -3,-3,0,-3,-3,  0,0,0.2,-0.5,0,1,
                          -6,6,6,6,-6,  1,0,0,0.2,-0.5,0,
                          -3,3,0,0,1, 0,0,0,1,0.2,-0.5),nrow=L,ncol=ptemp+q+1,byrow=TRUE)


### First part: transLASSO and transGLM
q=4
Data = Generate_Simulation_Data_contamination(L,p,q,Nlist,nlist,n0,mixture,delta,trues,MU_A,BETA,ptemp,
                                eps_A=0.5,eps_W=0.01,eps_Y = 0.1)
datte = Generate_test_data_contamination(L,p,q,n0,mixture,delta,trues,MU_A,BETA,ptemp,
                           eps_Y = 0.1,eps_A = 0.5,eps_W =0.01)


## Trans-LASSO

q = 199
ntar_arr = seq(from = 50,to = 400,by = 50)
ntar_num = length(ntar_arr)
df_TL <- data.frame(ntar = integer(), mse_TL = numeric(),mse_test = numeric())
beta_TL = matrix(0,nrow = ntar_num,ncol = q+1)

for (i in 1:(ntar_num)) {
  ntar = ntar_arr[i]
  TL = mytranslasso(Data$X0,Data$Xlist,Data$Y0,Data$Ylist,ntar,q,Inumd=3)
  mse_test = get_err(datte$X0,datte$Y0,TL$beta_TL,q,length(datte$Y0))
  df_TL <- rbind(df_TL, data.frame(ntar = ntar, mse_TL = TL$mse_TL,mse_test = mse_test))
  beta_TL[i,] = TL$beta_TL
}

df_TL
beta_TL

# Trans-GLM

df_TG <- data.frame(ntar = integer(), mse_TG = numeric(),mse_test = numeric())
beta_TG = matrix(0,nrow = ntar_num,ncol = q+1)

for (i in 1:(ntar_num)) {
  ntar = ntar_arr[i]
  
  D.training = List2TransGLM(Data$X0,Data$Xlist,Data$Y0,Data$Ylist,ntar,q)
  
  fit.gaussian <- glmtrans(D.training$target, D.training$source,family = 'gaussian')
  y.pred.glmtrans <- predict(fit.gaussian, Data$X0)
  
  mse_TG = mean((Data$Y0 - y.pred.glmtrans)^2)
  beta.glm = as.vector(fit.gaussian$beta)[-2]
  mse_test = get_err(datte$X0,datte$Y0,beta.glm,q,length(datte$Y0))
  
  df_TG <- rbind(df_TG, data.frame(ntar = ntar, mse_TG = mse_TG,mse_test = mse_test))
  beta_TG[i,] = beta.glm
  
}

df_TG
beta_TG

# PTL

df_PTL <- data.frame(ntar = integer(), mse_PTL = numeric(),mse_test = numeric())
beta_PTL = matrix(0,nrow = ntar_num,ncol = q+1)

for (i in 1:(ntar_num)) {
  ntar = ntar_arr[i]
  thePTL = PTL_algorithm(Data$Xlist,Data$Ylist,Data$X0,Data$Y0,ntar)
  mse_test = get_err(datte$X0,datte$Y0,thePTL$beta_PTL,q,length(datte$Y0))
  df_PTL <- rbind(df_PTL, data.frame(ntar = ntar, mse_PTL = thePTL$mse_PTL ,mse_test = mse_test))
  beta_PTL[i,] = thePTL$beta_PTL
}

df_PTL
beta_PTL


# Second part: our method

q = 199
## Implement parameter tuning. We need to put all beta(smax) in a matrix.
smax_array = seq(0.05,0.4,0.05)     # smax array to be tuned over.
nsmax = length(smax_array)
beta_array = matrix(0,nrow = nsmax, ncol = q+1)     # beta(smax), we put them in the corresponding rows.
result_list = list()
df_list1 = list()
df_list2 = list()

for(i in 1:nsmax){
  result = get_REMIX_beta(Data$Xlist,Data$Xtrainlist,Data$Ylist,Data$Ytrainlist,Data$X0,Data$X0train,nlist,q,smax_array[i],
                          penalty = TRUE,alpha = 0.5,rho_pseudo = 'pool', dr_type = "rf")
  beta_array[i,] = result$beta_star
  smax = smax_array[i]
  cat("smax =",smax_array[i],'\n')
  result_list[[as.character(smax)]] <- result
  
  df_list1[[as.character(smax)]] = benchmarks(Data$X0,Data$Y0,
                                              q,result$beta_star,result$beta_MI,result$rho,result$beta_RAP,result$DoublyR)
  df_list2[[as.character(smax)]] = benchmarks(datte$X0,datte$Y0,
                                              q,result$beta_star,result$beta_MI,result$rho,result$beta_RAP,result$DoublyR)
}

df_list1
df_list2
## Tuning: expect some of them choose same smax:
result_list
beta_array

ty_list = list()
#ts_list = list()
for(i in 1:ntar_num){
  nn = ntar_arr[i]
  ty_list[i] = tuning_Y(Data$X0[1:nn,],Data$Y0[1:nn],q,beta_array,smax_array)$best_smax
  #ts_list[i] = tuning_surrogate(Data$X0[1:nn,],Data$S0[1:nn],q,beta_array,smax_array)$best_smax
}

ty_list
#ts_list


save(df_TL, beta_TL, df_TG, beta_TG, df_PTL,beta_PTL,df_list1,
     df_list2,result_list,beta_array,ty_list,file = paste0(output_dir, "Contamination_r",r,"i",seed, ".RData"))

