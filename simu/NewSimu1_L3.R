## This R script runs over different L

source('DataGeneration3.R')
source('Functions3.R')
source('Tuning.R')

## Some parameter setting:

q = 4     # dim of low-d A
p = 195   # rest dim of covariate X. Total dim of X is (1+q+p).
ptemp = 6

trues = 0

# True setting of BETA. Y = X %*% BETA + eps.
smax_array = seq(0.05,0.5,0.05)  # CHANGE!
nsmax = length(smax_array)

# Setting for different L

L = 3     # number of source sites

delta = c(0.5,0.5,0)

mixture = c(0.5,0.5,0)  

MU_A = 1.0 * matrix(c(-2,0,1,0,
                      0,-1,2,1,
                      -1,2,0,1),nrow=L,ncol=q,byrow=TRUE) # Each row is the mean of A_l
# Each covariate X = (1,A,W). First few elements of W is correlated with A, while the rest is noise.

BETA = 0.8 * matrix(c(0,4,4,3,-3,  0.5,-0.2,0,0,1,0,
                     -3,3,-1,3,-3, 0,0.5,-0.2,0,0,1,
                     -5,5,5,5,-5,  1,0,0,0.2,-0.5,0),nrow=L,ncol=ptemp+q+1,byrow=TRUE)

Nlist = rep(1000,L)      # Source sample size of X. 
nlist = rep(500,L)      # Source sample size of Y.
n0 = 400                # Target sample size of X0.

beta_array1 = matrix(0,nrow = nsmax, ncol = q+1)
result_list1 = list()

for(i in 1:nsmax){
  
  # randomly choose trues and truedelta to generate data, which does not effect betas:
  # Because trues and truedelta only affect Y!
  dat = Generate_Simulation_Data(L,p,q,Nlist,nlist,n0,mixture,delta,trues,MU_A,BETA)
  smax = smax_array[i]
  result = get_REMIX_beta(dat$Xlist,dat$Xtrainlist,dat$Ylist,dat$Ytrainlist,
                          dat$X0,dat$X0train,nlist,q,smax,rho_pseudo = 'pool',dr_type = 'rf',normalize = FALSE)
  
  beta_array1[i,] = result$beta_star
  
  cat("smax =",smax,'\n')
  
  result_list1[[as.character(smax)]] <- result
}

result_list1
beta_array1

smaxn = names(result_list1)
smax_array = as.numeric(smaxn)
nsmax = length(smaxn)

trues_array = seq(0,1,0.05)     # CHANGE!!
ntrues = length(trues_array)

delta_num = 100

final_worst = data.frame(matrix(ncol = 15, nrow = 0))
final_ave = data.frame(matrix(ncol = 15, nrow = 0))

for(i in 1:nsmax){
  result = result_list1[[i]]
  smax = smax_array[i]
  # cat("smax =",smax_array[i],'\n')
  for(j in 1:ntrues){
    delta = c(0,1,0)
    trues = trues_array[j]
    datte = Generate_test_data(L,p,q,n0,mixture,delta,trues,MU_A,BETA,ptemp)
    
    resvec_worst = rep(0,6)
    single_worst = rep(0,L)
    
    resvec_ave = rep(0,6)
    single_ave = rep(0,L)
    # cat("trues =",trues_array[i],'\n')
    for(k in 1:delta_num){
      if(k == 1){
        delta = c(0,1,0)
      }else{
        delta = simplex_uniform(L)
      }
      X0 = datte$X0
      
      YS = Generate_test_Y(L,p,q,X0=X0,Assigns = datte$Assigns,delta,trues,BETA,ptemp)
      bb = benchmarks(X0,YS$Y0,q,result$beta_star,result$beta_MI,result$rho,result$beta_RAP,result$DoublyR)
      resvec = bb$errvec
      single = bb$single_site
      
      resvec_worst = pmax(resvec_worst,resvec)
      single_worst = pmax(single_worst,single)
      
      resvec_ave = resvec_ave + resvec
      single_ave = single_ave + single
    }
    
    resvec_ave = resvec_ave/delta_num
    single_ave = single_ave/delta_num
    
    new_roww =c(smax,trues,resvec_worst,result$deltas,single_worst)
    new_rowa =c(smax,trues,resvec_ave,result$deltas,single_ave)
    final_worst = rbind(final_worst, new_roww)
    final_ave = rbind(final_ave, new_rowa)
  }
}

colnames(final_worst) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                           "d1","d2","d3","s","SS1","SS2","SS3")
colnames(final_ave) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                         "d1","d2","d3","s","SS1","SS2","SS3")

final_worst
final_ave