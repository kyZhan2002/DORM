# This script processes beta.

library(here)

source(here('src', 'DataGeneration3.R'))
source(here('src', 'Functions3.R'))
source(here('src', 'Tuning.R'))

# fourth site is away from others.

evaluation_smax = function(L,p,q,ptemp,n0,mixture,MU_A,BETA,result_list,
                           trues_array = seq(0,1,0.05),delta_num = 100,
                           delta_best = c(0,0,0,0,1)){
  
  # only works for L=5
  smaxn = names(result_list)
  smax_array = as.numeric(smaxn)
  nsmax = length(smaxn)
  
  ntrues = length(trues_array)
  
  final_worst = data.frame(matrix(ncol = 19, nrow = 0))
  final_ave = data.frame(matrix(ncol = 19, nrow = 0))
  
  for(i in 1:nsmax){
    result = result_list[[i]]
    smax = smax_array[i]
    # cat("smax =",smax_array[i],'\n')
    for(j in 1:ntrues){
      delta = delta_best
      trues = trues_array[j]
      datte = Generate_test_data(L,p,q,n0,mixture,delta,trues,MU_A,BETA,ptemp)
      
      resvec_worst = rep(0,6)
      single_worst = rep(0,5)
      
      resvec_ave = rep(0,6)
      single_ave = rep(0,5)
      # cat("trues =",trues_array[i],'\n')
      for(k in 1:delta_num){
        if(k == 1){
          delta = delta_best
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
                             "d1","d2","d3","d4","d5","s","SS1","SS2","SS3","SS4","SS5")
  colnames(final_ave) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                           "d1","d2","d3","d4","d5","s","SS1","SS2","SS3","SS4","SS5")
  
  out = list(wst = final_worst,avg = final_ave)
  return(out)
}
