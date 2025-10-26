## Simulation script for Conditional A DORM method with L=5

## Check and install required packages
required_packages <- c("CVXR", "here")
for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE, quietly = TRUE)){
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

library(here)
source(here('src', 'DataGen_CondA.R'))
source(here('src', 'Functions3.R'))
source(here('src', 'Tuning.R'))

## Helper function to generate uniform random point on simplex
simplex_uniform = function(L) {
  # Generate L-1 uniform random variables, sort them, and compute differences
  u = sort(c(0, runif(L - 1), 1))
  delta = diff(u)
  return(delta)
}

## Get seed from command line argument
args = commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  seed = 123  # default seed for testing
  cat("No seed provided, using default seed = 123\n")
} else {
  seed = as.integer(args[1])
  cat("Using seed =", seed, "\n")
}

set.seed(seed)

## Parameter settings

L = 5     # number of source sites
q = 4     # dim of low-d A
p = 95    # rest dim of covariate X. Total dim of X is (1+q+p)
ptemp = 6

# Sample sizes
Nlist = rep(1000, L)
nlist = rep(500, L)
n0 = 500

# Mixture parameters
mixture = c(0.5, 0, 0.5, 0, 0)  # Target mixture weights
delta = c(0, 0, 1, 0, 0)    # Perturbation mixture

# Mean of A for each source
MU_A = 0.5 * matrix(c(-2, 0, 1, 0,
                       0, -1, 2, 1,
                      -1, 2, 0, 1,
                       1, 1, -1, 2,
                       2, 1, 1, 2), nrow = L, ncol = q, byrow = TRUE)

# True BETA coefficients
BETA = 1 * matrix(c( 0, 4, 4, 3, -3,  0.5, -0.2, 0, 0, 1, 0,
                    -3, 3, -1, 3, -3,  0, 0.5, -0.2, 0, 0, 1,
                    -3, -3, 0, -3, -3,  0, 0, 0.2, -0.5, 0, 1,
                    -6, 6, 6, 6, -6,   1, 0, 0, 0.2, -0.5, 0,
                    -3, 3, 0, 0, 1,    0, 0, 0, 1, 0.2, -0.5),
                  nrow = L, ncol = ptemp + q + 1, byrow = TRUE)

# Simulation parameters
smax_array = seq(0.05, 0.5, 0.05)
nsmax = length(smax_array)
trues_array = seq(0, 1, 0.05)
ntrues = length(trues_array)
delta_num = 100  # Number of random deltas to average over

## Step 1: Generate training data and fit DORM models for each smax

cat("=== Generating training data ===\n")
dat = Generate_Simulation_Data_CondA(L, p, q, Nlist, nlist, n0, mixture,
                                      delta, s = 0, MU_A, BETA, ptemp, 
                                      MU_Acoef = 1, eps_A = 1, eps_W = 0.1,
                                      eps_beta = 0.1, eps_Y = 0.25,
                                      sd_WA = 1, tarmixA = FALSE)

cat("=== Fitting DORM models for different smax values ===\n")
result_list = list()

for(i in 1:nsmax){
  smax = smax_array[i]
  cat(sprintf("Fitting smax = %.2f (%d/%d)\n", smax, i, nsmax))
  
  result = get_DORM_beta(dat$Xlist, dat$Xtrainlist, 
                         dat$Ylist, dat$Ytrainlist,
                         dat$X0, dat$X0train, nlist, q, smax,
                         penalty = FALSE, alpha = 0.5, 
                         rho_pseudo = 'pool', dr_type = 'logit',
                         normalize = FALSE, condA = TRUE)
  
  result_list[[as.character(smax)]] <- result
}

## Step 2: Evaluate on test data with varying trues and delta

cat("\n=== Evaluating on test data ===\n")

final_worst = data.frame(matrix(ncol = 2 + 6 + (L + 1) + L, nrow = 0))
final_ave = data.frame(matrix(ncol = 2 + 6 + (L + 1) + L, nrow = 0))

colnames(final_worst) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                           paste0("d", 1:L), "s", paste0("SS", 1:L))
colnames(final_ave) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                         paste0("d", 1:L), "s", paste0("SS", 1:L))

for(i in 1:nsmax){
  smax = smax_array[i]
  result = result_list[[as.character(smax)]]
  
  cat(sprintf("\nEvaluating smax = %.2f (%d/%d)\n", smax, i, nsmax))
  
  for(j in 1:ntrues){
    trues = trues_array[j]
    
    if(j %% 5 == 1) {
      cat(sprintf("  trues = %.2f (%d/%d)\n", trues, j, ntrues))
    }
    
    resvec_worst = rep(0, 6)
    single_worst = rep(0, L)
    
    resvec_ave = rep(0, 6)
    single_ave = rep(0, L)
    
    for(k in 1:delta_num){
      # Generate random delta (simplex uniform)
      if(k == 1){
        delta_test = rep(1/L, L)  # Start with uniform
      }else{
        delta_test = simplex_uniform(L)
      }
      
      # Generate test data
      datte = Generate_test_data_CondA(L, p, q, n0, mixture, delta_test, 
                                       trues, MU_A, BETA, 
                                       WtoA_list = dat$WtoA_list,
                                       ptemp = ptemp, MU_Acoef = 1,
                                       eps_A = 1, eps_W = 0.1,
                                       eps_beta = 0.1, eps_Y = 0.25,
                                       sd_WA = 1, tarmixA = FALSE)
      
      # Compute errors
      bb = benchmarks(datte$X0, datte$Y0, q, 
                     result$beta_star, result$beta_MI, result$rho,
                     result$beta_RAP, result$DoublyR)
      
      resvec = bb$errvec
      single = bb$single_site
      
      # Track worst case
      resvec_worst = pmax(resvec_worst, resvec)
      single_worst = pmax(single_worst, single)
      
      # Accumulate for average
      resvec_ave = resvec_ave + resvec
      single_ave = single_ave + single
    }
    
    # Compute averages
    resvec_ave = resvec_ave / delta_num
    single_ave = single_ave / delta_num
    
    # Store results
    new_roww = c(smax, trues, resvec_worst, result$deltas, single_worst)
    new_rowa = c(smax, trues, resvec_ave, result$deltas, single_ave)
    
    final_worst = rbind(final_worst, new_roww)
    final_ave = rbind(final_ave, new_rowa)
  }
}

# Assign column names again (rbind may strip them)
colnames(final_worst) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                           paste0("d", 1:L), "s", paste0("SS", 1:L))
colnames(final_ave) <- c("smax", "true-s", "Ours", "SS", "SA", "RA", "MI", "PA",
                         paste0("d", 1:L), "s", paste0("SS", 1:L))

## Step 3: Save results

output_dir = here('simu', 'simu_10242025', 'data')
if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

output_file_worst = file.path(output_dir, paste0("CondA_L5_worst_seed", seed, ".rds"))
output_file_ave = file.path(output_dir, paste0("CondA_L5_ave_seed", seed, ".rds"))

saveRDS(list(
  final_worst = final_worst,
  parameters = list(L = L, q = q, p = p, mixture = mixture,
                   smax_array = smax_array, trues_array = trues_array,
                   delta_num = delta_num, Nlist = Nlist, nlist = nlist, n0 = n0,
                   seed = seed)
), file = output_file_worst)

saveRDS(list(
  final_ave = final_ave,
  parameters = list(L = L, q = q, p = p, mixture = mixture,
                   smax_array = smax_array, trues_array = trues_array,
                   delta_num = delta_num, Nlist = Nlist, nlist = nlist, n0 = n0,
                   seed = seed)
), file = output_file_ave)

cat("\n=== Results saved ===\n")
cat("Worst case:", output_file_worst, "\n")
cat("Average case:", output_file_ave, "\n")

## Step 4: Print summary statistics

cat("\n=== Summary Statistics ===\n")
cat("Worst case performance:\n")
print(summary(final_worst[, c("Ours", "SS", "SA", "RA", "MI", "PA")]))

cat("\nAverage case performance:\n")
print(summary(final_ave[, c("Ours", "SS", "SA", "RA", "MI", "PA")]))

cat("\nSimulation completed successfully!\n")