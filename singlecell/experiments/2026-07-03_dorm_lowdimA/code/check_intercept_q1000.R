## Quick isolation of the intercept fix: re-fit DORM at the OLD exact config
## (q=1000, penalty=TRUE, smax=0.03). Only difference vs the old run is the
## intercept fix in PQCalculation. Reuses saved q=1000 data objects.
suppressPackageStartupMessages(library(data.table))
source("src/Functions3.R")
rd <- "singlecell/dorm_results"
combos <- list(
  list(t="Reticulocyte", y="CD71",  old=-0.010),
  list(t="Reticulocyte", y="CD49d", old=-0.134),
  list(t="MK_E_prog",    y="CD71",  old=-1.435))
out <- rbindlist(lapply(combos, function(c){
  d <- readRDS(file.path(rd, sprintf("Ver4_erythroid_subtype_sources_target_%s_Y_%s_dorm_data.rds", c$t, c$y)))
  r <- get_DORM_beta(d$Xlist,d$Xtrainlist,d$Ylist,d$Ytrainlist,d$X0,d$X0train,d$nlist,d$q,
                     smax=0.03, penalty=TRUE, alpha=0.5, rho_pseudo="NA", dr_type="logit",
                     normalize=FALSE, condA=FALSE)
  b <- as.vector(r$beta_star); q <- d$q
  mse <- get_err(d$X0,d$Y0,b,q,length(d$Y0)); r2 <- 1-mse/var(d$Y0)
  data.table(target=c$t, ADT=c$y, n0=nrow(d$X0), q_plus1=q+1,
             intercept_new=round(b[1],3), predmean=round(mean(d$X0[,1:(q+1)]%*%b),3),
             Ymean=round(mean(d$Y0),3), DORM_old_r2=c$old, DORM_new_r2=round(r2,3))
}))
cat("\n==== Intercept fix on OLD q=1000 penalty=TRUE config (DORM only) ====\n\n")
print(out, row.names=FALSE)
fwrite(out, "singlecell/experiments/2026-07-03_dorm_lowdimA/results/intercept_fix_q1000_check.csv")
