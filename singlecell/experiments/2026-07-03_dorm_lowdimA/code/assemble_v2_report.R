#!/usr/bin/env Rscript

## Assemble the DORM v2 (low-dim A + tuned smax) report:
##   1. improvement table: old DORM (q=1000) vs new DORM (best k) vs benchmarks
##   2. ablation: old  ->  new smax=0.03  ->  new tuned smax
##   3. per-source Q tracking highlights
suppressPackageStartupMessages(library(data.table))

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
res_dir  <- file.path(exp_dir, "results")
repo_dir <- local({ d <- code_dir; while (!dir.exists(file.path(d,"src")) && dirname(d)!=d) d <- dirname(d); d })

new <- fread(file.path(res_dir, "v2_all_summary.csv"))
setnames(new, "target", "target_group")

## ---- new DORM: best k per (ADT, target) by tuned-smax R^2 ----
dorm_new <- new[method == "DORM (tuned smax)"]
best_k <- dorm_new[, .SD[which.max(r2)], by = .(protein_y, target_group)]
best_k_lab <- best_k[, .(protein_y, target_group, DORM_new = round(r2, 3),
                         k = k, best_smax = best_smax)]

## new benchmarks (same A) at the best k, best over ntar already
combos <- best_k[, .(protein_y, target_group, k)]
new_at_bestk <- merge(new, combos, by = c("protein_y", "target_group", "k"))
wide_new <- dcast(new_at_bestk, protein_y + target_group ~ method, value.var = "r2")

## ---- old DORM (q=1000) matrix ----
old_path <- file.path(repo_dir, "singlecell", "dorm_results", "Ver4_erythroid_R2_matrix.csv")
if (file.exists(old_path)) {
  old <- fread(old_path)[, .(protein_y, target_group,
                             DORM_old = round(DORM, 3),
                             SimpleAvg_old = round(`Simple average`, 3))]
} else old <- unique(new[, .(protein_y, target_group)])[, `:=`(DORM_old = NA_real_, SimpleAvg_old = NA_real_)]

## ---- 1. improvement table ----
bench_cols <- intersect(c("TransLasso","TransGLM","PTL"), names(wide_new))
wb <- wide_new[, c("protein_y","target_group", bench_cols), with = FALSE]
for (cc in bench_cols) wb[[cc]] <- round(wb[[cc]], 3)
imp <- Reduce(function(a,b) merge(a,b, by=c("protein_y","target_group"), all=TRUE),
              list(old, best_k_lab, wb))
setcolorder(imp, intersect(c("protein_y","target_group","DORM_old","DORM_new","k","best_smax",
                   "TransLasso","TransGLM","PTL","SimpleAvg_old"), names(imp)))
setorder(imp, target_group, protein_y)
cat("\n================ IMPROVEMENT: old DORM (q=1000) -> new DORM (low-dim A, tuned smax) ================\n\n")
print(imp, row.names = FALSE)
cat(sprintf("\nMean DORM R^2:  old = %.3f   new = %.3f   (delta = %+.3f)\n",
            mean(imp$DORM_old, na.rm=TRUE), mean(imp$DORM_new, na.rm=TRUE),
            mean(imp$DORM_new, na.rm=TRUE) - mean(imp$DORM_old, na.rm=TRUE)))

## ---- 2. ablation at best k ----
abl <- dcast(new_at_bestk[method %in% c("DORM (smax=0.03)","DORM (tuned smax)")],
             protein_y + target_group ~ method, value.var = "r2")
abl <- merge(old[, .(protein_y, target_group, DORM_old)], abl, by=c("protein_y","target_group"))
setnames(abl, c("DORM (smax=0.03)","DORM (tuned smax)"), c("new_smax0.03","new_tuned"))
abl[, `:=`(new_smax0.03 = round(new_smax0.03,3), new_tuned = round(new_tuned,3))]
setorder(abl, target_group, protein_y)
cat("\n================ ABLATION: old(q=1000)  ->  low-dim A @smax0.03  ->  low-dim A @tuned smax ================\n\n")
print(abl, row.names = FALSE)

## ---- 3. Q tracking highlights ----
cat("\n================ Per-source Q tracking (best-k run per ADT x target) ================\n\n")
qfiles <- combos[, .(f = file.path(res_dir,
  sprintf("v2_%s_target_%s_k%d_Qtrack.csv",
          gsub("[^0-9A-Za-z_.-]+","_",protein_y),
          gsub("[^0-9A-Za-z_.-]+","_",target_group), k))), by=.(protein_y,target_group,k)]
qall <- rbindlist(lapply(seq_len(nrow(qfiles)), function(i){
  f <- qfiles$f[i]; if (!file.exists(f)) return(NULL)
  d <- fread(f); d[, `:=`(protein_y=qfiles$protein_y[i], target_group=qfiles$target_group[i])]; d
}), fill=TRUE)
if (nrow(qall)) {
  setcolorder(qall, c("protein_y","target_group","source","Q_L2norm","standalone_target_r2","rho"))
  print(qall, row.names=FALSE)
  fwrite(qall, file.path(res_dir, "v2_Qtracking_all.csv"))
}

fwrite(imp, file.path(res_dir, "v2_improvement_table.csv"))
fwrite(abl, file.path(res_dir, "v2_ablation_table.csv"))
cat("\nSaved: v2_improvement_table.csv, v2_ablation_table.csv, v2_Qtracking_all.csv\n")
