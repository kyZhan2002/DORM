#!/usr/bin/env Rscript

## Aggregate the erythroid-design DORM runs into a MIMAL-Table-5-style matrix:
##   rows    = ADT x held-out target cell type (8 combos)
##   columns = method (DORM family)
##   values  = out-of-sample R^2 on the held-out target labels
## Reads the per-target *_summary.csv files written by run_dorm_from_cite_nips_csv.R
## under the Ver4_erythroid_subtype_sources_* stub.

suppressPackageStartupMessages(library(data.table))

singlecell_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))[1])),
  error = function(e) getwd())
out_dir <- file.path(singlecell_dir, "dorm_results")

files <- list.files(out_dir, pattern = "^Ver4_erythroid_subtype_sources_target_.*_summary\\.csv$",
                    full.names = TRUE)
files <- files[!grepl("all_targets_summary", files)]
if (length(files) == 0) stop("No Ver4 erythroid summary files found in ", out_dir,
                             ". Run run_mimal_erythroid_grid.R first.")

dat <- rbindlist(lapply(files, fread), fill = TRUE)

## For transfer methods with several ntar rows, keep the best R^2 per method/combo.
dat <- dat[!is.na(r2)]
best <- dat[, .SD[which.max(r2)], by = .(target_group, protein_y, method)]

## Long -> wide: one row per (ADT, target), one column per method.
best[, combo := paste0(protein_y, " / ", target_group)]
wide <- dcast(best, protein_y + target_group ~ method, value.var = "r2")

## Order methods with DORM first if present.
method_cols <- setdiff(names(wide), c("protein_y", "target_group"))
preferred <- c("DORM", "MI", "RAP", "Best single source", "Simple average",
               "Rho average", "TransLasso", "TransGLM", "PTL")
method_cols <- c(intersect(preferred, method_cols), setdiff(method_cols, preferred))
setcolorder(wide, c("protein_y", "target_group", method_cols))
setorder(wide, target_group, protein_y)

## Round for display.
disp <- copy(wide)
for (m in method_cols) disp[[m]] <- round(disp[[m]], 3)

cat("\n==== Out-of-sample target R^2 by ADT x held-out cell type (DORM family) ====\n\n")
print(disp, row.names = FALSE)

## Averages per method (over the 8 combos) and per target cell type.
avg_overall <- best[, .(mean_r2 = round(mean(r2), 3)), by = method][order(-mean_r2)]
avg_by_target <- dcast(best[, .(mean_r2 = mean(r2)), by = .(method, target_group)],
                       method ~ target_group, value.var = "mean_r2")
num_cols <- setdiff(names(avg_by_target), "method")
for (m in num_cols) avg_by_target[[m]] <- round(avg_by_target[[m]], 3)

cat("\n==== Mean R^2 per method (averaged over all ADT x target combos) ====\n\n")
print(avg_overall, row.names = FALSE)
cat("\n==== Mean R^2 per method, split by held-out target cell type ====\n\n")
print(avg_by_target, row.names = FALSE)

fwrite(wide, file.path(out_dir, "Ver4_erythroid_R2_matrix.csv"))
fwrite(avg_overall, file.path(out_dir, "Ver4_erythroid_R2_mean_per_method.csv"))
fwrite(best, file.path(out_dir, "Ver4_erythroid_best_per_method_long.csv"))
cat("\nSaved: Ver4_erythroid_R2_matrix.csv, Ver4_erythroid_R2_mean_per_method.csv, ",
    "Ver4_erythroid_best_per_method_long.csv (in ", out_dir, ")\n", sep = "")
