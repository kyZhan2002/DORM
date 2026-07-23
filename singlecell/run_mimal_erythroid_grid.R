#!/usr/bin/env Rscript

## Driver: DORM on the erythroid developmental design borrowed from the MIMAL
## BMMC real-data analysis. We reuse ONLY their experimental design
##   sources  = {Proerythroblast, Erythroblast, Normoblast}   (nucleated stages)
##   targets  = {Reticulocyte (terminal), MK/E prog (progenitor)}   held-out, eval only
##   ADTs     = {CD71, CD36_1, CD105, CD49d}
##   features = 1000 HVG (elastic-net penalty)
## All methods/benchmarks are the DORM family (DORM, MI, RAP, best-single-source,
## simple/rho average, TransLasso, TransGLM, PTL). MIMAL's own method is NOT run.
##
## This wrapper loops the 4 ADTs, invoking run_dorm_from_cite_nips_csv.R once per
## ADT (each invocation runs both target cell types). Usage:
##   Rscript singlecell/run_mimal_erythroid_grid.R
##   DORM_PROTEINS=CD71,CD49d Rscript singlecell/run_mimal_erythroid_grid.R   # subset

this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))[1])
singlecell_dir <- dirname(this_file)
run_script <- file.path(singlecell_dir, "run_dorm_from_cite_nips_csv.R")
if (!file.exists(run_script)) stop("Cannot find run_dorm_from_cite_nips_csv.R next to this driver.")

proteins_env <- Sys.getenv("DORM_PROTEINS", unset = "CD71,CD36_1,CD105,CD49d")
proteins <- trimws(strsplit(proteins_env, ",", fixed = TRUE)[[1]])

rscript <- file.path(R.home("bin"), "Rscript")

for (p in proteins) {
  message("\n########################################################")
  message("### DORM erythroid grid | ADT = ", p)
  message("########################################################")
  env <- c(
    paste0("DORM_PARTITION=erythroid_subtype"),
    paste0("DORM_MAX_RNA_FEATURES=1000"),
    paste0("DORM_PROTEIN_Y=", p),
    paste0("DORM_RUN_TRANSFER=true")
  )
  status <- system2(rscript, args = run_script, env = env)
  if (status != 0) {
    message("!!! Run for ADT ", p, " exited with status ", status, "; continuing to next ADT.")
  }
}

message("\nAll ADTs done. Aggregate with: Rscript singlecell/assemble_mimal_table.R")
