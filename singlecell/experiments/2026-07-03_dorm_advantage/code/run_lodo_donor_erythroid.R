#!/usr/bin/env Rscript
## Leave-one-DONOR-out on erythroid cells (the well-motivated covariate/donor-shift design).
## Sources = the 3 nucleated erythroid stages (Proery, Ery, Normo) pooled across all OTHER
## donors -> 3 sources. Target = the held-out donor's {Proery,Ery,Normo} cells = a NATURAL
## mixture of exactly those source types (within-convex-hull covariate/donor shift, no concept
## shift). Predict CD71 & CD49d with no target labels; compare all benchmarks (ntar=50).
## DORM config: q=20, penalty=TRUE (regularised Q). NOTE: penalty=FALSE was tried first but
## its doubly-robust normal-equations estimate BLOWS UP under donor batch shift (RAP reached
## -58); penalty=TRUE bounds Q via elastic-net and is stable. Env: DORM_SEEDS=1,2,3
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

ADTS <- c("CD71", "CD49d"); STAGES <- c("Proerythroblast", "Erythroblast", "Normoblast")
SEEDS <- as.integer(strsplit(Sys.getenv("DORM_SEEDS", "1,2,3"), ",")[[1]])
MAXN_SRC <- 2500; MIN_TGT <- 200; K <- 20
D <- load_cite_data(csv_dir)

## usable target donors: >= MIN_TGT nucleated erythroid cells
ecnt <- D$meta[cell_type %in% STAGES, .N, by = batch][N >= MIN_TGT]
DONORS <- ecnt[order(-N), batch]
message("Target donors (", length(DONORS), "): ", paste(DONORS, collapse = ", "))

out <- file.path(res_dir, "lodo_donor_ery_results.csv"); if (file.exists(out)) file.remove(out)
for (dtar in DONORS) for (y in ADTS) for (seed in SEEDS) {
  set.seed(seed)
  ## sources = each stage pooled over OTHER donors
  src_rows <- lapply(STAGES, function(g) sample_at_most(which(D$meta$cell_type == g & D$meta$batch != dtar), MAXN_SRC))
  names(src_rows) <- STAGES
  if (any(lengths(src_rows) < 60)) { message("  skip ", dtar, "/", y, ": a source too small"); next }
  ## target = held-out donor's nucleated erythroid cells (natural mixture)
  tgt_rows <- which(D$meta$cell_type %in% STAGES & D$meta$batch == dtar)
  if (length(tgt_rows) < MIN_TGT) next

  A <- screen_A(D, unlist(src_rows, use.names = FALSE), y, K); q <- K
  ssp <- lapply(src_rows, split_half); tsp <- split_half(tgt_rows)
  Xl  <- lapply(ssp, function(s) make_X(D, s$a, A)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A))
  Yl  <- lapply(ssp, function(s) make_Y(D, s$a, y)); Yt <- lapply(ssp, function(s) make_Y(D, s$b, y))
  names(Xl) <- names(Xt) <- names(Yl) <- names(Yt) <- STAGES; nlist <- as.integer(lengths(Yl))
  X0 <- make_X(D, tsp$a, A); Y0 <- make_Y(D, tsp$a, y); X0t <- make_X(D, tsp$b, A); Y0t <- make_Y(D, tsp$b, y)

  fit <- tryCatch(dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nlist, q, penalty = TRUE),
                  error = function(e) { message("  fit err ", dtar, "/", y, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  ev <- eval_all_methods(fit, Xl, Yl, X0, Y0, X0t, Y0t, nlist, q)
  ## record the target donor's natural composition (for interpretation)
  comp <- D$meta[tsp$a][, .N, by = cell_type]
  ev[, `:=`(donor = dtar, ADT = y, seed = seed, n_target = length(Y0),
            wEry = round(comp[cell_type=="Erythroblast", N] / length(tsp$a), 2))]
  fwrite(ev, out, append = file.exists(out))
  message(sprintf("  %s | %s | seed%d : DORM(tuned)=%.3f  RAP=%.3f  TransGLM=%.3f",
                  dtar, y, seed, ev[method=="DORM (tuned smax)", r2],
                  ev[method=="RAP", r2], ev[method=="TransGLM", r2]))
}
message("DONE lodo_donor_ery -> ", out)
