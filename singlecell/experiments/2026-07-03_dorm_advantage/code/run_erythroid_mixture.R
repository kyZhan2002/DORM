#!/usr/bin/env Rscript
## Phase 1: adversarial mixture-shift stress test (erythroid).
## Sources = 3 nucleated stages (trained on split-A cells). Target = HELD-OUT (split-B)
## cells resampled to a mixture weight vector over the 3 stages, swept from balanced
## (centroid) to extreme (vertices). DORM uses no target labels. Shows DORM's worst-case
## R2 across the mixture spectrum vs label-free baselines (pooling collapses on extremes).
##
## Env: DORM_ADTS=CD71,CD49d  DORM_K=20  DORM_SEEDS=1,2  DORM_NTARGET=600

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

env_csv <- function(n, d) { v <- Sys.getenv(n); if (v == "") d else trimws(strsplit(v, ",")[[1]]) }
ADTS    <- env_csv("DORM_ADTS", c("CD71", "CD49d"))
K       <- as.integer(Sys.getenv("DORM_K", "20"))
SEEDS   <- as.integer(env_csv("DORM_SEEDS", c("1", "2")))
NTARGET <- as.integer(Sys.getenv("DORM_NTARGET", "600"))
STAGES  <- c("Proerythroblast", "Erythroblast", "Normoblast")

## mixture weight vectors over (Proery, Ery, Normo)
MIX <- list(
  vProery = c(1,0,0), vEry = c(0,1,0), vNormo = c(0,0,1),
  ePE = c(.5,.5,0), ePN = c(.5,0,.5), eEN = c(0,.5,.5),
  centroid = c(1,1,1)/3,
  nearProery = c(.7,.15,.15), nearEry = c(.15,.7,.15), nearNormo = c(.15,.15,.7))

D <- load_cite_data(csv_dir)
out_path <- file.path(res_dir, "mixture_results.csv")
w_path   <- file.path(res_dir, "mixture_weights.csv")
if (file.exists(out_path)) file.remove(out_path)
if (file.exists(w_path)) file.remove(w_path)

for (y_col in ADTS) for (seed in SEEDS) {
  set.seed(seed)
  ## reserve split-A (source training) vs split-B (target pool) per stage
  stage_rows <- lapply(STAGES, function(g) which(D$meta$cell_type == g)); names(stage_rows) <- STAGES
  srcA <- lapply(stage_rows, function(r) split_half(r)$a)   # sources
  tgtB <- lapply(stage_rows, function(r) split_half(r)$b)   # target pool

  A_genes <- screen_A(D, unlist(srcA, use.names = FALSE), y_col, K); q <- K
  ## build cross-fit source halves from split-A
  ssp <- lapply(srcA, split_half)
  Xlist      <- lapply(ssp, function(s) make_X(D, s$a, A_genes))
  Xtrainlist <- lapply(ssp, function(s) make_X(D, s$b, A_genes))
  Ylist      <- lapply(ssp, function(s) make_Y(D, s$a, y_col))
  Ytrainlist <- lapply(ssp, function(s) make_Y(D, s$b, y_col))
  names(Xlist) <- names(Xtrainlist) <- names(Ylist) <- names(Ytrainlist) <- STAGES
  nlist <- as.integer(lengths(Ylist))

  for (mn in names(MIX)) {
    w <- MIX[[mn]]
    ns <- pmin(round(w * NTARGET), lengths(tgtB))          # cells per stage for this mixture
    tgt <- unlist(lapply(seq_along(STAGES), function(i)
      if (ns[i] > 0) sample(tgtB[[i]], ns[i]) else integer(0)))
    if (length(tgt) < 100) next
    tsp <- split_half(tgt)
    X0 <- make_X(D, tsp$a, A_genes); Y0 <- make_Y(D, tsp$a, y_col)
    X0train <- make_X(D, tsp$b, A_genes); Y0train <- make_Y(D, tsp$b, y_col)

    fit <- tryCatch(dorm_fit(Xlist, Xtrainlist, Ylist, Ytrainlist, X0, X0train, nlist, q),
                    error = function(e) { message("fit err ", y_col, "/", mn, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    ev <- eval_all_methods(fit, Xlist, Ylist, X0, Y0, X0train, Y0train, nlist, q)
    ev[, `:=`(ADT = y_col, mixture = mn, seed = seed, n_target = length(tgt),
              wP = w[1], wE = w[2], wN = w[3])]
    fwrite(ev, out_path, append = file.exists(out_path))
    wt <- data.table(ADT = y_col, mixture = mn, seed = seed, source = STAGES,
                     rho = round(as.vector(fit$rho), 4))
    fwrite(wt, w_path, append = file.exists(w_path))
    message(sprintf("  %s | %s | seed %d : DORM(tuned) r2=%.3f  pooled=%.3f",
                    y_col, mn, seed,
                    ev[method == "DORM (tuned smax)", r2], ev[method == "Pooled elastic-net", r2]))
  }
}
message("DONE mixture -> ", out_path)
