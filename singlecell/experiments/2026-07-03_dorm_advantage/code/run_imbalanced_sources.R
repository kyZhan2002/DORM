#!/usr/bin/env Rscript
## Imbalanced-source demonstration (the regime DORM is built for).
## Sources = 3 erythroid stages, but Erythroblast is made DOMINANT (size = base*imb) while
## Proerythroblast/Normoblast stay small. Target = a MINORITY stage (Normoblast), held out.
## As the imbalance grows, naive pooling collapses toward the Erythroblast model (biased for
## the Normoblast target), while DORM reweights (rho -> Normoblast) and holds up.
## Money result: DORM - Pooled R2 gap grows with source imbalance.
##
## Env: DORM_ADTS=CD71,CD49d  DORM_IMB=1,2,4,6  DORM_SEEDS=1,2  DORM_TARGET_STAGE=Normoblast

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results")
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

env_csv <- function(n, d) { v <- Sys.getenv(n); if (v == "") d else trimws(strsplit(v, ",")[[1]]) }
ADTS   <- env_csv("DORM_ADTS", c("CD71", "CD49d"))
IMB    <- as.numeric(env_csv("DORM_IMB", c("1", "2", "4", "6")))
SEEDS  <- as.integer(env_csv("DORM_SEEDS", c("1", "2")))
TGT    <- Sys.getenv("DORM_TARGET_STAGE", "Normoblast")
K      <- 20; NBASE <- 300; NTGT <- 400
DOM    <- "Erythroblast"                    # the dominant source
MINOR  <- setdiff(c("Proerythroblast", "Erythroblast", "Normoblast"), c(DOM))
STAGES <- c("Proerythroblast", "Erythroblast", "Normoblast")

D <- load_cite_data(csv_dir)
out <- file.path(res_dir, "imbalanced_results.csv"); if (file.exists(out)) file.remove(out)

for (y_col in ADTS) for (seed in SEEDS) for (imb in IMB) {
  set.seed(seed * 100 + round(imb))
  ## per-stage: reserve split-A (sources) vs split-B (target pool)
  sa <- list(); tb <- list()
  for (g in STAGES) { s <- split_half(which(D$meta$cell_type == g)); sa[[g]] <- s$a; tb[[g]] <- s$b }
  ## imbalanced source sizes: dominant = base*imb, others = base
  size <- setNames(rep(NBASE, 3), STAGES); size[DOM] <- round(NBASE * imb)
  src_rows <- lapply(STAGES, function(g) sample(sa[[g]], min(size[g], length(sa[[g]]))))
  names(src_rows) <- STAGES

  A_genes <- screen_A(D, unlist(src_rows, use.names = FALSE), y_col, K); q <- K
  ssp <- lapply(src_rows, split_half)
  Xlist      <- lapply(ssp, function(s) make_X(D, s$a, A_genes))
  Xtrainlist <- lapply(ssp, function(s) make_X(D, s$b, A_genes))
  Ylist      <- lapply(ssp, function(s) make_Y(D, s$a, y_col))
  Ytrainlist <- lapply(ssp, function(s) make_Y(D, s$b, y_col))
  names(Xlist) <- names(Xtrainlist) <- names(Ylist) <- names(Ytrainlist) <- STAGES
  nlist <- as.integer(lengths(Ylist))

  ## target = pure MINORITY stage (held-out cells)
  tgt <- sample(tb[[TGT]], min(NTGT, length(tb[[TGT]])))
  tsp <- split_half(tgt)
  X0 <- make_X(D, tsp$a, A_genes); Y0 <- make_Y(D, tsp$a, y_col)
  X0train <- make_X(D, tsp$b, A_genes); Y0train <- make_Y(D, tsp$b, y_col)

  fit <- tryCatch(dorm_fit(Xlist, Xtrainlist, Ylist, Ytrainlist, X0, X0train, nlist, q),
                  error = function(e) { message("err ", y_col, "/imb", imb, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  yv <- var(Y0); R2 <- function(b) 1 - get_err(X0, Y0, as.vector(b), q, length(Y0)) / yv
  nt <- min(100, length(Y0train))
  tune <- tuning_Y(X0train[1:nt, , drop = FALSE], Y0train[1:nt], q, fit$beta_by_smax, SMAX_GRID)
  beta_dorm <- fit$beta_by_smax[which(SMAX_GRID == tune$best_smax)[1], ]
  src_mse <- apply(fit$Q, 2, function(b) get_err(X0, Y0, b, q, length(Y0)))
  row <- data.table(
    ADT = y_col, seed = seed, imb = imb, dom_frac = round(size[DOM]/sum(size), 2),
    target = TGT, rho_target = round(fit$rho[which(STAGES == TGT)], 3),
    DORM = round(R2(beta_dorm), 3), RAP = round(R2(fit$P), 3),
    Pooled = round(R2(pooled_en(Xlist, Ylist, q)), 3),
    BestSingle = round(1 - min(src_mse)/yv, 3),
    MI = round(R2(fit$beta_MI), 3))
  row[, gap := DORM - Pooled]
  fwrite(row, out, append = file.exists(out))
  message(sprintf("  %s seed%d imb=%g (Ery %.0f%%): DORM=%.3f Pooled=%.3f gap=%+.3f rho_%s=%.2f",
                  y_col, seed, imb, 100*size[DOM]/sum(size), row$DORM, row$Pooled, row$gap, TGT, row$rho_target))
}
message("DONE imbalanced -> ", out)
