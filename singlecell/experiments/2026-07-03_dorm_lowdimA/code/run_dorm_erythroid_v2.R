#!/usr/bin/env Rscript

## DORM v2 on the erythroid design, run in DORM's INTENDED regime.
##
## Diagnosis of the old (q=1000) runs (all R^2 <= 0):
##   (1) A = all 1000 genes, W empty -> high-dim degenerate maximin; the final
##       predictor X[,1:(q+1)] %*% beta_star is A-only but over-regularized, and
##       singular when target n0 < q+1 (MK/E prog).
##   (2) penalty=TRUE shrinks the INTERCEPT to 0 -> predictions not centered on
##       E[Y] -> negative R^2.
##   (3) smax=0.03 pins s* at the cap -> beta_star == RAP; smax never tuned.
##
## Fix (this script): use the contract from src/Example_new.R --
##   X = (1, A, W) with A = a SMALL screened set of predictive genes (low-dim,
##   q = k), W = the remaining HVG used only as nuisance controls (imputation +
##   density ratio). penalty = FALSE (so solve(Sigma) estimates a proper,
##   UNPENALIZED intercept), and smax is TUNED on a small labeled target split.
##
## A-genes are screened LABEL-FREE on pooled SOURCE data only.
##
## Env knobs (all optional):
##   DORM_ADTS=CD71,CD36_1,CD105,CD49d   DORM_TARGETS=Reticulocyte,MK/E prog
##   DORM_KS=10,20,30                     DORM_NTUNE=100
##   DORM_RUN_TRANSFER=true

suppressPackageStartupMessages({ library(data.table); library(glmnet) })

## ---- locate repo + IO ----
this_file <- normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))[1])
code_dir <- dirname(this_file)
exp_dir  <- dirname(code_dir)
repo_dir <- local({
  d <- code_dir
  while (!file.exists(file.path(d, "src", "Functions3.R")) && dirname(d) != d) d <- dirname(d)
  d
})
source(file.path(repo_dir, "src", "Functions3.R"))
source(file.path(repo_dir, "src", "TransLasso-functions.R"))
source(file.path(repo_dir, "src", "TransLG-functions.R"))
source(file.path(repo_dir, "src", "Tuning.R"))
csv_dir <- file.path(repo_dir, "singlecell", "processed", "csv_for_R")
res_dir <- file.path(exp_dir, "results"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(2026)
env_csv <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") return(default)
  trimws(strsplit(v, ",", fixed = TRUE)[[1]])
}
env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_); if (is.na(v) || v == "") default else as.integer(v)
}
env_flag <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") return(default); tolower(v) %in% c("1","true","t","yes","y")
}

ADTS      <- env_csv("DORM_ADTS", c("CD71", "CD36_1", "CD105", "CD49d"))
TARGETS   <- env_csv("DORM_TARGETS", c("Reticulocyte", "MK/E prog"))
KS        <- as.integer(env_csv("DORM_KS", c("10", "20", "30")))
NTUNE     <- env_int("DORM_NTUNE", 100)
RUN_TL    <- env_flag("DORM_RUN_TRANSFER", TRUE)
SMAX_GRID <- c(0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9)   # 0.03 kept = old fixed value
SOURCES   <- c("Proerythroblast", "Erythroblast", "Normoblast")
ERY5      <- c(SOURCES, "Reticulocyte", "MK/E prog")
MAX_ROWS  <- 12000

## ---- load data ----
meta <- fread(file.path(csv_dir, "cite_nips_obs_metadata_for_R.csv"))
fmap <- fread(file.path(csv_dir, "cite_nips_feature_name_map_for_R.csv"))
hvg  <- fread(file.path(csv_dir, "cite_nips_hvg_genes_for_R.csv"))
rna_features <- intersect(hvg[!is.na(csv_name), csv_name], fmap[modality == "RNA", csv_name])
adt_features <- fmap[modality == "ADT", csv_name]
metadata_cols <- c("obs_id", "batch", "cell_type", "cell_type_major")

message("Reading RNA (", length(rna_features), " HVG) + ADT (", length(adt_features), " proteins) ...")
rna_dt <- fread(file.path(csv_dir, "cite_nips_RNA_for_R.csv"),
                select = c(metadata_cols, rna_features))
adt_dt <- fread(file.path(csv_dir, "cite_nips_ADT_for_R.csv"),
                select = c(metadata_cols, adt_features))
stopifnot(identical(rna_dt$obs_id, adt_dt$obs_id), identical(rna_dt$obs_id, meta$obs_id))

## ---- helpers ----
sample_at_most <- function(idx, nmax) if (length(idx) <= nmax) as.integer(idx) else sort(sample(as.integer(idx), nmax))
## equal-size, disjoint halves (DORM requires Xlist and Xtrainlist to match in size)
split_half <- function(idx) { idx <- sample(as.integer(idx)); n <- floor(length(idx)/2); list(a = idx[seq_len(n)], b = idx[n + seq_len(n)]) }
safe <- function(x) gsub("[^0-9A-Za-z_.-]+", "_", x)

## Screen A genes: elastic-net of Y on ALL HVG using pooled SOURCE data only,
## rank by |coef|, take top-k. Label-free w.r.t. the target.
screen_A <- function(src_rows, y_col, kmax) {
  X <- as.matrix(rna_dt[src_rows, ..rna_features]); storage.mode(X) <- "double"
  y <- as.numeric(adt_dt[[y_col]][src_rows])
  fit <- cv.glmnet(X, y, alpha = 0.5, nfolds = 5, standardize = FALSE)
  b <- as.vector(coef(fit, s = "lambda.min"))[-1]
  ord <- order(abs(b), decreasing = TRUE)
  ord <- c(ord[abs(b[ord]) > 0], setdiff(order(abs(b), decreasing = TRUE), ord[abs(b[ord]) > 0]))
  ## greedy pick top-k by |coef|, skipping near-duplicates (|corr|>0.999) so that
  ## Sigma0 = X_A'X_A stays invertible in the penalty=FALSE solve().
  picked <- integer(0)
  for (j in ord) {
    if (length(picked) == 0 ||
        max(abs(cor(X[, j], X[, picked, drop = FALSE]))) < 0.999) {
      picked <- c(picked, j)
      if (length(picked) == kmax) break
    }
  }
  rna_features[picked]
}

## Build X = [intercept, A-genes, W-genes]; q = length(A).
make_X <- function(rows, A_genes) {
  W_genes <- setdiff(rna_features, A_genes)
  cols <- c(A_genes, W_genes)
  X <- as.matrix(rna_dt[rows, ..cols]); storage.mode(X) <- "double"
  cbind(Intercept = 1, X)
}
make_Y <- function(rows, y_col) as.numeric(adt_dt[[y_col]][rows])

## Combine (P,Q,Sigma0) at a given smax into a DORM beta.
combine_beta <- function(P, Q, Sigma0, smax) {
  L <- ncol(Q)
  dts <- opt_s_delta(P, Q, Sigma0, smax)
  s <- dts[L + 1]
  as.vector((1 - s) * P + s * Q %*% dts[1:L])
}

## ---- one (ADT, target, k) run ----
run_combo <- function(y_col, target, k) {
  keep_src <- lapply(SOURCES, function(g) which(meta$cell_type == g))
  names(keep_src) <- SOURCES
  tgt_rows <- sample_at_most(which(meta$cell_type == target), MAX_ROWS)
  src_rows <- lapply(keep_src, sample_at_most, nmax = MAX_ROWS)

  ## screen A on pooled sources (label-free wrt target)
  A_genes <- screen_A(unlist(src_rows, use.names = FALSE), y_col, k)
  q <- k

  ## source splits (equal halves for cross-fitting)
  ssp <- lapply(src_rows, split_half)
  Xlist      <- lapply(ssp, function(s) make_X(s$a, A_genes))
  Xtrainlist <- lapply(ssp, function(s) make_X(s$b, A_genes))
  Ylist      <- lapply(ssp, function(s) make_Y(s$a, y_col))
  Ytrainlist <- lapply(ssp, function(s) make_Y(s$b, y_col))
  names(Xlist) <- names(Xtrainlist) <- names(Ylist) <- names(Ytrainlist) <- SOURCES
  nlist <- as.integer(lengths(Ylist))

  ## target: eval split (X0,Y0) + aux split (X0train,Y0train); tuning uses aux LABELS only
  tsp <- split_half(tgt_rows)
  X0 <- make_X(tsp$a, A_genes); Y0 <- make_Y(tsp$a, y_col)             # evaluation (labels used only here)
  X0train <- make_X(tsp$b, A_genes); Y0train <- make_Y(tsp$b, y_col)   # DR auxiliary
  ntune <- min(NTUNE, length(Y0train))

  message(sprintf("  [%s | %s | k=%d]  q+1=%d  n0(eval)=%d  aux=%d  sources=%s",
                  y_col, target, k, q + 1, nrow(X0), nrow(X0train), paste(nlist, collapse = "/")))

  ## two cross-fit halves at penalty=FALSE; sweep smax EXACTLY from the returned P,Q.
  o1 <- maximin_s_beta(Xlist, Xtrainlist, Ylist, Ytrainlist, X0, X0train,
                       nlist, q, smax = 0.5, penalty = FALSE, alpha = 0.5,
                       rho_pseudo = "NA", dr_type = "logit", normalize = FALSE, condA = FALSE)
  o2 <- maximin_s_beta(Xtrainlist, Xlist, Ytrainlist, Ylist, X0train, X0,
                       nlist, q, smax = 0.5, penalty = FALSE, alpha = 0.5,
                       rho_pseudo = "NA", dr_type = "logit", normalize = FALSE, condA = FALSE)
  S1 <- (1/nrow(X0))      * t(X0[,1:(q+1)])      %*% X0[,1:(q+1)]
  S2 <- (1/nrow(X0train)) * t(X0train[,1:(q+1)]) %*% X0train[,1:(q+1)]

  beta_array <- t(sapply(SMAX_GRID, function(sm)
    0.5 * (combine_beta(o1$beta_RAP, o1$DoublyR, S1, sm) +
           combine_beta(o2$beta_RAP, o2$DoublyR, S2, sm))))

  ## tune smax on aux labels (disjoint from eval)
  tune <- tuning_Y(X0train[1:ntune, , drop = FALSE], Y0train[1:ntune], q, beta_array, SMAX_GRID)
  best_i <- which(SMAX_GRID == tune$best_smax)[1]
  beta_dorm_tuned <- beta_array[best_i, ]
  beta_dorm_s03   <- beta_array[1, ]   # smax = 0.03 (old fixed value) for ablation

  ## averaged nuisance for DORM-family baselines + Q tracking
  Q   <- 0.5 * (o1$DoublyR + o2$DoublyR)
  P   <- 0.5 * (o1$beta_RAP + o2$beta_RAP)
  bMI <- 0.5 * (o1$beta_MI  + o2$beta_MI)
  rho <- 0.5 * (o1$rho + o2$rho)
  src_mse <- apply(Q, 2, function(b) get_err(X0, Y0, b, q, length(Y0)))
  best_src <- which.min(src_mse)
  yvar <- var(Y0)
  ev <- function(method, beta, ntar = NA_real_) {
    m <- get_err(X0, Y0, as.vector(beta), q, length(Y0))
    data.table(method = method, mse = m, r2 = 1 - m / yvar, ntar = as.numeric(ntar))
  }

  fam <- rbindlist(list(
    ev("DORM (tuned smax)",  beta_dorm_tuned),
    ev("DORM (smax=0.03)",   beta_dorm_s03),
    ev("RAP",                P),
    ev("MI",                 bMI),
    ev("Best single source", Q[, best_src]),
    ev("Simple average",     rowMeans(Q)),
    ev("Rho average",        as.vector(Q %*% rho))
  ))

  ## transfer benchmarks on the SAME A (apples-to-apples), best over ntar
  bench <- data.table()
  if (RUN_TL) {
    ntar_grid <- c(20, 50, 100, 200); ntar_grid <- ntar_grid[ntar_grid <= length(Y0train)]
    run1 <- function(m, nt) tryCatch({
      if (m == "TransLasso") { f <- mytranslasso(X0train, Xlist, Y0train, Ylist, nt, q); b <- f$beta_TL }
      else if (m == "TransGLM") { tr <- List2TransGLM(X0train, Xlist, Y0train, Ylist, nt, q)
        fit <- glmtrans(tr$target, tr$source, family = "gaussian"); b <- as.vector(fit$beta)
        if (length(b) == q + 2) b <- b[-2] }
      else { f <- PTL_algorithm(Xlist, Ylist, X0train, Y0train, nt); b <- f$beta_PTL }
      ev(m, b, nt)
    }, error = function(e) data.table(method = m, mse = NA_real_, r2 = NA_real_, ntar = as.numeric(nt)))
    bench <- rbindlist(lapply(c("TransLasso", "TransGLM", "PTL"),
      function(m) rbindlist(lapply(ntar_grid, function(nt) run1(m, nt)))), fill = TRUE)
    bench <- bench[!is.na(r2), .SD[which.max(r2)], by = method]  # best ntar per method
  }

  summary <- rbindlist(list(fam, bench), fill = TRUE)
  summary[, `:=`(protein_y = y_col, target = target, k = k, q = q,
                 best_smax = tune$best_smax, n0_eval = nrow(X0))]

  ## ---- Q tracking (per-source betas over A) ----
  qtrack <- data.table(
    source = SOURCES,
    Q_L2norm = round(sqrt(colSums(Q^2)), 4),
    standalone_target_r2 = round(1 - src_mse / var(Y0), 4),
    rho = round(as.vector(rho), 4)
  )
  Qmat <- as.data.table(Q); setnames(Qmat, SOURCES)
  Qmat[, coef := c("Intercept", A_genes)]; setcolorder(Qmat, "coef")

  stub <- sprintf("v2_%s_target_%s_k%d", safe(y_col), safe(target), k)
  fwrite(summary, file.path(res_dir, paste0(stub, "_summary.csv")))
  fwrite(qtrack,  file.path(res_dir, paste0(stub, "_Qtrack.csv")))
  fwrite(Qmat,    file.path(res_dir, paste0(stub, "_Qmatrix.csv")))
  saveRDS(list(A_genes = A_genes, beta_array = beta_array, smax_grid = SMAX_GRID,
               best_smax = tune$best_smax, tune_losses = tune$losses, Q = Q, P = P,
               beta_MI = bMI, rho = rho, summary = summary, qtrack = qtrack),
          file.path(res_dir, paste0(stub, "_result.rds")))
  print(summary[order(-r2), .(method, mse = round(mse, 3), r2 = round(r2, 3), ntar)])
  message("  best_smax = ", tune$best_smax, "  -> saved ", stub)
  summary
}

## ---- run grid ----
all <- list()
for (y_col in ADTS) for (target in TARGETS) for (k in KS) {
  message(sprintf("\n=== %s | %s | k=%d ===", y_col, target, k))
  all[[length(all) + 1]] <- tryCatch(run_combo(y_col, target, k),
    error = function(e) { message("  ERROR: ", conditionMessage(e)); NULL })
}
alldt <- rbindlist(all, fill = TRUE)
fwrite(alldt, file.path(res_dir, "v2_all_summary.csv"))
message("\nDONE. Combined: ", file.path(res_dir, "v2_all_summary.csv"))
