#!/usr/bin/env Rscript
## Leave-one-donor-out (donors as sources), GEX->ADT, comparing standard DORM vs DORM-intrinsic
## (signed-affine prior) against the benchmark suite. See 2026-07-16 study for the base design.
## For each held-out donor: L sources = other donors (all cells); target = held-out donor.
## DORM q=20, penalty=TRUE, tuned smax; label-using benchmarks at ntar=20 (set in lib).
## Env: DORM_ADTS=CD72_1,CD19_1,CD11c,CD3  DORM_SEEDS=1,2  DORM_DONORS=<subset>  DORM_HVG=250
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

env_csv <- function(n, d) { v <- Sys.getenv(n); if (v == "") d else trimws(strsplit(v, ",")[[1]]) }
ADTS   <- env_csv("DORM_ADTS", c("CD72_1", "CD19_1", "CD11c", "CD3"))
SEEDS  <- as.integer(env_csv("DORM_SEEDS", c("1", "2")))
HVG    <- as.integer(Sys.getenv("DORM_HVG", "250"))
MAXN_SRC <- as.integer(Sys.getenv("DORM_MAXN_SRC", "500"))
MAXN_TGT <- as.integer(Sys.getenv("DORM_MAXN_TGT", "1500"))
K <- 20; site_of <- function(b) substr(b, 1, 2)
PRIOR_SPACE <- Sys.getenv("DORM_PRIOR_SPACE", "")   # "harmony" => estimate the signed-affine prior in the Harmony embedding

D <- load_cite_data(csv_dir)
D$rna_features <- head(D$rna_features, HVG)
if (PRIOR_SPACE == "harmony") D <- load_embedding(D, csv_dir)   # attaches D$emb (Harmony dims) for the prior space
ALL_DONORS <- sort(unique(D$meta$batch))
DONORS <- env_csv("DORM_DONORS", ALL_DONORS)
message("Donors (", length(DONORS), "): ", paste(DONORS, collapse = ", "))

tag  <- if (PRIOR_SPACE == "harmony") "_harmony" else ""
out  <- file.path(res_dir, sprintf("lodo_donor_intrinsic%s_results.csv", tag))
wout <- file.path(res_dir, sprintf("lodo_donor_intrinsic%s_rho.csv", tag))
if (Sys.getenv("DORM_APPEND", "") == "") for (f in c(out, wout)) if (file.exists(f)) file.remove(f)

## tuned-smax R^2 for a fit (same logic as eval_all_methods' DORM row)
dorm_tuned_r2 <- function(fit, X0, Y0, X0t, Y0t, q, ntune = 100) {
  nt <- min(ntune, length(Y0t))
  tn <- tuning_Y(X0t[1:nt, , drop = FALSE], Y0t[1:nt], q, fit$beta_by_smax, SMAX_GRID)
  beta <- fit$beta_by_smax[which(SMAX_GRID == tn$best_smax)[1], ]
  list(r2 = 1 - get_err(X0, Y0, as.vector(beta), q, length(Y0)) / var(Y0), best_smax = tn$best_smax)
}

for (dtar in DONORS) for (y in ADTS) for (seed in SEEDS) {
  set.seed(seed)
  src_donors <- setdiff(ALL_DONORS, dtar)
  src_rows <- lapply(src_donors, function(b) sample_at_most(which(D$meta$batch == b), MAXN_SRC))
  names(src_rows) <- src_donors
  src_rows <- src_rows[lengths(src_rows) >= 60]; src_donors <- names(src_rows)
  tgt_rows <- sample_at_most(which(D$meta$batch == dtar), MAXN_TGT)
  if (length(tgt_rows) < 150) { message("  skip ", dtar, ": too few target cells"); next }

  A <- screen_A(D, unlist(src_rows, use.names = FALSE), y, K); q <- K
  ssp <- lapply(src_rows, split_half); tsp <- split_half(tgt_rows)
  Xl <- lapply(ssp, function(s) make_X(D, s$a, A)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A))
  Yl <- lapply(ssp, function(s) make_Y(D, s$a, y)); Yt <- lapply(ssp, function(s) make_Y(D, s$b, y))
  names(Xl) <- names(Xt) <- names(Yl) <- names(Yt) <- src_donors; nlist <- as.integer(lengths(Yl))
  X0 <- make_X(D, tsp$a, A); Y0 <- make_Y(D, tsp$a, y); X0t <- make_X(D, tsp$b, A); Y0t <- make_Y(D, tsp$b, y)
  if (PRIOR_SPACE == "harmony") {                                  # parallel Harmony splits on the SAME rows
    pXl <- lapply(ssp, function(s) make_X_embed(D, s$a)); pXt <- lapply(ssp, function(s) make_X_embed(D, s$b))
    pX0 <- make_X_embed(D, tsp$a); pX0t <- make_X_embed(D, tsp$b)
  } else { pXl <- pXt <- NULL; pX0 <- pX0t <- NULL }

  fit <- tryCatch(dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nlist, q, penalty = TRUE, signed_prior = FALSE),
                  error = function(e) { message("  std fit err ", dtar, "/", y, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  ev <- eval_all_methods(fit, Xl, Yl, X0, Y0, X0t, Y0t, nlist, q)   # DORM + all benchmarks

  fit_i <- tryCatch(dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nlist, q, penalty = TRUE, signed_prior = TRUE,
                             prior_Xlist = pXl, prior_Xtrainlist = pXt, prior_X0 = pX0, prior_X0train = pX0t),
                    error = function(e) { message("  int fit err ", dtar, "/", y, ": ", conditionMessage(e)); NULL })
  if (!is.null(fit_i)) {
    di <- dorm_tuned_r2(fit_i, X0, Y0, X0t, Y0t, q)
    ev <- rbind(ev, data.table(method = "DORM_intrinsic", r2 = di$r2, tier = "label-free", best_smax = di$best_smax))
  }
  ev[, `:=`(donor = dtar, site = site_of(dtar), ADT = y, seed = seed, L = length(Xl), n_target = length(Y0))]
  fwrite(ev, out, append = file.exists(out))

  ss <- site_of(src_donors) == site_of(dtar)
  fwrite(data.table(target_donor = dtar, ADT = y, seed = seed, source_donor = src_donors, same_site = ss,
                    rho_std = round(as.vector(fit$rho), 4),
                    rho_int = if (!is.null(fit_i)) round(as.vector(fit_i$rho), 4) else NA_real_),
         wout, append = file.exists(wout))
  message(sprintf("  %s (%s) | %s | seed%d : DORM=%.3f DORM_int=%.3f Pooled=%.3f TransGLM=%.3f oracle=%.3f",
                  dtar, site_of(dtar), y, seed, ev[method=="DORM (tuned smax)", r2],
                  if (!is.null(fit_i)) di$r2 else NA, ev[method=="Pooled elastic-net", r2],
                  ev[method=="TransGLM", r2], ev[method=="Target-oracle EN", r2]))
}
message("DONE lodo_donor_intrinsic -> ", out)
