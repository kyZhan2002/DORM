#!/usr/bin/env Rscript
## Phase C: extrapolation validation of DORM-intrinsic on the Harmony-demixed structure.
## Idea: donor-LODO only tests targets INSIDE the source-donor hull, where the enlarged (signed-affine)
## set gives no benefit. Here we construct targets that LEAVE the hull and test whether DORM-intrinsic
## (signed-affine prior estimated in the Harmony embedding) beats standard DORM as the target
## extrapolates -- judged against the full benchmark suite. Fully annotation-free:
##   * "states" = unsupervised k-means clusters in the Harmony embedding (NOT cell-type labels)
##   * sources  = the other donors (natural annotation-free mixtures)
##   * target   = a held-out donor resampled to composition (1-a)*natural + a*(pure state c), a in ALPHAS
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

env_csv <- function(n, d) { v <- Sys.getenv(n); if (v == "") d else trimws(strsplit(v, ",")[[1]]) }
ADTS    <- env_csv("DORM_ADTS", c("CD72_1", "CD3"))
DONORS  <- env_csv("DORM_DONORS", c("s3d7", "s4d9"))       # held-out targets to shift
ALPHAS  <- as.numeric(env_csv("DORM_ALPHAS", c("0", "0.25", "0.5", "0.75", "1")))
KM      <- as.integer(Sys.getenv("DORM_KMEANS", "8"))      # unsupervised Harmony states
HVG     <- as.integer(Sys.getenv("DORM_HVG", "250"))
MAXN_SRC <- as.integer(Sys.getenv("DORM_MAXN_SRC", "500"))
N_TGT   <- as.integer(Sys.getenv("DORM_NTGT", "1200"))     # target size (split into two halves)
K <- 20; DRY <- Sys.getenv("DORM_DRY", "") != ""; site_of <- function(b) substr(b, 1, 2)

D <- load_cite_data(csv_dir); D$rna_features <- head(D$rna_features, HVG); D <- load_embedding(D, csv_dir)
ALL_DONORS <- sort(unique(D$meta$batch))

## Unsupervised states in the Harmony embedding (annotation-free proxies for cell states).
set.seed(1); km <- kmeans(D$emb, centers = KM, iter.max = 50, nstart = 5); klab <- km$cluster
cat(sprintf("k-means: %d states on Harmony embedding; sizes: %s\n", KM, paste(tabulate(klab, KM), collapse = ",")))

## Composition of a set of cells over the KM states.
comp_of <- function(rows) tabulate(klab[rows], KM) / length(rows)
## Build a target of size N with state-composition (1-alpha)*natural(donor) + alpha*e_c, drawing REAL
## held-out-donor cells (with replacement) so any composition is reachable. Returns row indices.
build_target <- function(d_idx, c, alpha, N) {
  p_nat <- comp_of(d_idx); p_tgt <- (1 - alpha) * p_nat; p_tgt[c] <- p_tgt[c] + alpha
  n_j <- round(N * p_tgt); rows <- integer(0)
  for (j in which(n_j > 0)) {
    pool <- d_idx[klab[d_idx] == j]
    if (length(pool)) rows <- c(rows, sample(pool, n_j[j], replace = TRUE))
  }
  rows
}

out <- file.path(res_dir, "extrapolation_shift_results.csv")
if (Sys.getenv("DORM_APPEND", "") == "" && file.exists(out)) file.remove(out)

for (dtar in DONORS) {
  d_idx <- which(D$meta$batch == dtar)
  c <- which.max(tabulate(klab[d_idx], KM))               # donor's dominant state (guaranteed present)
  src_donors <- setdiff(ALL_DONORS, dtar)
  src_comp_c <- max(sapply(src_donors, function(b) comp_of(which(D$meta$batch == b))[c]))  # max source frac in state c
  cat(sprintf("\n== target %s | dominant state c=%d (donor frac %.2f, max source frac %.2f) ==\n",
              dtar, c, comp_of(d_idx)[c], src_comp_c))

  for (y in ADTS) {
    set.seed(1)
    src_rows <- lapply(src_donors, function(b) sample_at_most(which(D$meta$batch == b), MAXN_SRC))
    names(src_rows) <- src_donors
    A <- screen_A(D, unlist(src_rows, use.names = FALSE), y, K); q <- K
    ssp <- lapply(src_rows, split_half)
    Xl <- lapply(ssp, function(s) make_X(D, s$a, A)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A))
    Yl <- lapply(ssp, function(s) make_Y(D, s$a, y)); Yt <- lapply(ssp, function(s) make_Y(D, s$b, y))
    nlist <- as.integer(lengths(Yl))
    pXl <- lapply(ssp, function(s) make_X_embed(D, s$a)); pXt <- lapply(ssp, function(s) make_X_embed(D, s$b))

    for (alpha in ALPHAS) {
      ra <- build_target(d_idx, c, alpha, N_TGT %/% 2); rb <- build_target(d_idx, c, alpha, N_TGT %/% 2)
      p_tgt_c <- comp_of(ra)[c]; excess_c <- p_tgt_c - src_comp_c   # >0 => outside source hull in state c
      X0  <- make_X(D, ra, A); Y0  <- make_Y(D, ra, y); X0t  <- make_X(D, rb, A); Y0t  <- make_Y(D, rb, y)
      pX0 <- make_X_embed(D, ra);                       pX0t <- make_X_embed(D, rb)
      if (DRY) { cat(sprintf("  %s %s a=%.2f: target state-c frac=%.2f excess=%+.2f (N=%d)\n",
                             dtar, y, alpha, p_tgt_c, excess_c, length(ra))); next }

      fit <- tryCatch(dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nlist, q, penalty = TRUE, signed_prior = FALSE),
                      error = function(e) { message("  std err: ", conditionMessage(e)); NULL })
      if (is.null(fit)) next
      ev <- eval_all_methods(fit, Xl, Yl, X0, Y0, X0t, Y0t, nlist, q)
      fit_i <- tryCatch(dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nlist, q, penalty = TRUE, signed_prior = TRUE,
                                 prior_Xlist = pXl, prior_Xtrainlist = pXt, prior_X0 = pX0, prior_X0train = pX0t),
                        error = function(e) { message("  int err: ", conditionMessage(e)); NULL })
      if (!is.null(fit_i)) {
        di <- dorm_tuned_r2(fit_i, X0, Y0, X0t, Y0t, q)
        ev <- rbind(ev, data.table(method = "DORM_intrinsic", r2 = di$r2, tier = "label-free", best_smax = di$best_smax))
      }
      ev[, `:=`(donor = dtar, ADT = y, state_c = c, alpha = alpha,
                tgt_frac_c = round(p_tgt_c, 3), excess_c = round(excess_c, 3), n_target = length(Y0))]
      fwrite(ev, out, append = file.exists(out))
      cat(sprintf("  %s %s a=%.2f (excess=%+.2f): DORM=%.3f DORM_int=%.3f Pooled=%.3f oracle=%.3f\n",
                  dtar, y, alpha, excess_c, ev[method == "DORM (tuned smax)", r2],
                  if (!is.null(fit_i)) di$r2 else NA_real_, ev[method == "Pooled elastic-net", r2],
                  ev[method == "Target-oracle EN", r2]))
    }
  }
}
if (!DRY) message("DONE extrapolation_shift -> ", out)
