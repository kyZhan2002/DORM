## Shared library for the DORM "advantage" analysis.
## Factored from experiments/2026-07-03_dorm_lowdimA/code/run_dorm_erythroid_v2.R,
## plus new label-free / label-using baselines. Sourced by the run scripts.

suppressPackageStartupMessages({ library(data.table); library(glmnet) })

## Full smax grid; label-free reports a fixed band, small-label tunes over all of it.
SMAX_GRID     <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9)
SMAX_FREE_BAND <- c(0.1, 0.25, 0.5)   # label-free fixed smax values reported as a band
NTAR_FIXED    <- 20                    # small labeled-target budget for the label-using tier

find_repo <- function(start) {
  d <- start; while (!file.exists(file.path(d, "src", "Functions3.R")) && dirname(d) != d) d <- dirname(d); d
}
## Ridge for the signed-affine (intrinsic) prior; larger => shrink signed weights toward uniform
## => more stable (default in IntrinsicSources.R is 0.1). Read dynamically by the wrapper below.
SIGNED_LAMBDA <- as.numeric(Sys.getenv("DORM_SIGNED_LAMBDA", "0.1"))

load_dorm_src <- function(repo) {
  source(file.path(repo, "src", "Functions3.R"))
  source(file.path(repo, "src", "IntrinsicSources.R"), local = FALSE)   # signed_affine_prior()
  source(file.path(repo, "src", "TransLasso-functions.R"))
  source(file.path(repo, "src", "TransLG-functions.R"))
  source(file.path(repo, "src", "Tuning.R"))
  ## Override signed_affine_prior's lambda default with SIGNED_LAMBDA (read at call time), so
  ## maximin_s_beta's plain call signed_affine_prior(dr, constr) uses our stabilizing ridge.
  .orig_sap <- signed_affine_prior
  assign("signed_affine_prior", function(dr_obj, dr_constr = NULL, lambda = NULL, eps = 1e-3) {
    if (is.null(lambda)) lambda <- get("SIGNED_LAMBDA", envir = .GlobalEnv)
    .orig_sap(dr_obj, dr_constr, lambda = lambda, eps = eps)
  }, envir = .GlobalEnv)
}

load_cite_data <- function(csv_dir) {
  meta <- fread(file.path(csv_dir, "cite_nips_obs_metadata_for_R.csv"))
  fmap <- fread(file.path(csv_dir, "cite_nips_feature_name_map_for_R.csv"))
  hvg  <- fread(file.path(csv_dir, "cite_nips_hvg_genes_for_R.csv"))
  rna_features <- intersect(hvg[!is.na(csv_name), csv_name], fmap[modality == "RNA", csv_name])
  adt_features <- fmap[modality == "ADT", csv_name]
  mc <- c("obs_id", "batch", "cell_type", "cell_type_major")
  rna_dt <- fread(file.path(csv_dir, "cite_nips_RNA_for_R.csv"), select = c(mc, rna_features))
  adt_dt <- fread(file.path(csv_dir, "cite_nips_ADT_for_R.csv"), select = c(mc, adt_features))
  stopifnot(identical(rna_dt$obs_id, adt_dt$obs_id), identical(rna_dt$obs_id, meta$obs_id))
  list(meta = meta, rna = rna_dt, adt = adt_dt, rna_features = rna_features, adt_features = adt_features)
}

sample_at_most <- function(idx, nmax) if (length(idx) <= nmax) as.integer(idx) else sort(sample(as.integer(idx), nmax))
split_half <- function(idx) { idx <- sample(as.integer(idx)); n <- floor(length(idx)/2); list(a = idx[seq_len(n)], b = idx[n + seq_len(n)]) }
safe <- function(x) gsub("[^0-9A-Za-z_.-]+", "_", x)

## Label-free gene screen: elastic-net of Y on all HVG using SOURCE rows only,
## rank by |coef|, greedy-drop near-duplicates so Sigma stays invertible. Returns k gene names.
screen_A <- function(D, src_rows, y_col, kmax) {
  rf <- D$rna_features
  X <- as.matrix(D$rna[src_rows, ..rf]); storage.mode(X) <- "double"
  y <- as.numeric(D$adt[[y_col]][src_rows])
  fit <- cv.glmnet(X, y, alpha = 0.5, nfolds = 5, standardize = FALSE)
  b <- as.vector(coef(fit, s = "lambda.min"))[-1]
  nz <- which(abs(b) > 0)
  ord <- c(nz[order(abs(b[nz]), decreasing = TRUE)], setdiff(order(abs(b), decreasing = TRUE), nz))
  picked <- integer(0)
  for (j in ord) {
    if (length(picked) == 0 || max(abs(cor(X[, j], X[, picked, drop = FALSE]))) < 0.999) {
      picked <- c(picked, j); if (length(picked) == kmax) break
    }
  }
  D$rna_features[picked]
}

## X = [intercept, A-genes, W-genes]; q = length(A).
make_X <- function(D, rows, A_genes) {
  cols <- c(A_genes, setdiff(D$rna_features, A_genes))
  X <- as.matrix(D$rna[rows, ..cols]); storage.mode(X) <- "double"
  cbind(Intercept = 1, X)
}
make_Y <- function(D, rows, y_col) as.numeric(D$adt[[y_col]][rows])

## Learned embedding (e.g. Harmony) for the intrinsic-source prior space. Attaches D$emb, a
## cell x k matrix aligned to D$rna/D$meta obs_id order (mirrors load_cite_data's obs_id join).
load_embedding <- function(D, csv_dir, prefix = "H", k = 50, file = "cite_nips_EMBED_for_R.csv") {
  emb <- fread(file.path(csv_dir, file))
  idx <- match(D$meta$obs_id, emb$obs_id); stopifnot(!any(is.na(idx)))
  cols <- paste0(prefix, seq_len(k)); D$emb <- as.matrix(emb[idx, ..cols]); D
}
## Embedding design matrix [intercept, embedding dims] for the given rows (analogue of make_X).
make_X_embed <- function(D, rows, k = ncol(D$emb)) cbind(Intercept = 1, D$emb[rows, seq_len(k), drop = FALSE])

combine_beta <- function(P, Q, Sigma0, smax) {
  L <- ncol(Q); dts <- opt_s_delta(P, Q, Sigma0, smax); s <- dts[L + 1]
  as.vector((1 - s) * P + s * Q %*% dts[1:L])
}

## Two cross-fit halves; exact beta_star for every smax in SMAX_GRID.
## penalty=TRUE regularises Q via elastic-net (bounded on rank-deficient / homogeneous
## targets); penalty=FALSE is the unregularised normal-equations path (can blow up on
## single-subtype targets even with the ridge). Default TRUE for robustness.
## prior_* (optional): parallel covariate splits for the intrinsic-source prior space (e.g. Harmony
## embedding), built on the SAME row halves as the gene splits. When supplied with signed_prior=TRUE,
## the signed-affine prior rho is estimated there while the predictor/posterior keep the gene space.
dorm_fit <- function(Xlist, Xtrainlist, Ylist, Ytrainlist, X0, X0train, nlist, q, penalty = TRUE,
                     signed_prior = FALSE,
                     prior_Xlist = NULL, prior_Xtrainlist = NULL, prior_X0 = NULL, prior_X0train = NULL) {
  o1 <- maximin_s_beta(Xlist, Xtrainlist, Ylist, Ytrainlist, X0, X0train, nlist, q,
                       smax = 0.5, penalty = penalty, alpha = 0.5, rho_pseudo = "NA",
                       dr_type = "logit", normalize = FALSE, condA = FALSE, signed_prior = signed_prior,
                       prior_Xlist = prior_Xlist, prior_Xtrainlist = prior_Xtrainlist,
                       prior_X0 = prior_X0, prior_X0train = prior_X0train)
  o2 <- maximin_s_beta(Xtrainlist, Xlist, Ytrainlist, Ylist, X0train, X0, nlist, q,
                       smax = 0.5, penalty = penalty, alpha = 0.5, rho_pseudo = "NA",
                       dr_type = "logit", normalize = FALSE, condA = FALSE, signed_prior = signed_prior,
                       prior_Xlist = prior_Xtrainlist, prior_Xtrainlist = prior_Xlist,
                       prior_X0 = prior_X0train, prior_X0train = prior_X0)
  S1 <- (1/nrow(X0))      * t(X0[, 1:(q+1)])      %*% X0[, 1:(q+1)]
  S2 <- (1/nrow(X0train)) * t(X0train[, 1:(q+1)]) %*% X0train[, 1:(q+1)]
  beta_by_smax <- t(sapply(SMAX_GRID, function(sm)
    0.5 * (combine_beta(o1$beta_RAP, o1$DoublyR, S1, sm) +
           combine_beta(o2$beta_RAP, o2$DoublyR, S2, sm))))
  list(beta_by_smax = beta_by_smax,
       Q = 0.5 * (o1$DoublyR + o2$DoublyR), P = 0.5 * (o1$beta_RAP + o2$beta_RAP),
       beta_MI = 0.5 * (o1$beta_MI + o2$beta_MI), rho = 0.5 * (o1$rho + o2$rho))
}

## Pooled elastic-net over A on pooled sources (label-free baseline).
pooled_en <- function(Xlist, Ylist, q) {
  Xp <- do.call(rbind, lapply(Xlist, function(X) X[, 2:(q+1), drop = FALSE]))
  yp <- unlist(Ylist, use.names = FALSE)
  fit <- cv.glmnet(Xp, yp, alpha = 0.5, nfolds = 5, standardize = FALSE)
  as.vector(coef(fit, s = "lambda.min"))   # length q+1: [intercept, A coefs]
}
## Target-oracle elastic-net (LABEL-USING upper bound / regret reference).
target_oracle_en <- function(X0, Y0, q) {
  fit <- cv.glmnet(X0[, 2:(q+1), drop = FALSE], Y0, alpha = 0.5, nfolds = 5, standardize = FALSE)
  as.vector(coef(fit, s = "lambda.min"))
}

## Tuned-smax R^2 for a dorm_fit (same logic as eval_all_methods' DORM row); reused for DORM_intrinsic.
dorm_tuned_r2 <- function(fit, X0, Y0, X0t, Y0t, q, ntune = 100) {
  nt <- min(ntune, length(Y0t))
  tn <- tuning_Y(X0t[1:nt, , drop = FALSE], Y0t[1:nt], q, fit$beta_by_smax, SMAX_GRID)
  beta <- fit$beta_by_smax[which(SMAX_GRID == tn$best_smax)[1], ]
  list(r2 = 1 - get_err(X0, Y0, as.vector(beta), q, length(Y0)) / var(Y0), best_smax = tn$best_smax)
}

## Evaluate all methods on (X0,Y0). Returns data.table(method, r2, tier).
## X0train/Y0train supply the ~100-label smax tuning split and the ntar for transfer methods.
eval_all_methods <- function(fit, Xlist, Ylist, X0, Y0, X0train, Y0train, nlist, q, ntune = 100) {
  yv <- var(Y0)
  R2 <- function(beta) 1 - get_err(X0, Y0, as.vector(beta), q, length(Y0)) / yv
  rows <- list()
  add <- function(method, beta, tier) rows[[length(rows)+1]] <<- data.table(method = method, r2 = R2(beta), tier = tier)

  ## ---- DORM: label-free fixed-smax band + small-label tuned ----
  for (sm in SMAX_FREE_BAND) add(sprintf("DORM (smax=%.2g)", sm),
                                 fit$beta_by_smax[which.min(abs(SMAX_GRID - sm)), ], "label-free")
  ntune <- min(ntune, length(Y0train))
  tune <- tuning_Y(X0train[1:ntune, , drop = FALSE], Y0train[1:ntune], q, fit$beta_by_smax, SMAX_GRID)
  best_i <- which(SMAX_GRID == tune$best_smax)[1]
  add("DORM (tuned smax)", fit$beta_by_smax[best_i, ], "label-free")

  ## ---- other label-free aggregations ----
  src_mse <- apply(fit$Q, 2, function(b) get_err(X0, Y0, b, q, length(Y0)))
  add("MI", fit$beta_MI, "label-free")
  add("RAP", fit$P, "label-free")
  add("Best single source", fit$Q[, which.min(src_mse)], "label-free")
  add("Simple average", rowMeans(fit$Q), "label-free")
  add("Rho average", as.vector(fit$Q %*% fit$rho), "label-free")
  add("Pooled elastic-net", pooled_en(Xlist, Ylist, q), "label-free")
  add("Target-mean (null)", c(mean(Y0), rep(0, q)), "label-free")

  ## ---- label-using tier (context): FIXED modest labeled-target budget ----
  nt <- min(NTAR_FIXED, length(Y0train))
  XlA <- lapply(Xlist, function(x) x[, 1:(q+1), drop = FALSE])   # A-only (PTL else uses full X)
  tl_r2 <- function(m) tryCatch({
    if (m == "TransLasso") b <- mytranslasso(X0train, Xlist, Y0train, Ylist, nt, q)$beta_TL
    else if (m == "TransGLM") { tr <- List2TransGLM(X0train, Xlist, Y0train, Ylist, nt, q)
      b <- as.vector(glmtrans(tr$target, tr$source, family = "gaussian")$beta); if (length(b)==q+2) b <- b[-2] }
    else b <- PTL_algorithm(XlA, Ylist, X0train[, 1:(q+1), drop = FALSE], Y0train, nt)$beta_PTL
    R2(b)
  }, error = function(e) NA_real_)
  for (m in c("TransLasso", "TransGLM", "PTL"))
    rows[[length(rows)+1]] <- data.table(method = m, r2 = tl_r2(m), tier = "label-using")
  add("Target-oracle EN", target_oracle_en(X0, Y0, q), "label-using")

  out <- rbindlist(rows); out[, best_smax := tune$best_smax]; out
}
