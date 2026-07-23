#!/usr/bin/env Rscript

## CHECK 1 diagnostic: why is the "Simple average" beta (rowMeans of the per-source
## doubly-robust matrix DoublyR = Q) so much worse than the other DORM-family betas?
## Hypothesis: one source's Q column is ill-conditioned / large-norm (a near-singular
## target design Sigma0 = X0'X0/n0 amplifies it), and equal 1/L weighting lets that
## single column dominate rowMeans, while MI/DORM down-weight it via the simplex delta.
##
## Usage: Rscript diagnose_simple_average.R <stub>
##   e.g. Rscript diagnose_simple_average.R Ver4_erythroid_subtype_sources_target_MK_E_prog_Y_CD71

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stub <- if (length(args) >= 1) args[1] else
  "Ver4_erythroid_subtype_sources_target_MK_E_prog_Y_CD71"

sc_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))[1])),
  error = function(e) getwd())
repo <- dirname(sc_dir)
source(file.path(repo, "src", "Functions3.R"))
res_dir <- file.path(sc_dir, "dorm_results")

data <- readRDS(file.path(res_dir, paste0(stub, "_dorm_data.rds")))
res  <- readRDS(file.path(res_dir, paste0(stub, "_dorm_result.rds")))

X0 <- data$X0; Y0 <- data$Y0; q <- data$q
n0 <- nrow(X0)
Q  <- res$DoublyR                 # (q+1) x L : per-source doubly-robust betas
L  <- ncol(Q)
Sigma0 <- (1 / n0) * t(X0[, 1:(q + 1)]) %*% X0[, 1:(q + 1)]

cat("\n==== ", stub, " ====\n", sep = "")
cat(sprintf("q+1 = %d,  target n0 = %d,  L sources = %d\n", q + 1, n0, L))

## Target design conditioning (the amplifier in the un-penalized path).
ev <- eigen(Sigma0, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("Sigma0 eigenvalues: max=%.3g  min=%.3g  condition=%.3g  (#<1e-8: %d)\n",
            max(ev), min(ev), max(ev) / max(min(ev), .Machine$double.eps),
            sum(ev < 1e-8)))

## Per-source column norms and their standalone target MSE.
col_norm <- sqrt(colSums(Q^2))
col_mse  <- apply(Q, 2, function(b) get_err(X0, Y0, as.vector(b), q, length(Y0)))
src_names <- names(data$Xlist); if (is.null(src_names)) src_names <- paste0("src", 1:L)

diag_tab <- data.table(source = src_names,
                       col_L2norm = round(col_norm, 3),
                       standalone_mse = round(col_mse, 3),
                       rho = round(as.vector(res$rho), 3),
                       delta_MI = round(as.vector(res$deltas)[1:L], 3))
cat("\nPer-source DoublyR columns:\n"); print(diag_tab, row.names = FALSE)

## Compare the aggregate betas.
beta_simple <- rowMeans(Q)
beta_rho    <- as.vector(Q %*% as.vector(res$rho))
compare <- data.table(
  beta = c("Simple average (rowMeans Q)", "MI (Q %*% delta_MI)",
           "DORM (beta_star)", "RAP", "Rho average (Q %*% rho)"),
  L2norm = round(c(sqrt(sum(beta_simple^2)), sqrt(sum(res$beta_MI^2)),
                   sqrt(sum(res$beta_star^2)), sqrt(sum(res$beta_RAP^2)),
                   sqrt(sum(beta_rho^2))), 3),
  target_mse = round(c(get_err(X0, Y0, beta_simple, q, n0),
                       get_err(X0, Y0, as.vector(res$beta_MI), q, n0),
                       get_err(X0, Y0, as.vector(res$beta_star), q, n0),
                       get_err(X0, Y0, as.vector(res$beta_RAP), q, n0),
                       get_err(X0, Y0, beta_rho, q, n0)), 3))
compare[, target_r2 := round(1 - target_mse / var(Y0), 3)]
cat("\nAggregate betas vs target:\n"); print(compare, row.names = FALSE)

## Is one column dominating rowMeans? Leave-one-source-out of the simple average.
cat("\nLeave-one-source-out simple average (drop each source, average the rest):\n")
loo <- data.table(dropped = src_names,
  mse = round(sapply(1:L, function(l)
    get_err(X0, Y0, rowMeans(Q[, -l, drop = FALSE]), q, n0)), 3))
loo[, r2 := round(1 - mse / var(Y0), 3)]
print(loo, row.names = FALSE)
cat("\n(If dropping ONE source restores a sane MSE, that source's column is the culprit.)\n")
