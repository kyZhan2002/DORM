#!/usr/bin/env Rscript

## Build DORM inputs from the CITE-NIPS RNA/ADT CSV exports.
##
## Goal:
##   X = processed RNA HVG features
##   Y = one ADT protein. Default is manually chosen as CD72_1; set
##       DORM_PROTEIN_Y=auto to instead screen for the most RNA-predictable ADT.
##
## DORM data contract, following src/Example_new.R:
##   Xlist[[l]]      source-l covariates, matrix with intercept in column 1
##   Xtrainlist[[l]] auxiliary source-l covariates for nuisance fitting
##   Ylist[[l]]      scalar ADT response aligned to Xlist[[l]]
##   Ytrainlist[[l]] scalar ADT response aligned to Xtrainlist[[l]]
##   X0, X0train     target covariates; target Y is kept only for evaluation
##   nlist           number of labeled source observations used per source
##   q               length(beta) - 1; here q = number of selected RNA features
##
## Example runs:
##   Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Default Ver1: sources are batches, target is one held-out batch.
##   DORM_PARTITION=major_cell_type Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Ver2: sources are major cell types, target is one held-out major type.
##   DORM_PARTITION=t_cell_subtype Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Ver3: restrict to T cells; sources/target are fine T-cell subtypes.
##   DORM_SOURCE_GROUPS=s1d1,s1d2,s2d1 DORM_TARGET_GROUPS=s4d8,s4d9 Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Batch Ver1 with only selected source batches; run DORM separately for each target batch.
##   DORM_PROTEIN_Y=CD72_1 Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Force one ADT protein as Y instead of the default CD72_1.
##   DORM_PROTEIN_Y=auto Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Choose Y by quick cross-validated elastic-net screening over all ADT proteins.
##   DORM_BENCHMARK_NTAR=50,100,200 Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Add TransLasso, TransGLM, and PTL rows using these numbers of labeled target cells.
##   DORM_RUN_TRANSFER=false Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Skip transfer benchmarks and report only DORM-family methods.
##   DORM_RUN=false DORM_PROTEIN_Y=CD3 Rscript singlecell/run_dorm_from_cite_nips_csv.R
##     Build/save the DORM input object only; useful for a quick smoke test.

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

candidate_roots <- unique(normalizePath(c(getwd(), file.path(getwd(), "..")), mustWork = FALSE))
repo_dir <- candidate_roots[file.exists(file.path(candidate_roots, "src", "Functions3.R"))][1]
if (is.na(repo_dir)) {
  stop("Could not locate repo root containing src/Functions3.R. Run from repo root or singlecell/.", call. = FALSE)
}

source(file.path(repo_dir, "src", "Functions3.R"))
source(file.path(repo_dir, "src", "TransLasso-functions.R"))
source(file.path(repo_dir, "src", "TransLG-functions.R"))

singlecell_dir <- file.path(repo_dir, "singlecell")
csv_dir <- file.path(singlecell_dir, "processed", "csv_for_R")
out_dir <- file.path(singlecell_dir, "dorm_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  rna = file.path(csv_dir, "cite_nips_RNA_for_R.csv"),
  adt = file.path(csv_dir, "cite_nips_ADT_for_R.csv"),
  metadata = file.path(csv_dir, "cite_nips_obs_metadata_for_R.csv"),
  feature_map = file.path(csv_dir, "cite_nips_feature_name_map_for_R.csv"),
  hvg = file.path(csv_dir, "cite_nips_hvg_genes_for_R.csv")
)
invisible(lapply(paths, function(x) {
  if (!file.exists(x)) stop("Missing required file: ", x, call. = FALSE)
}))

## -----------------------------
## Configuration
## -----------------------------

set.seed(2026)

active_partition <- Sys.getenv("DORM_PARTITION", unset = "batch")
protein_y <- Sys.getenv("DORM_PROTEIN_Y", unset = "CD72_1")
if (identical(protein_y, "")) protein_y <- "CD72_1"

cfg <- list(
  max_rna_features = 500,          # use NULL for all exported HVGs
  min_rows_per_group = 200,
  max_source_rows_per_group = 12000,
  max_target_rows = 12000,
  target_train_fraction = 0.5,
  source_train_fraction = 0.5,
  protein_screen_n = 5000,
  protein_screen_folds = 5,
  smax = 0.03,
  penalty = TRUE,
  alpha = 0.5,
  rho_pseudo = "NA",
  dr_type = "logit",
  normalize = FALSE,
  condA = FALSE,
  run_dorm = TRUE,
  run_transfer_benchmarks = TRUE,
  benchmark_ntar = c(10, 20, 50, 100, 200)
)

env_flag <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  as.integer(value)
}

env_csv <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

env_int_csv <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || value == "") return(default)
  as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

cfg$run_dorm <- env_flag("DORM_RUN", cfg$run_dorm)
cfg$run_transfer_benchmarks <- env_flag("DORM_RUN_TRANSFER", cfg$run_transfer_benchmarks)
cfg$max_rna_features <- env_int("DORM_MAX_RNA_FEATURES", cfg$max_rna_features)
cfg$max_source_rows_per_group <- env_int("DORM_MAX_SOURCE_ROWS", cfg$max_source_rows_per_group)
cfg$max_target_rows <- env_int("DORM_MAX_TARGET_ROWS", cfg$max_target_rows)
cfg$protein_screen_n <- env_int("DORM_PROTEIN_SCREEN_N", cfg$protein_screen_n)
cfg$benchmark_ntar <- env_int_csv("DORM_BENCHMARK_NTAR", cfg$benchmark_ntar)

partition_configs <- list(
  ## Ver1: treat experimental batches as domains.
  ## DORM_SOURCE_GROUPS may contain any subset of batches; DORM_TARGET_GROUPS
  ## may contain one or more held-out target batches.
  batch = list(
    label = "Ver1_batch_sources",
    group_col = "batch",
    subset_fun = function(meta) rep(TRUE, nrow(meta)),
    # source_groups = c("s1d1", "s1d2", "s1d3", "s2d1", "s2d4", "s2d5",
    #                  "s3d1", "s3d6", "s3d7", "s4d1", "s4d8"),
    source_groups = c("s1d1", "s1d3", "s2d1", "s2d4", "s2d5",
                      "s3d1", "s3d6", "s3d7", "s4d1"),
    target_groups = c("s4d8","s1d2")
  ),

  ## Ver2: treat major cell types as domains.
  ## Here B cells are the held-out target; all listed other major types are sources.
  major_cell_type = list(
    label = "Ver2_major_cell_type_sources",
    group_col = "cell_type_major",
    subset_fun = function(meta) rep(TRUE, nrow(meta)),
    source_groups = c("T cell", "monocyte", "NK/ILC", "DC", "Progenitor", "Erythroid"),
    target_groups = "B cell"
  ),

  ## Ver3: restrict to T cells, then use fine T-cell subtypes as domains.
  ## This is the immune-subtype version; edit target_group to hold out a different subtype.
  t_cell_subtype = list(
    label = "Ver3_T_cell_subtype_sources",
    group_col = "cell_type",
    subset_fun = function(meta) meta$cell_type_major == "T cell",
    source_groups = c("CD4+ T activated", "CD4+ T naive", "CD8+ T naive",
                      "CD8+ T CD57+ CD45RO+", "CD8+ T TIGIT+ CD45RO+",
                      "CD8+ T CD57+ CD45RA+", "CD4+ T activated integrinB7+"),
    target_groups = "T reg"
  ),

  ## Ver4: erythroid developmental subtypes as domains, mirroring the design of
  ## the MIMAL BMMC real-data analysis (only the source/target/ADT DESIGN is borrowed;
  ## methods and benchmarks remain the DORM family). Restrict to the 5 red-cell types
  ## (11,948 cells total). Sources are the 3 nucleated stages; the two held-out targets
  ## are the terminal state (Reticulocyte) and the upstream progenitor (MK/E prog).
  ## Response ADTs to pair with this: CD71, CD36_1, CD105, CD49d.
  erythroid_subtype = list(
    label = "Ver4_erythroid_subtype_sources",
    group_col = "cell_type",
    subset_fun = function(meta) meta$cell_type %in% c(
      "Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte", "MK/E prog"),
    source_groups = c("Proerythroblast", "Erythroblast", "Normoblast"),
    target_groups = c("Reticulocyte", "MK/E prog")
  )
)

if (!active_partition %in% names(partition_configs)) {
  stop("Unknown DORM_PARTITION='", active_partition, "'. Choose one of: ",
       paste(names(partition_configs), collapse = ", "), call. = FALSE)
}
part_cfg <- partition_configs[[active_partition]]
source_groups_was_set <- Sys.getenv("DORM_SOURCE_GROUPS", unset = "") != ""
part_cfg$source_groups <- env_csv("DORM_SOURCE_GROUPS", part_cfg$source_groups)
part_cfg$target_groups <- env_csv("DORM_TARGET_GROUPS", part_cfg$target_groups)
if (!source_groups_was_set) {
  part_cfg$source_groups <- setdiff(part_cfg$source_groups, part_cfg$target_groups)
}

if (length(intersect(part_cfg$source_groups, part_cfg$target_groups)) > 0) {
  stop("Source and target groups must be disjoint: ",
       paste(intersect(part_cfg$source_groups, part_cfg$target_groups), collapse = ", "),
       call. = FALSE)
}

## -----------------------------
## Read exported CSVs
## -----------------------------

metadata_cols <- c("obs_id", "batch", "cell_type", "cell_type_major")

meta <- fread(paths$metadata)
feature_map <- fread(paths$feature_map)
hvg <- fread(paths$hvg)

rna_features <- hvg[!is.na(csv_name), csv_name]
rna_features <- intersect(rna_features, feature_map[modality == "RNA", csv_name])
if (!is.null(cfg$max_rna_features)) {
  rna_features <- head(rna_features, cfg$max_rna_features)
}
if (length(rna_features) < 2) stop("Need at least two RNA features for DORM.", call. = FALSE)

adt_features <- feature_map[modality == "ADT", csv_name]

message("Reading RNA columns: metadata + ", length(rna_features), " selected HVG features")
rna_dt <- fread(paths$rna, select = c(metadata_cols, rna_features), showProgress = TRUE)

message("Reading ADT columns: metadata + ", length(adt_features), " proteins")
adt_dt <- fread(paths$adt, select = c(metadata_cols, adt_features), showProgress = TRUE)

if (!identical(rna_dt$obs_id, adt_dt$obs_id)) {
  stop("RNA and ADT CSV rows are not aligned by obs_id.", call. = FALSE)
}
if (!identical(rna_dt$obs_id, meta$obs_id)) {
  stop("RNA CSV and metadata CSV rows are not aligned by obs_id.", call. = FALSE)
}

## -----------------------------
## Partition, build DORM data, and run each target
## -----------------------------

sample_at_most <- function(idx, nmax) {
  idx <- as.integer(idx)
  if (is.infinite(nmax) || length(idx) <= nmax) return(idx)
  sort(sample(idx, nmax))
}

split_two_equal_parts <- function(idx, fraction = 0.5) {
  idx <- sample(as.integer(idx))
  n1 <- floor(length(idx) * fraction)
  n1 <- max(1, min(n1, length(idx) - 1))
  n2 <- length(idx) - n1
  n <- min(n1, n2)
  list(first = idx[seq_len(n)], second = idx[n1 + seq_len(n)])
}

safe_file_part <- function(x) {
  gsub("[^0-9A-Za-z_.-]+", "_", x)
}

make_partition <- function(meta, part_cfg, target_group, cfg) {
  keep <- part_cfg$subset_fun(meta)
  if (!is.logical(keep) || length(keep) != nrow(meta)) {
    stop("partition subset_fun must return a logical vector with nrow(meta) entries.", call. = FALSE)
  }

  group <- meta[[part_cfg$group_col]]
  available <- unique(group[keep])
  missing_sources <- setdiff(part_cfg$source_groups, available)
  if (length(missing_sources) > 0) {
    stop("Configured source groups not found: ", paste(missing_sources, collapse = ", "), call. = FALSE)
  }
  if (!target_group %in% available) {
    stop("Configured target group not found: ", target_group, call. = FALSE)
  }

  source_idx <- lapply(part_cfg$source_groups, function(g) which(keep & group == g))
  names(source_idx) <- part_cfg$source_groups
  source_idx <- source_idx[lengths(source_idx) >= cfg$min_rows_per_group]
  target_idx <- which(keep & group == target_group)

  if (length(source_idx) < 2) stop("DORM needs at least two source groups after filtering.", call. = FALSE)
  if (length(target_idx) < cfg$min_rows_per_group) {
    stop("Target group has too few rows after filtering: ", length(target_idx), call. = FALSE)
  }

  source_idx <- lapply(source_idx, sample_at_most, nmax = cfg$max_source_rows_per_group)
  target_idx <- sample_at_most(target_idx, cfg$max_target_rows)

  list(source_idx = source_idx, target_idx = target_idx)
}

choose_protein_y <- function(rna_dt, adt_dt, row_idx, x_cols, y_cols,
                             sample_n = 5000, nfolds = 5, seed = 1) {
  set.seed(seed)
  screen_idx <- sample_at_most(row_idx, sample_n)
  x <- as.matrix(rna_dt[screen_idx, ..x_cols])
  storage.mode(x) <- "double"

  score_one <- function(y_col) {
    y <- as.numeric(adt_dt[[y_col]][screen_idx])
    if (sd(y) <= 1e-8 || anyNA(y)) {
      return(data.table(protein = y_col, cv_mse = NA_real_, cv_r2 = NA_real_))
    }
    fit <- tryCatch(
      cv.glmnet(x = x, y = y, family = "gaussian", alpha = 0.5,
                nfolds = nfolds, standardize = FALSE, nlambda = 40),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.table(protein = y_col, cv_mse = NA_real_, cv_r2 = NA_real_))
    }
    cv_mse <- min(fit$cvm)
    cv_r2 <- 1 - cv_mse / var(y)
    data.table(protein = y_col, cv_mse = cv_mse, cv_r2 = cv_r2)
  }

  scores <- rbindlist(lapply(y_cols, score_one))
  scores <- scores[order(-cv_r2)]
  if (nrow(scores[!is.na(cv_r2)]) == 0) {
    stop("Could not select a protein: all protein screens failed.", call. = FALSE)
  }
  scores
}

make_x <- function(dt, idx, x_cols) {
  x <- as.matrix(dt[idx, ..x_cols])
  storage.mode(x) <- "double"
  x <- cbind(Intercept = 1, x)
  colnames(x)[1] <- "Intercept"
  x
}

make_y <- function(dt, idx, y_col) {
  as.numeric(dt[[y_col]][idx])
}

coerce_beta <- function(beta, q, method) {
  beta <- as.vector(beta)
  if (length(beta) == q + 1) return(beta)
  if (method == "TransGLM" && length(beta) == q + 2) return(beta[-2])
  stop(method, " returned beta length ", length(beta), "; expected ", q + 1, call. = FALSE)
}

method_metrics <- function(method, beta, X_eval, Y_eval, q, detail = NA_character_) {
  beta <- coerce_beta(beta, q, method)
  mse <- get_err(X_eval, Y_eval, beta, q, length(Y_eval))
  y_var <- var(Y_eval)
  data.table(
    method = method,
    mse = mse,
    r2 = if (y_var > 0) 1 - mse / y_var else NA_real_,
    detail = detail
  )
}

run_transfer_benchmarks <- function(Xlist, Ylist, X0_labeled, Y0_labeled,
                                    X_eval, Y_eval, q, ntar_values) {
  ntar_values <- sort(unique(ntar_values[!is.na(ntar_values)]))
  ntar_values <- ntar_values[ntar_values >= 5 & ntar_values <= length(Y0_labeled)]
  if (length(ntar_values) == 0) {
    warning("No valid ntar values for transfer benchmarks; skipping them.", call. = FALSE)
    return(data.table())
  }

  run_one <- function(method, ntar) {
    tryCatch({
      if (method == "TransLasso") {
        fit <- mytranslasso(X0_labeled, Xlist, Y0_labeled, Ylist, ntar, q)
        out <- method_metrics(method, fit$beta_TL, X_eval, Y_eval, q,
                              detail = paste0("ntar=", ntar, "; train_mse=", signif(fit$mse_TL, 4)))
      } else if (method == "TransGLM") {
        training <- List2TransGLM(X0_labeled, Xlist, Y0_labeled, Ylist, ntar, q)
        fit <- glmtrans(training$target, training$source, family = "gaussian")
        beta <- coerce_beta(fit$beta, q, method)
        train_pred <- as.vector(X0_labeled[seq_len(ntar), 1:(q + 1), drop = FALSE] %*% beta)
        train_mse <- mean((Y0_labeled[seq_len(ntar)] - train_pred)^2)
        out <- method_metrics(method, beta, X_eval, Y_eval, q,
                              detail = paste0("ntar=", ntar, "; train_mse=", signif(train_mse, 4)))
      } else if (method == "PTL") {
        fit <- PTL_algorithm(Xlist, Ylist, X0_labeled, Y0_labeled, ntar)
        out <- method_metrics(method, fit$beta_PTL, X_eval, Y_eval, q,
                              detail = paste0("ntar=", ntar, "; train_mse=", signif(fit$mse_PTL, 4)))
      } else {
        stop("Unknown transfer benchmark: ", method)
      }
      out[, ntar := ntar]
      out
    }, error = function(e) {
      data.table(method = method, mse = NA_real_, r2 = NA_real_,
                 detail = paste0("ntar=", ntar, "; error: ", conditionMessage(e)),
                 ntar = ntar)
    })
  }

  rbindlist(lapply(c("TransLasso", "TransGLM", "PTL"), function(method) {
    rbindlist(lapply(ntar_values, function(ntar) run_one(method, ntar)), fill = TRUE)
  }), fill = TRUE)
}

if (!identical(tolower(protein_y), "auto") && !protein_y %in% adt_features) {
  stop("Configured DORM_PROTEIN_Y='", protein_y, "' is not an ADT csv_name.", call. = FALSE)
}

run_one_target <- function(target_group) {
  message("\n=== Target group: ", target_group, " ===")
  partition <- make_partition(meta, part_cfg, target_group, cfg)
  allocation_summary <- rbind(
    data.table(role = "source", group = names(partition$source_idx),
               n = as.integer(lengths(partition$source_idx))),
    data.table(role = "target", group = target_group,
               n = length(partition$target_idx))
  )
  print(allocation_summary)

  task_rows <- sort(unique(c(unlist(partition$source_idx, use.names = FALSE), partition$target_idx)))
  this_protein_y <- protein_y
  if (identical(tolower(this_protein_y), "auto")) {
    message("Screening ADT proteins to choose the easiest RNA-predictable Y for target ", target_group, "...")
    protein_scores <- choose_protein_y(
      rna_dt = rna_dt,
      adt_dt = adt_dt,
      row_idx = task_rows,
      x_cols = rna_features,
      y_cols = adt_features,
      sample_n = cfg$protein_screen_n,
      nfolds = cfg$protein_screen_folds,
      seed = 1
    )
    score_path <- file.path(out_dir, paste0(part_cfg$label, "_target_",
                                            safe_file_part(target_group), "_protein_screen.csv"))
    fwrite(protein_scores, score_path)
    this_protein_y <- protein_scores[!is.na(cv_r2), protein][1]
    message("Selected protein_y = ", this_protein_y,
            " (screen cv_r2 = ", round(protein_scores[protein == this_protein_y, cv_r2], 3), ")")
  } else {
    message("Using manually selected protein_y = ", this_protein_y)
  }

  source_splits <- lapply(partition$source_idx, split_two_equal_parts,
                          fraction = cfg$source_train_fraction)
  target_split <- split_two_equal_parts(partition$target_idx,
                                        fraction = cfg$target_train_fraction)

  Xlist <- lapply(source_splits, function(s) make_x(rna_dt, s$first, rna_features))
  Xtrainlist <- lapply(source_splits, function(s) make_x(rna_dt, s$second, rna_features))
  Ylist <- lapply(source_splits, function(s) make_y(adt_dt, s$first, this_protein_y))
  Ytrainlist <- lapply(source_splits, function(s) make_y(adt_dt, s$second, this_protein_y))
  names(Xlist) <- names(Xtrainlist) <- names(Ylist) <- names(Ytrainlist) <- names(source_splits)

  X0 <- make_x(rna_dt, target_split$first, rna_features)
  X0train <- make_x(rna_dt, target_split$second, rna_features)
  Y0 <- make_y(adt_dt, target_split$first, this_protein_y)
  Y0train <- make_y(adt_dt, target_split$second, this_protein_y)

  nlist <- as.integer(lengths(Ylist))
  q <- length(rna_features)

  dorm_data <- list(
    Xlist = Xlist,
    Xtrainlist = Xtrainlist,
    Ylist = Ylist,
    Ytrainlist = Ytrainlist,
    X0 = X0,
    X0train = X0train,
    Y0 = Y0,
    Y0train = Y0train,
    nlist = nlist,
    q = q,
    protein_y = this_protein_y,
    rna_features = rna_features,
    partition = modifyList(part_cfg, list(target_group = target_group)),
    allocation_summary = allocation_summary
  )

  message("DORM input summary:")
  message("  L sources: ", length(Xlist))
  message("  q RNA features: ", q)
  message("  source nlist: ", paste(nlist, collapse = ", "))
  message("  target n0: ", nrow(X0), " (+ auxiliary ", nrow(X0train), ")")

  file_stub <- paste0(part_cfg$label, "_target_", safe_file_part(target_group),
                      "_Y_", safe_file_part(this_protein_y))
  data_path <- file.path(out_dir, paste0(file_stub, "_dorm_data.rds"))
  saveRDS(dorm_data, data_path)
  message("Saved DORM data object: ", data_path)

  if (!isTRUE(cfg$run_dorm)) return(NULL)

  message("Running DORM...")
  result <- get_DORM_beta(
    Xlist = Xlist,
    Xtrainlist = Xtrainlist,
    Ylist = Ylist,
    Ytrainlist = Ytrainlist,
    X0 = X0,
    X0train = X0train,
    nlist = nlist,
    q = q,
    smax = cfg$smax,
    penalty = cfg$penalty,
    alpha = cfg$alpha,
    rho_pseudo = cfg$rho_pseudo,
    dr_type = cfg$dr_type,
    normalize = cfg$normalize,
    condA = cfg$condA
  )

  mse_for_beta <- function(beta) {
    get_err(X0, Y0, as.vector(beta), q, length(Y0))
  }

  source_mse <- apply(result$DoublyR, 2, mse_for_beta)
  best_source_idx <- which.min(source_mse)
  beta_simple_average <- rowMeans(result$DoublyR)
  beta_rho_average <- as.vector(result$DoublyR %*% as.vector(result$rho))

  method_summary <- data.table(
    method = c("DORM", "Best single source", "Simple average",
               "Rho average", "MI", "RAP"),
    mse = c(
      mse_for_beta(result$beta_star),
      source_mse[best_source_idx],
      mse_for_beta(beta_simple_average),
      mse_for_beta(beta_rho_average),
      mse_for_beta(result$beta_MI),
      mse_for_beta(result$beta_RAP)
    ),
    detail = c(
      "beta_star",
      names(Xlist)[best_source_idx],
      "rowMeans(source-specific Q)",
      "DoublyR %*% rho",
      "beta_MI",
      "beta_RAP"
    )
  )

  target_y_var <- var(Y0)
  method_summary[, r2 := if (target_y_var > 0) 1 - mse / target_y_var else NA_real_]
  method_summary[, `:=`(
    partition = active_partition,
    target_group = target_group,
    protein_y = this_protein_y,
    q = q,
    L = length(Xlist),
    target_n = length(Y0),
    ntar = NA_integer_
  )]
  setcolorder(method_summary, c("partition", "target_group", "protein_y", "q",
                                "L", "target_n", "ntar", "method", "mse", "r2", "detail"))

  if (isTRUE(cfg$run_transfer_benchmarks)) {
    message("Running transfer benchmarks: TransLasso, TransGLM, PTL...")
    transfer_summary <- run_transfer_benchmarks(
      Xlist = Xlist,
      Ylist = Ylist,
      X0_labeled = X0train,
      Y0_labeled = Y0train,
      X_eval = X0,
      Y_eval = Y0,
      q = q,
      ntar_values = cfg$benchmark_ntar
    )
    if (nrow(transfer_summary) > 0) {
      transfer_summary[, `:=`(
        partition = active_partition,
        target_group = target_group,
        protein_y = this_protein_y,
        q = q,
        L = length(Xlist),
        target_n = length(Y0)
      )]
      setcolorder(transfer_summary, c("partition", "target_group", "protein_y", "q",
                                      "L", "target_n", "ntar", "method", "mse", "r2", "detail"))
      benchmark_path <- file.path(out_dir, paste0(file_stub, "_transfer_benchmarks.csv"))
      fwrite(transfer_summary, benchmark_path)
      message("Saved transfer benchmark summary: ", benchmark_path)
      method_summary <- rbindlist(list(method_summary, transfer_summary), fill = TRUE)
    }
  }

  result_path <- file.path(out_dir, paste0(file_stub, "_dorm_result.rds"))
  summary_path <- file.path(out_dir, paste0(file_stub, "_summary.csv"))
  saveRDS(result, result_path)
  fwrite(method_summary, summary_path)

  print(method_summary[order(-r2)])
  message("Saved DORM result: ", result_path)
  message("Saved method R^2 summary: ", summary_path)
  method_summary
}

all_method_summaries <- rbindlist(lapply(part_cfg$target_groups, run_one_target), fill = TRUE)
if (nrow(all_method_summaries) > 0) {
  summary_proteins <- unique(all_method_summaries$protein_y)
  combined_y_label <- if (length(summary_proteins) == 1) summary_proteins else "multiple_proteins"
  combined_path <- file.path(out_dir, paste0(part_cfg$label, "_Y_",
                                             safe_file_part(combined_y_label),
                                             "_all_targets_summary.csv"))
  fwrite(all_method_summaries, combined_path)
  message("\nCombined target summary saved: ", combined_path)
}
