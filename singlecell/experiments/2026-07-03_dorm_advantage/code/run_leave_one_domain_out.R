#!/usr/bin/env Rscript
## Phase 2: natural leave-one-domain-out worst-case robustness.
##   DORM_AXIS=erythroid : domains = erythroid stages (cell_type)
##   DORM_AXIS=donor     : cell type CD14+ Mono; domains = batches (donor covariate shift)
## For each domain held out as target (no target labels), sources = the rest.
## Reports per-target R2 for all methods -> worst-case / average / spread downstream.
##
## Env: DORM_AXIS=erythroid|donor  DORM_ADTS=...  DORM_K=20  DORM_SEEDS=1,2  DORM_MAXN=4000

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

env_csv <- function(n, d) { v <- Sys.getenv(n); if (v == "") d else trimws(strsplit(v, ",")[[1]]) }
AXIS  <- Sys.getenv("DORM_AXIS", "erythroid")
K     <- as.integer(Sys.getenv("DORM_K", "20"))
SEEDS <- as.integer(env_csv("DORM_SEEDS", c("1", "2")))
MAXN  <- as.integer(Sys.getenv("DORM_MAXN", "4000"))

if (AXIS == "erythroid") {
  GROUP_COL <- "cell_type"; SUBSET <- NULL
  DOMAINS <- c("Proerythroblast", "Erythroblast", "Normoblast", "Reticulocyte")  # covariate-shift set
  ADTS <- env_csv("DORM_ADTS", c("CD71", "CD49d"))
} else {                                        # donor axis
  GROUP_COL <- "batch"; SUBSET <- "CD14+ Mono"
  DOMAINS <- c("s3d7","s2d1","s3d6","s3d1","s2d5","s2d4","s1d1","s1d2")  # >=700 mono cells
  ADTS <- env_csv("DORM_ADTS", c("CD11c", "CD33_1"))
}

D <- load_cite_data(csv_dir)
keep <- if (is.null(SUBSET)) rep(TRUE, nrow(D$meta)) else D$meta$cell_type == SUBSET
grp  <- D$meta[[GROUP_COL]]
out_path <- file.path(res_dir, sprintf("lodo_%s_results.csv", AXIS))
w_path   <- file.path(res_dir, sprintf("lodo_%s_weights.csv", AXIS))
if (file.exists(out_path)) file.remove(out_path)
if (file.exists(w_path)) file.remove(w_path)

for (y_col in ADTS) for (seed in SEEDS) {
  set.seed(seed)
  dom_rows <- lapply(DOMAINS, function(g) sample_at_most(which(keep & grp == g), MAXN))
  names(dom_rows) <- DOMAINS
  for (tgt in DOMAINS) {
    src_names <- setdiff(DOMAINS, tgt)
    src_rows <- dom_rows[src_names]
    A_genes <- screen_A(D, unlist(src_rows, use.names = FALSE), y_col, K); q <- K
    ssp <- lapply(src_rows, split_half)
    Xlist      <- lapply(ssp, function(s) make_X(D, s$a, A_genes))
    Xtrainlist <- lapply(ssp, function(s) make_X(D, s$b, A_genes))
    Ylist      <- lapply(ssp, function(s) make_Y(D, s$a, y_col))
    Ytrainlist <- lapply(ssp, function(s) make_Y(D, s$b, y_col))
    names(Xlist) <- names(Xtrainlist) <- names(Ylist) <- names(Ytrainlist) <- src_names
    nlist <- as.integer(lengths(Ylist))
    tsp <- split_half(dom_rows[[tgt]])
    X0 <- make_X(D, tsp$a, A_genes); Y0 <- make_Y(D, tsp$a, y_col)
    X0train <- make_X(D, tsp$b, A_genes); Y0train <- make_Y(D, tsp$b, y_col)

    fit <- tryCatch(dorm_fit(Xlist, Xtrainlist, Ylist, Ytrainlist, X0, X0train, nlist, q),
                    error = function(e) { message("fit err ", y_col, "/", tgt, ": ", conditionMessage(e)); NULL })
    if (is.null(fit)) next
    ev <- eval_all_methods(fit, Xlist, Ylist, X0, Y0, X0train, Y0train, nlist, q)
    ev[, `:=`(axis = AXIS, ADT = y_col, target = tgt, seed = seed, L = length(Xlist), n_target = length(Y0))]
    fwrite(ev, out_path, append = file.exists(out_path))
    fwrite(data.table(axis = AXIS, ADT = y_col, target = tgt, seed = seed,
                      source = src_names, rho = round(as.vector(fit$rho), 4)),
           w_path, append = file.exists(w_path))
    message(sprintf("  %s | %s | target %s | seed %d : DORM(tuned)=%.3f pooled=%.3f best1=%.3f",
                    AXIS, y_col, tgt, seed, ev[method=="DORM (tuned smax)", r2],
                    ev[method=="Pooled elastic-net", r2], ev[method=="Best single source", r2]))
  }
}
message("DONE lodo ", AXIS, " -> ", out_path)
