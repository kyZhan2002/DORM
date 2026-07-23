#!/usr/bin/env Rscript
## "Test on a NEW cell type" benchmark run, for TWO DORM configurations, each vs the
## label-free + label-using benchmark suite. OOD targets = Reticulocyte, MK/E prog;
## sources = the 3 nucleated erythroid stages; ADTs = CD71, CD49d.
##   ver.a: q=20,   penalty=FALSE (low-dim A + high-dim W controls)
##   ver.b: q=1000, penalty=TRUE  (all genes as A, elastic-net regularised)
## All core fixes active (unpenalized intercept, ridge, PTL A-slice); label-using at ntar=50.
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

ADTS <- c("CD71", "CD49d"); TARGETS <- c("Reticulocyte", "MK/E prog")
SOURCES <- c("Proerythroblast", "Erythroblast", "Normoblast")
SEEDS <- as.integer(strsplit(Sys.getenv("DORM_SEEDS", "1"), ",")[[1]])
MAXN <- 4000
D <- load_cite_data(csv_dir)
out <- file.path(res_dir, "ood_bench_results.csv"); if (file.exists(out)) file.remove(out)

for (y in ADTS) for (tg in TARGETS) for (seed in SEEDS) {
  set.seed(seed)
  src_rows <- lapply(SOURCES, function(g) sample_at_most(which(D$meta$cell_type == g), MAXN)); names(src_rows) <- SOURCES
  tgt_rows <- sample_at_most(which(D$meta$cell_type == tg), MAXN)
  ssp <- lapply(src_rows, split_half); tsp <- split_half(tgt_rows)
  Ylist <- lapply(ssp, function(s) make_Y(D, s$a, y)); Ytrainlist <- lapply(ssp, function(s) make_Y(D, s$b, y))
  Y0 <- make_Y(D, tsp$a, y); Y0train <- make_Y(D, tsp$b, y)

  run_ver <- function(A_genes, q, penalty, tag) {
    Xl  <- lapply(ssp, function(s) make_X(D, s$a, A_genes)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A_genes))
    names(Xl) <- names(Xt) <- SOURCES
    X0 <- make_X(D, tsp$a, A_genes); X0t <- make_X(D, tsp$b, A_genes)
    nlist <- as.integer(lengths(Ylist))
    fit <- dorm_fit(Xl, Xt, Ylist, Ytrainlist, X0, X0t, nlist, q, penalty = penalty)
    ev <- eval_all_methods(fit, Xl, Ylist, X0, Y0, X0t, Y0train, nlist, q)
    ev[, `:=`(version = tag, ADT = y, target = tg, seed = seed)]
    fwrite(ev, out, append = file.exists(out))
    message(sprintf("  %s | %s | %s seed%d : DORM(tuned)=%.3f", tag, y, tg, seed,
                    ev[method == "DORM (tuned smax)", r2]))
  }
  message(sprintf("== %s | %s | seed %d ==", y, tg, seed))
  run_ver(screen_A(D, unlist(src_rows, use.names = FALSE), y, 20), 20, FALSE, "q=20, penalty=FALSE")
  if (tolower(Sys.getenv("DORM_RUN_VERB", "false")) %in% c("1","true","t","yes","y"))
    run_ver(D$rna_features, length(D$rna_features), TRUE, "q=1000, penalty=TRUE")
}
message("DONE ood_bench -> ", out)
