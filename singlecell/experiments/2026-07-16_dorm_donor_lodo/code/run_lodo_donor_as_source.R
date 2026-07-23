#!/usr/bin/env Rscript
## Leave-one-donor-out with DONORS AS SOURCES (annotation-free), GEX->ADT over ALL cell types.
## For each held-out donor d: L sources = the OTHER donors (each = its cells, a natural cell-type
## mixture); target = donor d's cells (also a mixture). Predict ADT from RNA, no target labels.
## The batch code is site+donor (e.g. s3d7 = site 3, donor 7); on leave-one-donor-out the target's
## same-site siblings are among the sources, so DORM's density ratio can up-weight them.
## DORM: q=20, penalty=TRUE. Label-using benchmarks at ntar=20 (set in lib). Env:
##   DORM_ADTS=CD71,CD72_1,CD19_1,CD11c  DORM_SEEDS=1,2  DORM_DONORS=<subset>  DORM_HVG=250
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
source(file.path(code_dir, "lib_dorm_common.R"))
repo <- find_repo(code_dir); load_dorm_src(repo)
res_dir <- file.path(exp_dir, "results"); csv_dir <- file.path(repo, "singlecell", "processed", "csv_for_R")

env_csv <- function(n, d) { v <- Sys.getenv(n); if (v == "") d else trimws(strsplit(v, ",")[[1]]) }
ADTS   <- env_csv("DORM_ADTS", c("CD71", "CD72_1", "CD19_1", "CD11c"))
SEEDS  <- as.integer(env_csv("DORM_SEEDS", c("1", "2")))
HVG    <- as.integer(Sys.getenv("DORM_HVG", "250"))
MAXN_SRC <- as.integer(Sys.getenv("DORM_MAXN_SRC", "500"))   # cells per source donor
MAXN_TGT <- as.integer(Sys.getenv("DORM_MAXN_TGT", "1500"))  # cells for the target donor
K <- 20; site_of <- function(b) substr(b, 1, 2)

D <- load_cite_data(csv_dir)
D$rna_features <- head(D$rna_features, HVG)            # subset HVG for speed (W controls)
ALL_DONORS <- sort(unique(D$meta$batch))
DONORS <- env_csv("DORM_DONORS", ALL_DONORS)
message("Donors (", length(DONORS), "): ", paste(DONORS, collapse = ", "))

out  <- file.path(res_dir, "lodo_donor_as_source_results.csv")
wout <- file.path(res_dir, "lodo_donor_as_source_rho.csv")
if (Sys.getenv("DORM_APPEND", "") == "" ) { for (f in c(out, wout)) if (file.exists(f)) file.remove(f) }

for (dtar in DONORS) for (y in ADTS) for (seed in SEEDS) {
  set.seed(seed)
  src_donors <- setdiff(ALL_DONORS, dtar)
  src_rows <- lapply(src_donors, function(b) sample_at_most(which(D$meta$batch == b), MAXN_SRC))
  names(src_rows) <- src_donors
  keep <- lengths(src_rows) >= 60
  src_rows <- src_rows[keep]; src_donors <- names(src_rows)
  tgt_rows <- sample_at_most(which(D$meta$batch == dtar), MAXN_TGT)
  if (length(tgt_rows) < 150) { message("  skip ", dtar, ": too few target cells"); next }

  A <- screen_A(D, unlist(src_rows, use.names = FALSE), y, K); q <- K
  ssp <- lapply(src_rows, split_half); tsp <- split_half(tgt_rows)
  Xl <- lapply(ssp, function(s) make_X(D, s$a, A)); Xt <- lapply(ssp, function(s) make_X(D, s$b, A))
  Yl <- lapply(ssp, function(s) make_Y(D, s$a, y)); Yt <- lapply(ssp, function(s) make_Y(D, s$b, y))
  names(Xl) <- names(Xt) <- names(Yl) <- names(Yt) <- src_donors; nlist <- as.integer(lengths(Yl))
  X0 <- make_X(D, tsp$a, A); Y0 <- make_Y(D, tsp$a, y); X0t <- make_X(D, tsp$b, A); Y0t <- make_Y(D, tsp$b, y)

  fit <- tryCatch(dorm_fit(Xl, Xt, Yl, Yt, X0, X0t, nlist, q, penalty = TRUE),
                  error = function(e) { message("  fit err ", dtar, "/", y, ": ", conditionMessage(e)); NULL })
  if (is.null(fit)) next
  ev <- eval_all_methods(fit, Xl, Yl, X0, Y0, X0t, Y0t, nlist, q)
  ev[, `:=`(donor = dtar, site = site_of(dtar), ADT = y, seed = seed, L = length(Xl), n_target = length(Y0))]
  fwrite(ev, out, append = file.exists(out))
  ## rho per source donor + same-site flag (mechanism check)
  fwrite(data.table(target_donor = dtar, ADT = y, seed = seed, source_donor = src_donors,
                    same_site = site_of(src_donors) == site_of(dtar), rho = round(as.vector(fit$rho), 4)),
         wout, append = file.exists(wout))
  message(sprintf("  %s (%s) | %s | seed%d : DORM=%.3f Pooled=%.3f TransGLM@20=%.3f oracle=%.3f  rho(same-site)=%.2f",
                  dtar, site_of(dtar), y, seed,
                  ev[method=="DORM (tuned smax)", r2], ev[method=="Pooled elastic-net", r2],
                  ev[method=="TransGLM", r2], ev[method=="Target-oracle EN", r2],
                  sum(fit$rho[site_of(src_donors) == site_of(dtar)])))
}
message("DONE lodo_donor_as_source -> ", out)
