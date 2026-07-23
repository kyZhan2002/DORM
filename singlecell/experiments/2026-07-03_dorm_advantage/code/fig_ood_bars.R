#!/usr/bin/env Rscript
## Two bar plots, "test on a NEW cell type" (OOD Reticulocyte, MK/E prog), CD71 & CD49d.
## Each compares ONE DORM version vs benchmarks (label-free + label-using at ntar=50):
##   fig6a: DORM q=20,  penalty=FALSE  (fresh 3-seed run, ood_bench_results.csv)
##   fig6b: DORM q=1000, penalty=TRUE   (Ver4 results, singlecell/dorm_results)
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); repo <- dirname(dirname(dirname(exp_dir)))
res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
v4_dir  <- file.path(repo, "singlecell", "dorm_results")

ADTS <- c("CD71", "CD49d"); TARGETS <- c("Reticulocyte", "MK/E prog"); tsafe <- c("Reticulocyte", "MK_E_prog")
METHODS <- c("DORM", "RAP", "MI", "Simple average", "Best single source", "TransLasso", "TransGLM", "PTL")
LABUSE  <- c("TransLasso", "TransGLM", "PTL")
COL <- c(DORM = "#1b7837", RAP = "#4575b4", MI = "#74add1", `Simple average` = "#b2182b",
         `Best single source` = "#999999", TransLasso = "#fee0b6", TransGLM = "#e08214", PTL = "#8073ac")
NTAR <- 50

## ver.a: fresh run (mean over seeds); label-using already at ntar=50
gather_a <- function() {
  R <- fread(file.path(res_dir, "ood_bench_results.csv"))
  R <- R[version == "q=20, penalty=FALSE"]
  R[method == "DORM (tuned smax)", method := "DORM"]
  ## median over seeds: robust to occasional per-seed blow-ups in the unregularised baselines
  R[method %in% METHODS, .(r2 = median(r2, na.rm = TRUE)), by = .(ADT, target, method)]
}
## ver.b: Ver4; label-free from NA-ntar rows, label-using from the ntar=50 rows
gather_b <- function() {
  rows <- list()
  for (i in seq_along(ADTS)) for (j in seq_along(TARGETS)) {
    f <- file.path(v4_dir, sprintf("Ver4_erythroid_subtype_sources_target_%s_Y_%s_summary.csv", tsafe[j], ADTS[i]))
    dt <- fread(f)
    for (m in METHODS) {
      r2 <- if (m %in% LABUSE) dt[method == m & ntar == NTAR, r2] else dt[method == m & is.na(ntar), r2]
      rows[[length(rows)+1]] <- data.table(ADT = ADTS[i], target = TARGETS[j], method = m,
                                           r2 = if (length(r2)) r2[1] else NA_real_)
    }
  }
  rbindlist(rows)
}

make_fig <- function(D, title, outfile) {
  YLO <- -1.5; YHI <- 0.65
  png(file.path(fig_dir, outfile), 1300, 780, res = 110)
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE), heights = c(0.8, 5))
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", METHODS, fill = COL[METHODS], ncol = 4, bty = "n", cex = 1.15,
         title = sprintf("DORM  %s   vs benchmarks    (label-using = TransLasso/TransGLM/PTL at ntar=50)", title))
  par(mar = c(4.5, 5, 3, 1))
  for (a in ADTS) {
    M <- sapply(TARGETS, function(tg) sapply(METHODS, function(m) {
      v <- D[ADT == a & target == tg & method == m, r2]; if (length(v)) v else NA_real_ }))
    Mc <- pmax(pmin(M, YHI), YLO)
    bp <- barplot(Mc, beside = TRUE, col = COL[METHODS], border = NA, ylim = c(YLO, YHI),
                  ylab = "target R^2 (clamped)", main = paste0(a, ":  DORM (", title, ") vs benchmarks"),
                  names.arg = TARGETS, cex.names = 1.4, cex.axis = 1.25, cex.lab = 1.4, cex.main = 1.4)
    abline(h = 0, col = "grey40")
    off <- which(!is.na(M) & M < YLO, arr.ind = TRUE)
    if (nrow(off)) for (r in seq_len(nrow(off)))
      text(bp[off[r,1], off[r,2]], YLO + 0.10, sprintf("%.1f", M[off[r,1], off[r,2]]),
           srt = 90, adj = 0, cex = 0.85, col = "red")
  }
  dev.off()
  cat("\n== DORM", title, "==\n"); print(dcast(D, method ~ ADT + target, value.var = "r2")[match(METHODS, method)], row.names = FALSE)
  cat("Saved", outfile, "\n")
}

make_fig(gather_b(), "q=1000, penalty=TRUE", "fig6b_ood_dorm_q1000_penT.png")  # Ver4 available now
if (file.exists(file.path(res_dir, "ood_bench_results.csv")))
  make_fig(gather_a(), "q=20, penalty=FALSE", "fig6a_ood_dorm_q20_penF.png")   # after fresh run
