#!/usr/bin/env Rscript
## Analyse leave-one-DONOR-out (erythroid) results: per-donor R^2 curves + worst-case/mean
## summary, DORM (q=20, penT) vs benchmarks (label-free + label-using ntar=50).
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
R <- fread(file.path(res_dir, "lodo_donor_ery_results.csv"))
R[method == "DORM (tuned smax)", method := "DORM"]

## median over seeds (robust to noisy folds)
M <- R[, .(r2 = median(r2, na.rm = TRUE), wEry = mean(wEry)), by = .(ADT, donor, method)]

## ---- summary table: worst-case (min over donors) + mean, per method ----
summ <- M[is.finite(r2), .(worstcase = round(min(r2), 3), mean = round(mean(r2), 3)), by = .(ADT, method)]
for (a in c("CD71", "CD49d")) {
  cat("\n==== leave-one-donor-out [", a, "] : worst-case & mean R^2 over donors ====\n")
  print(summ[ADT == a][order(-worstcase), .(method, worstcase, mean)], row.names = FALSE)
}
fwrite(summ, file.path(res_dir, "summary_lodo_donor_ery.csv"))

## ---- per-donor curves ----
STYLE <- data.table(
  method = c("DORM","RAP","Pooled elastic-net","TransGLM","Target-oracle EN"),
  col    = c("#1b7837","#4575b4","#000000","#e08214","#111111"), lty = c(1,1,1,2,4))
pch_of <- function(lty) c(`1`=19,`2`=17,`4`=18)[as.character(lty)]
YLIM <- c(-1.0, 0.6)
png(file.path(fig_dir, "fig7_lodo_donor_ery.png"), 1400, 620, res = 115)
par(mfrow = c(1, 2), mar = c(6.5, 5, 3.2, 1), oma = c(0, 0, 2.5, 0))
for (a in c("CD71", "CD49d")) {
  d <- M[ADT == a]
  ord <- unique(d[order(wEry), .(donor, wEry)]); labs <- sprintf("%s\n(%.0f%%Ery)", ord$donor, 100*ord$wEry)
  xs <- seq_len(nrow(ord))
  plot(NA, xlim = c(1, nrow(ord)), ylim = YLIM, xaxt = "n", xlab = "", ylab = "target R^2 (clamped)",
       main = sprintf("%s: leave-one-donor-out", a), cex.axis = 1.2, cex.lab = 1.35, cex.main = 1.4)
  axis(1, at = xs, labels = labs, las = 2, cex.axis = 0.95); abline(h = 0, col = "grey70", lty = 3)
  for (i in seq_len(nrow(STYLE))) {
    dm <- d[method == STYLE$method[i]][match(ord$donor, donor)]
    y <- pmax(pmin(dm$r2, YLIM[2]), YLIM[1])
    lines(xs, y, col = STYLE$col[i], lty = STYLE$lty[i], lwd = if (STYLE$lty[i]==4) 3 else 2.4)
    points(xs, y, col = STYLE$col[i], pch = pch_of(STYLE$lty[i]), cex = 1.1)
  }
  if (a == "CD71") legend("bottomright", STYLE$method, col = STYLE$col, lty = STYLE$lty, lwd = 2.4,
                          pch = pch_of(STYLE$lty), cex = 1.0, bty = "n")
}
mtext("Leave-one-donor-out on erythroid cells: 3 cell-type sources (11 donors) -> held-out donor mixture   (donors ordered by %Erythroblast)",
      side = 3, outer = TRUE, cex = 1.0)
dev.off(); cat("\nSaved fig7_lodo_donor_ery.png\n")
