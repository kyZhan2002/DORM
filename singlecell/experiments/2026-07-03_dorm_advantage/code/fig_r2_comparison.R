#!/usr/bin/env Rscript
## Clear R^2 comparison across methods on the erythroid mixture stress test.
## Fixes the unreadable y-range of the earlier figures by CLAMPING to a sensible
## window and annotating any method that falls off-scale with its true value.
## Methods: DORM(tuned), RAP, MI, Simple average, Pooled EN, Best single source (label-free)
##          + TransLasso, TransGLM, PTL, Target-oracle EN (label-using, best ntar<=200).
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")

m <- fread(file.path(res_dir, "mixture_results.csv"))
methods <- c("DORM (tuned smax)", "RAP", "MI", "Simple average", "Pooled elastic-net",
             "Best single source", "TransLasso", "TransGLM", "PTL", "Target-oracle EN")
lab_using <- c("TransLasso", "TransGLM", "PTL", "Target-oracle EN")

s <- m[method %in% methods & is.finite(r2),
       .(mean = mean(r2), sd = sd(r2), worst = min(r2)), by = .(ADT, method)]
YLO <- -0.6; YHI <- 0.75                          # clamp window

png(file.path(fig_dir, "fig5_r2_comparison.png"), 1250, 620, res = 115)
par(mfrow = c(1, 2), mar = c(9.5, 4, 3, 1))
for (a in c("CD71", "CD49d")) {
  d <- s[ADT == a]; d <- d[match(methods, method)]      # fixed method order
  col <- ifelse(d$method == "DORM (tuned smax)", "#1b7837",
         ifelse(d$method %in% lab_using, "#e08214", "#4575b4"))
  h <- pmax(pmin(d$mean, YHI), YLO)                       # clamped bar heights
  bp <- barplot(h, names.arg = rep("", nrow(d)), col = col, border = NA,
                ylim = c(YLO, YHI), ylab = "mean R^2 over mixtures (+/- sd)",
                main = paste0(a, " -- method comparison"))
  abline(h = 0, col = "grey40")
  ## error bars (clamped), and worst-case tick
  segs <- pmax(pmin(d$mean + d$sd, YHI), YLO); segl <- pmax(pmin(d$mean - d$sd, YHI), YLO)
  arrows(bp, segl, bp, segs, angle = 90, code = 3, length = 0.03, col = "grey30")
  ## annotate off-scale means with true value
  off <- which(d$mean < YLO)
  if (length(off)) text(bp[off], YLO + 0.03, sprintf("%.2f", d$mean[off]), srt = 90, adj = 0, cex = 0.7, col = "red")
  axis(1, at = bp, labels = d$method, las = 2, cex.axis = 0.75)
  if (a == "CD71") legend("topright", c("DORM (0 labels)", "other label-free", "label-using (ntar=50)"),
                          fill = c("#1b7837", "#4575b4", "#e08214"), bty = "n", cex = 0.75)
}
dev.off()
cat("Saved figures/fig5_r2_comparison.png\n")
cat("\nMean R^2 (clamp window [", YLO, ",", YHI, "]):\n", sep = "")
print(dcast(s, method ~ ADT, value.var = "mean")[match(methods, method)], row.names = FALSE)
