#!/usr/bin/env Rscript
## Aggregate the advantage analysis: worst-case / average / spread per method,
## split by tier, over seeds. Plus figures. Base R only (no extra deps).
suppressPackageStartupMessages(library(data.table))

this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
code_dir <- dirname(this_file); exp_dir <- dirname(code_dir)
res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## worst-case (min over shifts) / mean / spread per method, averaged over seeds
summarize_shifts <- function(dt) {
  ps <- dt[is.finite(r2), .(wc = min(r2), mean_r2 = mean(r2), sd = sd(r2)),
           by = .(ADT, method, tier, seed)]
  ps[, .(worstcase = round(mean(wc), 3), mean_r2 = round(mean(mean_r2), 3),
         spread = round(mean(sd), 3)), by = .(ADT, method, tier)]
}
show_tier <- function(s, title) {
  cat("\n==== ", title, " ====\n")
  for (a in unique(s$ADT)) {
    cat("\n-- ADT:", a, " (label-free tier, ranked by worst-case) --\n")
    print(s[ADT == a & tier == "label-free"][order(-worstcase),
          .(method, worstcase, mean_r2, spread)], row.names = FALSE)
    lu <- s[ADT == a & tier == "label-using"][order(-worstcase), .(method, worstcase, mean_r2)]
    if (nrow(lu)) { cat("   [label-using context]\n"); print(lu, row.names = FALSE) }
  }
}

## ---------- Phase 1: mixture ----------
mix_f <- file.path(res_dir, "mixture_results.csv")
if (file.exists(mix_f)) {
  mix <- fread(mix_f)
  s <- summarize_shifts(mix); fwrite(s, file.path(res_dir, "summary_mixture.csv"))
  show_tier(s, "PHASE 1: adversarial mixture stress test (worst-case over mixtures)")

  ## Figure 1: R2 vs mixture extremeness, DORM(tuned) vs Pooled vs Best-single
  mix[, extremeness := pmax(wP, wE, wN)]
  key <- c("DORM (tuned smax)", "Pooled elastic-net", "Best single source", "Simple average")
  agg <- mix[method %in% key, .(r2 = mean(r2)), by = .(ADT, method, mixture, extremeness)]
  png(file.path(fig_dir, "fig1_mixture_curves.png"), 1100, 500, res = 110)
  par(mfrow = c(1, length(unique(agg$ADT))), mar = c(7, 4, 3, 1))
  cols <- c("DORM (tuned smax)"="#1b7837","Pooled elastic-net"="#d73027",
            "Best single source"="#4575b4","Simple average"="#999999")
  for (a in unique(agg$ADT)) {
    d <- agg[ADT == a][order(extremeness)]
    ord <- unique(d[order(extremeness), mixture])
    d[, xi := match(mixture, ord)]
    plot(NA, xlim = c(1, length(ord)), ylim = range(d$r2, na.rm = TRUE),
         xaxt = "n", xlab = "", ylab = "target R^2", main = paste0(a, " (mixture shift)"))
    axis(1, at = seq_along(ord), labels = ord, las = 2, cex.axis = 0.7)
    for (m in key) { dm <- d[method == m][order(xi)]
      lines(dm$xi, dm$r2, col = cols[m], lwd = 2); points(dm$xi, dm$r2, col = cols[m], pch = 19) }
    if (a == unique(agg$ADT)[1]) legend("bottomleft", names(cols), col = cols, lwd = 2, cex = 0.7, bty = "n")
  }
  dev.off()

  ## Figure 2: worst-case bar chart (label-free), first ADT
  a1 <- unique(s$ADT)[1]
  sf <- s[ADT == a1 & tier == "label-free"][order(worstcase)]
  png(file.path(fig_dir, "fig2_worstcase_bar.png"), 800, 500, res = 110)
  par(mar = c(4, 12, 3, 1))
  barplot(sf$worstcase, names.arg = sf$method, horiz = TRUE, las = 1,
          col = ifelse(grepl("^DORM", sf$method), "#1b7837", "grey70"),
          xlab = "worst-case R^2 over mixtures", main = paste0("Worst-case robustness (", a1, ")"), cex.names = 0.7)
  abline(v = 0, lty = 2); dev.off()
  cat("\nFigures: fig1_mixture_curves.png, fig2_worstcase_bar.png\n")
}

## ---------- Phase 2: leave-one-domain-out ----------
for (ax in c("erythroid", "donor")) {
  f <- file.path(res_dir, sprintf("lodo_%s_results.csv", ax))
  if (!file.exists(f)) next
  d <- fread(f); s <- summarize_shifts(d); fwrite(s, file.path(res_dir, sprintf("summary_lodo_%s.csv", ax)))
  show_tier(s, sprintf("PHASE 2: leave-one-domain-out [%s] (worst-case over held-out domains)", ax))

  ## Figure: per-target R2 across domains, DORM(tuned) vs Pooled vs Best-single (first ADT)
  a1 <- unique(d$ADT)[1]
  key <- c("DORM (tuned smax)", "Pooled elastic-net", "Best single source")
  agg <- d[ADT == a1 & method %in% key, .(r2 = mean(r2)), by = .(method, target)]
  png(file.path(fig_dir, sprintf("fig3_lodo_%s.png", ax)), 900, 500, res = 110)
  par(mar = c(8, 4, 3, 1))
  tg <- unique(agg$target); cols <- c("#1b7837","#d73027","#4575b4"); names(cols) <- key
  plot(NA, xlim = c(1, length(tg)), ylim = range(agg$r2, na.rm = TRUE), xaxt = "n",
       xlab = "", ylab = "target R^2", main = sprintf("LODO %s: per-target R^2 (%s)", ax, a1))
  axis(1, at = seq_along(tg), labels = tg, las = 2, cex.axis = 0.7)
  for (m in key) { dm <- agg[method == m]; dm[, xi := match(target, tg)]
    lines(dm$xi, dm$r2, col = cols[m], lwd = 2); points(dm$xi, dm$r2, col = cols[m], pch = 19) }
  legend("bottomleft", key, col = cols, lwd = 2, cex = 0.7, bty = "n"); dev.off()
  cat(sprintf("Figure: fig3_lodo_%s.png\n", ax))
}
cat("\nDONE. Tables in results/summary_*.csv ; figures in figures/\n")
