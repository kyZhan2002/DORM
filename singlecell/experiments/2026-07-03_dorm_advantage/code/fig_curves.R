#!/usr/bin/env Rscript
## Clamped, legible curve figures with label-free AND label-using methods.
##   fig1_mixture_curves.png : R^2 vs mixture (ordered by distance from pooled centroid)
##   fig3_lodo_<axis>.png    : R^2 vs held-out target domain
## Usage: Rscript fig_curves.R [fig1|fig3|all]
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
which_fig <- { a <- commandArgs(TRUE); if (length(a)) a[1] else "all" }

## method -> (color, lty). Solid = label-free, dashed = label-using, dotted = unstable naive.
STYLE <- data.table(
  method = c("DORM (tuned smax)","RAP","Pooled elastic-net","TransGLM","TransLasso","PTL"),
  col    = c("#1b7837","#4575b4","#000000","#e08214","#d73027","#762a83"),
  lty    = c(1,1,1,2,2,2))
## fig1 also shows the naive aggregations MI and Simple average (often off-scale),
## plus the target-only ORACLE (elastic-net trained on target labels = upper bound).
STYLE1 <- rbind(STYLE, data.table(
  method = c("MI","Simple average","Target-oracle EN"),
  col    = c("#666666","#b2182b","#111111"), lty = c(3,3,4)))

pch_of <- function(lty) c(`1`=19,`2`=17,`3`=4,`4`=18)[as.character(lty)]

draw_panel <- function(d, xcol, xorder, xlabels, ylim, main, ylab, style) {
  xs <- seq_along(xorder)
  plot(NA, xlim = c(1, length(xorder)), ylim = ylim, xaxt = "n", xlab = "", ylab = ylab, main = main,
       cex.axis = 1.25, cex.lab = 1.4, cex.main = 1.35)
  axis(1, at = xs, labels = xlabels, las = 2, cex.axis = 1.05)
  abline(h = 0, col = "grey70", lty = 3)
  for (i in seq_len(nrow(style))) {
    m <- style$method[i]
    dm <- d[method == m]; if (!nrow(dm)) next
    dm <- dm[match(xorder, get(xcol))]
    y <- pmax(pmin(dm$r2, ylim[2]), ylim[1])            # clamp into view
    lwd <- if (style$lty[i] == 4) 3 else 2.2            # bolder oracle reference
    lines(xs, y, col = style$col[i], lty = style$lty[i], lwd = lwd)
    points(xs, y, col = style$col[i], pch = pch_of(style$lty[i]), cex = 1.05)
  }
}

## ---------- fig1: mixture curves ----------
if (which_fig %in% c("fig1", "all") && file.exists(file.path(res_dir, "mixture_results.csv"))) {
  m <- fread(file.path(res_dir, "mixture_results.csv"))
  nl <- c(756, 2019, 717); cen <- nl / sum(nl)
  m[, dist := sqrt((wP-cen[1])^2 + (wE-cen[2])^2 + (wN-cen[3])^2)]
  agg <- m[method %in% STYLE1$method, .(r2 = mean(r2), dist = mean(dist),
           wP = wP[1], wE = wE[1], wN = wN[1]), by = .(ADT, method, mixture)]
  ylim <- c(-0.5, 0.7)
  ## min R^2 of the naive methods (for the off-scale note)
  offnote <- m[method %in% c("MI","Simple average"), .(mn = min(r2)), by = method]
  png(file.path(fig_dir, "fig1_mixture_curves.png"), 1400, 860, res = 115)
  layout(matrix(c(1, 1, 2, 3, 4, 4), nrow = 3, byrow = TRUE), heights = c(1.0, 5, 1.0))
  ## row 1: legend strip
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", STYLE1$method, col = STYLE1$col, lty = STYLE1$lty, lwd = 2.4,
         pch = pch_of(STYLE1$lty), cex = 1.1, ncol = 5, bty = "n",
         title = "solid = label-free    dashed = label-using (ntar=50)    dotted = naive aggregation    dot-dash = target-only ORACLE (uses labels)")
  ## row 2: the two panels
  par(mar = c(6.5, 5, 2.6, 1))
  for (a in c("CD71", "CD49d")) {
    d <- agg[ADT == a]; meta <- unique(d[, .(mixture, dist, wP, wE, wN)])[order(dist)]
    ord <- meta$mixture
    xlab <- sprintf("%d/%d/%d", round(100*meta$wP), round(100*meta$wE), round(100*meta$wN))
    draw_panel(d, "mixture", ord, xlab, ylim, paste0(a, ": R^2 vs target composition"),
               "target R^2 (clamped)", STYLE1)
  }
  ## row 3: notes
  par(mar = c(0, 0, 0, 0)); plot.new()
  text(0.5, 0.68, "target composition  %Proerythroblast / %Erythroblast / %Normoblast   (LEFT near pooled mix -> RIGHT far/adversarial)", cex = 1.05)
  text(0.5, 0.22, sprintf("Note: MI and Simple average fall below the axis on several mixtures (min R^2: MI %.1f, Simple avg %.1f)",
       offnote[method=="MI", mn], offnote[method=="Simple average", mn]), cex = 1.0, col = "#b2182b")
  dev.off(); cat("Saved fig1_mixture_curves.png\n")
}

## ---------- fig3: leave-one-domain-out curves ----------
if (which_fig %in% c("fig3", "all")) {
  for (ax in c("erythroid", "donor")) {
    f <- file.path(res_dir, sprintf("lodo_%s_results.csv", ax)); if (!file.exists(f)) next
    d0 <- fread(f)
    ylim <- if (ax == "erythroid") c(-1.2, 0.6) else c(-1.0, 0.4)
    adts <- unique(d0$ADT)
    png(file.path(fig_dir, sprintf("fig3_lodo_%s.png", ax)), 620 * length(adts), 560, res = 115)
    par(mfrow = c(1, length(adts)), mar = c(9, 4, 3, 1))
    for (a in adts) {
      d <- d0[ADT == a & method %in% STYLE$method, .(r2 = mean(r2)), by = .(method, target)]
      ord <- unique(d0[ADT == a, target])
      draw_panel(d, "target", ord, ord, ylim, sprintf("LODO %s (%s)", ax, a), "target R^2 (clamped)", STYLE)
      if (a == adts[1]) legend("bottomleft", STYLE$method, col = STYLE$col, lty = STYLE$lty, lwd = 2,
                               pch = ifelse(STYLE$lty == 1, 19, 17), cex = 0.7, bty = "n",
                               title = "solid=label-free, dashed=label-using(ntar=50)")
    }
    dev.off(); cat(sprintf("Saved fig3_lodo_%s.png\n", ax))
  }
}
