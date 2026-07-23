#!/usr/bin/env Rscript
## Paper-style (ggplot2 + theme_bw) figures for the donor-LODO study, in TWO sets:
##   *_clean  : DORM + benchmarks (no intrinsic)          -> figs 9a/10a
##   *_full   : adds DORM_intrinsic (may be unstable)     -> figs 9b/10b
## Palette matches simu/graph_main.R.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
R <- fread(file.path(res_dir, "lodo_donor_intrinsic_results.csv"))
R[method == "DORM (tuned smax)", method := "DORM"]

BENCH <- c("DORM", "Pooled elastic-net", "Simple average", "TransLasso", "TransGLM", "PTL", "Target-oracle EN")
PAL <- c("DORM"="#FF0000", "DORM_intrinsic"="#8B0000", "Pooled elastic-net"="lightblue3",
         "Simple average"="blue2", "TransLasso"="orange", "TransGLM"="purple",
         "PTL"="green3", "Target-oracle EN"="grey50")
LTY <- c("DORM"="solid", "DORM_intrinsic"="longdash", "Pooled elastic-net"="dashed",
         "Simple average"="dashed", "TransLasso"="dashed", "TransGLM"="dashed",
         "PTL"="dashed", "Target-oracle EN"="dotted")
ADT_ORDER <- intersect(c("CD72_1","CD19_1","CD11c","CD3"), unique(R$ADT))
base_theme <- theme_bw() + theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.title = element_text(size = 13),
  legend.text = element_text(size = 12), strip.text = element_text(size = 13), plot.title = element_text(size = 14))

make_set <- function(methods, tag) {
  M <- R[method %in% methods, .(r2 = median(r2, na.rm = TRUE)), by = .(ADT, donor, site, method)]
  M[, r2 := pmax(r2, -5)]                              # clamp extreme blow-ups (view is [-0.6,0.85] anyway)
  M[, method := factor(method, levels = methods)]; M[, ADT := factor(ADT, levels = ADT_ORDER)]
  don_ord <- unique(M[order(site, donor), donor]); M[, donor := factor(donor, levels = don_ord)]

  ## per-donor, faceted by ADT
  p9 <- ggplot(M, aes(donor, r2, color = method, linetype = method, group = method)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
    facet_wrap(~ ADT, ncol = 2) + scale_color_manual(values = PAL) + scale_linetype_manual(values = LTY) +
    coord_cartesian(ylim = c(-0.6, 0.85)) +
    labs(x = "held-out donor (grouped by site)", y = expression("target " * R^2 * " (clamped)"),
         color = "Method", linetype = "Method",
         title = sprintf("Leave-one-donor-out (donors as sources) [%s]", tag)) + base_theme
  ggsave(file.path(fig_dir, sprintf("fig9_%s_per_donor.png", tag)), p9, width = 12, height = 8, dpi = 200)

  ## summary: median R^2 over all donor x ADT (single number per method)
  FLOOR <- -0.5                                          # truncate display at -0.5; off-scale value kept in label
  S <- M[is.finite(r2), .(median_r2 = median(r2)), by = method][order(-median_r2)]
  S[, disp := pmax(median_r2, FLOOR)]                    # clamp bar height to the floor
  S[, method := factor(method, levels = S$method)]
  p10 <- ggplot(S, aes(method, disp, fill = method)) + geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", median_r2)), vjust = ifelse(S$disp >= 0, -0.3, 1.2), size = 4) +
    scale_fill_manual(values = PAL, guide = "none") + coord_cartesian(ylim = c(FLOOR, max(S$median_r2)*1.15)) +
    labs(x = NULL, y = expression("median " * R^2 * " (all donor " %*% " ADT, truncated at -0.5)"),
         title = sprintf("Summary median R^2 [%s]", tag)) +
    theme_bw() + theme(axis.title = element_text(size = 14), axis.text.x = element_text(size = 11, angle = 30, hjust = 1),
                       axis.text.y = element_text(size = 12), plot.title = element_text(size = 13))
  ggsave(file.path(fig_dir, sprintf("fig10_%s_summary.png", tag)), p10, width = 9, height = 6, dpi = 200)
  cat(sprintf("\n[%s] overall median R^2 per method:\n", tag)); print(S, row.names = FALSE)
}

make_set(BENCH, "clean")                                   # figs 9_clean / 10_clean (no intrinsic)
if ("DORM_intrinsic" %in% R$method) make_set(c("DORM", "DORM_intrinsic", setdiff(BENCH, "DORM")), "full")
cat("\nSaved fig9_{clean,full}_per_donor.png, fig10_{clean,full}_summary.png\n")
