#!/usr/bin/env Rscript
## Phase C figures: does DORM-intrinsic (Harmony signed-affine prior) beat standard DORM as the target
## leaves the source-donor hull? Palette matches fig_intrinsic.R / simu/graph_main.R.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
R <- fread(file.path(res_dir, "extrapolation_shift_results.csv"))
R[method == "DORM (tuned smax)", method := "DORM"]

PAL <- c("DORM"="#FF0000", "DORM_intrinsic"="#8B0000", "Pooled elastic-net"="lightblue3",
         "Simple average"="blue2", "TransLasso"="orange", "TransGLM"="purple",
         "PTL"="green3", "Target-oracle EN"="grey50")
LTY <- c("DORM"="solid", "DORM_intrinsic"="solid", "Pooled elastic-net"="dashed",
         "Simple average"="dashed", "TransLasso"="dashed", "TransGLM"="dashed",
         "PTL"="dashed", "Target-oracle EN"="dotted")
METHODS <- names(PAL)
M <- R[method %in% METHODS]
M[, r2 := pmax(r2, -1)]                                   # clamp extreme negatives for display
M[, method := factor(method, levels = METHODS)]
M[, panel := paste0(donor, " | ", ADT, "  (state ", state_c, ")")]

## Panel A: R^2 vs alpha (composition shift), one line per method, facet donor x ADT
pA <- ggplot(M, aes(alpha, r2, color = method, linetype = method, group = method)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  facet_wrap(~ panel) + scale_color_manual(values = PAL) + scale_linetype_manual(values = LTY) +
  coord_cartesian(ylim = c(-1, 0.9)) +
  labs(x = expression("composition shift " * alpha * "  (0 = natural donor, 1 = pure state)"),
       y = expression("target " * R^2 * " (clamped at -1)"), color = "Method", linetype = "Method",
       title = "Extrapolation: prediction as the target leaves the source-donor hull") +
  theme_bw() + theme(strip.text = element_text(size = 11), legend.position = "right")
ggsave(file.path(fig_dir, "fig_extrap_r2.png"), pA, width = 12, height = 7, dpi = 200)

## Panel B: the decisive plot -- gap (DORM_intrinsic - DORM) vs how far outside the hull (excess_c)
W <- dcast(R[method %in% c("DORM", "DORM_intrinsic")], donor + ADT + alpha + excess_c ~ method, value.var = "r2")
if (all(c("DORM", "DORM_intrinsic") %in% names(W))) {
  W[, gap := DORM_intrinsic - DORM]
  pB <- ggplot(W, aes(excess_c, gap, color = paste0(donor, " | ", ADT))) +
    geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    labs(x = "excess state-c fraction vs sources  (>0 = target outside the donor hull)",
         y = expression(R^2 * "(DORM-intrinsic) - " * R^2 * "(DORM)"), color = "target | ADT",
         title = "Does the intrinsic prior help more the further the target extrapolates?") +
    theme_bw()
  ggsave(file.path(fig_dir, "fig_extrap_gap.png"), pB, width = 9, height = 5.5, dpi = 200)
}

## console summary
S <- R[method %in% c("DORM", "DORM_intrinsic", "Pooled elastic-net", "Target-oracle EN"),
       .(r2 = round(median(r2), 3)), by = .(method, interior = excess_c <= 0)]
cat("median R^2 by regime (interior = target inside hull):\n"); print(dcast(S, method ~ interior, value.var = "r2"))
cat("\nSaved fig_extrap_r2.png, fig_extrap_gap.png\n")
