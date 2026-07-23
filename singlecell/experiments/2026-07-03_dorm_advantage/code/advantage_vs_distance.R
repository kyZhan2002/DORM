#!/usr/bin/env Rscript
## Key finding figure: DORM's advantage over naive pooling grows with how far the
## target is from the pooled (sample-weighted) training distribution -- the regime
## transfer learning is designed for, and diagnosable from covariates alone.
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")

m <- fread(file.path(res_dir, "mixture_results.csv"))
nl <- c(756, 2019, 717); cen <- nl / sum(nl)                 # pooled centroid = source sample proportions
w <- dcast(m[method %in% c("DORM (tuned smax)", "Pooled elastic-net")],
           ADT + mixture + seed + wP + wE + wN ~ method, value.var = "r2")
w[, dist := sqrt((wP-cen[1])^2 + (wE-cen[2])^2 + (wN-cen[3])^2)]
w[, gap := `DORM (tuned smax)` - `Pooled elastic-net`]
agg <- w[, .(dist = mean(dist), gap = mean(gap)), by = .(ADT, mixture)]

png(file.path(fig_dir, "fig4_advantage_vs_distance.png"), 780, 560, res = 110)
par(mar = c(4.5, 4.5, 3, 1))
plot(agg$dist, agg$gap, pch = 19, col = ifelse(agg$ADT == "CD71", "#1b7837", "#762a83"),
     xlab = "target distance from pooled training distribution",
     ylab = "R^2 gap:  DORM (tuned) - Pooled elastic-net",
     main = "DORM's advantage grows with target shift")
abline(h = 0, lty = 2, col = "grey50")
fit <- lm(gap ~ dist, data = agg); abline(fit, col = "#d73027", lwd = 2)
legend("topleft", c("CD71", "CD49d", sprintf("trend (r=%.2f)", cor(agg$dist, agg$gap))),
       pch = c(19, 19, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
       col = c("#1b7837", "#762a83", "#d73027"), bty = "n", cex = 0.85)
dev.off()

tab <- w[, .(mean_dist = round(mean(dist),2),
             DORM = round(mean(`DORM (tuned smax)`),3),
             Pooled = round(mean(`Pooled elastic-net`),3),
             gap = round(mean(gap),3)), by = mixture][order(mean_dist)]
fwrite(tab, file.path(res_dir, "advantage_vs_distance.csv"))
cat("cor(distance, gap) =", round(cor(agg$dist, agg$gap), 2), "\n")
cat("far targets (dist>median): mean gap =", round(w[dist>median(dist), mean(gap)], 3), "\n")
cat("near targets:               mean gap =", round(w[dist<=median(dist), mean(gap)], 3), "\n")
print(tab, row.names = FALSE)
cat("\nSaved fig4_advantage_vs_distance.png, advantage_vs_distance.csv\n")
