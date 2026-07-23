#!/usr/bin/env Rscript
## Analyse leave-one-donor-out (donors as sources) GEX->ADT: per-donor R^2 curves + worst-case/
## median/mean summary (DORM vs benchmarks), and the same-site rho mechanism check.
suppressPackageStartupMessages(library(data.table))
this_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1])
exp_dir <- dirname(dirname(this_file)); res_dir <- file.path(exp_dir, "results"); fig_dir <- file.path(exp_dir, "figures")
R <- fread(file.path(res_dir, "lodo_donor_as_source_results.csv"))
R[method == "DORM (tuned smax)", method := "DORM"]
ADTS <- intersect(c("CD71","CD72_1","CD19_1","CD11c"), unique(R$ADT))

## median over seeds
M <- R[, .(r2 = median(r2, na.rm = TRUE)), by = .(ADT, donor, site, method)]

## ---- summary: worst-case / median / mean over donors, per method ----
summ <- M[is.finite(r2), .(worstcase = round(min(r2),3), median = round(median(r2),3),
                          mean = round(mean(r2),3)), by = .(ADT, method)]
keymeth <- c("DORM","RAP","MI","Simple average","Best single source","Pooled elastic-net",
             "TransLasso","TransGLM","PTL","Target-oracle EN")
for (a in ADTS) { cat("\n==== donor-LODO [",a,"] worst-case / median / mean over donors ====\n")
  print(summ[ADT==a][match(intersect(keymeth,summ[ADT==a,method]),method)][order(-median)], row.names=FALSE) }
fwrite(summ, file.path(res_dir, "summary_lodo_donor_as_source.csv"))

## ---- same-site rho mechanism ----
wf <- file.path(res_dir, "lodo_donor_as_source_rho.csv")
if (file.exists(wf)) {
  W <- fread(wf)
  mech <- W[, .(same_site_rho = sum(rho[same_site]), n_same = sum(same_site), L = .N), by = .(target_donor, ADT, seed)]
  mech[, expected_uniform := n_same / L]
  cat(sprintf("\n== same-site rho concentration ==\n DORM puts mean %.2f of rho on same-site source donors; uniform baseline would be %.2f (%.1fx)\n",
              mech[, mean(same_site_rho)], mech[, mean(expected_uniform)],
              mech[, mean(same_site_rho)] / mech[, mean(expected_uniform)]))
}

## ---- per-donor R^2 curves ----
STYLE <- data.table(method = c("DORM","Pooled elastic-net","TransGLM","Target-oracle EN"),
                    col = c("#1b7837","#000000","#e08214","#777777"), lty = c(1,1,2,4))
pch_of <- function(lty) c(`1`=19,`2`=17,`4`=18)[as.character(lty)]
YLIM <- c(-0.6, 0.85)
np <- length(ADTS)
png(file.path(fig_dir, "fig8_lodo_donor_as_source.png"), 720*min(np,2), 560*ceiling(np/2), res = 115)
par(mfrow = c(ceiling(np/2), min(np,2)), mar = c(6.5, 5, 3, 1))
for (a in ADTS) {
  d <- M[ADT == a]; ord <- unique(d[order(site, donor), .(donor, site)])
  xs <- seq_len(nrow(ord)); sitecol <- c(s1="#1b9e77",s2="#d95f02",s3="#7570b3",s4="#e7298a")
  plot(NA, xlim=c(1,nrow(ord)), ylim=YLIM, xaxt="n", xlab="", ylab="target R^2 (clamped)",
       main=sprintf("%s: leave-one-donor-out (donors as sources)", a), cex.axis=1.2, cex.lab=1.35, cex.main=1.3)
  axis(1, at=xs, labels=ord$donor, las=2, cex.axis=1.0,
       col.axis="black"); abline(h=0, col="grey70", lty=3)
  ## color-code donor tick by site
  for (i in xs) axis(1, at=i, labels=ord$donor[i], las=2, cex.axis=1.0, col.axis=sitecol[ord$site[i]], tick=FALSE)
  for (i in seq_len(nrow(STYLE))) {
    dm <- d[method==STYLE$method[i]][match(ord$donor, donor)]
    y <- pmax(pmin(dm$r2, YLIM[2]), YLIM[1])
    lines(xs, y, col=STYLE$col[i], lty=STYLE$lty[i], lwd=if(STYLE$lty[i]==4) 3 else 2.4)
    points(xs, y, col=STYLE$col[i], pch=pch_of(STYLE$lty[i]), cex=1.1)
  }
  if (a==ADTS[1]) legend("bottomright", STYLE$method, col=STYLE$col, lty=STYLE$lty, lwd=2.4,
                         pch=pch_of(STYLE$lty), cex=0.95, bty="n")
}
dev.off(); cat("\nSaved fig8_lodo_donor_as_source.png (donor ticks colored by site)\n")
