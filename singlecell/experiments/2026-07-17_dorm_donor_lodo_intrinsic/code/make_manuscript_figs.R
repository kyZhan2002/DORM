## Manuscript figures for Appendix (intrinsic sources): (A) enlarged feasible set R vs simplex,
## two examples (L=m=3 bounded; L=3,m=2 unbounded); (B) target mixing-weight recovery in the
## L=4 DORM-intrinsic simulation. Saved as vector PDFs to docs/Pictures/.
suppressPackageStartupMessages({ library(ggplot2); library(MASS); library(data.table) })
set.seed(1)
out <- "docs/Pictures"
oracle_rbar <- function(Xev, Pi, mu, sig) {
  ph <- sapply(seq_len(nrow(mu)), function(j) exp(-rowSums(sweep(Xev, 2, mu[j, ])^2) / (2 * sig^2)))
  pl <- ph %*% t(Pi); t(pl / rowMeans(pl))
}
g <- seq(-3, 4, length.out = 240); G <- as.matrix(expand.grid(u = g, v = g)); rho <- cbind(G[,1], G[,2], 1-G[,1]-G[,2])
cell <- (g[2]-g[1])^2; aS <- sum((G[,1]>=0)&(G[,2]>=0)&(G[,1]+G[,2]<=1))*cell
run <- function(Pi, mu, label) {
  L <- nrow(Pi); m <- ncol(Pi)
  X <- do.call(rbind, lapply(1:L, function(l){ comp <- sample(1:m, 700, TRUE, Pi[l,]); mu[comp,] + matrix(rnorm(700*ncol(mu), sd=1), 700, ncol(mu)) }))
  R <- oracle_rbar(X, Pi, mu, 1); Rs <- R[, sort(sample(ncol(R), 700))]
  inR <- apply(rho %*% Rs, 1, min) >= 1e-3; at <- t(ginv(Pi))
  list(dt = data.table(u=G[,1], v=G[,2], inR=inR, cfg=label), vt = data.table(u=at[1,], v=at[2,], cfg=label),
       ratio = sum(inR)*cell/aS)
}
mu3 <- 4*rbind(c(1,0,0,0),c(0,1,0,0),c(0,0,1,0)); Pi3 <- rbind(c(.6,.3,.1),c(.2,.6,.2),c(.1,.3,.6))
mu2 <- 4*rbind(c(1,0,0,0),c(0,1,0,0));            Pi2 <- rbind(c(.8,.2),c(.5,.5),c(.2,.8))
eA <- run(Pi3, mu3, sprintf("(a) L=3, m=3 (independent): R bounded, %.1fx simplex", run(Pi3,mu3,"")$ratio))
eB <- run(Pi2, mu2, "(b) L=3, m=2 (dependent): R unbounded")
tri <- data.frame(u=c(0,1,0), v=c(0,0,1))
dt <- rbind(eA$dt, eB$dt); vt <- rbind(eA$vt, eB$vt)
pA <- ggplot() + geom_raster(data=dt[inR==TRUE], aes(u,v), fill="#8B0000", alpha=0.5) +
  geom_polygon(data=tri, aes(u,v), fill=NA, color="red", linewidth=1) +
  geom_point(data=vt, aes(u,v), shape=17, size=2.6) +
  facet_wrap(~cfg) + coord_equal(xlim=c(-2.5,3.5), ylim=c(-2.5,3.5)) +
  labs(x=expression(bar(rho)[1]), y=expression(bar(rho)[2])) +
  theme_bw() + theme(strip.text=element_text(size=10))
ggsave(file.path(out, "intrinsic_signedset.pdf"), pA, width=8, height=4.2)

## Figure B: L=4 target mixing-weight recovery (numbers from Table intrinsic-dorm text)
W <- data.table(source=factor(paste0("P",1:4), levels=paste0("P",1:4)),
                Truth=c(.60,.39,.13,-.13), `Standard DORM`=c(.94,0,0,.06), `Signed-affine`=c(.58,.38,.06,-.02))
Wl <- melt(W, id.vars="source", variable.name="method", value.name="w")
Wl[, method := factor(method, levels=c("Truth","Standard DORM","Signed-affine"))]
pB <- ggplot(Wl, aes(source, w, fill=method)) +
  geom_col(position=position_dodge(0.75), width=0.7) + geom_hline(yintercept=0, color="grey40") +
  scale_fill_manual(values=c("Truth"="grey55","Standard DORM"="#4575b4","Signed-affine"="#FF0000")) +
  labs(x="observed source", y=expression("target mixing weight "*bar(rho)[l]), fill=NULL) +
  theme_bw() + theme(legend.position="bottom", axis.title=element_text(size=12), axis.text=element_text(size=11))
ggsave(file.path(out, "intrinsic_dorm_weights.pdf"), pB, width=6, height=4)
cat("wrote docs/Pictures/intrinsic_signedset.pdf and intrinsic_dorm_weights.pdf\n")
