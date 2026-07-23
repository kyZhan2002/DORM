# Intrinsic-source (signed-affine) machinery — handover summary

Illustrations only; DORM prediction is not involved. eps (validity floor) = 1e-3.

## 1. Enlarged rho set + demixing recovery (two synthetic examples with KNOWN mixing)

rho = weights over the L=3 observed sources; lives on plane sum(rho)=1. Signed-affine set
R = {sum(rho)=1, induced mixture sum_l rho_l*rbar_l(x) >= eps for all x}. Figure: `figures/intrinsic_signedset_L3.png`.

### Example A: L=3, m=3 (independent)

How observed P are mixed from intrinsic Q (TRUE Pi):  P1 = 0.6 Q1 + 0.3 Q2 + 0.1 Q3; P2 = 0.2 Q1 + 0.6 Q2 + 0.2 Q3; P3 = 0.1 Q1 + 0.3 Q2 + 0.6 Q3
- Enlarged rho set: area(R) = 3.32, **R/simplex = 6.5x**, rho range [-0.98, 2.30], bounded = **yes**.
- Demixing kappa* (max pairwise): oracle-ratios = 0.0000, estimated-ratios = 0.0000  (~0 => recovered, irreducible).
- alpha (Q_j = sum_l alpha_lj P_l), columns = Q1..Q3:
    TRUE  t(pinv(Pi)) : (2,-1,0)  (-0.67,2.33,-0.67)  (0,-1,2)
    ORACLE-recovered : (-0.67,2.33,-0.67)  (0,-1,2)  (2,-1,0)   (Frobenius err vs true, best-perm = 0.00)
    EST-recovered    : (-0.3,1.6,-0.29)  (-0.05,-0.35,1.4)  (1.25,-0.64,0.39)   (err = 1.57)

### Example B: L=3, m=2 (dependent)

How observed P are mixed from intrinsic Q (TRUE Pi):  P1 = 0.8 Q1 + 0.2 Q2; P2 = 0.5 Q1 + 0.5 Q2; P3 = 0.2 Q1 + 0.8 Q2
- Enlarged rho set: area(R) = 11.74, **R/simplex = 23.1x**, rho range [-3.00, 4.00], bounded = **NO (unbounded)**.
- Demixing kappa* (max pairwise): oracle-ratios = 0.0000, estimated-ratios = 0.0000  (~0 => recovered, irreducible).
- alpha (Q_j = sum_l alpha_lj P_l), columns = Q1..Q2:
    TRUE  t(pinv(Pi)) : (1.17,0.33,-0.5)  (-0.5,0.33,1.17)
    ORACLE-recovered : (1.33,0,-0.33)  (-0.33,0,1.33)   (Frobenius err vs true, best-perm = 0.58)
    EST-recovered    : (1.15,0,-0.15)  (-0.17,0,1.17)   (err = 0.67)

**Takeaways.** (A) m=3=L: sources affinely INDEPENDENT => R is a BOUNDED triangle ~6.5x the simplex,
with the 3 true latent Q_j at its corners (outside the observed simplex = 'purer than any source');
demixing recovers Pi (oracle err 0.00). (B) m=2 < L=3: sources affinely DEPENDENT => R is UNBOUNDED
(an infinite strip; a null direction over the 3 weights leaves the 2-latent mixture unchanged) — this
non-compactness is the geometric root of the signed-prior instability at large L. Demixing still
recovers the 2 latent (oracle err 0.58); note alpha is only unique up to that null direction, so the
recovered coefs (which use the two extreme sources) may differ from t(pinv(Pi)) while giving the same Q.

## 2. Demixing the real data, L = 12 (all donors), m = 3..12

Sources = all 12 donors. kappa* ~ 0 => recovered components mutually irreducible.

| m | max pairwise kappa* | donors chosen as vertices | max alpha |
|---|---|---|---|
| 3 | 0.0000 | s1d2, s3d1, s3d7 | 1.03 |
| 4 | 0.0000 | s1d2, s3d1, s3d7, s4d9 | 1.06 |
| 5 | 0.0000 | s1d2, s2d5, s3d1, s3d7, s4d9 | 1.08 |
| 6 | 0.0000 | s1d2, s2d4, s2d5, s3d1, s3d7, s4d9 | 1.08 |
| 7 | 0.0000 | s1d2, s2d1, s2d4, s2d5, s3d1, s3d7, s4d9 | 1.12 |
| 8 | 0.0000 | s1d2, s2d1, s2d4, s2d5, s3d1, s3d7, s4d8, s4d9 | 1.12 |
| 9 | 0.0000 | s1d2, s1d3, s2d1, s2d4, s2d5, s3d1, s3d7, s4d8, s4d9 | 1.13 |
| 10 | 0.0000 | s1d2, s1d3, s2d1, s2d4, s2d5, s3d1, s3d6, s3d7, s4d8, s4d9 | 1.13 |
| 11 | 0.0000 | s1d1, s1d2, s1d3, s2d1, s2d4, s2d5, s3d1, s3d6, s3d7, s4d8, s4d9 | 1.18 |
| 12 | 0.0000 | s1d1, s1d2, s1d3, s2d1, s2d4, s2d5, s3d1, s3d6, s3d7, s4d1, s4d8, s4d9 | 1.26 |

**kappa*=0 is uninformative here** (donors are pairwise mutually irreducible; mean pairwise kappa* 0.016,
affine dim 9.2/10 -- see diag_pairwise_kappa.R), so the recovered 'intrinsic sources' are individual donors,
NOT cell types, and the count does not reduce below ~L.

## 3. Forward reconstruction (P_l = sum_j Pi_lj Q_j) and affine dimension

### What "affine effective dimension = 9.2 / 10" means

Stack the L donor density-ratio vectors rbar_l (evaluated on the pooled sample) as rows, center each
column, and take the SVD. The **participation ratio** PR = (sum_k s_k)^2 / sum_k s_k^2 of the singular
values s_k counts how many directions carry appreciable variance: PR ~ 1 if a single direction dominates,
PR ~ d if d directions contribute about equally. It is a soft, continuous rank. Because the observed
sources lie on the affine slice {prior^T r = 1}, the maximum possible is **L-1** (the sum constraint
removes one degree of freedom).
- diag_pairwise_kappa.R (L=11 sources, target s3d7 held out): **PR = 9.2 out of a max of 10**.
- All L=12 donors here: **PR = 10.2 out of a max of 11**; singular values rel. to s1: 1.00 0.81 0.80 0.73 0.62 0.57 0.55 0.53 0.47 0.41 0.40 0.00.
PR ~ L-1 means the donors span **almost the full affine space**: there is essentially no low-dimensional
latent structure to compress into a handful of Q_j. Together with pairwise kappa* ~ 0 (every donor is
mutually irreducible), the genuine intrinsic-source count is ~L -- the Q_j the demixer returns are the
donors themselves, not a smaller set of shared cell-type programs.

### Forward-reconstruction fidelity by m

For each m we express every observed source P_l as an affine mixture of the m estimated Q_j,
Pi_l = argmin_w || rbar_l - sum_j w_j Qbar_j ||  s.t. sum_j w_j = 1, and report the relative L2 residual.
The m *vertex* donors (the ones the demixer picked as Q_j) reconstruct ~perfectly by construction; the
L-m non-vertex donors are the diagnostic.

| m | median resid (non-vertex) | max resid (non-vertex) | # non-vertex recon well (<0.10) | most-negative Pi |
|---|---|---|---|---|
| 3 | 0.952 | 1.052 | 0 / 9 | 0.00 |
| 4 | 0.874 | 0.924 | 0 / 8 | 0.00 |
| 5 | 0.868 | 0.891 | 0 / 7 | -0.07 |
| 6 | 0.818 | 0.868 | 0 / 6 | -0.68 |
| 7 | 0.800 | 0.867 | 0 / 5 | -0.64 |
| 8 | 0.753 | 0.824 | 0 / 4 | -0.59 |
| 9 | 0.704 | 0.745 | 0 / 3 | -1.38 |
| 10 | 0.717 | 0.745 | 0 / 2 | -0.75 |
| 11 | 0.690 | 0.690 | 0 / 1 | -2.03 |
| 12 | - | - | 0 / 0 (all donors are vertices) | -1.30 |

Residuals stay high and the Pi weights go **negative** for the non-vertex donors until m approaches L,
confirming that no small set of Q_j reconstructs the left-out donors -- exactly what PR ~ L-1 predicts.

### m = 3  (backward alpha: Q_j = sum_l alpha_lj P_l  |  forward Pi: P_l = sum_j Pi_lj Q_j)

| donor | vertex? | a:Q1 | a:Q2 | a:Q3 | Pi:Q1 | Pi:Q2 | Pi:Q3 | resid |
|---|---|---|---|---|---|---|---|---|
| s1d1 |  | 0.00 | 0.00 | 0.00 | 0.23 | 0.15 | 0.63 | 0.90 |
| s1d2 | yes | -0.00 | -0.00 | 1.01 | 0.00 | 0.00 | 0.99 | 0.00 |
| s1d3 |  | 0.00 | 0.00 | 0.00 | 0.23 | 0.15 | 0.62 | 0.94 |
| s2d1 |  | 0.00 | 0.00 | 0.00 | 0.50 | 0.20 | 0.30 | 0.90 |
| s2d4 |  | 0.00 | 0.00 | 0.00 | 0.33 | 0.25 | 0.42 | 0.98 |
| s2d5 |  | 0.00 | 0.00 | 0.00 | 0.37 | 0.25 | 0.38 | 0.95 |
| s3d1 | yes | 1.00 | -0.03 | -0.00 | 1.00 | 0.00 | 0.00 | 0.00 |
| s3d6 |  | 0.00 | 0.00 | 0.00 | 0.34 | 0.30 | 0.36 | 0.86 |
| s3d7 | yes | -0.00 | 1.03 | -0.00 | 0.03 | 0.97 | 0.00 | 0.00 |
| s4d1 |  | 0.00 | 0.00 | 0.00 | 0.30 | 0.21 | 0.50 | 0.99 |
| s4d8 |  | 0.00 | 0.00 | 0.00 | 0.23 | 0.24 | 0.53 | 1.02 |
| s4d9 |  | 0.00 | 0.00 | 0.00 | 0.30 | 0.23 | 0.48 | 1.05 |

### m = 6  (backward alpha: Q_j = sum_l alpha_lj P_l  |  forward Pi: P_l = sum_j Pi_lj Q_j)

| donor | vertex? | a:Q1 | a:Q2 | a:Q3 | a:Q4 | a:Q5 | a:Q6 | Pi:Q1 | Pi:Q2 | Pi:Q3 | Pi:Q4 | Pi:Q5 | Pi:Q6 | resid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| s1d1 |  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.04 | 0.16 | 0.06 | 0.02 | 0.10 | 0.60 | 0.83 |
| s1d2 | yes | -0.00 | -0.00 | -0.00 | -0.00 | -0.01 | 0.73 | 0.02 | 0.03 | 0.13 | 0.13 | -0.68 | 1.38 | 0.04 |
| s1d3 |  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 | 0.24 | 0.01 | 0.05 | 0.17 | 0.53 | 0.83 |
| s2d1 |  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.08 | 0.07 | 0.26 | 0.21 | 0.27 | 0.11 | 0.81 |
| s2d4 | yes | -0.03 | -0.07 | -0.03 | -0.01 | 0.83 | 0.43 | 0.02 | 0.14 | -0.26 | -0.12 | 1.21 | 0.01 | 0.03 |
| s2d5 | yes | -0.02 | -0.01 | -0.01 | 0.86 | 0.09 | -0.03 | 0.02 | 0.06 | -0.25 | 1.16 | 0.01 | 0.01 | 0.01 |
| s3d1 | yes | -0.03 | -0.02 | 0.86 | 0.19 | 0.21 | 0.01 | 0.05 | -0.24 | 1.15 | 0.01 | 0.02 | 0.00 | 0.01 |
| s3d6 |  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.17 | 0.16 | 0.13 | 0.12 | 0.21 | 0.20 | 0.74 |
| s3d7 | yes | 1.08 | 0.11 | -0.03 | -0.04 | -0.04 | -0.04 | 0.93 | 0.00 | 0.02 | 0.02 | 0.03 | 0.01 | 0.01 |
| s4d1 |  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 | 0.39 | 0.00 | 0.17 | 0.17 | 0.25 | 0.73 |
| s4d8 |  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.04 | 0.47 | -0.03 | 0.12 | 0.10 | 0.30 | 0.87 |
| s4d9 | yes | -0.00 | 0.99 | 0.21 | -0.00 | -0.08 | -0.09 | -0.10 | 1.02 | 0.01 | 0.00 | 0.07 | 0.01 | 0.02 |

### m = 12  (backward alpha: Q_j = sum_l alpha_lj P_l  |  forward Pi: P_l = sum_j Pi_lj Q_j)

| donor | vertex? | a:Q1 | a:Q2 | a:Q3 | a:Q4 | a:Q5 | a:Q6 | a:Q7 | a:Q8 | a:Q9 | a:Q10 | a:Q11 | a:Q12 | Pi:Q1 | Pi:Q2 | Pi:Q3 | Pi:Q4 | Pi:Q5 | Pi:Q6 | Pi:Q7 | Pi:Q8 | Pi:Q9 | Pi:Q10 | Pi:Q11 | Pi:Q12 | resid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| s1d1 | yes | -0.03 | -0.03 | -0.04 | -0.03 | -0.02 | -0.02 | -0.02 | -0.05 | -0.04 | -0.09 | -0.10 | 0.56 | 0.05 | 0.06 | -0.13 | 0.07 | 0.11 | -0.24 | -0.06 | 0.01 | 0.84 | -1.09 | -0.25 | 1.62 | 0.09 |
| s1d2 | yes | -0.00 | 0.00 | -0.01 | -0.00 | -0.01 | -0.00 | 0.00 | -0.00 | 0.00 | -0.03 | 0.82 | 0.11 | 0.04 | 0.03 | 0.09 | -0.02 | 0.03 | 0.09 | -0.22 | 0.12 | 0.02 | -0.51 | 1.19 | 0.15 | 0.05 |
| s1d3 | yes | -0.08 | -0.04 | -0.12 | -0.03 | -0.03 | -0.01 | -0.02 | -0.05 | -0.02 | 0.52 | 0.19 | 0.39 | 0.06 | 0.03 | 0.11 | 0.00 | -0.15 | 0.03 | 0.16 | 0.05 | -1.30 | 1.75 | 0.03 | 0.23 | 0.10 |
| s2d1 | yes | -0.01 | -0.05 | -0.06 | -0.07 | -0.15 | -0.16 | -0.15 | -0.09 | 0.87 | 0.66 | 0.27 | 0.01 | -0.00 | -0.00 | 0.03 | 0.00 | 0.12 | 0.01 | 0.11 | -0.46 | 1.16 | -0.00 | -0.01 | 0.04 | 0.04 |
| s2d4 | yes | -0.00 | -0.05 | -0.08 | -0.04 | -0.10 | -0.04 | -0.02 | 0.93 | 0.40 | 0.28 | 0.03 | -0.02 | 0.00 | 0.08 | 0.04 | -0.07 | -0.05 | 0.30 | -0.51 | 1.05 | 0.05 | 0.03 | -0.01 | 0.09 | 0.04 |
| s2d5 | yes | -0.02 | -0.00 | -0.03 | -0.03 | -0.02 | -0.02 | 0.58 | 0.29 | 0.07 | -0.01 | 0.09 | -0.01 | 0.03 | -0.05 | -0.04 | 0.19 | 0.02 | -0.97 | 1.72 | -0.04 | 0.09 | 0.04 | 0.01 | 0.01 | 0.04 |
| s3d1 | yes | -0.00 | -0.01 | -0.02 | -0.02 | -0.05 | 0.79 | 0.46 | -0.01 | -0.05 | -0.13 | -0.03 | 0.09 | -0.03 | 0.02 | 0.10 | -0.12 | -0.36 | 1.24 | 0.04 | -0.05 | 0.15 | 0.00 | -0.00 | 0.02 | 0.04 |
| s3d6 | yes | -0.08 | -0.09 | -0.06 | -0.16 | 1.21 | 0.36 | 0.22 | 0.07 | -0.14 | 0.01 | -0.02 | 0.06 | 0.05 | -0.02 | -0.05 | -0.06 | 0.82 | 0.06 | 0.00 | 0.02 | 0.12 | 0.02 | 0.01 | 0.03 | 0.05 |
| s3d7 | yes | -0.00 | 0.13 | 0.02 | 0.99 | 0.08 | 0.12 | -0.04 | 0.01 | -0.00 | 0.01 | 0.00 | -0.05 | 0.04 | -0.05 | -0.29 | 1.02 | 0.11 | 0.02 | 0.03 | -0.00 | 0.08 | -0.01 | 0.00 | 0.03 | 0.03 |
| s4d1 | yes | -0.04 | -0.06 | 1.26 | 0.35 | 0.13 | -0.02 | -0.01 | -0.00 | -0.05 | -0.11 | -0.14 | 0.04 | 0.00 | -0.07 | 0.81 | -0.01 | 0.02 | 0.01 | 0.03 | 0.05 | -0.02 | 0.09 | 0.01 | 0.08 | 0.06 |
| s4d8 | yes | 1.17 | 0.05 | 0.02 | -0.04 | -0.09 | -0.00 | -0.01 | -0.02 | 0.01 | -0.04 | -0.06 | -0.08 | 0.87 | -0.08 | 0.03 | 0.01 | 0.05 | -0.01 | 0.03 | 0.01 | -0.04 | 0.07 | 0.00 | 0.05 | 0.04 |
| s4d9 | yes | 0.11 | 1.14 | 0.12 | 0.09 | 0.06 | -0.00 | 0.03 | -0.08 | -0.04 | -0.05 | -0.06 | -0.09 | -0.03 | 0.89 | 0.08 | -0.12 | 0.04 | 0.02 | -0.00 | 0.03 | 0.02 | 0.03 | -0.01 | 0.05 | 0.04 |

### Interpretation (m = 3, 6, 12)

- **alpha (backward)** is supported only on the m vertex donors: each Q_j ~ one donor (its alpha column is
~ a unit vector) and the other donors contribute 0. The demixer does not blend donors into shared latent
programs; it simply *re-labels* m of them as the 'pure' sources Q_j.
- **Pi (forward)** shows the consequence. A vertex donor reconstructs as Pi_l ~ e_j (all weight on its own
Q_j, resid ~ 0). A non-vertex donor needs a spread, often **signed**, combination of the Q_j and *still*
leaves a sizeable residual -- it carries structure that none of the m chosen vertices span.
- At **m = 12** every donor is its own vertex, so Pi ~ identity and residuals ~ 0: perfect but with zero
compression. There is no intermediate m at which a few Q_j reconstruct all 12 donors. That is the concrete
meaning of affine dim ~ L-1 and pairwise kappa* ~ 0: **on this data the intrinsic sources are the donors,**
**so the signed-affine extension enlarges rho only marginally and buys no real demixing here.**

