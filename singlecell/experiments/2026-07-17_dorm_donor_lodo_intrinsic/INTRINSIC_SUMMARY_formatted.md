# Intrinsic-Source (Signed-Affine) Machinery

## Handover summary

> **Scope.** These results are illustrations of the intrinsic-source and signed-affine machinery only. DORM prediction is not involved.  
> **Validity floor.** \(\varepsilon = 10^{-3}\).

---

## Notation

- \(L\): number of observed source distributions.
- \(m\): number of intrinsic distributions.
- \(P_1,\ldots,P_L\): observed source distributions.
- \(Q_1,\ldots,Q_m\): intrinsic distributions.
- \(\rho=(\rho_1,\ldots,\rho_L)^\top\): signed-affine weights over the observed sources.
- \(\bar r_\ell(x)\): density-ratio representation of \(P_\ell\).

The weights lie on the affine hyperplane

\[
\mathbf 1^\top \rho = 1.
\]

The signed-affine feasible set is

\[
\mathcal R
=
\left\{
\rho\in\mathbb R^L:
\mathbf 1^\top\rho=1,\quad
\sum_{\ell=1}^L \rho_\ell \bar r_\ell(x)\ge \varepsilon
\ \text{for every }x
\right\}.
\]

Unlike the probability simplex, \(\mathcal R\) permits negative coordinates as long as the induced mixture remains valid.

---

## 1. Enlarged signed-affine set and demixing recovery

This section considers two synthetic examples with known mixing matrices.

![Signed-affine feasible set for \(L=3\)](figures/intrinsic_signedset_L3.png)

### 1.1 Example A: \(L=3\), \(m=3\) — affinely independent sources

The observed distributions satisfy

\[
\begin{bmatrix}
P_1\\
P_2\\
P_3
\end{bmatrix}
=
\underbrace{
\begin{bmatrix}
0.6 & 0.3 & 0.1\\
0.2 & 0.6 & 0.2\\
0.1 & 0.3 & 0.6
\end{bmatrix}
}_{\Pi}
\begin{bmatrix}
Q_1\\
Q_2\\
Q_3
\end{bmatrix}.
\]

Equivalently,

\[
\begin{aligned}
P_1 &= 0.6Q_1+0.3Q_2+0.1Q_3,\\
P_2 &= 0.2Q_1+0.6Q_2+0.2Q_3,\\
P_3 &= 0.1Q_1+0.3Q_2+0.6Q_3.
\end{aligned}
\]

#### Geometry and recovery

| Quantity | Result |
|---|---:|
| Area of \(\mathcal R\) | \(3.32\) |
| Enlargement relative to the simplex | \(6.5\times\) |
| Coordinate range of \(\rho\) | \([-0.98,\,2.30]\) |
| Is \(\mathcal R\) bounded? | **Yes** |
| Maximum pairwise \(\kappa^*\), oracle ratios | \(0.0000\) |
| Maximum pairwise \(\kappa^*\), estimated ratios | \(0.0000\) |

Values near zero for the maximum pairwise \(\kappa^*\) indicate that the recovered components are mutually irreducible.

Define the backward coefficients by

\[
Q_j=\sum_{\ell=1}^L \alpha_{\ell j}P_\ell,
\qquad
\alpha=(\Pi^+)^\top,
\]

where \(\Pi^+\) denotes the Moore–Penrose pseudoinverse. The coefficient columns are:

| Method | Column 1 | Column 2 | Column 3 | Best-permutation Frobenius error |
|---|---|---|---|---:|
| True \((\Pi^+)^\top\) | \((2,-1,0)^\top\) | \((-0.67,2.33,-0.67)^\top\) | \((0,-1,2)^\top\) | — |
| Oracle recovery | \((-0.67,2.33,-0.67)^\top\) | \((0,-1,2)^\top\) | \((2,-1,0)^\top\) | \(0.00\) |
| Estimated recovery | \((-0.30,1.60,-0.29)^\top\) | \((-0.05,-0.35,1.40)^\top\) | \((1.25,-0.64,0.39)^\top\) | \(1.57\) |

The oracle procedure recovers the true columns exactly up to permutation.

### 1.2 Example B: \(L=3\), \(m=2\) — affinely dependent sources

The observed distributions satisfy

\[
\begin{bmatrix}
P_1\\
P_2\\
P_3
\end{bmatrix}
=
\underbrace{
\begin{bmatrix}
0.8 & 0.2\\
0.5 & 0.5\\
0.2 & 0.8
\end{bmatrix}
}_{\Pi}
\begin{bmatrix}
Q_1\\
Q_2
\end{bmatrix}.
\]

Equivalently,

\[
\begin{aligned}
P_1 &= 0.8Q_1+0.2Q_2,\\
P_2 &= 0.5Q_1+0.5Q_2,\\
P_3 &= 0.2Q_1+0.8Q_2.
\end{aligned}
\]

#### Geometry and recovery

| Quantity | Result |
|---|---:|
| Area of the displayed portion of \(\mathcal R\) | \(11.74\) |
| Displayed enlargement relative to the simplex | \(23.1\times\) |
| Displayed coordinate range of \(\rho\) | \([-3.00,\,4.00]\) |
| Is \(\mathcal R\) bounded? | **No** |
| Maximum pairwise \(\kappa^*\), oracle ratios | \(0.0000\) |
| Maximum pairwise \(\kappa^*\), estimated ratios | \(0.0000\) |

The backward coefficient columns are:

| Method | Column 1 | Column 2 | Best-permutation Frobenius error |
|---|---|---|---:|
| True \((\Pi^+)^\top\) | \((1.17,0.33,-0.50)^\top\) | \((-0.50,0.33,1.17)^\top\) | — |
| Oracle recovery | \((1.33,0,-0.33)^\top\) | \((-0.33,0,1.33)^\top\) | \(0.58\) |
| Estimated recovery | \((1.15,0,-0.15)^\top\) | \((-0.17,0,1.17)^\top\) | \(0.67\) |

Because \(m<L\), the representation of each \(Q_j\) in terms of the \(P_\ell\)'s is not unique: coefficients may differ by a null-space direction while inducing the same distribution.

### 1.3 Main synthetic-example conclusions

1. **When \(m=L=3\)**, the observed sources are affinely independent. The set \(\mathcal R\) is a bounded triangle approximately \(6.5\) times as large as the simplex. Its corners correspond to the three latent \(Q_j\)'s, which lie outside the observed simplex and are therefore “purer” than any observed source. Oracle demixing recovers them exactly up to permutation.

2. **When \(m=2<L=3\)**, the observed sources are affinely dependent. The set \(\mathcal R\) is unbounded: a null direction in the three-dimensional weight vector leaves the two-component mixture unchanged. This noncompactness is the geometric source of signed-prior instability when \(L\) is large.

---

## 2. Demixing the real data: \(L=12\) donors and \(m=3,\ldots,12\)

All 12 donors are used as observed sources. **The density ratio is now estimated in the Harmony
embedding** (batch-corrected RNA PCA, batch = donor; RNA-only, no ADT leakage) instead of the full
gene space, so shared cell-type structure is no longer masked by technical batch effects.

| $m$ | $\max_{j\neq k}\kappa^*(Q_j\mid Q_k)$ | Vertex donors | $\max_{\ell,j}\lvert\alpha_{\ell j}\rvert$ |
| :---: | ---: | --- | ---: |
| $3$ | 0.0000 | `s1d3`, `s3d7`, `s4d8` | 2.95 |
| $4$ | 0.0000 | `s1d3`, `s3d7`, `s4d8`, `s4d9` | 3.33 |
| $5$ | 0.0000 | `s1d3`, `s2d5`, `s3d7`, `s4d8`, `s4d9` | 3.79 |
| $6$ | 0.0000 | `s1d3`, `s2d1`, `s2d5`, `s3d7`, `s4d8`, `s4d9` | 5.35 |
| $7$ | 0.0000 | `s1d3`, `s2d1`, `s2d4`, `s2d5`, `s3d7`, `s4d8`, `s4d9` | 4.40 |
| $8$ | 0.0000 | `s1d3`, `s2d1`, `s2d4`, `s2d5`, `s3d7`, `s4d1`, `s4d8`, `s4d9` | 4.83 |
| $9$ | 0.0000 | `s1d1`, `s1d3`, `s2d1`, `s2d4`, `s2d5`, `s3d7`, `s4d1`, `s4d8`, `s4d9` | 6.22 |
| $10$ | 0.0000 | `s1d1`, `s1d2`, `s1d3`, `s2d1`, `s2d4`, `s2d5`, `s3d7`, `s4d1`, `s4d8`, `s4d9` | 7.35 |
| $11$ | 0.0000 | `s1d1`, `s1d2`, `s1d3`, `s2d1`, `s2d4`, `s2d5`, `s3d6`, `s3d7`, `s4d1`, `s4d8`, `s4d9` | 8.32 |
| $12$ | 0.0000 | `s1d1`, `s1d2`, `s1d3`, `s2d1`, `s2d4`, `s2d5`, `s3d1`, `s3d6`, `s3d7`, `s4d1`, `s4d8`, `s4d9` | 7.30 |

Unlike the full-gene-space analysis (where every pair had \(\kappa^*\approx0\) and 4% of pairs were
reducible), in the embedding **100% of donor pairs are now mutually reducible** (mean pairwise
\(\kappa^*=0.518\)): the sources genuinely share structure. The demixer, however, still selects individual
donors as vertices \(Q_j\) — it operates on the 12 donor-mixtures, so it relabels donors rather than
isolating cell-type programs (which would require demixing over cells or clusters).

---

## 3. Affine effective dimension and forward reconstruction

### 3.1 Meaning of “affine effective dimension”

Stack the \(L\) donor density-ratio vectors \(\bar r_\ell\), evaluated on the pooled sample, as rows of a matrix. Center each column and compute its singular values \(s_1,s_2,\ldots\).

The participation ratio is

\[
\operatorname{PR}
=
\frac{\left(\sum_k s_k\right)^2}
{\sum_k s_k^2}.
\]

It acts as a continuous or “soft” rank:

- \(\operatorname{PR}\approx 1\) when one direction dominates;
- \(\operatorname{PR}\approx d\) when approximately \(d\) directions contribute comparably.

Because the observed sources lie on an affine slice of the form \(\pi^\top r=1\), the maximum possible affine dimension is \(L-1\).

#### Observed values (Harmony embedding)

- All \(L=12\) donors: \(\operatorname{PR}=5.9\) out of a maximum of 11.
  The singular values relative to \(s_1\) are
  \[
  1.00,\ 0.93,\ 0.39,\ 0.34,\ 0.25,\ 0.19,\ 0.17,\ 0.12,\ 0.11,\ 0.10,\ 0.08,\ 0.00.
  \]
- With donor `s3d7` held out (\(L=11\)): \(\operatorname{PR}=5.6\) out of a maximum of 10.

For comparison, in the **full gene space** these were \(\operatorname{PR}=10.2/11\) and \(9.2/10\) — i.e. nearly
full rank. The embedding **roughly halves the affine dimension to \(\approx6\)**: a genuine low-dimensional
latent structure now exists (consistent with the high reducibility above), whereas in the gene space there
was essentially none.

### 3.2 Forward-reconstruction criterion

For each \(m\), every observed source is approximated by an affine combination of the estimated intrinsic sources:

\[
\widehat\Pi_{\ell\cdot}
\in
\underset{w\in\mathbb R^m}{\arg\min}
\left\|
\bar r_\ell-\sum_{j=1}^m w_j\bar Q_j
\right\|_2
\quad\text{subject to}\quad
\sum_{j=1}^m w_j=1.
\]

The reported diagnostic is the relative \(\ell_2\) residual. The \(m\) vertex donors reconstruct almost perfectly by construction; the \(L-m\) non-vertex donors are the meaningful diagnostic cases.

| $m$ | Median residual<br>(non-vertex) | Maximum residual<br>(non-vertex) | Non-vertex donors with<br>residual $<0.10$ | $\min_{\ell,j}\Pi_{\ell j}$ |
| :---: | ---: | ---: | :---: | ---: |
| $3$ | 0.163 | 0.252 | 0 / 9 | -0.12 |
| $4$ | 0.127 | 0.253 | 1 / 8 | -0.21 |
| $5$ | 0.109 | 0.178 | 0 / 7 | -0.32 |
| $6$ | 0.104 | 0.173 | 3 / 6 | -0.28 |
| $7$ | 0.087 | 0.169 | 3 / 5 | -0.28 |
| $8$ | 0.098 | 0.113 | 2 / 4 | -0.27 |
| $9$ | 0.081 | 0.087 | 3 / 3 | -0.27 |
| $10$ | 0.090 | 0.091 | 2 / 2 | -0.37 |
| $11$ | 0.088 | 0.088 | 1 / 1 | -0.30 |
| $12$ | — | — | 0 / 0 (all vertices) | -0.30 |

The non-vertex residuals are **lower than in the full gene space** (there the median stayed \(0.7–1.0\);
here it is \(\approx0.10\) at \(m=6\)), reflecting the shared structure the embedding exposes. They are
nonetheless far from zero and the forward weights still turn negative, so a small set of \(Q_j\) still does not
reconstruct all donors — the compression is partial, not clean.

---

## 4. Detailed backward and forward coefficients

The backward and forward representations are

\[
Q_j=\sum_{\ell=1}^L\alpha_{\ell j}P_\ell,
\qquad
P_\ell\approx\sum_{j=1}^m\Pi_{\ell j}Q_j.
\]

All coefficients below are computed in the Harmony embedding space. To keep the tables readable, \(\alpha\) and \(\Pi\) are separated for \(m=6\) and \(m=12\).

### 4.1 Results for \(m=3\)  (vertices: `s1d3`, `s3d7`, `s4d8`)

| Donor | Vertex | $\alpha_{\ell 1}$ | $\alpha_{\ell 2}$ | $\alpha_{\ell 3}$ | $\Pi_{\ell 1}$ | $\Pi_{\ell 2}$ | $\Pi_{\ell 3}$ | Residual |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | — | 0.00 | 0.00 | 0.00 | 0.24 | -0.06 | 0.81 | 0.14 |
| `s1d2` | — | 0.00 | 0.00 | 0.00 | 0.35 | -0.12 | 0.77 | 0.16 |
| `s1d3` | Yes | -0.76 | -1.41 | 1.11 | 0.16 | -0.08 | 0.92 | 0.03 |
| `s2d1` | — | 0.00 | 0.00 | 0.00 | 0.42 | -0.06 | 0.64 | 0.15 |
| `s2d4` | — | 0.00 | 0.00 | 0.00 | 0.31 | -0.06 | 0.75 | 0.17 |
| `s2d5` | — | 0.00 | 0.00 | 0.00 | 0.14 | 0.17 | 0.70 | 0.25 |
| `s3d1` | — | 0.00 | 0.00 | 0.00 | 0.47 | -0.04 | 0.57 | 0.12 |
| `s3d6` | — | 0.00 | 0.00 | 0.00 | 0.32 | -0.07 | 0.75 | 0.13 |
| `s3d7` | Yes | 2.02 | -0.53 | -0.42 | 0.58 | 0.00 | 0.41 | 0.03 |
| `s4d1` | — | 0.00 | 0.00 | 0.00 | 0.24 | -0.07 | 0.84 | 0.18 |
| `s4d8` | Yes | -0.26 | 2.95 | 0.31 | 0.18 | 0.30 | 0.51 | 0.02 |
| `s4d9` | — | 0.00 | 0.00 | 0.00 | 0.37 | -0.12 | 0.75 | 0.24 |

### 4.2 Results for \(m=6\)  (vertices: `s1d3`, `s2d1`, `s2d5`, `s3d7`, `s4d8`, `s4d9`)

#### Backward coefficients \(\alpha\)

| Donor | Vertex | $\alpha_{\ell 1}$ | $\alpha_{\ell 2}$ | $\alpha_{\ell 3}$ | $\alpha_{\ell 4}$ | $\alpha_{\ell 5}$ | $\alpha_{\ell 6}$ |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | — | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `s1d2` | — | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `s1d3` | Yes | -1.63 | -0.51 | -0.81 | -1.77 | -1.64 | 1.75 |
| `s2d1` | Yes | 0.10 | -5.35 | -2.13 | -2.80 | 3.21 | 0.04 |
| `s2d4` | — | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `s2d5` | Yes | -1.65 | 0.87 | 0.46 | 3.48 | 0.10 | 0.26 |
| `s3d1` | — | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `s3d6` | — | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `s3d7` | Yes | -1.99 | 1.50 | 2.81 | 2.35 | -0.51 | 0.57 |
| `s4d1` | — | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| `s4d8` | Yes | 4.63 | 0.17 | -0.12 | -1.39 | -0.05 | -0.69 |
| `s4d9` | Yes | 1.55 | 4.32 | 0.78 | 1.14 | -0.11 | -0.93 |

#### Forward weights \(\Pi\) and reconstruction residuals

| Donor | Vertex | $\Pi_{\ell 1}$ | $\Pi_{\ell 2}$ | $\Pi_{\ell 3}$ | $\Pi_{\ell 4}$ | $\Pi_{\ell 5}$ | $\Pi_{\ell 6}$ | Residual |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | — | 0.08 | 0.14 | -0.12 | -0.04 | 0.31 | 0.62 | 0.11 |
| `s1d2` | — | 0.01 | 0.26 | -0.04 | -0.13 | 0.39 | 0.51 | 0.10 |
| `s1d3` | Yes | 0.07 | 0.17 | -0.11 | -0.04 | 0.15 | 0.76 | 0.06 |
| `s2d1` | Yes | 0.04 | 0.11 | 0.01 | -0.05 | 0.45 | 0.45 | 0.04 |
| `s2d4` | — | 0.05 | 0.16 | -0.07 | -0.03 | 0.35 | 0.54 | 0.13 |
| `s2d5` | Yes | 0.17 | 0.02 | -0.28 | 0.36 | 0.17 | 0.56 | 0.05 |
| `s3d1` | — | 0.04 | 0.08 | 0.13 | -0.07 | 0.39 | 0.43 | 0.08 |
| `s3d6` | — | 0.05 | 0.15 | -0.05 | -0.04 | 0.33 | 0.56 | 0.09 |
| `s3d7` | Yes | 0.03 | 0.02 | 0.40 | -0.07 | 0.24 | 0.38 | 0.04 |
| `s4d1` | — | 0.06 | 0.19 | -0.11 | -0.05 | 0.28 | 0.63 | 0.17 |
| `s4d8` | Yes | 0.32 | -0.05 | 0.05 | 0.11 | 0.07 | 0.49 | 0.05 |
| `s4d9` | Yes | -0.01 | 0.38 | -0.08 | -0.11 | 0.44 | 0.38 | 0.07 |

### 4.3 Results for \(m=12\)

Because the \(m=12\) matrices are wide, each is divided into two column blocks.

#### Backward coefficients \(\alpha\): columns 1–6

| Donor | Vertex | $\alpha_{\ell 1}$ | $\alpha_{\ell 2}$ | $\alpha_{\ell 3}$ | $\alpha_{\ell 4}$ | $\alpha_{\ell 5}$ | $\alpha_{\ell 6}$ |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | Yes | -2.28 | -1.92 | -3.29 | -0.01 | -2.97 | -1.73 |
| `s1d2` | Yes | -1.22 | -5.07 | -2.97 | -2.71 | -2.53 | -0.66 |
| `s1d3` | Yes | -0.74 | -0.12 | -0.99 | -0.97 | -0.89 | 0.10 |
| `s2d1` | Yes | -2.52 | -6.72 | -4.62 | 0.05 | -5.21 | -3.27 |
| `s2d4` | Yes | -0.05 | -1.19 | -1.32 | -1.07 | -0.76 | -0.80 |
| `s2d5` | Yes | -1.18 | 0.14 | 0.59 | 0.54 | 0.21 | 0.59 |
| `s3d1` | Yes | -2.26 | 2.55 | 1.23 | -7.30 | -1.99 | 3.23 |
| `s3d6` | Yes | 2.87 | 3.99 | 2.53 | 0.46 | 7.23 | 0.18 |
| `s3d7` | Yes | 0.80 | 0.35 | 1.12 | 5.17 | 3.29 | 1.94 |
| `s4d1` | Yes | 1.48 | 1.34 | 5.24 | 2.74 | 1.84 | 0.95 |
| `s4d8` | Yes | 4.19 | 1.25 | 0.78 | 0.76 | 0.77 | -0.20 |
| `s4d9` | Yes | 1.92 | 6.40 | 2.69 | 3.33 | 2.01 | 0.67 |

#### Backward coefficients \(\alpha\): columns 7–12

| Donor | Vertex | $\alpha_{\ell 7}$ | $\alpha_{\ell 8}$ | $\alpha_{\ell 9}$ | $\alpha_{\ell 10}$ | $\alpha_{\ell 11}$ | $\alpha_{\ell 12}$ |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | Yes | -2.73 | -2.38 | -4.30 | -5.65 | -4.02 | 2.35 |
| `s1d2` | Yes | 1.27 | -0.60 | -2.41 | -3.34 | 2.56 | 0.69 |
| `s1d3` | Yes | -0.08 | -0.18 | 0.27 | 3.95 | 0.95 | -0.51 |
| `s2d1` | Yes | -6.33 | -5.58 | 3.11 | 3.67 | 1.45 | -0.06 |
| `s2d4` | Yes | -0.26 | 4.20 | 2.08 | 1.39 | 1.28 | -0.33 |
| `s2d5` | Yes | 2.93 | 1.17 | 0.65 | 1.08 | -0.04 | 0.06 |
| `s3d1` | Yes | 1.83 | 1.32 | 0.25 | 0.90 | -0.14 | 0.04 |
| `s3d6` | Yes | 1.90 | 0.74 | 0.63 | -0.48 | -0.10 | -0.37 |
| `s3d7` | Yes | 2.74 | 2.67 | -0.50 | -0.36 | -0.51 | 0.11 |
| `s4d1` | Yes | 0.22 | 0.19 | 0.34 | 0.10 | -0.22 | -0.14 |
| `s4d8` | Yes | -0.15 | -0.41 | 0.48 | -0.57 | 0.27 | -0.43 |
| `s4d9` | Yes | -0.35 | -0.15 | 0.39 | 0.32 | -0.46 | -0.41 |

#### Forward weights \(\Pi\): columns 1–6

| Donor | Vertex | $\Pi_{\ell 1}$ | $\Pi_{\ell 2}$ | $\Pi_{\ell 3}$ | $\Pi_{\ell 4}$ | $\Pi_{\ell 5}$ | $\Pi_{\ell 6}$ |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | Yes | 0.07 | 0.04 | 0.01 | -0.04 | 0.03 | -0.05 |
| `s1d2` | Yes | 0.02 | 0.05 | 0.01 | -0.02 | 0.03 | 0.02 |
| `s1d3` | Yes | 0.12 | 0.02 | 0.02 | -0.06 | 0.07 | -0.08 |
| `s2d1` | Yes | 0.02 | 0.01 | -0.03 | -0.00 | 0.05 | 0.11 |
| `s2d4` | Yes | 0.04 | 0.04 | 0.02 | -0.04 | 0.02 | -0.09 |
| `s2d5` | Yes | 0.02 | 0.06 | 0.05 | 0.01 | -0.09 | -0.26 |
| `s3d1` | Yes | 0.04 | 0.03 | -0.02 | -0.07 | 0.02 | 0.25 |
| `s3d6` | Yes | 0.02 | 0.03 | -0.01 | -0.09 | 0.19 | -0.03 |
| `s3d7` | Yes | 0.05 | -0.05 | -0.11 | 0.05 | 0.10 | 0.35 |
| `s4d1` | Yes | 0.04 | -0.05 | 0.25 | -0.04 | 0.02 | -0.05 |
| `s4d8` | Yes | 0.32 | -0.02 | -0.04 | -0.05 | -0.06 | 0.04 |
| `s4d9` | Yes | -0.05 | 0.22 | -0.04 | 0.06 | -0.02 | 0.02 |

#### Forward weights \(\Pi\): columns 7–12 and residuals

| Donor | Vertex | $\Pi_{\ell 7}$ | $\Pi_{\ell 8}$ | $\Pi_{\ell 9}$ | $\Pi_{\ell 10}$ | $\Pi_{\ell 11}$ | $\Pi_{\ell 12}$ | Residual |
| --- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `s1d1` | Yes | -0.02 | 0.03 | -0.01 | 0.11 | 0.03 | 0.81 | 0.05 |
| `s1d2` | Yes | -0.01 | -0.02 | -0.04 | 0.01 | 0.36 | 0.59 | 0.06 |
| `s1d3` | Yes | -0.03 | 0.04 | -0.30 | 0.36 | 0.15 | 0.69 | 0.07 |
| `s2d1` | Yes | -0.05 | -0.07 | 0.18 | 0.03 | 0.14 | 0.60 | 0.05 |
| `s2d4` | Yes | -0.09 | 0.19 | 0.07 | 0.06 | 0.16 | 0.62 | 0.06 |
| `s2d5` | Yes | 0.34 | -0.09 | 0.23 | 0.08 | -0.01 | 0.67 | 0.06 |
| `s3d1` | Yes | -0.05 | -0.02 | 0.11 | 0.03 | 0.10 | 0.58 | 0.05 |
| `s3d6` | Yes | -0.02 | -0.03 | 0.06 | 0.08 | 0.14 | 0.66 | 0.06 |
| `s3d7` | Yes | -0.08 | 0.00 | 0.07 | 0.02 | 0.11 | 0.46 | 0.05 |
| `s4d1` | Yes | -0.07 | 0.03 | -0.06 | 0.12 | 0.15 | 0.66 | 0.07 |
| `s4d8` | Yes | 0.10 | -0.05 | 0.11 | 0.08 | -0.01 | 0.57 | 0.06 |
| `s4d9` | Yes | -0.05 | -0.02 | 0.06 | -0.01 | 0.32 | 0.53 | 0.08 |

---

## 5. Interpretation

- **Batch-robust embedding restores shared structure.** Moving the density ratio into the Harmony space
  drops the (soft) affine dimension from \(\approx L-1\) (full genes) to \(\approx6\) and makes 100% of donor
  pairs reducible (up from 4%). A genuine low-dimensional latent structure now exists.
- **Forward reconstruction improves sharply.** Non-vertex donors, essentially unreconstructable in the gene
  space (residual \(0.7–1.0\)), now reconstruct with median residual \(\approx0.10\) at \(m=6\): about six
  intrinsic sources already capture the donors well.
- **But the demixer still returns donor-vertices, not cell types.** \(\alpha\) stays supported on the selected
  donors, and its entries are now large (up to \(\approx8\)) because several embedding directions are weak —
  mild ill-conditioning, not clean cell-type separation (which would require demixing over cells/clusters).

\[
\boxed{
\begin{array}{c}
\text{In the batch-robust embedding the donor sources share a genuine low-dimensional structure}\\
\text{(soft affine dim }\approx6\text{; about six intrinsic sources reconstruct the donors), unlike the gene space.}\\
\text{The demixer still selects donor-vertices rather than cell-type programs.}\\
\text{Whether this improves the DORM-intrinsic predictor is a separate, not-yet-run question (Phase 4).}
\end{array}
}
\]
