# Validating DORM on the Harmony-demixed intrinsic sources

Do the batch-corrected, demixed intrinsic sources let the **signed-affine (intrinsic) prior improve
vanilla DORM**? Judged against the full benchmark suite (pooled EN, Simple avg, TransLasso, TransGLM,
PTL, target-oracle EN). Short answer: **no — the intrinsic prior does not improve DORM and hurts it
under composition shift; but vanilla DORM itself is validated as the robust label-free method.**

## Design
- **Intrinsic prior wiring** (`src/Functions3.R::maximin_s_beta`, `lib_dorm_common.R::dorm_fit`, new
  `prior_*` args): the signed-affine *objective* (which target-mixture to fit) is estimated in the
  **Harmony embedding** (batch-corrected RNA PCA, batch=donor), while the *validity constraint* stays
  on the gene-space ratio where the posterior applies ρ. Predictor A = screened genes, unchanged.
- This fixed the earlier catastrophic blow-up (full-gene intrinsic prior → R² ≈ −1e11): ρ is now finite
  and sums to 1. **Stability: solved.**

## Phase B — interior regime (donor-LODO, partial 6/48 combos)
DORM-intrinsic is stable but **below vanilla DORM** on every interior combo (e.g. s1d1/CD72_1 0.46 vs
0.75; s1d1/CD19_1 **−0.23** vs 0.71). Expected bias–variance cost of the enlarged set when the target is
inside the donor hull — but larger than benign.

## Phase C — extrapolation regime (the decisive test)
Unsupervised Harmony k-means states (annotation-free); a held-out donor is resampled to composition
`(1−α)·natural + α·(pure state)`, α: 0→1, leaving the source-donor hull. Figure: `figures/fig_extrap_r2.png`.

Median target R² by regime (label-free methods + oracle):

| method | target **inside** hull | target **outside** hull |
|---|---:|---:|
| **DORM (vanilla)** | **0.58** | **−0.08** |
| Pooled elastic-net | 0.45 | −1.34 |
| DORM_intrinsic (Harmony) | 0.27 | **−1.35** |
| Target-oracle EN (label-using) | 0.70 | 0.30 |

Representative sweeps:
- **s4d9/CD72_1:** α=0 (interior) DORM 0.58, intrinsic 0.27; α=0.5 (**just outside** hull, excess +0.14)
  DORM **0.61**, intrinsic **−1.44**; α=1 DORM −0.61, intrinsic −12.7.
- **s3d7/CD72_1:** DORM ~0.23–0.26 across α; intrinsic 0.12 → −0.19 → −2.31 → −4.53.

**The intrinsic prior collapses precisely where it was meant to help** (target outside the hull), while
**vanilla DORM degrades gracefully and stays the best label-free method** — at moderate extrapolation
(α=0.5) it even holds ~0.6 while pooling, simple-average, TransGLM/Lasso/PTL, and the intrinsic prior
all collapse. Only the label-using oracle beats it.

## Verdict
1. **Vanilla DORM is validated** as a robust label-free predictor under donor/composition shift — it
   dominates the benchmark suite both inside and (especially) outside the source hull.
2. **The intrinsic (signed-affine) extension does not add value on this data.** The Harmony embedding
   fixes its numerical stability and exposes real low-dimensional demixing structure (see
   `INTRINSIC_SUMMARY_formatted.md`), but routing that structure through the prior degrades prediction,
   badly under extrapolation. On this mixture data the observed-donor simplex is the better uncertainty
   set; the enlarged signed-affine set over-extrapolates.

## Caveats / what was not exhausted
- Hybrid prior design (Harmony objective + gene-space validity constraint); a fully-Harmony posterior or
  a larger ridge λ was not swept — but the collapse already appears at α=0.5 (a non-degenerate,
  ADT-variance-preserving shift), so the negative result is not merely an α=1 artifact (at α=1 the pure-
  state target has near-zero ADT variance, so all R², incl. oracle→0, are unreliable).
- Phase B stopped at 6/48 (background jobs were being killed by an environment limit); Phase C used the
  full s3d7 sweep + s4d9 extrapolation points (15 combos).

## Files
- `code/{make_harmony_embedding.py, run_lodo_donor_intrinsic.R (DORM_PRIOR_SPACE=harmony), run_extrapolation_shift.R, fig_extrapolation.R}`
- `results/{lodo_donor_intrinsic_harmony_results.csv (partial), extrapolation_shift_results.csv}`
- `figures/{fig_extrap_r2.png, fig_extrap_gap.png}`
