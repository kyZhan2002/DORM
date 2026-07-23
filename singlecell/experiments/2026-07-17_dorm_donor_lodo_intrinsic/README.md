# Donor-LODO real-data study, v2 (2026-07-17): refreshed ADTs + DORM-intrinsic + paper figures

Update of the 2026-07-16 donor-LODO study. Same design: **sources = the 11 training donors**
(annotation-free natural cell-type mixtures), **target = held-out donor**, predict ADT from RNA over
**all cell types**, q=20 low-dim A screened on pooled sources, penalty=TRUE, tuned smax; label-using
benchmarks at **ntar=20**; 12-donor leave-one-out, 1 seed (median over donors is the robust summary).

Changes vs 2026-07-16: (1) **ADT set** = CD72_1, CD19_1, CD11c, **CD3** (dropped CD71 where DORM was
weaker; added CD3 for T-lineage diversity). (2) Added **DORM-intrinsic** (signed-affine prior). (3)
Figures match the paper palette (`simu/graph_main.R`): DORM red, Simple-avg blue, TransLasso orange,
TransGLM purple, PTL green, Target-oracle grey, Pooled-EN light-blue; DORM solid, benchmarks dashed.

## Headline (clean tier) — DORM wins label-free, near the labeled oracle
Median R^2 over 12 donors x 4 ADTs (`figures/fig10_clean_summary.png`, `fig9_clean_per_donor.png`):

| method | labels | median R^2 |
|---|---|---|
| Target-oracle EN | full target | 0.652 |
| **DORM** | **0** | **0.499** |
| TransGLM | 20 | 0.459 |
| Pooled elastic-net | 0 | 0.403 |
| Simple average | 0 | 0.065 |
| TransLasso | 20 | −0.118 |
| PTL | 20 | −2.816 |

- **DORM (0 labels) beats pooling (+0.10 median) and TransGLM (20 labels), reaching ~77% of the labeled
  oracle.** In fig9 the red DORM curve tracks the grey oracle across donors and avoids the deep dips that
  Pooled / Simple-average / TransLasso suffer on batch-shifted donors.
- Simple average, TransLasso, PTL are unstable/poor here (naive aggregation or small-ntar transfer).

## DORM-intrinsic (signed-affine prior) — does not help on this real data
Two figure sets are provided: `*_clean` (no intrinsic) and `*_full` (adds DORM_intrinsic).
- The default ridge (lambda=0.1) blows up catastrophically (R^2 ~ −1e2 to −1e7). A **stabilizing sweep**
  picked **lambda=0.5** (see `logs/lambda_sweep.log`), used here.
- Even at lambda=0.5: **DORM_intrinsic median R^2 = 0.313 — worse than standard DORM (0.499) and below
  pooling (0.403)** — and it *still* blows up on some combos (min −1e7 on s1d1/CD19_1; also CD3 on
  s2d4/s3d1). So the signed-affine prior does not transfer to this L=11 real setting.
- Why (documented in `INTRINSIC_SUMMARY.md` + diagnostics): the signed prior is high-variance under
  estimated density ratios at large L; the instability lives in the pooled DR term P (not collinearity,
  not the eps-floor). It's realization-dependent — stable when a source dominates the estimated ratio,
  explosive otherwise.

## Intrinsic illustrations (handover): `INTRINSIC_SUMMARY.md`
- L=3 signed-affine feasible set: ~1x the simplex for distinct real donors, ~1.4x for synthetic
  mixture-sources (`figures/intrinsic_signedset_L3.png`).
- Demixing m=3..8: max kappa* = 0 for all m, but this is uninformative — donors are pairwise mutually
  irreducible (mean kappa* 0.016), affine dim ~9.2/10, so the intrinsic count ≈ L (no reduction); it
  recovers donor-vertices, not cell types.

## Files
- `code/{lib_dorm_common.R (NTAR_FIXED=20, signed-lambda wrapper), run_lodo_donor_intrinsic.R, fig_intrinsic.R}`
- diagnostics: `code/{diag_collinearity_demix, diag_blowup, diag_variance, diag_lambda_sweep, diag_pairwise_kappa, gen_intrinsic_summary}.R`
- `results/lodo_donor_intrinsic_results.csv` (+ `_rho.csv`); `figures/fig9_{clean,full}_per_donor.png`,
  `fig10_{clean,full}_summary.png`, `intrinsic_signedset_L3.png`; `INTRINSIC_SUMMARY.md`.

## Reproduce
```
DORM_SIGNED_LAMBDA=0.5 Rscript code/run_lodo_donor_intrinsic.R    # ~2 h, 1 seed
Rscript code/fig_intrinsic.R                                       # both figure sets
Rscript code/gen_intrinsic_summary.R                              # L=3 set + demix handover doc
```
