# Leave-one-donor-out, donors as sources — where DORM shines (2026-07-16)

## Design
Annotation-free multi-source transfer, GEX→ADT, on the NeurIPS-2021 BMMC CITE-seq data.
- **Sources = the 11 training donors** (each = all its cells, a natural cell-type mixture). No
  cell-type annotations used anywhere (avoids the circularity of RNA/protein-derived labels).
- **Target = held-out 12th donor** (all cells); leave-one-donor-out over all 12. Batch codes are
  site+donor (s3d7 = site 3, donor 7) with nested donor+site batch effects.
- **Task:** predict a surface protein from RNA over ALL cell types (strong signal — erythroid-only
  was near-noise). **ADTs:** CD71, CD72_1, CD19_1, CD11c (screened as strongly RNA-predictable,
  pooled R² 0.66–0.76).
- **DORM:** q=20 low-dim A screened label-free on the pooled sources, `penalty=TRUE` (regularised Q;
  penF blew up under donor shift), tuned smax. Benchmarks: label-free (RAP, MI, Simple avg, Pooled EN)
  and label-using at **ntar=20** (TransLasso, TransGLM, PTL) + Target-oracle EN. 2 seeds; median over seeds.

## Result: DORM beats pooling and matches label-using methods, label-free
Median / mean R² over 12 donors (`results/summary_lodo_donor_as_source.csv`, `figures/fig8_*.png`):

| ADT | DORM median | Pooled median | TransGLM@20 | Oracle | DORM mean | Pooled mean |
|---|---|---|---|---|---|---|
| CD71 | 0.635 | 0.664 | 0.668 | 0.784 | 0.585 | 0.549 |
| CD72_1 | 0.546 | 0.448 | 0.446 | 0.673 | 0.423 | 0.230 |
| CD19_1 | 0.467 | 0.386 | 0.412 | 0.685 | 0.379 | 0.189 |
| CD11c | 0.626 | 0.583 | 0.594 | 0.755 | 0.599 | 0.421 |

- **DORM > Pooled on mean R² for all 4 ADTs; on median for 3/4** (ties/slightly trails on CD71).
  The edge is **robustness** — pooling collapses to negative R² on batch-shifted donors (e.g. CD72_1
  on s3d6: Pooled ≈ −0.6, DORM stays positive), which fig8 shows as DORM (green) avoiding the pooled
  (black) dips and tracking the oracle (grey).
- **DORM (0 target labels) ≥ TransGLM (20 labels) on median for 3/4 ADTs.**
- **DORM reaches ~70–85% of the labeled-oracle ceiling, label-free.**

## Mechanism
DORM's `rho` puts **0.51 of its weight on the target donor's same-site source donors, vs 0.18 under
uniform (2.8×)** — it discovers the nested batch structure from RNA alone and up-weights the most
target-like donors. Per-donor, high same-site `rho` tracks DORM winning; when it can't find similar
donors (low same-site `rho`, e.g. s3d1) DORM can lag pooling. See `results/*_rho.csv`.

## Honest caveats
- "Best single source" in the summary selects the best source by target MSE (peeks at target labels →
  oracle-ish); excluded from fig8.
- MI / Simple average / TransLasso / PTL are numerically unstable here (some −1e6 R² from donors with
  near-zero ADT variance + blow-ups); DORM is far more stable than these naive aggregations.
- s3d1 is a hard fold for every method (DORM dips there too).

## Files
`code/run_lodo_donor_as_source.R`, `code/fig_lodo_donor_as_source.R`, `code/lib_dorm_common.R`
(NTAR_FIXED=20); `results/lodo_donor_as_source_results.csv`, `..._rho.csv`,
`summary_lodo_donor_as_source.csv`; `figures/fig8_lodo_donor_as_source.png`.

## Reproduce
```
Rscript code/run_lodo_donor_as_source.R      # DORM_ADTS=... DORM_SEEDS=1,2 (~1.5 h)
Rscript code/fig_lodo_donor_as_source.R      # tables + fig8
```
