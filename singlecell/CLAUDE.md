# Single-cell application (`singlecell/`)

_Loaded when working under `singlecell/`. The core estimator, data contract, and
simulations live in the repo-root `CLAUDE.md`._

Two-stage pipeline: Python preprocessing (notebooks) → R DORM run.

**Preprocessing (Python/scanpy, run once):**
1. `nips_preprocess_clean.ipynb` — from `fullbatch.h5ad` (paired CITE, 12 batches). RNA path:
   counts → HVG (Seurat v3) → normalize_total → log1p → scale (clip ±10). ADT path:
   counts → **CLR** (compositional, no clipping). Writes `processed/*_proc.h5ad`.
2. `nips_export_processed_to_R_csv.ipynb` — exports processed objects to
   `processed/csv_for_R/` as CSVs with R-safe, stable feature names (plus a
   `feature_name_map` back to original gene/protein names). This is the boundary the R
   script reads from.

**Running DORM** — `singlecell/run_dorm_from_cite_nips_csv.R`:
```
Rscript singlecell/run_dorm_from_cite_nips_csv.R                       # defaults
DORM_SOURCE_GROUPS=s2d5,s3d1 DORM_TARGET_GROUPS=s3d7 Rscript singlecell/run_dorm_from_cite_nips_csv.R
```
This script builds the DORM data contract from the CSVs, runs `get_DORM_beta`, runs the
transfer baselines, and writes results to `singlecell/dorm_results/`
(`*_dorm_data.rds`, `*_dorm_result.rds`, `*_summary.csv`). It is the reference example of
wiring real data into the estimator; mirror its `make_x`/`make_y`/split logic when adding data.

Configured entirely via env vars (all optional):
- `DORM_PARTITION` — what defines a domain: `batch` (Ver1, default), `major_cell_type`
  (Ver2), `t_cell_subtype` (Ver3, T cells only). See `partition_configs` in the script.
- `DORM_SOURCE_GROUPS` / `DORM_TARGET_GROUPS` — comma-separated group ids; run is repeated
  per target group. Source/target must be disjoint.
- `DORM_PROTEIN_Y` — which ADT protein is the response (default `CD72_1`); `auto` screens all
  proteins by CV elastic-net R² and picks the most RNA-predictable.
- `DORM_RUN=false` — build/save the input object only (smoke test, no fit).
- `DORM_RUN_TRANSFER=false` — skip TransLasso/TransGLM/PTL baselines.
- `DORM_MAX_RNA_FEATURES`, `DORM_MAX_SOURCE_ROWS`, `DORM_MAX_TARGET_ROWS`,
  `DORM_BENCHMARK_NTAR` — size/subsampling knobs (see `cfg` block).

The script auto-locates the repo root by finding `src/Functions3.R`, so it runs from either
the repo root or `singlecell/`.

---

# Handover: intrinsic-source (DORM-intrinsic) study — continue on cluster

_Status as of 2026-07-19. Work was done on a laptop where **background jobs kept getting killed
by an environment limit** (both a ~2h donor-LODO run and the extrapolation sweep were stopped
mid-way). That is why several runs below are **partial** — rerun them to completion on the
cluster, where long jobs survive. All scripts are deterministic given their seeds and are
idempotent (each writes/overwrites its own results CSV; pass `DORM_APPEND=1` to add to an
existing one)._

## Where everything lives
`singlecell/experiments/2026-07-17_dorm_donor_lodo_intrinsic/` — the active study.
- `code/lib_dorm_common.R` — shared lib: `load_cite_data`, `screen_A` (label-free elastic-net
  gene screen → the q=20 low-dim A), `make_X`, `make_X_embed`, `load_embedding`, `dorm_fit`
  (+ optional `prior_*` args), `eval_all_methods` (DORM + **all benchmarks**), `dorm_tuned_r2`,
  `pooled_en`, `target_oracle_en`. Constants: `NTAR_FIXED=20`, `SMAX_GRID`, `SIGNED_LAMBDA`.
- `code/make_harmony_embedding.py` — builds the batch-robust embedding (see below).
- `code/run_lodo_donor_intrinsic.R` — donor-LODO: 12 donors as annotation-free sources, one
  held out as target; standard DORM + benchmarks + DORM_intrinsic. `DORM_PRIOR_SPACE=harmony`
  switches the intrinsic prior into the Harmony space.
- `code/run_extrapolation_shift.R` — the decisive extrapolation test (composition-shifted targets).
- `code/fig_intrinsic.R`, `code/fig_extrapolation.R` — figures (paper palette).
- `code/{diag_embed_space.R, diag_drdim.R, diag_pairwise_kappa.R, gen_intrinsic_summary.R,
  update_formatted_harmony.R}` — diagnostics + the two handover write-ups.
- Write-ups (read these first): `VALIDATION_SUMMARY.md` (headline result),
  `INTRINSIC_SUMMARY_formatted.md` (demixing, Harmony space), `README.md`.
- Core estimator changes (additive, backward-compatible): `src/Functions3.R::maximin_s_beta`
  gained `prior_Xlist/prior_Xtrainlist/prior_X0/prior_X0train`; `src/IntrinsicSources.R` holds
  `signed_affine_prior()` + Katz-Samuels demixing (`recover_intrinsic_sources`).

## Design recap (the vanilla DORM real-data setup — this is the headline)
- **Sources = donors/batches** (annotation-free natural cell-type mixtures), **never cell-type
  annotations** (circular: they're RNA/protein-derived). Target = held-out donor. Predict an ADT
  from RNA over **all** cell types.
- **q=20** low-dim shared component **A** = genes chosen by `screen_A` (label-free EN on pooled
  sources); **p=250** HVG total (`DORM_HVG=250`) = A(20) + W(230). `penalty=TRUE`, tuned `smax`.
- **Density ratio uses `condA=FALSE` → the full 250-dim X** (A+W), not just A. The predictor
  `beta` uses only A (q+1 coefs). Benchmarks (`eval_all_methods`): Pooled EN, Simple avg,
  TransLasso, TransGLM, PTL, Target-oracle EN.

## The Harmony embedding (intrinsic prior space)
`make_harmony_embedding.py` runs Harmony (`harmonypy`, `pip install harmonypy` into the repo
`.venv`) on the dataset's 50-dim `GEX_X_pca` with **batch=donor**, RNA-only (no ADT → no
leakage), and writes `processed/csv_for_R/cite_nips_EMBED_for_R.csv` = `obs_id` + `H1..H50`
(Harmony) + `PC1..PC50` (raw PCA). Regenerate if the CSV is missing:
`.venv/bin/python .../code/make_harmony_embedding.py`.

## How the intrinsic prior is wired (`prior_*`)
`dorm_fit(..., signed_prior=TRUE, prior_Xlist=, prior_Xtrainlist=, prior_X0=, prior_X0train=)`
builds those from `make_X_embed` on the **same row splits** as the gene splits. In
`maximin_s_beta`, when `prior_*` is given, the signed-affine **objective** (which target-mixture
to fit) is estimated in the embedding, but the **validity constraint stays on the gene-space
ratio** where the posterior applies ρ (`Functions3.R` ~line 785 block + the `rho=` line). This
fixed a catastrophic blow-up (full-gene intrinsic prior gave R²≈−1e11). Default (no `prior_*`)
reproduces standard DORM exactly.

## Findings so far (do NOT need to re-derive)
1. **Vanilla DORM is the robust label-free winner** on donor-LODO: median R²≈0.50 over 12
   donors × 4 ADTs, beats pooling/TransGLM/TransLasso/PTL, near the label-using oracle
   (fig9/fig10; `README.md`).
2. **Intrinsic prior in the full 250-dim gene space blows up** (R² −1e2…−1e11); demixing finds
   no reduction (every donor pair mutually irreducible, affine dim ≈ L−1). Root cause = the
   density-ratio **feature dimension**: batch effects make donors trivially separable
   (`diag_drdim.R`, `diag_pairwise_kappa.R`).
3. **Harmony embedding fixes stability + reveals structure** (affine dim 10.2→5.9, 100% donor
   pairs reducible, forward-recon residual 0.10 at m=6) — `INTRINSIC_SUMMARY_formatted.md`,
   `diag_embed_space.R`.
4. **But the intrinsic prior does NOT improve prediction.** Interior donor-LODO: below vanilla
   DORM. Extrapolation test (`run_extrapolation_shift.R`): as the target leaves the source-donor
   hull, DORM_intrinsic **collapses** (R² −1.4 to −12 at moderate shift) while **vanilla DORM
   degrades gracefully and stays best label-free**. Conclusion: on this mixture data the
   observed-donor simplex is the better uncertainty set; the signed-affine enlargement
   over-extrapolates. **`VALIDATION_SUMMARY.md` is the authoritative write-up.**
   (Guiding principle, from memory: the intrinsic extension exists to *improve* vanilla DORM and
   must always be judged against all benchmarks — here it does not, and we report that honestly.)

## To finish on the cluster (rerun to completion — no kill limit there)
1. **Phase B — full interior donor-LODO with the Harmony prior** (laptop stopped at 6/48):
   ```
   DORM_PRIOR_SPACE=harmony DORM_SEEDS=1 Rscript code/run_lodo_donor_intrinsic.R
   ```
   → `results/lodo_donor_intrinsic_harmony_results.csv` (+ `_rho.csv`). ~2 h. Gives the full
   stability + benchmark table (DORM vs DORM_intrinsic vs all benchmarks, 12×4).
2. **Phase C — extrapolation, all donors/ADTs/α** (laptop did 15 combos: s3d7 full + s4d9 partial):
   ```
   DORM_DONORS=<all or subset> DORM_ADTS=CD72_1,CD3,CD19_1,CD11c \
   DORM_ALPHAS=0,0.25,0.5,0.75,1 DORM_KMEANS=8 Rscript code/run_extrapolation_shift.R
   Rscript code/fig_extrapolation.R
   ```
   → `results/extrapolation_shift_results.csv`, `figures/fig_extrap_{r2,gap}.png`. Use `DORM_DRY=1`
   for a no-fit sanity check of the states/target construction.
3. **Optional — confirm the negative result is design-robust** (only if reviewers ask): sweep
   `DORM_SIGNED_LAMBDA` (0.1/0.5/2) and/or try a fully-Harmony posterior (would require routing
   `prior_*` ratios through `posterior()`/`PQCalculation`, currently gene-space only). The
   collapse is consistent across donors/ADTs, so this is unlikely to reverse it.

## Cluster gotchas
- R 4.4.1, deps via `renv` (`renv::restore()`); Python side is the separate repo `.venv`
  (scanpy/anndata/sklearn + `harmonypy`). Only `harmonypy` was added — reinstall on the cluster.
- Runners auto-locate the repo root via `find_repo`; run from anywhere. Long jobs are safe on the
  cluster — no need for the laptop's foreground-chunking workarounds.
- Large data (`.h5ad`, `csv_for_R/*.csv`, incl. the new `cite_nips_EMBED_for_R.csv`) is
  gitignored; regenerate the embedding CSV with `make_harmony_embedding.py` if it's not synced.
