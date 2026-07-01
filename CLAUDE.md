# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

DORM (Domain adaptation Optimized for Robustness in mixture population) — an R implementation of a distributionally-robust
transfer-learning estimator that aggregates a linear model from multiple **source**
sites so it performs well on a shifted **target** distribution, without needing target
labels. `src/` holds the estimator and simulation code; `singlecell/` is the real-data
application on NeurIPS 2021 paired CITE-seq data (predict a surface protein/ADT from RNA).

R 4.4.1, dependencies pinned via `renv` (`renv.lock`). Run `renv::restore()` once to install.

## Core estimator (src/)

`get_DORM_beta()` in `src/Functions3.R` is the single entry point. Everything else in that
file is machinery it calls. Read `src/Example_new.R` first — it documents the **data contract**
that every caller (including the singlecell app) must satisfy.

Data contract (all covariate matrices have an intercept in column 1; `q+1` = dim of the
low-dimensional shared component `A`, remaining columns are `W`):
- `Xlist[[l]]` / `Ylist[[l]]` — source-l covariates and scalar response
- `Xtrainlist[[l]]` / `Ytrainlist[[l]]` — auxiliary source-l data, **must be equal size** to
  `Xlist`/`Ylist`; used for cross-fitting the nuisance (density-ratio) models
- `X0` / `X0train` — target covariates (two splits); target `Y0` is used only for evaluation
- `nlist` — labeled obs per source; `q` — `length(beta) - 1`

Algorithm flow: `get_DORM_beta` runs `maximin_s_beta` twice with the two sample splits swapped
and averages (cross-fitting). `maximin_s_beta` fits a density-ratio model per source
(`or_estimation_*` + `dratio_*`, selected by `dr_type` ∈ `rf`/`logit`/`XGB`/`nnet`), builds the
`P`/`Q` matrices (`PQCalculation`), then solves the maximin problem over the source simplex
(`opt_s_delta`, CVXR). Returns `beta_star` (the DORM estimate), plus baselines `beta_MI`,
`beta_RAP`, per-source estimates `DoublyR`, and weights `rho`.

Key `get_DORM_beta` args: `penalty` (TRUE → elastic-net for high-dim `A`, `alpha` is the
mixing param), `rho_pseudo` ∈ `max`/`pool`/`NA`, `dr_type`, `condA` (density ratio conditional
on `A`), `normalize` (truncate/normalize density ratio — keep FALSE). `smax` bounds the
adversarial shift; tune it with `src/Tuning.R` (`tuning_Y`, `tuning_surrogate`).

Baselines for comparison live in `src/TransLasso-functions.R` and `src/TransLG-functions.R`
(TransLasso, TransGLM via `glmtrans`, PTL). `get_err(X,Y,beta,q,nn)` computes MSE.

## Single-cell application (singlecell/)

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

## Simulations (simu/)

`simu/Simu_*.R` are the paper's synthetic experiments (each sources data generators from
`src/DataGeneration*.R`); `simu/graph_main.R` produces the figures. Reference these for how
the estimator behaves under controlled shift, not for the real-data workflow.

## Notes

- Large data (`.h5ad`, `csv_for_R/*.csv`) and outputs (`.rds`, `.job`, `err`) are gitignored;
  the CSV inputs the R script needs are regenerated by the export notebook.
- `get_DORM_beta` takes several minutes on real data (fits an ML nuisance model per source).
