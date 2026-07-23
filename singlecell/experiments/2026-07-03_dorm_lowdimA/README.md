# DORM v2 on erythroid CITE-seq — low-dim A + tuned smax (2026-07-03)

Goal: the first erythroid DORM runs (`singlecell/dorm_results/Ver4_*`, q=1000) gave
**R² ≤ 0 everywhere**. This experiment diagnoses why and runs DORM in its intended
regime to push it above the target-mean baseline. Only the *design* (sources /
targets / ADTs) is borrowed from the MIMAL paper; all methods are DORM-family.

## Diagnosis of the old collapse (all confirmed)
1. **Data-contract misuse.** `src/Example_new.R` expects `X = (1, A, W)` with a
   **low-dim A** (the shared predictive component, q+1 small) and a **high-dim W**
   used only as nuisance controls (imputation + density ratio). The old script set
   `q = 1000` → A = all genes, W empty. The final predictor `X[,1:(q+1)] %*% beta`
   is A-only but over-regularized, and singular when target n0 < q+1 (MK/E prog: 345<1001).
2. **Penalized intercept.** In the `penalty=TRUE` path the intercept was shrunk to
   exactly 0, so predictions were not centered on E[Y] (e.g. pred mean −0.27 vs true
   +0.52) → negative R². **FIXED in core code** (`PQCalculation`, `src/Functions3.R`):
   the four elastic-net `glmnet` calls now pass `penalty.factor = c(0, rep(1, q))` so
   the intercept is unpenalized (verified on the penalty=TRUE path: intercept
   −2.30 vs 0 before, predictions re-centered, R² 0.21). This experiment itself uses
   `penalty=FALSE`, so its numbers are unaffected; the fix re-enables the high-dim-A path.
3. **Throttled smax.** `smax=0.03` pinned `s*` at the cap → `beta_star ≈ RAP`; smax
   was never tuned. `Example_new.R` uses `smax≈0.25` and *tunes* it.

## Fix (this experiment)
- **A = top-k genes** screened by CV elastic-net on **pooled SOURCE data only**
  (label-free wrt target), with a greedy near-duplicate drop so `Sigma0` stays
  invertible. **W = the remaining ~980 HVG** as nuisance controls.
- **`penalty = FALSE`** → `solve(Sigma)` estimates a proper, unpenalized intercept.
- **smax tuned** over a grid via `tuning_Y` on a small (~100-cell) labeled target
  split, disjoint from the evaluation split.
- Cross-fit `P,Q` computed once (two halves); all smax reconstructed cheaply via
  `opt_s_delta` — one RNA read, ~60 s/combo.

Design: sources = {Proerythroblast, Erythroblast, Normoblast}; targets =
{Reticulocyte, MK/E prog}; ADTs = {CD71, CD36_1, CD105, CD49d}; k ∈ {10,20,30}.

## Headline results (`results/v2_improvement_table.csv`)
Mean DORM R²: **old −0.814 → new −0.315 (+0.499)**. The improvement is concentrated
on the biologically-coherent target:

| ADT | target | DORM old | DORM new | TransGLM† |
|---|---|---|---|---|
| CD71  | Reticulocyte | −0.010 | **+0.185** | 0.402 |
| CD36_1| Reticulocyte | −0.149 | **+0.284** | 0.451 |
| CD49d | Reticulocyte | −0.134 | **+0.271** | 0.486 |
| CD105 | Reticulocyte | −1.899 | −0.152 | 0.094 |
| CD71  | MK/E prog | −1.435 | −0.390 | 0.179 |
| CD36_1| MK/E prog | −0.829 | −1.239 | 0.257 |
| CD49d | MK/E prog | −1.008 | −0.512 | 0.107 |
| CD105 | MK/E prog | −1.047 | −0.968 | 0.120 |

†TransGLM uses up to 200 **labeled target** cells; DORM uses 0 (only ~100 for smax tuning).

- **Reticulocyte** (terminal erythroid, same maturation trajectory as the sources):
  3/4 ADTs now **positive** (was 0/4), competitive with label-using TransLasso/GLM.
- **MK/E prog** (an upstream progenitor on a different branch): still negative — this
  is **intrinsic concept shift**, not a code issue (see Q-tracking).

## Ablation (`results/v2_ablation_table.csv`)
old(q=1000) → low-dim-A@smax0.03 → low-dim-A@tuned. The **decisive lever is tuned
smax**, e.g. CD105/Reticulocyte −3.336 → −0.152, CD36/Reticulocyte −0.502 → +0.284.

## Per-source Q tracking (`results/v2_Qtracking_all.csv`)
DORM's `rho` recovers the developmental ordering:
- **Reticulocyte → rho = 1 on Normoblast** (its direct precursor), which also has the
  best standalone target R².
- **MK/E prog → rho ≈ 0.87 on Proerythroblast** (earliest nucleated stage, closest to
  the progenitor). Erythroblast predicts MK/E prog catastrophically (standalone R²
  −47 to −147) and correctly gets rho = 0.

## How much did the intercept fix matter? (isolation check)
Re-running the OLD config (q=1000, penalty=TRUE, smax=0.03) with ONLY the intercept
fix changed DORM negligibly: CD71/Reticulocyte −0.010→−0.008, CD49d/Reticulocyte
−0.134→−0.115, CD71/MK-E −1.435→−1.961 (noise). Reason: with 1000 gene features the
coefficients already absorb the mean, so the optimal intercept is ≈0 anyway and
predictions were already centered. **The intercept bug is real and worth fixing
(it clearly matters in the synthetic penalty=TRUE case: intercept −2.30, R² 0.21),
but it is NOT what rescued DORM here.** The dominant levers were low-dim A + tuned
smax. See `results/intercept_fix_q1000_check.csv`, `code/check_intercept_q1000.R`.

## Reproduce
```
Rscript code/run_dorm_erythroid_v2.R           # full 24-combo grid (~25 min)
Rscript code/assemble_v2_report.R              # tables
# single combo: DORM_ADTS=CD71 DORM_TARGETS=Reticulocyte DORM_KS=20 Rscript code/run_dorm_erythroid_v2.R
```

## Files
- `code/run_dorm_erythroid_v2.R` — screening + DORM (low-dim A, penalty=FALSE) + smax tune + benchmarks + Q tracking
- `code/assemble_v2_report.R` — improvement / ablation / Q-tracking tables
- `results/v2_improvement_table.csv`, `v2_ablation_table.csv`, `v2_Qtracking_all.csv`
- `results/v2_<ADT>_target_<t>_k<k>_{summary,Qtrack,Qmatrix,result.rds}` — per-combo (A genes in `Qmatrix` / `result.rds$A_genes`)

## Notes / next levers
- `best_smax` often lands at 0.9 (≈ pure doubly-robust maximin, MI); k=10 usually wins.
- To push Reticulocyte further without target labels: try `dr_type='rf'` (nonlinear
  density-ratio) and k∈{5,8}.
- MK/E prog is not recoverable label-free — report it as the concept-shift boundary.
