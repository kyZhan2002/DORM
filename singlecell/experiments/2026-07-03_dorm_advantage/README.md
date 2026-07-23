# Does DORM show an advantage on CITE-seq? (2026-07-03)

Goal: design a real-data analysis matching DORM's actual guarantee — label-free predictive
robustness under target distribution shift — and test whether DORM beats its fair peer group
(label-free methods, esp. naive pooling). Two arms: erythroid trajectory (mixture stress test +
leave-one-stage-out) and cross-donor within CD14+ Mono. Uses the fixed `src/Functions3.R`
(unpenalized intercept + ridge-stabilised normal equations) and the low-dim-A / tuned-smax pipeline.

## Headline verdict (honest)
On this data, **DORM does not deliver a blanket worst-case win over naive pooling**, because the
clean covariate-shift regime is benign and the strong-shift regimes are concept-shift / weak-signal.
BUT there is a real, theory-consistent advantage in the regime transfer learning actually targets:

**DORM (and its density-ratio-reweighted variant RAP) beat pooling precisely when the target is
far from the pooled training distribution, and the advantage grows with that distance
(cor ≈ 0.6).** Far targets: mean gap +0.03 to +0.15; near targets: pooling suffices. This is
diagnosable from covariates alone (no labels). Figure: `figures/fig4_advantage_vs_distance.png`.

Second real advantage: **DORM is the stable/safe aggregation.** DORM (tuned smax) is positive on
every mixture, while the naive maximin (MI), simple-average, best-single, and rho-average all fail
catastrophically (worst-case −2 to −7e5). DORM's smax tuning correctly avoids MI's blow-ups.

## Per-analysis results
1. **Phase 1 — mixture stress test (erythroid, the clean covariate-shift test).** Worst-case over
   mixtures: RAP 0.291 ≈ Pooled EN 0.278 ≈ DORM(tuned) 0.156 (CD71); DORM wins on minority/off-
   dominant targets (vNormo 0.44 vs 0.32, vProery 0.56 vs 0.51, eEN 0.43 vs 0.32), loses on
   Erythroblast-dominated targets (pooling is already Ery-dominated by sample size). Advantage
   grows with target distance from the pooled centroid (fig4). MI catastrophic on balanced mixtures.
2. **Phase 2 — leave-one-stage-out (erythroid).** Concept shift: ALL label-free methods fail
   (worst-case ≤ −1.4; target-mean null ≈ 0 wins). Only label-using methods (oracle 0.42, TransGLM
   0.34) survive. Not an informative regime for label-free transfer.
3. **Phase 2 — cross-donor (CD14+ Mono, CD11c).** Weak signal: even the label-using oracle gets
   mean R² 0.19. Uninformative; DORM-family baselines (MI/simple-avg) blow up on some donors.

## Why DORM doesn't dominate here
DORM's edge over pooling needs source heterogeneity that biases pooling, in a covariate-shift
(shared Y|X) regime with real signal. This data gives either (a) benign shared-Y|X mixtures where
pooling is near-optimal, or (b) concept-shift / weak-signal settings where nothing label-free works.
The narrow window (moderate heterogeneity + strong signal) isn't cleanly hit by the natural splits.

## Two bugs fixed along the way (in core `src/Functions3.R`)
- Unpenalized intercept in `PQCalculation` penalty path (`penalty.factor = c(0, rep(1,q))`).
- Ridge-stabilised normal equations in the penalty=FALSE path (rank-deficient targets).
- Analysis switched to `penalty=TRUE` (regularised, bounded Q) after a single pure-subtype target
  produced a −8766 blow-up under the unregularised path — the worst-case metric surfaced it.

## Files
- `code/lib_dorm_common.R` — shared pipeline + baselines (pooled-EN, target-oracle-EN, two tiers).
- `code/run_erythroid_mixture.R`, `code/run_leave_one_domain_out.R` — the runs.
- `code/assemble_advantage_report.R`, `code/advantage_vs_distance.R` — tables + figures.
- `results/summary_*.csv`, `results/advantage_vs_distance.csv`; `figures/fig1..fig4`.

## Leave-one-donor-out on erythroid (added later; the best-motivated design)
Sources = 3 nucleated erythroid stages pooled over 11 donors; target = held-out donor's
{Proery,Ery,Normo} cells (natural mixture of the source types = covariate/donor shift, no
concept shift). q=20; penalty=FALSE BLEW UP (RAP -58) so switched to penalty=TRUE (stable).
Result (`fig7_lodo_donor_ery.png`, `summary_lodo_donor_ery.csv`): **all label-free methods fail
(DORM mean R^2 -4.8, no better than Pooled -2.3); only the label-using oracle (mean 0.28) and
TransGLM (0.11) are positive.** Even the oracle is modest -> weak cross-donor RNA->ADT signal.
DORM extrapolates catastrophically on some donors (s1d1 CD71 = -29.6 at normal ADT variance =
real failure; CD49d's -22M values are partly an R^2/near-zero-variance artifact). Conclusion:
label-free cross-donor ADT transfer does not work on this data; DORM = pooling here.

## Overall verdict across designs
DORM's only real win is the shifted-MIXTURE-of-related-sources regime (fig1/fig4): it approaches
the labeled oracle and beats pooling as target shift grows. On a genuinely new cell type (fig6)
or new donor (fig7), DORM ~ pooling and both fail; only label-using methods (weakly) transfer.

## Recommended next step
The natural splits are benign/adverse; to demonstrate DORM's advantage cleanly, run it in the
regime it is built for: **imbalanced, heterogeneous sources with a minority-dominated target**
(the realistic "abundant common populations, deploy on a rare one" scenario), and/or lean on the
controlled `simu/` study. See conversation for options.
