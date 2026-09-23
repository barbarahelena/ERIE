# Problems and fixes

This document records the problems found during model development and how each was resolved. It is the history behind the current design. The current behavior is described in `docs/pk-model.md`. Entries are grouped by topic and carry the commit that made the change where one is identifiable (`git show <hash>`). Subject-visit cases are written as subject ID plus visit (`baseline` or `intervention`).

## Model structure

### A one-compartment absorption curve cannot produce a post-peak dip (02eaaf9, 2026-09-18)

Many curves show a dip after an initial peak followed by a second, often higher, rise, most clearly on ER01 and ER09. A one-compartment absorption curve is monotonic after its single peak for any `ka`, `kel` and `F`. A "biphasic gastric emptying" variant with two gastric pools that both start at t=0 produced no second hump in thousands of random parameter draws and several hand-engineered attempts, because pathways that all start at t=0 and share one elimination rate converge to the same single-exponential tail (`simulate_biphasic_emptying()` remains in `pk_curves.R` as a record of this variant). The fix is a time-lagged second dose (`simulate_lagged_dose()`): a fraction `f_delayed` of the dose starts `t_lag` minutes late. A pilot on ER01, ER03 and ER09 changed the two dip cases from one smooth hump to a double-humped fit that follows the dip.

### Model names and selection (95c5998, 2026-09-18)

The two candidates were called "baseline" and "lagged" and were reported side by side. They are now `single_wave` and `two_wave`, and one is selected per subject x visit by AIC on 12C's residuals. Both candidates' fits are still written out.

### A lag informed by both curves misattributes 13C6's slow onset to a second wave

The delay of the second wave (`f_delayed`, `t_lag`) was first estimated jointly from both curves' residuals. 13C6 has its own slow onset from the capsule dissolution (`k_release`), and that onset was read as a second wave and imposed on 12C. ER32 baseline has a 12C t=30 sample already near its peak, so no onset delay is plausible, yet the joint fit pinned `f_delayed` near 1 and `t_lag` near 19 min and lowered both curves' R² below the `single_wave` fit. ER35 showed the same pattern. The lag is now estimated from 12C alone (stage 1), and 13C6 is fitted afterward with `ka`, `kel` and `t_lag` fixed (stage 2).

### 13C6 inherited 12C's delayed fraction (4be99cf, 2026-09-18)

Forcing 12C's `f_delayed` onto 13C6 failed on ER03 baseline. 12C has two waves (R² = 0.99), while 13C6 has one early peak that decays monotonically, and the inherited second wave gave `r2_13C6` = 0.075. 13C6 now has its own `f_delayed_13C6`, fitted in stage 2, and shares only `t_lag`, `ka` and `kel` with 12C.

### `two_wave` on every curve, and a gate that only saw dips (4be99cf, 0afa4fc, fa83913)

Attempting `two_wave` on every curve leaves AIC to discover afterwards that nothing needed explaining, so the attempt requires evidence in the raw 12C data. The first gate looked for a post-peak dip and missed plateaus. ER02 baseline (flat top, R² 0.955 to 0.997), ER06 intervention (peak neighbor at 86% of the peak, R² 0.981 to 0.9985) and ER30 baseline (irregular pre-peak rise, R² 0.859 to 0.957) all improved with `two_wave` and had no dip. `has_near_peak_neighbor()` and a dip check based on the chord between neighbors now cover these cases. t=30 is excluded as a near neighbor because closeness to the peak at the first sample reflects fast absorption: ER06 intervention has t=30 as its only near-peak point, and its `two_wave` fit pinned `t_lag` at the lower bound.

A rule that skipped `two_wave` when `single_wave` already had R² above 0.95 was dropped (fa83913). It was not the rule blocking the missed cases: ER30 baseline had a `single_wave` R² of 0.846 and was blocked by the dip-only gate.

### A smoothly accelerating rise with no dip and no plateau (ER21 baseline, 477a932)

ER21 baseline has 12C at t=30/60/90 of 32/61/132 mg/L, a `single_wave` R² of 0.31, and no signal for either detector. An unconstrained `two_wave` search found R² = 0.99 at `t_lag` = 57 min and `f_delayed` = 0.80. A `single_wave` R² below `R2_ALWAYS_TRY_TWO_WAVE` (0.50) now also triggers a `two_wave` attempt.

### A marginal AIC edge overruling clear evidence (e9c553c, 2026-09-19)

An unlucky local optimum can give `single_wave` a small AIC advantage on a curve with an unmistakable second wave. A dip excess of 0.90 or more (`EXCESS_DEFINITE_TWO_WAVE`) now selects `two_wave` without the AIC comparison. Only ER08 baseline (1.57) and ER11 baseline (1.01) exceed it, and both are confirmed by independent evidence.

### AICc tried and reverted

With 8 to 9 points per curve and `k` rising from 3 to 5, the AIC penalty looked too weak once `two_wave` was attempted more broadly. AICc (`+2k(k+1)/(n-k-1)`) raises the complexity gap between k=3 and k=5 from 4 to 28 at n=8, and it excluded all three confirmed plateau cases (ER02, ER06 intervention, ER30) despite independently verified R² gains. Plain AIC is used, and the raw-data gates and parameter bounds provide the protection against overfitting.

### `t_lag` at the edges of the timeline (4be99cf, 0afa4fc)

Upper end: sampling gaps widen to 60 min (180 to 240) and 120 min (240 to 360), and a `t_lag` in one of them places a second wave where no sample can check it. ER12 intervention reached `t_lag` = 179 min and `r2_13C6` fell from 0.95 (`single_wave`) to 0.75 with a large peak and no data around t=200 to 220. The upper bound is 150 min.

Lower end: a flat lower bound of 5 min let the optimizer place an early onset delay in the unsampled 0 to 30 min window. ER06 intervention reached `t_lag` = 19.7 min and `f_delayed` = 0.999 although its t=30 sample (about 80 mg/L) was already near its peak. The first fix interpolated the floor between 5 and 30 min as t=30 rose between two thresholds (`EARLY_LAG_OK_MGL` = 5 mg/L and `EARLY_LAG_BAD_MGL` = 75 mg/L). ER04 baseline (t=30 = 51.4 mg/L) then received a floor of 21.6 min, and the fit landed exactly there with a flat-then-rise artifact, since every value between 5 and 30 min lies in the same unsampled window. The floor is now a step (5 min when t=30 is at or below 5 mg/L, otherwise 30 min) and `EARLY_LAG_BAD_MGL` is removed. ER06 intervention then converged to `t_lag` = 35.3 min and `f_delayed` = 0.320.

### The search missed the evidence that justified `two_wave` (ER04 baseline, 0afa4fc)

ER04 baseline passed the gate through a plateau at t=150, but the generic seeds converged to `t_lag` = 21.6 (an early near-total-delay solution with a marginally better AIC). Seeds anchored on the trigger time of the dip check (at it, 20 and 40 min before it) now cover the evidence.

### `MIN_TMAX` constrained only the taller wave

The combined-curve check on the peak time constrains whichever wave is taller. ER09 baseline and intervention have a taller wave 1 (combined peak at t=38 and t=33) and a wave 2 that peaks alone at t=151 and t=147, and ER16 intervention has a wave 2 peaking at t=162. `MIN_TMAX` now applies to each wave's own peak (`tmax_wave1`, `tmax_wave2`).

### 13C6 needed its own onset lag (84a6075, 2026-09-19)

The first-order dissolution of delayed release is fastest at t=0 and cannot represent a 13C6 signal that stays near zero for about 30 min and then rises sharply. The optimizer compromised by releasing some dose early. A separate onset-delay candidate (`t_lag1_13C6`) is now available, attempted when the first observed point is below 25% of the curve's peak. Validation: ER25 baseline (R² 0.71 to 0.94, 25 min), ER06 intervention (0.64 to 0.94, 28 min) and ER25 intervention (R² 0.75, 27 min).

### `k_release` in the onset-lag candidate sat at its bound (ae8eded, 2026-09-19)

With a delay in the model, `k_release` reached its upper bound in all three validation cases (ER06 intervention, ER25 baseline, ER25 intervention), because 30-min sampling cannot separate fast dissolution from instant release once the delay explains the flat start. The onset-lag candidate is now an instant bolus after the delay (K=2 instead of 3). The same commit fixed a plotting bug for this mechanism.

### 13C6 needed its own second wave, with tighter bounds (b6170be, 2026-09-20)

13C6 is noisier than 12C, and the detectors flag 47/68 of its curves. ER21 baseline (dip excess 3.416, the strongest in the cohort) reached R² 0.508 to 0.741 only by pinning `f_delayed2_13C6` at 0.999, driven by a near-zero crash and rebound (0.0117 mg/L at t=150, 0.0354 at t=180) that looks like measurement noise. `f_delayed2_13C6` is bounded to [0.05, 0.95] and `t_lag2_13C6` to [30, 150]. ER11 baseline (R² 0.787 to 0.846, `f_delayed2_13C6` = 0.751) and ER18 intervention (0.860 to 0.879, 0.095) are accepted fits away from the extremes.

### 13C6's mechanisms were unavailable under `two_wave` (477a932)

With 12C selected as `two_wave`, 13C6 could only use delayed release with 12C's `t_lag`. ER23 intervention has 13C6 at 0.009 mg/L at t=30 and 0.165 at t=60 (5.6% of its peak), and delayed release collapsed to `f_delayed_13C6` = 0.001 and `r2_13C6` = 0.35. All three 13C6 candidates are now compared by AIC under both models, and ER23 intervention picks up an onset lag of about 27 min with `r2_13C6` = 0.78.

### A second-wave result at the lower bound (477a932)

ER18 intervention won its AIC comparison at `f_delayed2_13C6` = 0.05 and `t_lag2_13C6` = 150, both at their bounds, although its dip evidence lies at t=120. A result within 1e-3 of the lower bound of `f_delayed2_13C6` is now rejected without an AIC comparison. ER18 intervention falls back to delayed release with `r2_13C6` = 0.66 (0.71 with the spurious second wave).

## Fitting

### Bound-width `parscale` hid optima (7b84cff, 935b9d7)

`fit_multistart()` first used `parscale = upper - lower`. For `ka` and `F_12C`, bounded to [1e-4, 1] with fitted values around 0.01 to 0.05, that scale is 20 to 100 times too large, and L-BFGS-B stayed near its start. ER11 baseline scored 0.208 at `t_lag` = 47.3 under this scaling, while a search scaled by each start's own magnitude reached 0.188 at `t_lag` = 120, which matches the dip evidence (`trigger_time` = 150). `parscale` is now set per start from that start's own magnitude, floored at 1% of the bound range (935b9d7, 2026-09-20).

### Proportional weighting failed on curves with an early feature (a8ad159, 02eaaf9)

Residuals of an unweighted fit showed heteroscedasticity (the squared-residual scale was about 24 times higher in the top than in the bottom quartile of predicted concentration, for both curves), so a proportional weighting `1 / max(pred, floor)^2` was added (a8ad159, 2026-09-17). It discounts the region where an informative early feature lives. The 13C6 curve of ER03 baseline (a spike at t=30 that decays by t=90) gave R² = -2.18 weighted, and loosening `MIN_TMAX` from 30 to 10 moved it only from -0.86 to -0.67. The unweighted fit with the same `MIN_TMAX`, bounds and an ordinary multistart gave R² = 0.37, and the former model (unweighted) reached 0.71 on the same subject. The fits use the unweighted objective (02eaaf9, 2026-09-18), and the weighting stays available as an option.

### A rough penalty surface (5830762, 2026-09-17)

The Tmax penalty had a hard boundary, which makes the surface rough for L-BFGS-B's finite-difference gradients. It is now a squared shortfall that is zero at and below the threshold.

### Fits stuck at a bound or in a worse optimum (ffe3f6b, eaaf059, d0b93dd, 2026-09-17)

The former joint model had a subject stuck at the `kel` bound because the standard search kept landing on a worse local optimum on the boundary, while a better solution existed inside it. An adaptive retry pass (ffe3f6b) re-searches flagged fits with a denser search and keeps the result only when its objective is strictly lower. Boundary flags alone miss a worse optimum inside the bounds, so the triggers include `r2_12C_low` and `converged = FALSE`. `r2_13C6_low` alone is excluded because a poor 13C6 fit usually reflects a structural mismatch that a denser search cannot fix. `retried` and `retry_improved` columns record what happened (eaaf059), and progress is printed per item (d0b93dd).

### Missed optimum when `ka` and `kel` are close (ER06 intervention, 0afa4fc)

ER06 intervention (`single_wave`) was reported at `ka` = 0.016, `kel` = 0.051, R² = 0.945, while multiple `optim()` starts converging to one point found `ka` = 0.023, `kel` = 0.025, R² = 0.986. This is the "flip-flop" region of one-compartment models, where the `ka/(ka-kel)` term makes the surface narrow. Explicit `ka`=`kel` diagonal seeds are added in the independent pre-fit, the joint `single_wave` fit and `two_wave` stage 1. Cohort-wide R² values have not been re-checked without these seeds.

### Reproducibility under `mclapply` (353ed45)

A single `set.seed()` before `parallel::mclapply()` gave different random draws on separate runs of identical code, and `RNGkind("L'Ecuyer-CMRG")` fixes reproducibility only for a fixed task order and core count. Each subject x visit seeds itself from its own ID (`string_seed()`), and a regression test checks that the seeds are deterministic and distinct.

## Volume of distribution and dose

Vd was first a flat 0.15 L/kg (an extracellular fluid volume estimate from glucose literature), then the individual Nadler blood volume from weight, height and sex, then the individual extracellular fluid volume (Watson total body water / 3) from weight, height, age and sex (29ebdf0), and is now 0.4 x the measured total body water of each visit (TBW / 2.5). Blood volume (about 65 to 75 mL/kg) is much smaller than the ECFV (about 150 to 200 mL/kg), so the move from the flat 0.15 L/kg to Nadler made `F_12C` and `F_13C6` smaller. In an earlier run under the Nadler Vd the median `F_12C` was about 3.2% and the median `F_13C6` about 3.0%; the last Nadler run gave medians of about 1.9%. The move to the Watson ECF goes the other way: the mean ECF is about 2.7 times the mean Nadler blood volume in this cohort, and `F` scales with `Vd`, so `F` from the different versions cannot be compared with each other. `ka` and `kel` are unaffected by any of these changes. The argument for ECFV (Vd of glucose equals the extracellular space, van der Crabben et al. 2011) was not refuted by the Nadler version, and the choice of the ECF has not been checked against a fructose-specific Vd. The move to the measured TBW (2026-09-22) followed a comparison of Watson's TBW with the measurements in `data/TBW_ERIE.csv`: r = 0.90 over 66 visits, measured mean 43.9 L against 42.6 L, and individual differences of about 4 L (SD). Measured TBW / 2.5 averages 17.6 L against 14.2 L for Watson TBW / 3, so the median `F_12C` went from 4.7% to 5.9% and the median `F_13C6` from 4.7% to 5.8%. The selected model stayed the same for all 68 fits and the median R² barely moved (12C 0.975 in both runs, 13C6 0.929 to 0.927). `ka` and `kel` do not depend on Vd in principle, yet 13 fits ended at different values, four of them (ER14, ER15, ER19 and ER27 baseline) by swapping `ka` and `kel` at the same R², the flip-flop ambiguity described above, and ER12 intervention's 12C R² fell from 0.83 to 0.76. The measured TBW changes more between the two visits of a subject than the Watson estimate does (SD of the change 2.2 L against 0.5 L), and ER033's baseline value (57.1 L) is implausibly high.

The 13C6 dose is 120 mg / 644.78 µmol. The raw constants files implied about 100 mg, and the administered dose was confirmed as 120 mg (2026-09-17), which made the molar amount in the files wrong as well. A larger assumed dose gives a somewhat smaller `F_13C6` and leaves `ka` and `kel` unchanged.

## Diet summary

- **`converged` in the reliability gate.** It only records whether L-BFGS-B reached its stopping criterion before `maxit`. In the first check, `r2_12C` for `converged = FALSE` rows (mean 0.902, n=22) matched or exceeded `converged = TRUE` rows (mean 0.890, n=46), so the gate uses 12C's `kel` bound flag and R² only. The current run shows the same picture (mean 0.957 in both groups).
- **Blank 13C6 mean curves (cb5b0fb, 2026-09-20).** The 13C6 onset-lag and second-wave mechanisms set `k_release` to NA, and passing NA to `simulate_delayed_release()` returns NA for every t>0. One such subject blanked a whole diet x visit group's mean curve, and 14 of 21 `single_wave` rows had `k_release` = NA. The mean curves now dispatch on the mechanism exactly as `simulate_fit()` does.
- **Delta plot test (c9035fe, d74334e).** A Welch t-test for the between-arm comparison was tried and reverted to the Wilcoxon rank-sum test, which suits the small n per arm.
- **Mixed models (76ad4f5, aa84538, bf4d3e7, ecab4bd, 2026-09-20).** Log and raw scale were run side by side, and log was better behaved. Sex and BMI were added as covariates, and BMI was dropped. The final model is log scale with sex.
- **Inclusion cutoff (bd4dbcd).** `R2_INCLUDE_MIN` was raised to 0.85 and is now 0.90.
- **Capsule half-life (4614c09).** `capsule_dissolution_halflife_min` was excluded from the diet summary because `k_release` exists for 9 of 68 fits in the run at that time (14 of 68 now). The column was later removed from `fit_results.csv`.
- **Paired t-test block from #22 (1affb40).** It read `reliable_shared`, `reliable_12C` and `reliable_13C6` columns that no longer exist and included the capsule half-life. The paired Wilcoxon test and the mixed models cover the same question.
- **Layout (ea0c8c8, ee42666).** Curves are faceted by diet arm and colored by visit, and parameter and before/after boxplots have one PDF per parameter.
- **LMMs and the delta plot removed (2026-09-22).** The mixed models above were briefly extended to add AUC, Cmax and Tmax (per isotope) alongside `ka`/`kel`/`F_12C`/`F_13C6`, and age and body weight alongside sex as covariates. The extended LMMs and the delta plot (between-arm Wilcoxon rank-sum on each subject's log fold-change) were then removed, judged unnecessary next to the curve-derived metrics section, which has its own paired Wilcoxon test on AUC/Cmax/Tmax. `library(lmerTest)` and the `r-lme4`/`r-lmertest` pixi dependencies were dropped along with the LMMs. `R2_INCLUDE_MIN` is back to 0.90 (see "Inclusion cutoff" above; it had drifted to 0.75 in between).
- **Parameter and curve-metric boxplots unified (2026-09-22).** The separate `plot_before_after()` boxplot (its own PDF per parameter, `diet_before_after_boxplot_*.pdf`) is gone; its paired Wilcoxon test (baseline vs intervention within each diet arm) is now annotated directly on the plain parameter boxplot instead, one PDF per parameter as before. This matches the curve-derived metrics boxplots exactly: same annotation style (`geom_text`, not `stat_pvalue_manual`), same plain black jitter (previously colored by visit/diet, which blended into the box fill) on every boxplot in the script. The two-wave characteristics plot keeps its unpaired rank-sum test, since it compares different subjects across diet arms and cannot be paired. All of this script's outputs (previously loose in `results/`) now go in `results/diet_summary/`, alongside `plots_individual/` from `02_fit_erie_model.R`.
- **Vd plot split by visit (2026-09-22).** It previously showed one Vd per subject (mean over visits) by diet arm only. It now shows both visits, dodged and filled by visit within each diet arm (matching the R2 plot's dodge pattern), jitter still shaped by sex. The title was also shortened from "Extracellular fluid volume (Vd = 0.4 x measured TBW) by diet arm" to "Extracellular fluid volume".

## Pipeline and data

- **Visit labels (d7b46f8, 2026-09-17).** `visit` is recoded from `FCT1`/`FCT2` to `baseline`/`intervention` in `01_clean_data.R`, and every downstream script and doc uses the coded values. The visit label is part of each fit's random seed.
- **Raw label typo.** Raw column 34 of `ERIE_fructose_13C.csv` had a label typo (`"FCT -  34"`, missing the "1"), which the position fallback for the visit recovered from. The raw file has since been corrected.
- **Output names.** The selected-model table is `fit_results.csv` (formerly `fit_results_joint.csv`), and `04_compare_to_former_model.R` reads it (0a3df8a).
- **Data findings.** A decimal-mark mismatch between raw files (comma in some, period in others) and an undocumented column-to-subject mapping (reverse-engineered from `former_models/`) were found and handled in `01_clean_data.R`; see the 13C6 dose discrepancy above for the same kind of issue.

### R² compared with the former model (`04_compare_to_former_model.R`, removed 2026-09-22)

A one-off script compared each subject x visit's classical R² with the former model's saved results (`former_models/MixedModel/Results/fit_results_joint.csv`), to check that the current joint model was actually an improvement before relying on it. In its last run the median 12C R² was 0.975 against the former model's 0.933 (higher for 57 of 68 subject x visits, lower for 11), and the median 13C6 R² was 0.927 against 0.793 (higher for 58, lower for 10); the correlation between former and current R² across subject x visits was modest (0.36 for 12C, 0.38 for 13C6). R² was a fair yardstick between the two: both objectives are unweighted and normalized by each curve's own total variance. The comparison script and its outputs (`r2_comparison_*.csv/pdf`) were removed once this was established.

## Documentation corrections

- The dip prevalence figures of an earlier version (13/68 dips, about 31% of curves with a smaller-amplitude version, 85% of dips at t=60 or t=90) came from a stricter trough-then-rise check. They are replaced by counts from the current detectors (see `docs/pk-model.md`).
- A statement that 45 of 68 13C6 fits fall below R² 0.70 described an earlier model. The current run has 4 of 68.
- A section explaining a lower R² than the former model by the proportional weighting described the earlier weighted objective. The current R² is higher for 57 of 68 (12C) and 58 of 68 (13C6) subject x visits.
- `t_lag`'s lower bound was documented as an interpolation with `EARLY_LAG_BAD_MGL`, and the `ka` and `kel` issue as an open TODO. The step function and the diagonal seeds are documented instead.
