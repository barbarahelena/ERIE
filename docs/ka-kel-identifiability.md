# Why ka and kel cannot be separated, and what we tested

This explains why the absorption rate `ka` and the elimination rate `kel` are hard to tell apart in the ERIE fructose curves, which tests were run to show it, and what can still be said about a diet effect on the rates. The model itself is in `docs/pk-model.md`.

## Summary

With one plasma curve per visit, `ka` and `kel` are interchangeable: swapping them gives the same curve once `F` is rescaled. A fit can report which rate is the fast one and which is the slow one. It cannot report which of the two is the absorption rate. The data support statements about "the fast rate" and "the slow rate", and a statement about `ka` alone needs information from outside the curve.

In the single-compartment fits, 44 of 68 curves end with `ka` and `kel` equal to within 0.1%. In simulations of this study's design, a `kel` change shows up as a `ka` change in 74% of the replicates, and the fitted labels are swapped in about half of all fits.

## Why the two rates are interchangeable

The model sends the dose through two first-order steps in series. The dose leaves the gut at rate `ka`, then leaves the blood at rate `kel`. The plasma curve is the dose passed through both delays.

The result of two first-order delays in series does not depend on their order. A drug that waits in a compartment with half-life A and then in one with half-life B produces the same output curve as one that waits in B and then in A. The sampled compartment is the second one, so a rapid step followed by a slow step looks the same as a slow step followed by a rapid one, apart from the height of the curve.

In formulas, the Bateman curve satisfies this identity:

```
C(t; ka, kel, F) = C(t; kel, ka, F * ka / kel)
```

Swapping the two rates and multiplying `F` by `ka / kel` gives exactly the same concentrations at every time. With `ka = 0.04` and `kel = 0.025` per minute, the swapped solution has `ka = 0.025`, `kel = 0.04` and an `F` 1.6 times larger. Typical fitted `F` values are 0.05 to 0.09 (`results/fructose_independent/diet_summary/diet_parameter_summary.csv`), so the swapped `F` stays far below its upper bound of 1.

A numerical check with illustrative values, `ka = 0.0251`, `kel = 0.0124`, `F = 0.08`, a dose of 80000 mg and `Vd` = 16 L, at the eight sampling times 30 to 360 min, gives concentrations of 172.7, 200.3, 176.4, 139.6, 104.8, 76.2, 38.4 and 9.0 mg/L. The swapped parameters (`ka = 0.0124`, `kel = 0.0251`, `F = 0.162`) give the same eight values, with a largest difference of 6e-14 mg/L, which is floating-point noise. The `two_wave` model, which reuses `ka`, `kel` and `F` in both waves, has the same symmetry and the same difference for a second wave of 50% of the dose at 90 min.

Identical predictions mean identical residuals, R² and objective value. The Tmax and Cmax penalties in the objective depend only on the predicted curve, so they are identical as well.

`F` is free in every fit, so nothing in the objective prefers one of the two solutions. `F` and `Vd` also enter only as `F * D / Vd`, so a known `Vd` does not break the symmetry either.

The parameter bounds can make one of the two solutions unreachable. `kel` is limited to [0.005, 0.1], `ka` to [1e-4, 1] and `F` to [1e-4, 1]. ER26 intervention has `ka` = 0.31 and `kel` = 0.010, and its mirror would need `kel` = 0.31, above the ceiling of 0.1. In that case the bounds exclude the swapped solution and the data have no say. When both versions lie inside the bounds, as for ER17 baseline (0.025 / 0.012) and ER24 baseline (0.021 / 0.042), the data cannot prefer either, and the one the optimizer returns depends on its starting points.

Joining 13C6 to the fit does not help. The 13C6 capsule adds a third first-order step (`k_release`), and a chain of first-order steps is symmetric in all of its rates.

## What the curve does show

A difference of two exponentials has a fast part and a slow part. The rising edge and the peak are shaped by both rates, and the tail after the peak is set by the slower one. The curve therefore identifies `k_fast = max(ka, kel)` and `k_slow = min(ka, kel)` without ambiguity, and nothing in the curve says which is which.

If the true `ka` is the slower rate (flip-flop kinetics, where absorption is the slow step), the tail carries `ka` and the rise carries `kel`. Gastric emptying of a 1 g/kg liquid load is slow, so this case is plausible for fructose. The labels in `fit_results.csv` therefore say little about physiology.

## Why the data make it worse

When `ka` and `kel` are close, the curve becomes `t * exp(-k * t)`, a single-rate curve in which the two rates have merged. Of the 68 fits, 44 end within 0.1% of this diagonal: 17 of 22 `single_wave` fits and 27 of 46 `two_wave` fits.

This is where the data put the optimum. A grid search over `ka` and `kel` (with `F` solved exactly) puts the minimum within 10% of the diagonal for 17 of the 22 `single_wave` fits, and the best diagonal point has an R² within 0.0007 of the grid minimum. The other five have their grid minimum off the diagonal: ER17 baseline, ER24 baseline and ER27 baseline (mirrored pairs such as `ka` 0.026 / `kel` 0.012, with R² gains of 0.001 to 0.008 over the diagonal), ER26 intervention (R² 0.961 against 0.859 on the diagonal) and ER12 intervention (grid minimum at `kel` = 0.3, outside the allowed range).

The sampling grid adds to this. Each curve has about 8 samples (0, 30, 60, 90, 120, 150, 180, 240 and 360 min), the first after baseline is at 30 min, and the residual noise is about 4.3 mg/L (median over the real fits) on peaks of roughly 60 to 100 mg/L. At that resolution a curve with `ka = 0.04, kel = 0.025` and one with `ka = kel = 0.03` differ by less than the noise. The curve is broad and symmetric, and only its overall width and tail are measured.

## Why `ka` looks noisier than `kel`

The bounds are asymmetric: `ka` is limited to [1e-4, 1] and `kel` to [0.005, 0.1]. Where the optimizer leaves the diagonal, the mirror image needs one rate to rise. `ka` has room to rise (12 of 68 fits have `ka` above 0.1, up to 0.33, and all of them have `kel` at 0.008 to 0.023), while `kel` hits its ceiling at 0.1. The scatter lands in `ka` and `kel` stays near 0.02. A `ka` effect in a diet comparison can therefore come from which fits slid off the diagonal.

The group means show this. Mean `ka` is 0.058 against a median of 0.026, with SDs as large as the means.

## Tests that were run

All scripts run from the project root. Each model setup has its own folder under `results/`: `fructose_independent/` holds the independent fits (scripts `02_fit_erie_model_fructose_independent.R` and `03_diet_summary_fructose_independent.R`), `fructose_shared_kel/` the shared-`kel` fit, and `diagnostics_ka_kel/` the tests of the independent-fit model (scripts 04 to 07).

| Question | Script | Output |
|---|---|---|
| How many fits sit on the `ka = kel` diagonal, and does it affect fit quality? | inline summary of `fit_results.csv` | counts above |
| Is the diagonal a true minimum or an optimizer artifact? | `scripts/diagnostics_ka_kel/04_profile_ka_kel.R` | `diagnostics_ka_kel/profile_ka_kel.csv` |
| Does a diet effect show in quantities the curve identifies? | `scripts/diagnostics_ka_kel/05_identified_quantities.R` | `diagnostics_ka_kel/identified_quantities/` |
| Could this design tell a `ka` change from a `kel` change? | `scripts/diagnostics_ka_kel/06_simulate_ka_kel_recovery.R` | `diagnostics_ka_kel/simulation/simulation_summary.csv` |
| What does assuming a stable `kel` buy, and what does it cost? | `scripts/diagnostics_ka_kel/07_simulate_shared_kel.R` | `diagnostics_ka_kel/simulation/simulation_shared_kel_summary.csv` |
| What does the shared-`kel` fit say about the real cohort? (diet summary: `03_diet_summary_fructose_shared_kel.R`) | `scripts/fructose_shared_kel/02_fit_erie_model_fructose_shared_kel.R` | `fructose_shared_kel/` |

`bateman_conc()` in `scripts/assets/pk_curves.R` uses the `expm1` form of the Bateman curve, which has no cancellation as `ka` approaches `kel`. Two tests in `scripts/assets/tests.R` check it against the textbook form and across `ka - kel` from 1e-12 to 1e-6.

### Profile of the objective (04)

For each of the 22 `single_wave` curves, the residual sum of squares is computed on a 120 x 120 log-spaced grid of (`ka`, `kel`) with `F` solved by least squares at each point. The profile is symmetric about the diagonal by construction, so an off-diagonal minimum would come as a mirrored pair. For 17 curves the grid minimum is on the diagonal within one grid step.

### Identified quantities (05)

The paired Wilcoxon signed-rank test (intervention against baseline, within each diet arm, no multiplicity correction) is applied to quantities the curve identifies. The fit-derived rows use subjects with a reliable fit at both visits (`R² >= 0.9`, `kel` not at a bound). The observed-data rows use every completer. Rates and ratios are tested on the log scale, Tmax in minutes.

| Quantity | High-fructose | Low-fructose |
|---|---|---|
| `ka` (fit) | p = 0.093 (n = 16) | p = 0.038 (n = 14) |
| `kel` (fit) | p = 0.82 | p = 0.12 |
| `max(ka, kel)` (fit) | p = 0.29 | p = 0.020 |
| `min(ka, kel)` (fit) | p = 0.37 | p = 0.069 |
| Tmax, observed | p = 0.18 (n = 16) | p = 0.71 (n = 17) |
| Time to 50% of Cmax on the rise, observed | p = 0.33 | p = 0.80 |
| Cmax, observed | p = 0.003 (up) | p = 0.74 |
| AUC 0-360, observed | p = 0.13 | p = 0.48 |
| Terminal rate from t >= 180 min, observed (n = 8 per arm) | p = 0.73 | p = 0.44 |

The observed Tmax and time-to-half-Cmax medians are zero because the sampling grid is 30 minutes wide. The terminal rate needs three positive samples at or after 180 min, which only 8 subjects per arm have, so it carries little weight.

In the low-fructose arm the fast rate falls while the slow rate rises. At baseline the median gap between the two (log of `k_fast / k_slow`) is 0.17 and 50% of fits are on the diagonal. At the intervention visit the gap is about 0 and 86% are on the diagonal. The gap change is not significant (p = 0.15). The pattern fits curves becoming more symmetric, and the fit cannot say whether that is physiology or a change in which solution the optimizer reached.

### Simulation with independent fits (06)

The study design is simulated from `single_wave` truth: the 33 completers' actual sampling times, `Vd` and doses at both visits, between-subject variation of 20% in `ka`, `kel` and Cmax, and additive noise of 4.3 mg/L (the median residual SD of the real fits). Each curve is refitted with the same objective, penalties and bounds as `02b`, but with 15 random starts instead of 40. The paired Wilcoxon is applied to `ka`, `kel`, `max` and `min`. Each scenario has 100 replicates. A true change is a factor of 1.5 at the intervention visit.

Share of replicates with p < 0.05:

| Truth (baseline `ka` / `kel`) | `ka` | `kel` | `max(ka, kel)` | `min(ka, kel)` |
|---|---|---|---|---|
| No change (0.040 / 0.025) | 0.01 | 0.02 | 0.03 | 0.04 |
| No change (0.025 / 0.025) | 0.05 | 0.07 | 0.03 | 0.04 |
| `ka` x 1.5 (0.040 / 0.025) | 0.71 | 0.12 | 0.92 | 0.21 |
| `kel` x 1.5 (0.040 / 0.025) | 0.74 | 0.47 | 0.47 | 0.99 |
| `ka` x 1.5 (0.025 / 0.025) | 0.74 | 0.49 | 0.92 | 0.70 |
| `kel` x 1.5 (0.025 / 0.025) | 0.69 | 0.49 | 0.88 | 0.71 |
| `ka` x 1.5 (0.020 / 0.050) | 0.48 | 0.24 | 0.12 | 1.00 |

With no change, no test rejects more often than 7%. A `kel` change is flagged as a `ka` change in 74% of replicates, more often than as a `kel` change (47%). In the two diagonal scenarios the `ka`-up and `kel`-up rows are the same. Across scenarios, 45 to 59% of the fits had the fitted `ka` and `kel` in the wrong order relative to truth, and 32 to 47% ended on the diagonal. In the separated scenarios the median fitted `ka` is 0.09 to 0.26 too low (log) and `kel` 0.08 to 0.22 too high, because the fit pulls the two together.

The sorted pair works. A change in the fast rate gives `max` power 0.92 and `min` power 0.21 (well separated), and a change in the slow rate gives the reverse (0.47 and 0.99). So the fast or the slow rate can be shown to have changed with about 90% power at n = 33 for a 1.5-fold change, and the label is the missing piece.

### Assuming a stable kel (07 and the shared-`kel` fit)

The hypothesis for this study is that `ka` changes with the diet and `kel` stays the same. Fitting each subject's two visits together with one shared `kel` and a separate `ka` and `F` per visit uses that assumption. The swap no longer works, because a swapped solution would need the baseline `ka` and the intervention `ka` to be equal.

The simulation in `07` repeats the scenarios of `06` with the shared fit (40 replicates per scenario):

| Truth (baseline `ka` / `kel`) | `ka` rejection, shared `kel` | Bias in fitted `ka` change (log) |
|---|---|---|
| No change (0.040 / 0.025) | 0.08 | -0.01 |
| No change (0.025 / 0.025) | 0.15 | 0.00 |
| `ka` x 1.5 (0.040 / 0.025) | 1.00 | -0.03 |
| `ka` x 1.5 (0.025 / 0.025) | 1.00 | -0.01 |
| `ka` x 1.5 (0.020 / 0.050) | 1.00 | +0.05 |
| `kel` x 1.5 (0.040 / 0.025), assumption violated | 1.00 | +0.44 |
| `kel` x 1.5 (0.025 / 0.025), assumption violated | 1.00 | +0.41 |

When `kel` is stable, the shared fit detects a 1.5-fold `ka` change in every replicate and estimates its size without bias (the true log change is 0.41). When `kel` changes, the shared fit moves the entire `kel` change into `ka` (bias +0.44 against a true `kel` change of 0.41 in log units), and the `ka` test rejects in every replicate. The shared fit has a false-positive rate of 8% to 15% under no change, higher than the 1% to 7% of the independent fits. Only 0 to 3% of shared fits end on the diagonal.

The shared fit therefore detects a rate change with high power. It cannot tell a `ka` change from a `kel` change, because the assumption forces any change into `ka`. A significant `ka` result under this fit has to be reported as conditional on a stable `kel`.

On the real cohort (`fructose_shared_kel/02_fit_erie_model_fructose_shared_kel.R`), each visit keeps the mechanism that the independent fit selected, with its `f_delayed` and `t_lag` held at the independent estimates. The shared `kel` has the lower AIC in 26 of 33 subjects and is clearly worse (AIC higher by more than 2) in 5. Since the shared model has one parameter fewer, and since the simulation shows that it absorbs a real `kel` change without a large loss of fit, this agreement is weak evidence for a stable `kel`.

| Paired Wilcoxon on log(`ka` intervention / `ka` baseline) | High-fructose (n = 16) | Low-fructose (n = 17) |
|---|---|---|
| Independent fits | median +0.30, p = 0.093 | median -0.19, p = 0.11 |
| Shared `kel` | median +0.50, p = 0.053 | median -0.09, p = 0.81 |

## Statistics in the diet summaries

Each diet summary reports two tests side by side, and the plots show both. The paired Wilcoxon signed-rank test compares baseline with intervention within each arm. The linear mixed model (LMM, `scripts/assets/pk_stats.R`) is `log(value) ~ diet * visit + (1 | subject)` fitted with `nlme::lme` (raw minutes for Tmax). It gives the within-arm change with a 95% CI and p-value, and the diet x visit interaction, which tests whether the change differs between the arms. The interaction is the comparison that answers "did the diets act differently", and it appears as the subtitle of each plot. Tables are in `diet_parameter_tests.csv` and `diet_curve_metrics_tests.csv` in each `diet_summary/` folder. The shared-`kel` summary does not test `kel`, which is assumed.

Interaction p-values in the independent-fit summary (low-fructose arm n = 14, high-fructose arm n = 16): `ka` 0.010, Cmax 0.030, AUC 0.052, Tmax 0.16, `kel` 0.39, `F` 0.92. In the shared-`kel` summary `ka` has p = 0.038 and `F` 0.17. The `ka` interaction inherits the identifiability limits above, because a `kel` change shows up as a `ka` change, and the shared-`kel` one is conditional on a stable `kel`. No multiplicity correction is applied, and these summaries run 12 or more tests each.

## What can be claimed

- The fast and the slow rate can be reported for each curve and compared between visits (`k_fast`, `k_slow`). These are identified, and the tests above give about 90% power for a 1.5-fold change.
- `ka` and `kel` separately cannot be reported from these curves. The labels in `fit_results.csv` are the optimizer's choice and were wrong about half the time in the simulations.
- A `ka` effect can be reported as conditional on a stable `kel` when the shared-`kel` fit is used. The text has to say so, and the conditional result cannot be used to argue that `kel` was stable.
- The model-free curve metrics (Cmax, AUC, Tmax) do not depend on the split. In the high-fructose arm, observed Cmax rose (p = 0.003), which points at the amount absorbed or `Vd`. The rates are unaffected by that finding.

## What would separate the two rates

- A direct measure of `kel` at both visits, for example an intravenous dose or a literature value that applies to the cohort. This turns the shared-`kel` assumption into a measurement.
- Denser early sampling. The first sample after dosing is at 30 min and the rise is invisible below that. Time points at 5, 10, 15 and 20 min would show the rising edge directly.
- More late samples with positive concentrations after 180 min, which would give an observed terminal slope for more than 8 subjects per arm.
- A population model with a prior on `kel`, which borrows strength across subjects. The labels are still symmetric without an outside prior.

## Limits of the tests

- The simulations have one wave in the truth. 46 of the 68 real curves use `two_wave`, so real recovery is probably worse than the simulated one.
- Only 1.5-fold changes and one noise level (4.3 mg/L, additive) were simulated.
- The simulated refits used 15 random starts for the independent fit and 25 for the shared fit, fewer than the 40 used in the independent fit. The shared simulation has 40 replicates per scenario, so its rejection rates carry a standard error of about 5 to 8 percentage points.
- Wave structure in the real shared fit is fixed at the independent estimates, so the shared fit tests only `ka`, `kel` and `F`.
- The Wilcoxon tests use no multiplicity correction.
