# Fructose PK model specification

This documents the pharmacokinetic (PK) model fit by `scripts/02_fit_erie_model.R`, the biological and statistical reasoning behind it, and - most importantly - what its output can and cannot be used to claim. The reasoning draws on the documents in `former_models/*/*.docx` and was re-derived and re-validated for this pipeline (see "Origin of the model"). The history of problems and fixes is in `docs/problems-and-fixes.md`.

## Summary

**What is modelled.** For every subject x visit, two plasma curves come from the same blood draws: 12C-fructose (1 g/kg, drunk as a liquid) and 13C6-fructose (a fixed 120 mg tracer in an enteric capsule). Both are described by one-compartment first-order oral kinetics with the same absorption rate `ka` and elimination rate `kel`, a separate bioavailable fraction `F` per curve, and a volume of distribution `Vd` computed from the subject's weight, height and sex.

```
12C    dose D12 --(ka, F_12C)---------------------> blood --kel--> cleared
13C6   dose D13 --(k_release)--> gut --(ka, F_13C6)--> blood --kel--> cleared

concentration = amount in blood / Vd
```

**Formulas.** `t` is in minutes, doses in mg, `Vd` in litres and concentrations in mg/L. Two building blocks are used:

```
bateman(t; D, F)   one dose D absorbed with ka, cleared with kel (closed form)
  = F * ka * D / (Vd * (ka - kel)) * (exp(-kel * t) - exp(-ka * t))      0 for t <= 0
  (when ka = kel: F * ka * D * t * exp(-ka * t) / Vd)

capsule(t; D, F)   the dose first dissolves (rate k_release), then behaves as above
  dR/dt = -k_release * R              R(0) = D
  dG/dt =  k_release * R - ka * G     G(0) = 0
  dA/dt =  F * ka * G - kel * A       A(0) = 0
  = A(t) / Vd                          (solved numerically)
```

The plain model (`single_wave`) is:

```
C_12C(t)  = bateman(t; D12, F_12C)
C_13C6(t) = capsule(t; D13, F_13C6)
```

Many curves have a post-peak dip that one absorption phase cannot produce, so `two_wave` splits the 12C dose into an early part and a part that starts `t_lag` minutes later. `f_delayed` is the late fraction:

```
C_12C(t) = bateman(t; (1 - f_delayed) * D12, F_12C) + bateman(t - t_lag; f_delayed * D12, F_12C)
```

13C6 then uses whichever of three mechanisms fits it best (chosen by AIC on its own residuals, independently of 12C's model):

```
delayed release   capsule(t; D13, F_13C6)
                  under two_wave, split like 12C: (1 - f_delayed_13C6) at t and f_delayed_13C6 at t - t_lag
onset lag         bateman(t - t_lag1_13C6; D13, F_13C6)
second wave       bateman(t; (1 - f_delayed2_13C6) * D13, F_13C6) + bateman(t - t_lag2_13C6; f_delayed2_13C6 * D13, F_13C6)
```

`ka` and `kel` are shared by every term in both curves.

**Inputs that are not fitted.** `D12` = 1000 mg x body weight in kg. `D13` = 120 mg. `Vd` is Nadler's blood volume, with `h` in m and `w` in kg: men `0.3669 h^3 + 0.03219 w + 0.6041`, women `0.3561 h^3 + 0.03308 w + 0.1833` (see "Volume of distribution"). 12C concentrations are baseline-corrected (the fasted level is subtracted).

**Parameters.** Which ones exist in a fit depends on the model and mechanism chosen: `single_wave` has `ka`, `kel`, `F_12C`, `F_13C6` and, for 13C6's delayed release, `k_release`. The rest only appear when their mechanism is selected.

| Parameter | Meaning | Applies to | Bounds | Explained in |
|---|---|---|---|---|
| `ka` | absorption rate constant (1/min) | both curves, shared | [1e-4, 1] | "Design choices", "Fitting bounds" |
| `kel` | elimination rate constant (1/min) | both curves, shared | [0.005, 0.1] | "Fitting bounds" |
| `F_12C`, `F_13C6` | fraction of the dose that reaches the blood | one per curve | [1e-4, 1] | "Design choices", "Fit quality and what to trust" |
| `k_release` | capsule dissolution rate constant (1/min) | 13C6 delayed release only | [0.001, 1] | "Design choices" |
| `f_delayed` | fraction of the 12C dose in the late wave | 12C, `two_wave` | [0.001, 0.999] | "The post-peak dip and the lagged-dose model" |
| `t_lag` | start of the late wave (min) | 12C `two_wave`; shared with 13C6's delayed release under `two_wave` | [5, 150]; 30 to 150 if 12C at t=30 is above 5 mg/L | "The post-peak dip and the lagged-dose model" |
| `f_delayed_13C6` | fraction of the 13C6 dose in the late wave | 13C6 delayed release, under `two_wave` | [0.001, 0.999] | "13C6's own wave choice" |
| `t_lag1_13C6` | onset delay of the 13C6 curve (min) | 13C6 onset lag | [0, 90] | "13C6's onset-lag and second-wave candidates" |
| `f_delayed2_13C6`, `t_lag2_13C6` | fraction and start (min) of 13C6's own second wave | 13C6 second wave | [0.05, 0.95] and [30, 150] | "13C6's onset-lag and second-wave candidates" |

**How a fit works.** Each subject x visit is fitted with `single_wave`, and with `two_wave` when the raw 12C data show evidence of a second wave (or `single_wave` fits very poorly); the choice between them is made by AIC on 12C's residuals, except that unambiguous dip evidence forces `two_wave`. Every fit uses multi-start bounded optimization (`optim`, L-BFGS-B). The objective of the main fits is, per curve, the squared error divided by that curve's total variance (unweighted), plus two soft penalties: predicted Tmax at least 30 min, and predicted Cmax within +/-10% of the observed Cmax. 13C6's onset-lag and second-wave candidates are refitted afterwards with `ka` and `kel` fixed, using plain squared error. Fits flagged for a bound, a low R² or non-convergence get a denser retry. See "Fitting procedure" and "Choosing between single_wave and two_wave".

**What to trust.** `ka` and `kel` are reasonably well identified. `F` is only interpretable within a subject (for example baseline versus intervention) because `F` and `Vd` enter the model only as `F * D / Vd`; report it as conditional on the Nadler `Vd`. Check `r2_12C`, `r2_13C6`, `converged` and the bound flags before using an individual fit. See "Fit quality and what to trust".

## Background: what's being measured

Each fructose challenge test (FCT) gives two simultaneous plasma curves from the same blood draws, at t = 30, 60, 90, 120, 150, 180, 240, 360 min after a combined dose:

- **12C-fructose**: 1 g/kg body weight of unlabeled fructose, drunk in water (an individualized, large dose).
- **13C6-fructose**: a fixed tracer dose in an enteric capsule (120 mg / 644.78 µmol, confirmed administered dose; see `docs/problems-and-fixes.md` for the dose discrepancy this project resolved).

Every subject has this pair of curves at both FCT1 (coded `baseline`) and FCT2 (coded `intervention`, post 4-week diet intervention) - `visit` in the cleaned data and in every downstream script is `baseline`/`intervention`, and `FCT1`/`FCT2` appear only as raw-file labels.

### What the model is fit to

`02_fit_erie_model.R` prepares the cleaned data as follows before any fit:
- Concentrations are converted from µmol/L to mg/L with each isotope's own molecular weight (`MW_12C`, `MW_13C6` in `data/processed/erie_constants.csv`), so they match the mass units of the doses.
- 12C is baseline-corrected: the subject's own fasted (t=0) 12C level is subtracted from every 12C observation, since the dose is added on top of endogenous fructose. 13C6 needs no correction (the tracer is absent at baseline).
- t=0 is excluded from fitting - it carries no information (12C is 0 by construction after correction, and 13C6 hasn't been dosed yet; both models predict exactly 0 there for any parameters) - but is kept for plotting so the observed point still shows.
- Rows without a body weight (the intervention visits of the two study dropouts) are dropped.

## Design choices: shared and separate parameters

The equations are in the Summary. This section explains why the two curves share `ka` and `kel` while each has its own `F`, and why the capsule has a dissolution step.

### Why `ka` and `kel` are shared

Once fructose is in the gut, the same transporters absorb it and the same machinery clears it from blood, whatever the route (liquid or capsule), so `ka` and `kel` should be the same for both curves of a subject and visit. Fitting the two curves independently gave poor cross-curve agreement in `kel` for many subjects, which led to testing the shared-parameter model directly. In the former joint model, sharing raised the 13C6 median R² from about 0.70 (independent fits) to about 0.78 (joint fit) and left far fewer subjects below R² 0.3 (see "Origin of the model").

### Why each curve has its own `F`

The two doses differ roughly 1000-fold (1 g/kg liquid, 120 mg capsule). First-pass clearance of fructose in the small intestine depends on dose: at low doses the small intestine metabolizes most of the fructose before it reaches the liver, and this capacity saturates above roughly 1 g/kg (Jang et al. 2018), which is the 12C dose used here. The two curves are therefore expected to differ in intestinal saturation and in `F` even with shared `ka` and `kel`, so `F_12C` and `F_13C6` are fitted separately and their difference is interpreted.

### Why the capsule has a dissolution step

Instant availability of the full capsule dose at t=0 assumes infinitely fast dissolution. The dissolution compartment (`k_release`) lets the data determine how fast the capsule releases its contents. It is the delayed-release mechanism of 13C6; in the current results 14 of 68 fits use it, and the others use the onset-lag or second-wave mechanisms described below.

## The post-peak dip and the lagged-dose model

Many curves show a dip after an initial peak, then a second, often higher, rise before the final decline, most clearly on `ER01` and `ER09`, whose 12C curves are flagged at both visits. How common this is depends on the detector. With the current detectors (see "How the dip is detected" and "Choosing between single_wave and two_wave"), counted over the 68 subject x visit curves of each isotope from `data/processed/` (t=0 excluded and 12C baseline-corrected, as in the fit):

| Check | 12C | 13C6 |
|---|---|---|
| `has_peak_dip_rise()` | 43/68 (63%) | 37/68 (54%) |
| `has_near_peak_neighbor()` | 29/68 (43%) | 21/68 (31%) |
| either | 53/68 (78%) | 47/68 (69%) |

These checks are permissive gates that decide whether `two_wave` is attempted, and the true number of curves with a second wave is unknown. Where the dip check fires on 12C, the point that triggers it falls at 90 min (11 curves), 120 (10), 150 (17) or 180 (5), spread across the sampling grid.

**No parameter choice in the base model (or any variant with a single shared elimination compartment) can produce this shape.** A one-compartment absorption curve is monotonic after its single peak for any `ka`, `kel` and `F`: it is a difference of two exponentials and has no second local maximum. Any number of absorption pathways that all start at t=0 and share one terminal elimination rate converge to the same single-exponential decay tail, so summing them can reshape or delay a single peak and cannot create a local minimum followed by a second rise. A "biphasic gastric emptying" variant (two gastric pools emptying at different rates from t=0 into a shared absorption compartment) confirmed this: thousands of random parameter draws and several hand-engineered attempts never produced a second hump.

**A time-lagged second dose produces the shape** (`simulate_lagged_dose()` in `scripts/assets/pk_curves.R`). A fraction `f_delayed` of a dose contributes nothing until `t_lag`, then behaves like a fresh dose given at that later time, which is mathematically two separately timed doses. The lag creates a discontinuity in the delivery rate and a gap between the two absorption phases, which pathways that all start at t=0 cannot have. Analytically and numerically it yields a local minimum followed by a second, higher peak. In a pilot on `ER01`, `ER03` and `ER09` (`results/model-post-peak-dip-pilot/`), the two dip cases (`ER01`, `ER09`) changed from a single smooth hump to a double-humped fit that follows the dip, and `f_delayed` falls to near 0 for subjects that need no second wave.

`02_fit_erie_model.R` fits two candidates per subject x visit, `fit_subject_visit_single_wave()` (the plain joint delayed-release model) and `fit_subject_visit_two_wave()` (the lagged-dose extension), and selects one (see "Choosing between single_wave and two_wave"). Both candidates' fits are written out in full (`results/fit_results_single_wave.csv`, `results/fit_results_two_wave.csv`) alongside the selected combination (`results/fit_results.csv`).

`t_lag` is bounded above by 150 min. After 180 min the sampling gaps widen to 60 min (180-240) and 120 min (240-360), and a `t_lag` in one of them lets the optimizer place a second wave where no sample can confirm or refute it. `ER12` intervention illustrates this: at `t_lag`=179 min, `r2_13C6` falls from single_wave's 0.95 to 0.75, with a large peak and no supporting data around t=200-220.

The lower bound of `t_lag` depends on the subject's own 12C concentration at t=30: 5 min when it is at or below `EARLY_LAG_OK_MGL` (5 mg/L), 30 min otherwise (`t_lag_lower` in `fit_subject_visit_two_wave()`). The 0-30 min window contains no sample. A low bound for everyone lets the optimizer place an early onset delay there (small `t_lag`, `f_delayed` near 1, so that almost the whole dose is late) whatever the subject's own data say: `ER06` intervention reaches `t_lag`=19.7 min and `f_delayed`=0.999 although its t=30 sample (~80 mg/L) is already near its peak. A floor of 30 min for everyone would exclude a real early delay in a subject whose t=30 concentration is near zero. With the per-subject bound, `ER06` intervention converges to `t_lag`=35.3 min and `f_delayed`=0.320, and `ER01` and `ER09` (t=30 values well above 5 mg/L, `t_lag` 85-114 min) are unaffected.

The bound is a step because every value between 5 and 30 min lies in the same unsampled window. For `ER04` baseline (t=30 = 51.4 mg/L) a floor interpolated between the two values would be 21.6 min, and the fit lands exactly there with a flat-then-rise artifact. Either t=30 is near zero, a delay somewhere in the window is plausible and 5 min is the permissive floor (the data cannot locate the delay within the window), or t=30 is elevated, no point in the window is supported and the floor is 30 min. The `BOUNDS_LAGGED` list in `02_fit_erie_model.R` keeps the permissive floor (5) for the retry driver's seed generation; the subject-specific floor is what the optimizer receives, and seeds are clipped to it.

### How the dip is detected

The raw-data checks in `scripts/assets/pk_diagnostics.R` involve no curve fitting. They use only a curve's own observed points and gate whether the more flexible `two_wave` model is attempted. Fitting first and relying on AIC or R² to notice that there was nothing to explain lets a flexible model explain ordinary sampling noise around a single peak.

**`peak_dip_rise_info()`** flags a post-peak point that lies at least 15% (`min_rise_frac`) above the straight line joining its two neighbors. A single-exponential elimination phase is convex (it decays ever more slowly), so a real point on it lies at or below the chord between its neighbors, and a point clearly above the chord shows non-monotonic decay whether or not it becomes a new local maximum. This catches a plateau or slowed decline that never rises above what came before, which a test of whether a later point exceeds an earlier local minimum misses. Design details:
- The first-wave candidate is the first point after which concentration turns down and that is at least half (`min_first_peak_frac`) of the curve's overall peak. A second wave is often taller than the first (`ER01`, `ER09`), so anchoring on the global maximum would miss it, and the half-of-peak floor keeps a small early wobble, while concentration is still climbing, from counting as a first wave.
- At least one point must separate the first peak from any candidate second peak, because the point right after the peak still describes how sharply the first wave turns over.
- When several points deviate, the one with the largest deviation is reported as `trigger_time`, with its magnitude as `excess`. `trigger_time` anchors the search seeds (see "Two-stage fit") and `excess` feeds the "definite two_wave" override.

**`has_near_peak_neighbor()`** covers curves with no visible dip (see "Choosing between single_wave and two_wave").

**`has_onset_lag_evidence()`** addresses a delay before a single wave starts (a standard PK absorption lag time). It fires when the curve's first observed point is below 25% of its own peak, which indicates that absorption had barely started. The 25% threshold lies above both validation cases (`ER25` baseline: 18.4%, `ER06` intervention: 2.2%).

### Choosing between single_wave and two_wave

Selection is by AIC, `n*ln(RSS/n) + 2k`, computed on 12C's own raw residuals (`K_12C_SINGLE_WAVE = 3`: `ka`, `kel`, `F_12C`; `K_12C_TWO_WAVE = 5`: plus `f_delayed`, `t_lag`). 13C6 is fitted after 12C's structure is decided, so basing the selection on 12C alone keeps 13C6's fit quality out of the question whether 12C's curve earns a second wave, and avoids mixing the two concentration scales that a combined AIC would need.

`two_wave` is attempted when 12C's t=30 sample exists and at least one of three conditions holds:
- The raw 12C data pass `has_peak_dip_rise()`.
- The raw 12C data pass `has_near_peak_neighbor()`: the point immediately before or after the curve's peak is at least 85% of the peak's height. This catches two waves that overlap closely enough to blend into a flat top or an irregular rise without a visible dip: `ER02` baseline (flat top, R² 0.955 to 0.997) and `ER30` baseline (irregular pre-peak rise, R² 0.859 to 0.957). t=30 is excluded as a near neighbor, since closeness to the peak at the first sample reflects fast absorption. `ER06` intervention illustrates the reason: its only near-peak point was t=30, and its `two_wave` fit pinned `t_lag` at the lower bound, where the fits of `ER02` (51.5 min) and `ER30` (75.7 min) settle in the interior.
- `single_wave` fits 12C very poorly (R² below `R2_ALWAYS_TRY_TWO_WAVE` = 0.50). Two waves close together in time can add up to one smoothly accelerating rise with no down-turn and no plateau, which neither detector sees, and the duodenal-brake feedback that paces gastric emptying of a large caloric and osmotic load such as the 1 g/kg dose pauses and resumes on a timescale of tens of minutes, so such overlap is plausible. `ER21` baseline shows it: 12C at t=30/60/90 is 32/61/132 mg/L (more than doubling with no earlier down-turn), `single_wave` gives R²=0.31, and an unconstrained `two_wave` search found R²=0.99 at `t_lag`=57 min and `f_delayed`=0.80. The threshold lies well below `R2_RELIABLE_MIN` (0.70), so this path opens only when `single_wave` fails to fit.

A t=30 sample is required because without it the 0-60 min window has only its t=0 and t=60 endpoints to anchor a 5-parameter fit. That fit is under-identified, and AIC accepts any shape that 5 sparse points fit almost exactly (`ER05` baseline: t=30, 90, 150 and 240 are all missing, leaving 5 points for 12C). Checking the raw data before fitting keeps a flexible model from producing a plausible-looking fit for a curve without a second wave.

AIC is used without the small-sample correction. AICc (`+2k(k+1)/(n-k-1)`) at n=8 raises the complexity gap between k=3 and k=5 from 4 points to 28 and would exclude all three confirmed plateau cases (`ER02`, `ER06` intervention, `ER30`), which have independently verified R² gains. Protection against overfitting comes from the raw-data gates and the bounds (the data-dependent `t_lag` floor, the per-wave `MIN_TMAX`, the `ka`=`kel` diagonal seeding).

### Definite two_wave override

Where the raw data make a second wave unambiguous, `two_wave` is selected without the AIC comparison. The trigger is `peak_dip_rise_info()`'s `excess` at or above `EXCESS_DEFINITE_TWO_WAVE` (0.90 in `02_fit_erie_model.R`): a marginal AIC edge for `single_wave` (for example from an unlucky local optimum) should not overrule direct evidence in the data. Only a failed `two_wave` fit (no AIC) falls back to `single_wave`.

The threshold lies well above the 15% detection floor and above `ER01`'s 0.40, and is set from the cohort distribution: only `ER08` baseline (1.57) and `ER11` baseline (1.01) exceed it, and both are confirmed by independent evidence (`ER08`: 13C6 peaks at the same timepoint; `ER11`: a sustained dip followed by a peak over several points). Most dip cases, including `ER01` and `ER09`, go through the AIC comparison.

### Two-stage fit: 12C first, then 13C6

`fit_subject_visit_two_wave()` fits the two curves in two stages:
1. Stage 1 fits `ka`, `kel`, `F_12C`, `f_delayed`, `t_lag` to 12C alone.
2. Stage 2 fixes those and fits `F_13C6`, `k_release` and 13C6's own `f_delayed_13C6` to 13C6 alone (see "13C6's own wave choice").

The shared lag is estimated from 12C alone. 13C6 has its own onset mechanism (`k_release`, the dissolution of the enteric capsule), and a lag informed by 13C6's residuals attributes 13C6's slow, `k_release`-driven onset to a second wave that then constrains 12C even when 12C's own data give no support for one. `ER32` baseline shows this: its 12C t=30 sample is already near its peak, so no onset delay is plausible, yet a joint fit pins `f_delayed`~1 and `t_lag`~19 min and lowers both curves' R² below the single_wave fit; `ER35` behaves the same way. 12C is also the better source for the lag: its signal is much larger and cleaner, and the 1 g/kg osmotic and caloric liquid load is the physiological trigger of biphasic emptying (the tracer dose is tiny), so a real dip shows up in 12C's data if it exists.

Stage 1 mixes several seed families, because passing the raw-data gate does not place the optimizer's generic seeds near the evidence that justified the attempt: the independent 12C pre-fit (with `f_delayed` near 0, which recovers a plain Bateman fit); a few generic no-lag and dip-like seeds (real dip cases converge near `t_lag`~90); random log-uniform seeds; the `ka`=`kel` diagonal (see "Known limitations and open questions"); and evidence-anchored seeds at `trigger_time` and 20 and 40 min before it, since the trigger point is where the second wave becomes visible, which can lie after its start. The evidence-anchored seeds matter for `ER04` baseline: without them the search converges to `t_lag`=21.6 (an early near-total-delay solution with a marginally better AIC) and misses the t=150 plateau that triggered the attempt.

### 13C6's own wave choice

13C6 has its own free `f_delayed_13C6`, fitted in stage 2 after 12C's `ka`, `kel` and `t_lag` are fixed. 13C6 then decides independently whether it rode only the first wave (`f_delayed_13C6` near 0), only the second (near 1) or both, and shares only `t_lag` (the timing of gastric emptying's second wave, a systemic event) with 12C. The capsule contents need not split between the waves in the proportion of the much larger liquid dose, and the capsule may have emptied in one wave. `ER03` baseline shows the need: 12C has two waves (R² = 0.99), while 13C6 has one early peak that decays monotonically with nothing where an inherited second wave would put one, and inheriting 12C's fraction gives `r2_13C6` = 0.075.

`MIN_TMAX` applies to each wave's own peak (`tmax_wave1`, `tmax_wave2`, computed by simulating each wave alone) as well as to the combined curve, because the combined maximum only constrains the taller wave. `ER09` baseline and intervention illustrate this: wave 1 is taller (combined peak at t=38 and t=33) while wave 2 alone peaks at t=151 and t=147, and `ER16` intervention has a wave 2 peaking at t=162. This bounds the peak time from below. The upper end is limited only by `t_lag`'s 150 min cap plus the rise time set by `ka` and `kel` (see "Known limitations and open questions").

### 13C6's onset-lag and second-wave candidates

Besides delayed release, 13C6 has two instant-bolus candidates. All three are available under both `single_wave` and `two_wave` and are compared by AIC on 13C6's own residuals, so 13C6's mechanism is independent of 12C's model. The candidates use `ka` and `kel` fixed at the values from the 12C-anchored fit and leave 12C's reported parameters untouched. The AIC uses K=2 for delayed release under `single_wave` (`F_13C6`, `k_release`) and K=3 under `two_wave` (plus `f_delayed_13C6`), K=2 for the onset lag (`F_13C6`, `t_lag1_13C6`) and K=3 for the second wave (`F_13C6`, `f_delayed2_13C6`, `t_lag2_13C6`).

**Onset lag.** The candidate gives 13C6 a delay `t_lag1_13C6` (bounds 0-90 min) before absorption starts and then simulates the whole dose as an instant bolus (`simulate_two_lag_dose()` with `f_delayed` fixed at 0). The smooth first-order dissolution of delayed release is fastest at t=0 and cannot represent a 13C6 signal that stays near zero for about 30 min and then rises sharply. An instant bolus after a delay has no second undelayed portion, so the near-total-delay solution that `t_lag` and `f_delayed` allow for 12C cannot occur. `k_release` is omitted: with a delay in the model it sat at its upper bound in all three validation cases (`ER06` intervention, `ER25` baseline, `ER25` intervention), because 30-min sampling cannot separate fast dissolution from instant release once the delay explains the flat start. The candidate is attempted when `has_onset_lag_evidence()` fires. Validation: `ER25` baseline (R² 0.71 to 0.94, `t_lag1_13C6` = 25 min), `ER06` intervention (0.64 to 0.94, 28 min) and `ER25` intervention (R² 0.75, 27 min). `ER23` intervention has 13C6 at 0.009 mg/L at t=30 and 0.165 at t=60 (5.6% of its peak); with 12C as `two_wave`, delayed release collapses to `f_delayed_13C6`=0.001 and `r2_13C6`=0.35, while the onset-lag candidate gives about 27 min and `r2_13C6`=0.78.

**Second wave.** The candidate (`f_delayed2_13C6`, `t_lag2_13C6`) has the shape of 12C's `two_wave` model (an instant bolus per wave) and is attempted from 13C6's own raw data with the same two detectors. It is mutually exclusive with the onset-lag candidate, since both explain a flat or slow start. 13C6 is noisier than 12C, and the detectors flag 47/68 of its curves (69%; 12C: 53/68, 78%), a high rate for a low-signal curve, so the candidate carries two safeguards. `f_delayed2_13C6` is bounded to [0.05, 0.95] (12C: [0.001, 0.999]), which excludes the near-total-delay solution: `ER21` baseline (the strongest dip excess in the cohort, 3.416) reaches R² 0.508 to 0.741 only with `f_delayed2_13C6`=0.999, driven by a near-zero crash and rebound (0.0117 mg/L at t=150, 0.0354 at t=180) that looks like measurement noise. `t_lag2_13C6` is bounded to [30, 150], since the onset-lag candidate covers delays before the wave starts and this candidate models a later, separate hump. `ER11` baseline (R² 0.787 to 0.846, `f_delayed2_13C6`=0.751) and `ER18` intervention (0.860 to 0.879, `f_delayed2_13C6`=0.095) are accepted fits away from the extremes. A result within 1e-3 of the lower bound of `f_delayed2_13C6` is rejected without an AIC comparison, because it indicates no second wave: `ER18` intervention wins its AIC comparison at `f_delayed2_13C6`=0.05 and `t_lag2_13C6`=150 simultaneously (both bounds) although its dip evidence lies at t=120, and with the rejection it falls back to delayed release (`r2_13C6`=0.66).

### Known limitations and open questions

- **The multistart search can miss a much better optimum when `ka` and `kel` are close to each other.** This is a known hard region for one-compartment models ("flip-flop" kinetics): the `ka/(ka-kel)` term of `bateman_conc` makes the objective surface narrow for gradient-based search, and generic seeds rarely land near it. `ER06` intervention (`single_wave`) showed it: the pipeline reported `ka`=0.016, `kel`=0.051, R²=0.945 (neither at a bound), while multiple `optim()` starts converging to one point found `ka`=0.023, `kel`=0.025, R²=0.986. Explicit `ka`=`kel` diagonal seeds (`v` from 0.008 to 0.07, `F` fixed low) are added in `fit_curve_independent()`, the joint `single_wave` fit and `two_wave` stage 1. Cohort-wide R² values have not been re-checked against a search without these seeds, and reparameterizing to avoid the `ka/(ka-kel)` singularity is untried.
- **A wave's peak can land late.** Nothing limits the peak time of a wave beyond `t_lag`'s 150 min cap, and a slow `ka` can place the peak well after `t_lag` (`ER16` intervention: wave 2 with `t_lag` = 111 min peaks at t=162 min). t=162 lies inside the densely sampled 0-180 min window, so it has data support, unlike a peak inside the 180-240 or 240-360 min gaps. An upper bound is worth adding if a fitted wave peaks inside one of those gaps.
- **`MIN_TMAX` is 30 min.** A sweep of 30, 20 and 10 min on `ER03` baseline under a proportionally weighted objective moved R² only from -0.86 to -0.67, so the floor was not the limiting constraint there. The effect under the unweighted objective and the lagged-dose model is untested.
- **A two-lag model for 12C is untested at cohort scale.** `simulate_two_lag_dose()` can give 12C's first wave its own onset delay `t_lag1` in addition to the second-wave delay. On `ER03` baseline it matched the former model's R² (0.71 vs 0.71) with `f_delayed` near 0, so it mostly used the onset lag. That is one subject, the fit did not converge, and its 12C R² (0.70) was below the simpler unweighted fit (0.85).

## Volume of distribution (Vd)

`Vd` is the total blood volume of each subject x visit, computed from weight, height and sex with Nadler's equation:

```
Men:   Vd (L) = 0.3669 x height(m)^3 + 0.03219 x weight(kg) + 0.6041
Women: Vd (L) = 0.3561 x height(m)^3 + 0.03308 x weight(kg) + 0.1833
```

> Nadler DA, Hidalgo JU, Bloch T. *Prediction of blood volume in normal human adults.* Surgery. 1962;51(2):224-232.

It is implemented as `nadler_blood_volume()` in `scripts/assets/pk_curves.R`, which is project-agnostic like the rest of that file.

Blood volume (about 65-75 mL/kg for a typical adult) is an assumption about fructose's distribution volume. Fructose is small, freely water-soluble and unbound to protein, so it is expected to equilibrate into interstitial fluid as well as the vascular compartment. For glucose, van der Crabben et al. measured a `Vd` of 191-206 mL/kg across three tracers and showed that it equals the extracellular fluid space (ECFV, about 150-200 mL/kg) and exceeds blood or plasma volume. An ECFV-based `Vd` would be 2-3 times larger than the Nadler blood volume and would give correspondingly larger `F` values. This project uses the individual weight, height and sex estimate over a flat per-kg ratio, and the choice of blood volume has not been checked against a fructose-specific `Vd`.

`F` enters the model only as `F * dose / Vd`, so `F` scales with the assumed `Vd`, while `ka` and `kel` do not depend on it. Within-subject comparisons of `F` (baseline vs. intervention) stay valid, because each subject's own Nadler `Vd` applies to both visits and any systematic bias in it cancels in a paired comparison. `Vd` differs between the two visits only when the subject's weight changed.

> van der Crabben SN et al. *Relationship between glucose volume of distribution and the extracellular space: a multiple tracer study.* Metabolism. 2011.

## Fitting bounds (why these numbers)

The Summary lists the bounds of every parameter. The rationale for the four base parameters is here, and the bounds of the wave and lag parameters are explained where those parameters are introduced.

| Parameter | Bounds | Rationale |
|---|---|---|
| `kel` | [0.005, 0.1] /min (t½ ≈ 7-140 min) | Anchored to Hannou et al. 2018's review of fructose metabolism. Unconstrained fits give physiologically implausible and inconsistent `kel` across subjects. |
| `ka` | [1e-4, 1] /min | Wide, because absorption-rate differences (baseline vs. intervention, liquid vs. capsule) are part of the research question. The upper bound of 1/min is a generous ceiling. |
| `F_12C`, `F_13C6` | [1e-4, 1] | A fraction of a dose. |
| `k_release` | [0.001, 1] /min (t½ ≈ 0.7-700 min) | Wide enough to cover anything from near-instant to very slow capsule dissolution. |

Two soft penalties are added to the objective. The first pulls the predicted Tmax above 30 min (`MIN_TMAX`, weight `TMAX_LAMBDA`), which rules out a spurious early-spike solution in which the model absorbs and clears almost instantly and produces a sharp early peak between the sparse observed timepoints. The second pulls the predicted Cmax toward the observed Cmax within +/-10%, which discourages systematic over- or undershoot of the peak and leaves the rest of the curve free. Both are evaluated on a dense time grid, because degenerate solutions look fine only at the sparse observed points. Both are smooth (a squared shortfall or excess that is zero at and below the threshold), because `optim()`'s L-BFGS-B relies on finite-difference gradients and a hard penalty boundary makes the surface rough to search.

## Fitting procedure

PK objective surfaces are multimodal (a fast-absorption, fast-clearance solution can fit nearly as well as a slow/slow one), so every fit uses multi-start optimization; a single run from one starting point is insufficient. Each subject x visit is fitted in two stages:

1. An independent pre-fit per curve (closed form for 12C, numerical for 13C6, whose lone curve follows the same ODE) generates informed starting points for stage 2.
2. The joint fit is seeded from the pre-fit estimates, a systematic grid and random starting points (log-uniform for rate constants). Its objective is the sum over curves of each curve's `sse / ss_tot` (roughly 1 - R² of that curve), plus the two penalties above. Dividing by each curve's total variance keeps the 12C curve, whose concentrations are about 1000 times those of the 13C6 tracer, from dominating the objective and the optimizer from ignoring 13C6.

Both stages use the same generic multi-start engine (`scripts/assets/pk_fit.R`).

**Objective weighting.** `build_joint_objective()` supports an optional proportional weighting `1 / max(pred, floor)^2` (floor: 1% of the curve's Cmax, `weight_floor_frac`). The fits use `proportional_weighting = FALSE`. Residuals of an unweighted fit show heteroscedasticity (the squared-residual scale is about 24 times higher in the top than in the bottom quartile of predicted concentration, for both curves), so unweighted SSE lets each curve's peak region dominate at the expense of the tail, which carries most of the information about `kel`. Weighting discounts exactly the region where an informative early feature lives, because a large predicted value receives a small weight. The 13C6 curve of `ER03` baseline shows the effect: it spikes at t=30 and decays by t=90, and its weighted fit gives R² = -2.18, which loosening `MIN_TMAX` from 30 to 10 barely changes (-0.86 to -0.67), while the unweighted fit with the same `MIN_TMAX`, bounds and an ordinary multistart gives R² = 0.37. The former model's fit of the same subject reaches R² = 0.71, and its objective is unweighted (it normalizes each curve by its own `ssTot` only, in `objective_joint()` of `former_models/MixedModel/Scripts/fructose_joint_model_final.R`), with cohort medians of about 0.9+ for 12C and 0.7-0.9 for 13C6 from the same objective. Weighting remains available for a sensitivity check. The heteroscedasticity persists in the unweighted fit, so `r2_12C` and `r2_13C6` should be checked across the cohort.

**Parameter scaling.** `fit_multistart()` sets `parscale` for L-BFGS-B per start from that start's own magnitude (floored at 1% of the bound range). L-BFGS-B's step sizing assumes roughly unit-scaled parameters, and bounds such as [1e-4, 1] for `ka` and `F_12C` are wide relative to typical fitted values (0.01-0.05), so a bound-width `parscale` is 20-100 times too large and leaves the search stuck near its start. `ER11` baseline shows the effect: the objective (with the `MIN_TMAX` and Cmax penalties) scores a fit at `t_lag`=47.3 as 0.208, a search scaled by each start's own magnitude reaches 0.188 at `t_lag`=120, which matches the dip evidence (`trigger_time`=150) that justified the `two_wave` attempt, and bound-width scaling stays at 0.26-0.28 like a search without scaling. A multistart seed sits in the right neighborhood (pre-fit, evidence anchor or systematic grid), so its magnitude is a better estimate of the local scale than the global bound width.

**Reproducibility.** Fitting is parallelized across subject x visit pairs with `parallel::mclapply()`, and each pair seeds its random starting points from its own subject and visit ID (`string_seed()` in `pk_fit.R`). A single `set.seed()` before the parallel loop does not give reproducible results under `mclapply()` (identical code and seed drew differently on separate runs), and the `RNGkind("L'Ecuyer-CMRG")` workaround fixes reproducibility only for a fixed task order and core count. Per-ID seeding reproduces each subject's fit for any `N_CORES` and any row order of `subject_visits`.

**Adaptive retry.** Fits flagged after the standard pass get a denser retry: the prior estimate as a seed (so the retry cannot end worse), a systematic grid over the `kel` and `k_release` bound ranges when either bound flag is set, and 80 more random starts at a higher `maxit`. A solution pinned at a bound can hide a better one inside it, and a worse local optimum inside the bounds shows no boundary flag and appears as a low R² or a non-converged fit, so a subject x visit is retried when `kel_at_bound`, `k_release_at_bound` or `r2_12C_low` is TRUE or `converged` is FALSE. `r2_13C6_low` alone does not trigger a retry, in either model: a poor 13C6 fit usually reflects a structural mismatch, and 13C6 has its own `f_delayed_13C6`, so a persistently poor `r2_13C6` means the extra freedom found no good shape, and a denser search of the same structure is unlikely to help. The retry replaces the original only when its `objective_value` is strictly lower (comparable only within the same subject x visit), so retrying cannot make a result worse.

A flagged fit can reflect a weak curve or a search failure, and the per-model results tables (`fit_results_single_wave.csv`, `fit_results_two_wave.csv`) record which. `retried` is TRUE when the fit was flagged and a denser retry was attempted. `retry_improved` is TRUE when the retry found something strictly better, FALSE when it could not beat the original, and NA when the fit was never flagged. `retried = TRUE, retry_improved = FALSE` points to a low-information curve: for example, slow capsule opening makes the 13C6 shape depend on `k_release` more than on `ka` and `kel`, which leaves a broad, nearly flat region of the objective surface. `retried = TRUE, retry_improved = TRUE` means the original result was under-searched.

## Fit quality and what to trust

In the current run (`results/fit_results.csv`, 68 subject x visit fits) the median R² is 0.975 for the 12C curve (10th percentile 0.916, minimum 0.784) and 0.930 for the 13C6 curve (10th percentile 0.769): 41 of the 13C6 fits are at or above 0.90, 23 between 0.70 and 0.90, 4 below 0.70, and none below 0.30. Per-subject R² is in `results/fit_results.csv`. Check it before trusting an individual subject's parameters, and inspect that subject's plot in `results/plots_individual/` if R² is low.

`results/fit_results.csv` also carries `r2_12C_low` and `r2_13C6_low`, TRUE when that curve's R² falls below `R2_RELIABLE_MIN` (0.70, set in `02_fit_erie_model.R`). They are per-curve flags. A low flag on one curve says nothing by itself about the other curve or about the `ka` and `kel` shared by both, but because `ka` and `kel` are fitted jointly, a poor fit on one curve can still bias them, so a low flag is a prompt to inspect the subject's plot. In the current run the flag is set for 0/68 fits on `r2_12C` and 4/68 on `r2_13C6`.

**Interpretable as approximately absolute numbers:**
- `kel` and `ka` are reasonably well identified given the bounds and the multi-start search.
- Within-subject, paired comparisons (for example baseline vs. intervention `F`): a systematic `Vd` bias applies equally to both visits of a subject and cancels in a paired comparison.

**Interpretable only with caveats:**
- `F` on its own. `F` and `Vd` enter the model as the product `F * dose / Vd`, so the data identify only their ratio and cannot distinguish a small `F` with a small `Vd` from a large `F` with a large `Vd`. This is a structural limit of oral-only concentration data. An IV tracer dose in the same subjects or a measured individual `Vd` would resolve it, and this protocol has neither. Report `F` conditional on the Nadler `Vd`. Under that `Vd` the median `F_12C` and the median `F_13C6` are both about 1.9% in the current run (ranges 0.8-6.8% and 0.4-8.0%). These values are low and have not been compared with a literature fructose or glucose bioavailability at a comparable dose. Substantially higher literature values would indicate that the Nadler blood volume is too small for fructose's distribution volume (see "Volume of distribution").
- Any boundary-flagged parameter (`kel_at_bound` or `k_release_at_bound` = TRUE). After the retry the data may still leave the parameter at the bound, for example `k_release` at its ceiling because the first post-dose sample already shows near-peak tracer concentration and nothing argues for slower dissolution, which is a sampling-resolution limit. A flag that survives the retry (`retried = TRUE, retry_improved = FALSE`) is more informative than one from a single pass, but it cannot separate "unconstrained by the data" from "the search missed it", so inspect the plot either way.
- Any fit with `converged = FALSE`. The winning multi-start result did not meet `optim()`'s convergence criterion (it had the lowest objective among the seeds tried), typically because it reached the `maxit` cap.
- Extrapolated quantities (for example AUC beyond 360 min or the time to full clearance) are projections beyond the last observation.

## Output files

Generated by `02_fit_erie_model.R` into `results/` (git-ignored):

- `fit_results.csv`: one row per subject x visit from the model that won the selection (`model` = `single_wave` or `two_wave`), with both candidates' AIC (`aic_single_wave`, `aic_two_wave`). `f_delayed`, `t_lag` and `f_delayed_13C6` are filled only when `two_wave` was selected, and a `two_wave` fit that lost the comparison is set to NA. 13C6's mechanism columns (`t_lag1_13C6`, `f_delayed2_13C6`, `t_lag2_13C6`) can be filled under either 12C model. `k_release` is NA whenever 13C6 uses the onset-lag or second-wave mechanism, because both model the capsule contents as an instant bolus with no dissolution step.
- `fit_results_single_wave.csv`, `fit_results_two_wave.csv`: each candidate's full fit, including `retried`, `retry_improved`, `objective_value`, `rss_12C` and `n_12C` (what the AIC is computed from) and `t30_present`. `two_wave` rows are NA where it was never attempted (no t=30 sample, or no evidence of a second wave).
- `plots_individual/ER##_joint.pdf`: one figure per subject with both isotopes (rows) and both visits (columns), the observed points (including t=0), the selected model's fitted curve, each panel annotated with that curve's R², and the 12C panel also with both candidates' AIC. Each curve is drawn until it has cleared to 1% of its own peak (`CLEARANCE_FRAC`), which sets the x-axis limit.
- `r2_comparison_vs_former_model.csv`, `r2_comparison_summary.csv`, `r2_comparison_plot.pdf`: written by `04_compare_to_former_model.R` (see "R² compared with the former model").

`t30_present` describes the data and is independent of the selected model: without the t=30 sample `ka` is poorly anchored even under `single_wave`, so downstream reliability filtering should consider it alongside R², `converged` and the boundary flags. `scripts/03_diet_summary.R` reads `fit_results.csv`; see `docs/diet-summary.md`.

## Origin of the model

The model follows the final joint model of `former_models/MixedModel/Scripts/fructose_joint_model_final.R` (Melany's), and the joint version with shared `ka` and `kel` was chosen over fitting either curve alone and over Bas's original 3-compartment Julia model. In her analysis the joint fit raised the 13C6 median R² from about 0.70 (independent fits) to about 0.78 and left far fewer subjects below R² 0.3, which supports the shared-clearance hypothesis. The engine is split into generic files (`scripts/assets/pk_curves.R`, `pk_fit.R`, `pk_diagnostics.R`, meant to be copied into other PK projects) and the ERIE-specific `02_fit_erie_model.R`. The 13C6 dose is the confirmed administered 120 mg / 644.78 µmol (see `docs/problems-and-fixes.md`).

### R² compared with the former model

`scripts/04_compare_to_former_model.R` (`pixi run compare-former`) compares each subject x visit's classical R² (`1 - SSE/SS_tot` on the observed points) with the former model's saved results (`former_models/MixedModel/Results/fit_results_joint.csv`). In the current run the median 12C R² is 0.975 (former: 0.933), higher for 57 of 68 subject x visits and lower for 11; the median 13C6 R² is 0.930 (former: 0.793), higher for 58 and lower for 10. Per-subject numbers are in `results/r2_comparison_vs_former_model.csv`, the summary in `results/r2_comparison_summary.csv`, and a scatter plot in `results/r2_comparison_plot.pdf`. The correlation between former and current R² across subject x visits is modest (0.37 for both isotopes).

R² is a fair yardstick between the two: the objective is unweighted and normalized by each curve's own total variance, like the former model's, so both optimize essentially the quantity R² measures. R² does not penalize the extra parameters of `two_wave` and of 13C6's mechanisms, so a higher R² alone does not show a better model; for 12C the AIC-based selection accounts for them.

## References

- Hannou SA, Haslam DE, McKeown NM, Herman MA. *Fructose metabolism and metabolic disease.* J Clin Invest. 2018;128(2):545-555.
- Jang C, Hui S, Litchfield B, et al. *The small intestine converts dietary fructose into glucose and organic acids.* Cell Metabolism. 2018;27(2):351-361.
- van der Crabben SN, et al. *Relationship between glucose volume of distribution and the extracellular space: a multiple tracer study.* Metabolism. 2011.
- Nadler DA, Hidalgo JU, Bloch T. *Prediction of blood volume in normal human adults.* Surgery. 1962;51(2):224-232.
