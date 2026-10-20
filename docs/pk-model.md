# Fructose PK model specification

This documents the pharmacokinetic (PK) model of plasma fructose fitted by `scripts/02_fit_shared_kel.R`, the reasoning behind it, and what its output can and cannot be used to claim. The script starts from the cleaned data and fits every parameter itself, including the wave structure of each curve. The manuscript text is in `docs/methods.md`, and the history of problems and fixes is in `docs/problems-and-fixes.md`.

## Summary

**What is modelled.** For every subject x visit there is one plasma curve of unlabeled (12C) fructose after a drink of 1 g/kg body weight. It is described by one-compartment first-order oral kinetics with an absorption rate `ka`, an elimination rate `kel`, a bioavailable fraction `F` and a volume of distribution `Vd` equal to the measured extracellular water of that visit. The two visits of a subject are fitted together with one `kel` per subject and a separate `ka` and `F` per visit.

```
dose D --(ka, F)--> blood --kel--> cleared

concentration = amount in blood / Vd
```

**Formulas.** `t` is in minutes, doses in mg, `Vd` in litres and concentrations in mg/L.

```
bateman(t; D) = F * ka * D / (Vd * (ka - kel)) * (exp(-kel * t) - exp(-ka * t))      0 for t <= 0
  (when ka = kel: F * ka * D * t * exp(-ka * t) / Vd)
```

The plain model (`single_wave`) is `C(t) = bateman(t; D)`. Many curves have a post-peak dip that one absorption phase cannot produce, so `two_wave` splits the dose into an early part and a part that starts `t_lag` minutes later, with `f_delayed` the late fraction:

```
C(t) = bateman(t; (1 - f_delayed) * D) + bateman(t - t_lag; f_delayed * D)
```

Both waves share `ka`, `kel` and `F`.

**How the fit works.** `02_fit_shared_kel.R` fits the 33 subjects who completed both visits in two stages. Stage 1 fits each visit alone, with its own `kel`, and chooses its structure: `single_wave` always, `two_wave` when the raw data point to a second wave or `single_wave` fits less than well, and the lower AIC on that curve's residuals decides. Stage 2 fits both visits together with one `kel` per subject, keeps each visit's structure from stage 1 and re-estimates `ka`, `F`, `f_delayed` and `t_lag`. A subject has 5 parameters plus 2 for each `two_wave` visit. In the current fit 17 of the 66 curves use `single_wave` and 49 `two_wave`.

**Inputs that are not fitted.** `D` = 1000 mg x body weight in kg. `Vd` is the measured extracellular water of that visit in litres (see "Volume of distribution"). Concentrations are baseline-corrected (the fasted level is subtracted).

**Parameters.**

| Parameter | Meaning | Fitted in | Bounds | Explained in |
|---|---|---|---|---|
| `kel` | elimination rate constant (1/min) | shared fit, one per subject | [0.0099, 0.099] (t½ 7-70 min) | "Fitting bounds" |
| `ka` | absorption rate constant (1/min) | shared fit, one per visit | [1e-4, 1] | "Fitting bounds", "Why one kel per subject" |
| `F` | fraction of the dose that reaches the blood | shared fit, one per visit | [1e-4, 1] | "Fit quality and what to trust" |
| `f_delayed` | fraction of the dose in the late wave | stage 1 chooses the structure, stage 2 the value, per `two_wave` visit | [0.001, 0.999] | "The post-peak dip and the lagged-dose model" |
| `t_lag` | start of the late wave (min) | stage 1 chooses the structure, stage 2 the value, per `two_wave` visit | [5, 150]; 30 to 150 if t=30 is above 5 mg/L | "The post-peak dip and the lagged-dose model" |

**What to trust.** Every `ka` and `F` result is conditional on `kel` being the same at both visits: the fit attributes any change in elimination to `ka` or `F`. `F` is only interpretable within a subject because `F` and `Vd` enter the model only as `F * D / Vd`. Check `r2_b`, `r2_i` and `kel_at_bound` before using an individual fit. See "Why one kel per subject" and "Fit quality and what to trust".

## Background: what's being measured

Each fructose challenge test (FCT) samples plasma at t = 0, 30, 60, 90, 120, 150, 180, 240 and 360 min after a drink of 1 g/kg body weight of unlabeled fructose. Every subject has a test at FCT1 (coded `baseline`) and FCT2 (coded `intervention`, after the 4-week diet intervention). `visit` in the cleaned data and in every downstream script is `baseline`/`intervention`, and `FCT1`/`FCT2` appear only as raw-file labels. A capsule with 120 mg of 13C6-fructose tracer was taken with the drink; its curve is fitted afterwards with the rates of the 12C fit held fixed (see "The 13C6 tracer with fixed rates").

### What the model is fit to

Before any fit:
- Concentrations are converted from µmol/L to mg/L with the molecular weight `MW_12C` in `data/processed/erie_constants.csv`, so they match the mass unit of the dose.
- The subject's own fasted (t=0) level is subtracted from every observation of that test, since the dose is added on top of endogenous fructose.
- t=0 is excluded from fitting, because the corrected value is 0 by construction and the model predicts exactly 0 there for any parameters.
- Rows without a body weight (the intervention visits of the two study dropouts, ER33 and ER34) are dropped, and the shared fit uses only the 33 subjects with both visits.

## Why one kel per subject

With one curve per visit, `ka` and `kel` are interchangeable. Swapping them and multiplying `F` by `ka / kel` gives exactly the same curve:

```
C(t; ka, kel, F) = C(t; kel, ka, F * ka / kel)
```

The dose passes through two first-order steps in series, and the output of two first-order delays does not depend on their order. Identical predictions give identical residuals, R² and penalties, and `F` is free, so nothing in the objective prefers one of the two solutions. A single curve identifies a fast rate `max(ka, kel)` and a slow rate `min(ka, kel)` without saying which one is absorption. When absorption is the slower step (flip-flop kinetics), the tail of the curve carries `ka`. Slow gastric emptying of a large liquid load makes this plausible for fructose. The `two_wave` model reuses `ka`, `kel` and `F` in both waves and has the same symmetry.

The data add to this. Each curve has about eight samples after the dose, the first at 30 min, and when `ka` and `kel` are close the curve becomes `t * exp(-k * t)`, a single-rate curve. When each visit is fitted alone (stage 1 of `02_fit_shared_kel.R`), 42 of the 66 curves end with `ka` and `kel` equal within 0.1%, so for most curves only one rate is visible.

The study hypothesis is that the diet changes absorption and leaves elimination unchanged. The shared fit uses this assumption: with one `kel` for both visits, an exact swap would need the same `ka` at both visits, so the change in `ka` between the visits becomes identifiable. Which of the two rates is absorption stays open: in the simulations below, the shared fit still puts `ka` and `kel` in the wrong order in 26% to 69% of fits. The price is that any real change in `kel` is attributed to `ka` or `F`. A `ka` result from this fit has to be reported as conditional on a stable `kel`, and it cannot be used to argue that `kel` was stable.

The AIC compares the shared fit (stage 2) with the stage-1 fits of the same subject, with a residual variance per visit and one parameter fewer for the shared fit. The shared `kel` has the lower AIC in 25 of 33 subjects and is worse by more than 2 in 3 (ER15, ER17 and ER26). A model with one parameter fewer that absorbs a real `kel` change without a large loss of fit would give the same picture, so this agreement is weak evidence for a stable `kel`.

The simulations in `scripts/simulations/01_ka_kel_recovery.R` and `scripts/simulations/02_shared_kel.R` repeat the study design with known `ka` and `kel`: the 33 completers' sampling times, `Vd` and doses at both visits, `single_wave` truth, 20% between-subject variation in `ka`, `kel` and Cmax, and additive noise equal to the median residual SD of the `02_fit_shared_kel.R` fit. The results below were run at 4.5 mg/L; the current fit gives 4.3 mg/L. The intervention visit has a 1.5-fold change in `ka` or in `kel`, or no change. Every curve is refitted with each approach: each visit alone (01, 15 random starts, 100 replicates per scenario) and one `kel` per subject (02, 25 random starts, 40 replicates). `scripts/simulations/03_plot.R` draws the comparison in `results/02_diagnostics_ka_kel/simulation/`: `simulation_figure_ka.pdf` shows the scenarios the shared-`kel` model is built for (no change and a `ka` change), and `simulation_figure.pdf` adds the `kel` change that violates its assumption.

- **Each visit fitted alone.** A 1.5-fold `ka` change is detected in 71% to 84% of replicates, and its size is underestimated (median error in the log change -0.07 to -0.20, against a true change of 0.41). A 1.5-fold `kel` change is reported as a `ka` change in 80% to 87% of replicates. The fitted `ka` and `kel` are in the wrong order in 35% to 70% of fits. The sorted rates do work: a change in the fast rate gives `max(ka, kel)` a power of 0.87 to 0.91, and a change in the slow rate gives `min(ka, kel)` a power of 0.98 to 1.00.
- **One `kel` per subject.** A 1.5-fold `ka` change is detected in every replicate, with a median error in the log change of -0.03 to +0.05. With no change the `ka` test rejects in 7.5% and 15% of replicates (95% CI up to 29% with 40 replicates), against 5% to 6% for the fits of each visit alone.
- **When the assumption is violated.** With a 1.5-fold `kel` change, the shared fit reports the whole change as a `ka` change (median error +0.41 and +0.44 in the log change) and the `ka` test rejects in every replicate.

## The post-peak dip and the lagged-dose model

Many curves show a dip after an initial peak, then a second, often higher, rise before the final decline, most clearly in ER01 and ER09 at both visits. Over the 68 subject x visit curves in `data/processed/` (t=0 excluded, baseline-corrected), the raw-data checks described below fire as follows:

| Check | Curves |
|---|---|
| `has_peak_dip_rise()` | 43/68 (63%) |
| `has_near_peak_neighbor()` | 17/68 (25%) |
| either | 49/68 (72%) |

These checks are permissive gates that decide whether `two_wave` is attempted, and the true number of curves with a second wave is unknown. With the tail limit described below the dip check fires for the same 43 curves.

**A one-compartment curve cannot produce this shape.** A difference of two exponentials has one maximum and decreases monotonically after it, for any `ka`, `kel` and `F`. Absorption pathways that all start at t=0 and share one elimination rate converge to the same single-exponential tail, so summing them can reshape or delay a single peak and cannot create a dip followed by a second rise.

**A time-lagged second dose produces it** (`simulate_lagged_dose()` in `scripts/assets/pk_equations.R`). A fraction `f_delayed` of the dose contributes nothing until `t_lag` and then behaves like a fresh dose given at that time. The lag creates a gap between the two absorption phases, which yields a local minimum followed by a second peak. Physiologically this matches a second phase of gastric emptying: the duodenal feedback that paces emptying of a large caloric and osmotic load pauses and resumes on a timescale of tens of minutes.

`t_lag` is bounded above by 150 min. After 180 min the sampling gaps widen to 60 min (180-240) and 120 min (240-360), and a `t_lag` there would let the optimizer place a second wave where no sample can confirm or refute it.

The lower bound of `t_lag` depends on the subject's own concentration at t=30: 5 min when it is at or below `EARLY_LAG_OK_MGL` (5 mg/L), 30 min otherwise. The 0-30 min window contains no sample. A low floor for everyone lets the optimizer place an early onset delay there (small `t_lag`, `f_delayed` near 1, almost the whole dose late) even when the t=30 sample is already high. A floor of 30 min for everyone would exclude a real early delay in a subject whose t=30 concentration is near zero. The bound is a step because every value between 5 and 30 min lies in the same unsampled window, so an interpolated floor (21.6 min for ER04 baseline, t=30 = 51.4 mg/L) has no more support than none.

### How the dip is detected

The raw-data checks in `scripts/assets/curve_checks.R` use only a curve's observed points and involve no fitting. Fitting first and relying on AIC or R² lets a flexible model explain ordinary sampling noise around a single peak.

**`peak_dip_rise_info()`** flags a post-peak point that lies at least 15% (`min_rise_frac`) above the straight line joining its two neighbors. A single-exponential decline is convex, so a point on it lies at or below the chord between its neighbors, and a point clearly above the chord shows non-monotonic decay even when it never exceeds an earlier point. Design details:
- The first-wave candidate is the first point after which concentration turns down and that is at least half (`min_first_peak_frac`) of the curve's overall peak. A second wave is often taller than the first (ER01, ER09), so anchoring on the global maximum would miss it, and the half-of-peak floor keeps a small early wobble from counting as a first wave.
- At least one point must separate the first peak from any candidate second peak.
- When several points deviate, the one with the largest deviation is reported as `trigger_time`, with its size as `excess`. `trigger_time` anchors search seeds and `excess` feeds the "definite two_wave" override.
- Only points up to 180 min count (`DIP_MAX_TIME`), and only rises of at least 5% of the curve's peak (`DIP_MIN_PEAK_FRAC`). A second wave starts by `t_lag`'s 150-min bound and is visible by the 180-min sample, so a rise later than that cannot belong to an allowed wave. Near zero in the tail the line between two neighbors is close to 0 mg/L, and a few mg/L of noise would count as a rise of 100% or more without the 5% floor. In the simulations without these limits most false triggers sat in the tail (median 180 min, 8 mg/L); in the real data the dip check fires for the same 43 curves with and without them.

**`has_near_peak_neighbor()`** flags a curve where the point immediately before or after the peak is at least 85% of the peak's height. This catches two waves that overlap into a flat top or an irregular rise without a visible dip, such as ER02 baseline (`single_wave` R² 0.956, `two_wave` 0.998) and ER30 baseline (0.859 to 0.957). t=30 is excluded as a neighbor, since closeness to the peak at the first sample reflects fast absorption.

### Choosing between single_wave and two_wave

The choice is made in stage 1, for each visit on its own curve, so it does not depend on the shared-`kel` assumption: a second wave cannot be selected to absorb a `kel` change between the visits. `two_wave` is attempted when the t=30 sample exists and at least one of three conditions holds:
- `has_peak_dip_rise()` fires (with the tail limit above).
- `has_near_peak_neighbor()` fires.
- `single_wave` fits less than well (R² below `R2_ALWAYS_TRY_TWO_WAVE` = 0.95). Two waves close together can add up to one smoothly accelerating rise or a broad shoulder that neither detector sees. ER21 baseline shows it: the concentration at t=30/60/90 is 32/61/132 mg/L, `single_wave` gives R² 0.33 and `two_wave` 0.99. A curve that `single_wave` already fits well is left alone.

Without a t=30 sample, the 0-60 min window has only its two endpoints, and a 5-parameter fit would be under-identified. One curve lacks it. In the current fit `two_wave` was attempted on 53 of the 66 curves and selected on 49.

Selection is by AIC, `n ln(RSS/n) + 2k`, on the raw residuals of that curve (`k` = 3 for `single_wave`, 5 for `two_wave`). AIC is used without the small-sample correction: at n = 8, AICc raises the penalty gap between k = 3 and k = 5 from 4 to 28 and would exclude plateau cases such as ER02 and ER30, whose R² gains are large. Protection against overfitting comes from the raw-data gates and the bounds.

**Definite two_wave override.** When `excess` is at or above `EXCESS_DEFINITE_TWO_WAVE` (0.90), `two_wave` is selected without the AIC comparison, so that an unlucky local optimum cannot overrule unambiguous evidence. Only ER08 baseline (1.57, at 120 min) and ER11 baseline (1.01, at 150 min) reach it.

**Other rules.** `STRUCTURE_RULE` in `scripts/assets/model_12C.R` selects the rule: `"tail_limit"` (the one above), `"gated"` (no tail limit, `single_wave` R² threshold 0.50) and `"always_try"` (`two_wave` attempted on every curve with a t=30 sample, AIC only, no override). Compared with `gated`, `tail_limit` changes 5 curves to `two_wave`: ER17 baseline, ER21 intervention, ER28 at both visits and ER32 intervention, whose `single_wave` R² is 0.84 to 0.94. `always_try` adds `two_wave` to 14 curves, mostly ones `single_wave` already fits with R² 0.98 to 0.99, where the second wave may follow noise.

### Testing the choice by simulation

`scripts/simulations/04_structure_selection.R` simulates the study with a mix of `single_wave` and `two_wave` curves: the completers' sampling times, `Vd` and doses, `(f_delayed, t_lag)` drawn from 02's `two_wave` fits, `two_wave` in the share of curves 02 selects (and in 35% in one scenario), noise at 02's residual SD and the same `kel` at both visits. Each visit is fitted with both structures as in stage 1, the rules decide which one to keep, and stage 2 follows. `scripts/simulations/05_plot_structure_examples.R` draws example curves for each outcome (`structure_examples.pdf`). Results are in `results/02_diagnostics_ka_kel/structure_selection/`.

- **The structure choice is the main source of error in the `ka` change.** With the true structure given, the log `ka` change of a subject is off by more than 1 in 1% to 2% of subjects; with the `gated` rule in about 7%, and nearly all of these subjects have a wrong structure at one visit. A missed wave costs more than an added one (median error 0.63 against 0.31 in the first run).
- **Tail noise caused most added waves.** Without the tail limit, the definite override fired in 26% of simulated `single_wave` curves and 7% of `two_wave` curves, nearly always on a point after 150 min close to zero. The tail limit halves the added waves (34% to 18% of `single_wave` curves).
- **The gate missed waves.** Under `gated`, 21% of true second waves were missed, most because `two_wave` was never tried: neither check fired for 17% of simulated `two_wave` curves. The `single_wave` R² threshold of 0.95 lowers this to 17%; trying `two_wave` everywhere lowers it to 8% but adds waves to 30% to 38% of `single_wave` curves.
- **Overall.** In the quick run (5 replicates, 165 subjects per cell), `tail_limit` picks the right structure for 82% of curves (`gated` 75%, `always_try` 85%) and errs by more than 1 in the log `ka` change in 1% of subjects (`gated` 7%, `always_try` 4%); with second waves in 35% of curves, 80% and 4% (`always_try` 75% and 4%). The `ka` test across subjects behaves under every rule (median bias in the log `ka` change below 0.1).
- **Who is not fitted by the right structure in the real data is unknown.** About one in five curves gets the wrong structure in the simulation, so individual large `ka` changes deserve a look at the plot.

The simulated curves decay faster than many real curves and the simulated noise is the same at every concentration, so the tail problem is larger in the simulation than in the real data, where 5 of 43 dip triggers lie after 150 min.

## Volume of distribution (Vd)

`Vd` is the measured extracellular water (ECW, in litres) of each subject x visit, read from `data/ECW_ERIE.csv` in `scripts/01_clean_data.R`. It reaches the fits as the `vd_L` column of `erie_covariates.csv`. ER01 and ER03 have no baseline ECW measurement and use their intervention value (`ecw_source` = `other_visit`). Two entries look suspicious and are used as recorded: ER033 baseline (28.55 L in an 89.5 kg woman, 319 mL/kg) and ER009 (20.70 L at both visits). ER033 was lost to follow-up and is left out of the diet summary.

The ECF is an assumption about fructose's distribution volume. Fructose is small, freely water-soluble and unbound to protein, so it is expected to equilibrate into interstitial fluid as well as the vascular compartment. For glucose, van der Crabben et al. measured a `Vd` of 191-206 mL/kg across three tracers and showed that it equals the extracellular fluid space (ECFV, about 150-200 mL/kg) and exceeds blood or plasma volume. In this cohort `Vd` averages 228 mL/kg (range 159-319), above the ECFV range quoted for glucose. The choice of ECF has not been checked against a fructose-specific `Vd`, and the measurement method of the ECW is not recorded in the data.

`F` enters the model only as `F * D / Vd`, so `F` scales with the assumed `Vd`, while `ka` and `kel` do not depend on it. Within-subject comparisons of `F` cancel any systematic bias in `Vd`. The measured ECW differs between the two visits of a subject by a median of 4.2% (SD 1.2 L, 31 subjects with both values), so a paired `F` comparison also carries the visit-to-visit variation of the ECW measurement.

## Fitting bounds (why these numbers)

| Parameter | Bounds | Rationale |
|---|---|---|
| `kel` | [ln 2 / 70, ln 2 / 7] = [0.0099, 0.099] /min (t½ 7-70 min) | Smith et al. 1953 report a fructose half-life of about 20-25 min in normal subjects, longer in hepatic disease. The range is centered on that value on a log scale (geometric midpoint 22 min). |
| `ka` | [1e-4, 1] /min | Wide, because absorption-rate differences between visits are the research question. |
| `F` | [1e-4, 1] | A fraction of a dose. |

Two soft penalties are added to the objective. The first pulls the predicted Tmax above 10 min (`MIN_TMAX`, weight `TMAX_LAMBDA` = 50), which rules out a degenerate solution that absorbs and clears almost instantly and produces a sharp spike between the fasting sample and the first sample at 30 min. At a `two_wave` visit it also applies to each wave's own peak, because the combined maximum constrains only the taller wave. Both waves have the shape of the plain Bateman curve, so the floor on that shape covers each wave. The second pulls the predicted Cmax to within ±10% of the observed maximum (`CMAX_TOL`, weight `CMAX_LAMBDA` = 20), which discourages systematic over- or undershoot of the peak and leaves the rest of the curve free. Both are evaluated on a fine time grid, because degenerate solutions look fine at the sparse observed points, and both are smooth (zero within the limit, squared beyond it), because L-BFGS-B relies on finite-difference gradients.

## Fitting procedure

PK objective surfaces have several local minima, so every fit uses multi-start bounded optimization (`fit_multistart()` in `scripts/assets/optimizer.R`, `optim` with L-BFGS-B) and keeps the best solution.

**Objective.** Per curve, the residual sum of squares divided by the curve's total sum of squares (roughly 1 - R²), plus the two penalties. The shared fit sums this over the two visits of a subject. Residuals are unweighted: a proportional weighting `1 / pred²` is available in `build_objective()` and is switched off, because it discounts the peak region where early features of the curve live.

**Starting points.** Every fit starts from four fixed rate seeds and nine `ka` = `kel` diagonal seeds (`v` from 0.008 to 0.07), each paired with one wave seed for every `two_wave` visit, plus 30 random seeds (log-uniform for the rates). The wave seeds include `f_delayed` near 0 (which recovers a plain Bateman fit), dip-like seeds, and `t_lag` at the dip's `trigger_time` and 20 and 40 min before it, since the trigger point is where the second wave becomes visible, which can be after it starts. The diagonal seeds matter because the objective near `ka` = `kel` is narrow and generic seeds rarely land there. In stage 1, `two_wave` also starts from the `single_wave` fit combined with every wave seed. Stage 2 also starts from the stage-1 estimates: with the geometric mean of the two `kel` values, and with each visit's slower rate taken as `kel` and the faster rates as `ka`. Seeds outside the bounds are moved just inside them. `maxit` is 80.

**Retry.** A fit with `kel` within 1e-4 of a bound, a visit with R² below 0.70 or no convergence may sit in a worse local optimum. It is searched again from its estimate, a grid over `kel`'s range when `kel` is at a bound and 80 random seeds, with `maxit` 100, and replaced only by a fit with a lower objective. In the current fit 5 stage-1 fits were retried (4 improved) and 8 stage-2 fits (6 improved); only ER04's `ka` change moved noticeably.

**Parameter scaling.** `fit_multistart()` sets L-BFGS-B's `parscale` per start from that start's own magnitude (floored at 1% of the bound range). Bounds such as [1e-4, 1] for `ka` are wide relative to typical fitted values (0.01-0.05), and a scale set from the bound width leaves the search stuck near its start.

**Reproducibility.** Subjects run in parallel with `parallel::mclapply()`. Each fit seeds its random starts from the subject ID, the visit and the structure (`string_seed()`), so results do not depend on the number of cores or the row order.

## Fit quality and what to trust

In the current fit (`results/01_shared_kel/shared_kel_fits.csv`, 33 subjects, 66 curves) the median R² is 0.982 (10th percentile 0.938, minimum 0.794). No curve is below 0.70 and 2 are below 0.90. The shared `kel` has a median of 0.025 /min (t½ 28 min, interquartile range 0.020-0.031). ER14 and ER15 have their `kel` at the upper bound (0.099), so 31 subjects have a reliable fit at both visits (R² ≥ 0.70 and `kel` away from both bounds), which is the inclusion rule of the diet summary. In 33 of the 66 curves the fitted `ka` is below `kel`. The median `F` is 0.067 (range 0.017-0.42).

**Interpretable within a subject:**
- The change in `ka` between visits, conditional on a stable `kel`.
- The change in `F` between visits: a systematic `Vd` bias applies equally to both visits and cancels.
- Curve metrics derived from the fitted curve (Cmax, Tmax, AUC), which do not depend on how the rates are labeled.

**Interpretable only with caveats:**
- `ka` and `kel` as absolute values. Without the shared-`kel` assumption the curve identifies only the fast and the slow rate (see "Why one kel per subject").
- `F` on its own. The data identify only `F / Vd`, so a small `F` with a small `Vd` cannot be told from a large `F` with a large `Vd`. An IV tracer dose or a measured fructose `Vd` would resolve it, and this protocol has neither. Report `F` conditional on `Vd` = measured ECW.
- A fit with `kel_at_bound` = TRUE. The data may push `kel` past the allowed half-life range, or the bound may hide a better solution, and the plot is needed to tell which.
- Extrapolated quantities, such as AUC beyond 360 min, which are projections beyond the last observation.

### Known limitations and open questions

- **Two shared fits stop at `maxit` after the retry.** ER11 and ER22 (`converged` = FALSE); their retry found no better optimum.
- **About one curve in five may have the wrong structure** (see "Testing the choice by simulation"). A large `ka` change in a single subject is more likely when the structure of one visit is wrong.
- **A wave's peak can land late.** Nothing limits the peak time of a wave beyond `t_lag`'s 150 min cap, and a slow `ka` can place the peak well after `t_lag` (ER16 intervention: `t_lag` 111 min). An upper bound is worth adding if a fitted wave peaks inside the 180-240 or 240-360 min gaps.

## The 13C6 tracer with fixed rates

`scripts/04_fit_13C6.R` fits the 13C6 curve of each completer x visit with `ka` (that visit) and `kel` (shared per subject) held at the estimates of `02_fit_shared_kel.R`. The tracer and the unlabeled fructose are the same molecule and are cleared by the same pathways, so they share `kel`, and the tracer is assumed to be absorbed at the 12C rate once it is released. When 12C is `two_wave`, its `t_lag` from `02_fit_shared_kel.R` is available to 13C6. 13C6 is not baseline-corrected (there is no endogenous tracer), t=0 is excluded, and the dose is 120 mg.

Only the parameters specific to the capsule are fitted. Three mechanisms are compared, as in the archived joint model (`scripts/archive/joint_12C_13C6/02_fit_erie_model.R`), whose history is in `docs/problems-and-fixes.md`:

| Mechanism | Fitted parameters | K | Tried when |
|---|---|---|---|
| `release` | `F_13C6`, `k_release` [0.001, 1] (capsule dissolution, first order); with a `two_wave` 12C also `f_delayed_13C6` [0.001, 0.999] at 12C's `t_lag` | 2 or 3 | always |
| `onset_lag` | `F_13C6`, `t_lag1_13C6` [0, 90] (instant bolus that starts late) | 2 | the first sample is below 25% of the peak (`has_onset_lag_evidence()`) |
| `second_wave` | `F_13C6`, `f_delayed2_13C6` [0.05, 0.95], `t_lag2_13C6` [30, 150] (its own late wave, instant bolus) | 3 | `has_peak_dip_rise()` or `has_near_peak_neighbor()` fires on 13C6's own data |

`second_wave` is rejected without comparison when `f_delayed2_13C6` ends within 1e-3 of its lower bound. The mechanism with the lowest AIC on 13C6's raw residuals is kept. The objective, penalties and multi-start settings are those of the 12C fit (Tmax floor 10 min, Cmax band ±10%, `maxit` 60).

In the current run (66 curves) 13C6 uses `release` in 26 curves, `onset_lag` in 21 and `second_wave` in 19. The median R² is 0.93 (minimum 0.07), and 7 curves are below 0.70. ER22 baseline (R² 0.07) peaks at 30 min and halves by 60 min, faster than the fixed 12C rates of that visit allow. `k_release` sits at its upper bound in 7 curves, which makes the release effectively instant. The median `F_13C6` is 0.074 against 0.067 for `F_12C`, and the ratio `F_13C6 / F_12C` has a median of 1.11 (10th-90th percentile 0.67-1.60).

A 13C6 curve that the fixed rates fit poorly is evidence against the assumption that the tracer follows the 12C kinetics in that visit, or against the 12C rates themselves. The fit cannot tell these apart.

## Output files

Generated by `02_fit_shared_kel.R` into `results/01_shared_kel/` (git-ignored):

- `visit_fits.csv`: stage 1, one row per subject x visit: the selected structure, `ka`, `kel`, `F`, `f_delayed`, `t_lag`, R², `kel_at_bound`, `converged`, `retried`, `retry_improved`, whether `two_wave` was tried and the definite override fired, the AIC of both structures, RSS, number of points, `diet` and `structure_rule`.
- `shared_kel_fits.csv`: stage 2, one row per subject: the structure of each visit (`model_b`, `model_i`), `kel_shared`, `ka_b`, `ka_i`, `F_b`, `F_i`, the wave parameters of each `two_wave` visit (`f_delayed_b`, `t_lag_b`, `f_delayed_i`, `t_lag_i`), the R² of each visit (`r2_b`, `r2_i`), `kel_at_bound` (TRUE at either bound), `converged`, `retried`, `retry_improved`, the stage-1 estimates (`ka_b_indep`, `kel_b_indep`, ...), the RSS and number of points per visit, the summed RSS of the shared and the stage-1 fits, `delta_aic` (shared minus stage 1; negative favors the shared `kel`), `ka_change_log`, `diet` and `structure_rule`.
- `plots_individual/`: one figure per subject with both visits, the shared fit and the stage-1 fit of each visit, and the stage-1 AIC of both structures.
- `shared_kel_wilcoxon.csv`: paired Wilcoxon on the log ratio of intervention to baseline, per diet arm, for `ka` from the shared fit and `ka` and `kel` from the stage-1 fits.

`03_diet_summary_shared_kel.R` reads these and writes the diet-arm tables and plots to `results/01_shared_kel/diet_summary/`. It leaves `kel` out of every table, test and boxplot, because it is assumed equal at both visits.

Generated by `04_fit_13C6.R` into `results/03_13C6_fixed_rates/`: `fits_13C6.csv` (one row per subject x visit with the fixed `ka`, `kel` and 12C `t_lag`, `F_12C`, the selected `mechanism_13C6`, its parameters, `r2_13C6`, `k_release_at_bound` and the AIC of each candidate) and `plots_individual/`, one figure per subject with both isotopes at both visits.

`05_diet_summary_13C6.R` reads `fits_13C6.csv` and writes the diet-arm tables and plots to `results/03_13C6_fixed_rates/diet_summary/`, with the statistics of 03 (paired Wilcoxon and the random-intercept mixed model). A curve is included when its 13C6 R² is at least 0.70 and the 12C fit that supplies its rates passes the inclusion rule of 03 (55 of 66 curves; 7 fail on 13C6, 4 on 12C). The outcomes are `F_13C6`, the ratio `F_13C6 / F_12C` within a visit (in which `Vd` cancels), and Cmax, Tmax and AUC of the fitted 13C6 curve. `ka` and `kel` are inputs from 12C and are not tested again. The 13C6 mechanism is counted per diet arm and visit, with Fisher's exact test on diet x mechanism. Before plotting, the script recomputes every R² from the curve it draws and stops if one differs from `fits_13C6.csv`.

## Origin of the model

The model follows the final joint model of `former_models/MixedModel/Scripts/fructose_joint_model_final.R` (Melany's), which fitted the 12C curve and the 13C6 tracer curve together with shared `ka` and `kel`. That joint 12C/13C6 version and the independent fructose fit are archived in `scripts/archive/`, and the current analysis fits the 12C curve first and the 13C6 curve afterwards with the 12C rates held fixed. The model code is in `scripts/assets/`: the PK equations, the optimizer, the raw-data checks, one file per model (`model_12C.R`, `model_13C6.R`) and the statistics. The numbered scripts read the data, run a model and write the output.

## References

- Smith LH Jr, Ettinger RH, Seligson D, Lightcap PL. *A comparison of the metabolism of fructose and glucose in hepatic disease and diabetes mellitus.* J Clin Invest. 1953;32(4):273-282. doi:10.1172/JCI102736
- van der Crabben SN, et al. *Relationship between glucose volume of distribution and the extracellular space: a multiple tracer study.* Metabolism. 2011.
