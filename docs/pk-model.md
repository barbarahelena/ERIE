# Fructose PK model specification

This documents the pharmacokinetic (PK) model fit by `scripts/02_fit_erie_model.R`,
the biological and statistical reasoning behind it, and - most importantly -
what its output can and cannot be used to claim. It consolidates and updates
the reasoning previously scattered across `former_models/*/*.docx`, after
independently re-deriving and re-validating it (see "Continuity with prior
work" below).

## Background: what's being measured

Each fructose challenge test (FCT) gives two simultaneous plasma curves from
the same blood draws, at t = 30, 60, 90, 120, 150, 180, 240, 360 min after a
combined dose:

- **12C-fructose**: 1 g/kg body weight of unlabeled fructose, drunk in water
  (an individualized, large dose).
- **13C6-fructose**: a fixed tracer dose in an enteric capsule (120 mg /
  644.78 µmol, confirmed administered dose; see `docs/data-cleaning-notes.md`
  for the dose discrepancy this project resolved).

Every subject has this pair of curves at both FCT1 (coded `baseline`) and
FCT2 (coded `intervention`, post 4-week diet intervention) - `visit` in the
cleaned data and everywhere downstream is `baseline`/`intervention`, not the
raw `FCT1`/`FCT2` labels.

## Model structure

Both curves are modeled as one-compartment first-order oral absorption
(gut/absorption pool -> sampled central compartment), but the capsule route
gets an extra upstream compartment for gradual dissolution:

**12C (liquid dose)**
```
dA1_12C/dt = -ka * A1_12C                          A1_12C(0) = dose_12C
dA2_12C/dt =  F_12C * ka * A1_12C - kel * A2_12C    A2_12C(0) = 0
C_12C(t)   =  A2_12C(t) / Vd
```

**13C6 (capsule dose)**
```
dRelease/dt <- -k_release * Release                    Release(0) = dose_13C6
dA1_13C6/dt <-  k_release * Release - ka * A1_13C6      A1_13C6(0) = 0
dA2_13C6/dt <-  F_13C6 * ka * A1_13C6 - kel * A2_13C6   A2_13C6(0) = 0
C_13C6(t)   =  A2_13C6(t) / Vd
```

| Parameter | Meaning | Shared or separate across the two curves? |
|---|---|---|
| `ka` | absorption rate constant | **Shared** |
| `kel` | peripheral clearance rate constant | **Shared** |
| `F_12C`, `F_13C6` | bioavailable fraction of the dose | Separate |
| `k_release` | capsule dissolution rate constant | 13C6 only (no liquid-dose equivalent) |
| `Vd` | volume of distribution | Not fit - computed per subject x visit from weight/height/sex (see below) |

### Why share `ka`/`kel` across the two curves

The working hypothesis: once fructose is in the gut, it's absorbed by the
same transporters and cleared from blood by the same physiological
machinery regardless of how it got there (liquid vs. capsule) - so `ka` and
`kel` should be the same for both curves in a given subject/visit. Fitting
the two curves *independently* (no sharing) produced poor cross-curve
agreement in `kel` for many subjects, which motivated testing this
hypothesis directly rather than assuming either the sharing or the
disagreement.

### Why *not* share `F` across the two curves

The two doses differ by roughly 1000-fold (1 g/kg liquid vs. 120 mg fixed
capsule). Intestinal first-pass fructose clearance is dose-dependent: at low
doses the small intestine can metabolize the large majority of fructose
before it reaches the liver, but this capacity saturates once intake
exceeds roughly 1 g/kg (Jang et al. 2018) - which is exactly the 12C dose
used here. So the two curves are expected to see *different* degrees of
intestinal saturation, and therefore different `F`, even under the
shared-`ka`/`kel` hypothesis. Treating a `F_12C` vs `F_13C6` difference as
informative (not noise to reconcile) is a deliberate modeling choice, not an
oversight.

### Why the capsule gets a dissolution step

Assuming the full capsule dose is instantly available at t=0 implicitly
treats capsule dissolution as infinitely fast. Adding an explicit
dissolution compartment (`k_release`) relaxes that assumption and lets the
data inform how fast the capsule actually releases its contents - this
consistently improved 13C6 fit quality over assuming instant availability
(see "Fit quality" below), which is itself evidence that gradual dissolution
is real, not just an extra free parameter absorbing noise.

## The post-peak dip and the lagged-dose model

A meaningful minority of curves (13/68 = 19.1%, by `has_peak_dip_rise()` in
`pk_diagnostics.R` - a threshold-free check for a post-peak local trough
followed by a rise of >=15% of that curve's own peak, requiring the
candidate first peak itself to be a substantial fraction of the overall
peak; see "Choosing between single_wave and two_wave" below, where this
same function gates whether `two_wave` is even attempted) show a dip after
an initial peak, then a second, often higher, rise before the final decline
- most clearly on `ER01` and `ER09` (every one of their four curves), and
at lower relative amplitude (a flattening rather than a full second rise)
on several more, including `ER03`. A smaller-amplitude
version of the same signature - a step in the decline much flatter than its
neighbors - shows up on ~31% of curves. Both patterns cluster tightly in
time (85% of the larger dips fall at exactly t=60 or t=90), which is too
systematic to be per-subject noise and argues for a real, shared
physiological trigger rather than assay artifact.

**No parameter choice in the base model (or any variant that keeps a single
shared elimination compartment) can produce this shape.** A one-compartment
absorption curve is monotonic after its single peak for any `ka`/`kel`/`F` -
provably so, not just empirically: it's a two-term exponential difference
with no room for a second local maximum. This was tested exhaustively before
concluding it's a structural limit, not a search failure: thousands of
random parameter draws (and several hand-engineered attempts) across a
"biphasic gastric emptying" variant - two gastric pools emptying at
different *rates* but both from t=0, feeding a shared absorption
compartment - never produced a second hump, and the reason generalizes:
*any* number of absorption pathways that all still start at t=0 and share
one terminal elimination rate converge to the same single-exponential decay
tail, so summing them can reshape or delay a single peak but never create a
genuine local minimum followed by a second rise.

**What does work: a genuinely time-*lagged* second dose**
(`simulate_lagged_dose()` in `scripts/assets/pk_curves.R`), not a second
pool active from t=0 at a different rate. A fraction `f_delayed` of a dose
contributes nothing until `t_lag`, then behaves exactly like a fresh dose
given at that later time - mathematically equivalent to two separately-timed
doses. This breaks the shared-tail constraint above because the delayed
portion has a genuine gap, not just a different decay rate, and was
confirmed both analytically (a lag creates an actual discontinuity in the
delivery rate) and numerically (produces a real local minimum followed by a
second, higher peak) before being wired into the fitting pipeline. Piloted
on `ER01`/`ER03`/`ER09` (`results/model-post-peak-dip-pilot/`) with a clear,
visually-confirmed win on the genuine dip cases (`ER01`, `ER09`) - both
curves transformed from a single smooth hump to a real double-humped fit
tracking the dip - and degrades gracefully (`f_delayed` near 0) on subjects
that don't need it.

`02_fit_erie_model.R` fits two candidate models per subject x visit,
`fit_subject_visit_single_wave()` (the plain joint delayed-release model,
formerly called "baseline") and `fit_subject_visit_two_wave()` (this
lagged-dose extension, formerly "lagged"), and picks one - see "Choosing
between single_wave and two_wave" below - rather than reporting both side
by side. Both candidates' own fits are still written out in full
(`results/fit_results_single_wave.csv`, `results/fit_results_two_wave.csv`)
alongside the selected combination (`results/fit_results.csv`).

`t_lag`'s upper bound is 150 min - past 180, sampling widens to 60-
(180->240) and 120-min (240->360) gaps, and a `t_lag` landing in one lets
the optimizer place an entire invented second wave where no sample can
confirm or refute it: found on `ER12` FCT2, `t_lag`=179min, `r2_13C6`
dropped from single_wave's 0.95 to 0.75, a large invented peak with zero
supporting data around t=200-220.

`t_lag`'s **lower** bound is not a fixed constant - it depends on that
subject's own observed 12C concentration at t=30 (`EARLY_LAG_OK_MGL = 5`,
`EARLY_LAG_BAD_MGL = 75` in `02_fit_erie_model.R`, linearly interpolated: 5
min when t=30 is at or below 5 mg/L, rising to a 30 min floor once t=30 is
at or above 75 mg/L). A flat low bound (5 min for everyone) let the
optimizer fake an early *onset* delay (small `t_lag`, `f_delayed` near 1 -
"almost the whole dose is late") in the same always-unsampled 0-30min gap,
regardless of whether the subject's own data gave any reason to believe
something was delayed that early: found on `ER06` FCT2, which converged to
`t_lag`=19.7min, `f_delayed`=0.999 even though its own t=30 sample was
already near its peak (~80 mg/L, well above `EARLY_LAG_BAD_MGL`). But a
flat *high* bound (30 min for everyone) would equally wrongly forbid a
genuine early delay for a subject whose t=30 really is near 0 - a real,
separately-observed phenomenon - so the bound is computed per subject
instead of set once for the whole cohort. Refit with this bound, `ER06`
FCT2 converged to a qualitatively different, more moderate solution
(`t_lag`=35.3min, `f_delayed`=0.320, not just clamped to the new floor with
the same extreme `f_delayed`), while the genuine-dip subjects (`ER01`,
`ER09`, whose own t=30 values are also well above `EARLY_LAG_OK_MGL`) were
unaffected - their real `t_lag` (85-114min) was never near the floor to
begin with.

### Choosing between single_wave and two_wave

Selection is by AIC, `n*ln(RSS/n) + 2k`, computed on **12C's own raw
residuals only** (`K_12C_SINGLE_WAVE = 3`: `ka`, `kel`, `F_12C`; vs.
`K_12C_TWO_WAVE = 5`: + `f_delayed`, `t_lag`) - not a combined-curve AIC.
13C6 is fit separately once 12C's structure is already decided (see below),
so basing selection on 12C alone keeps the question "does 12C's own curve
earn a second wave" from being answered by 13C6's fit quality instead, and
avoids the cross-curve concentration-scale mixing problem a combined AIC
would have.

Two rules override a pure AIC comparison:
- **`two_wave` is never attempted when 12C's t=30 sample is missing.**
  Without it, the 0-60min window has only its t=0/t=60 endpoints to anchor a
  5-parameter fit - under-identified, and AIC would otherwise rubber-stamp
  whatever shape 5 sparse points can fit almost exactly. Found via `ER05`
  FCT1 (t=30, 90, 150, 240 all missing - only 5 real points for 12C).
- **`two_wave` is never attempted unless 12C's own raw observations show
  EITHER a peak-dip-rise pattern OR a plateau/near-equal peak neighbor.**
  Fitting `two_wave` first and relying on AIC to notice afterward that
  there was nothing to explain lets a flexible model produce a
  statistically-plausible-*looking* fit for a curve that never had a
  second wave at all - checking the raw data directly catches this before
  it happens. Two complementary checks, both in `pk_diagnostics.R`:
  - `has_peak_dip_rise()`: a post-peak local trough followed by a
    subsequent rise of >=15% of that curve's own peak (the same
    threshold-free check used to characterize dip prevalence across the
    cohort, at the top of this section).
  - `has_near_peak_neighbor()`: the point immediately before or after the
    curve's own peak is at least 85% of the peak's height - catching two
    waves that overlap closely enough in time to blend into a flat top or
    an irregular rise WITHOUT ever producing a visible dip, which
    `has_peak_dip_rise()` structurally cannot see. Found via `ER02` FCT1 (a
    genuine flat top, R² 0.955->0.997), `ER30` FCT1 (an irregular
    pre-peak rise, R² 0.859->0.957) - both real, well-bounded improvements
    the dip-only check missed entirely. Deliberately excludes t=30 as a
    valid "near" neighbor: being close to the eventual peak at the very
    first sample just reflects ordinary fast absorption (there's no
    earlier sample to show a genuinely different pre-peak trajectory), not
    evidence of a second wave - found on `ER06` FCT2, whose only "near-peak"
    point was t=30, and whose `two_wave` fit consistently pinned `t_lag`
    exactly at its own lower bound rather than settling in the interior
    the way genuine cases do (`ER02`'s 51.5min, `ER30`'s 75.7min) - a sign
    the fit still wanted to exploit the unsampled 0-30min gap, just
    clipped at the boundary, not a freely-preferred timing.

**Added: a third path in, alongside the two raw-data checks above -
`single_wave`'s own R² being very poor (`< R2_ALWAYS_TRY_TWO_WAVE = 0.50`)
also triggers a `two_wave` attempt, even with neither dip nor plateau
evidence.** Two waves close enough together in time - physiologically
plausible here: the duodenal-brake/enterogastric feedback that paces
gastric emptying of a large caloric/osmotic load (this study's 1 g/kg
fructose dose) typically pauses and resumes on a tens-of-minutes
timescale - can compound into a single, smoothly-ACCELERATING rise with no
down-turn and no near-peak plateau at all, a shape neither detector is
built to catch (both require the data to already show some sign of
non-monotonic or blended behavior). Confirmed concretely on `ER21` FCT1:
raw 12C at t=30/60/90 is 32/61/132 mg/L, a sudden more-than-doubling with
no down-turn beforehand - `single_wave` R²=0.31, and neither gate fires.
An unconstrained `two_wave` search (bypassing the gate entirely, for
diagnosis) found R²=0.99 at `t_lag`=57min, `f_delayed`=0.80 - essentially
a perfect fit. Threshold set well below `R2_RELIABLE_MIN` (0.70) so this
only fires when `single_wave` doesn't fit AT ALL, not as a general
substitute for the two raw-data checks.

**Fixed: a `parscale` bug in `fit_multistart()` (`pk_fit.R`) that was
causing genuinely missed optima, not just failed-to-attempt cases.**
`optim(method = "L-BFGS-B")`'s internal step sizing assumes roughly
unit-scaled parameters; the previous default (`parscale = upper - lower`)
is a poor proxy for that whenever a parameter's bounds are deliberately
wide open relative to its typical fitted value - `ka`/`F_12C` are bounded
`[1e-4, 1]` (to not bias the research question) even though real values
sit around 0.01-0.05, making their parscale ~20-100x too large. Confirmed
concretely on `ER11` FCT1: production's own objective (including its
`MIN_TMAX`/`CMAX` penalties) scored the *reported* fit (`t_lag`=47.3) at
0.208, while a search using each start's OWN magnitude as parscale found
0.188 at `t_lag`=120 - which lines up with the actual post-peak dip
evidence (`trigger_time`=150) that justified attempting `two_wave` in the
first place; the old bound-width parscale reproduced the same stuck
behavior as no parscale at all. Fixed by computing parscale PER START from
that start's own magnitude (floored at 1% of the bound range) rather than
once globally from the bounds - a multistart seed is already meant to be
in the right neighborhood (informed by a pre-fit, evidence anchor, or a
systematic grid), so its own magnitude is a far better local-scale
estimate than the global bound width.

**Removed:** the earlier third rule (never attempt `two_wave` when
`single_wave` already fits with R² > 0.95) turned out not to be doing the
actual protective work its rationale implied. It wasn't even the thing
blocking the two cases above - `ER30` FCT1's `single_wave` R² was 0.846,
nowhere near the 0.95 threshold, yet it was still being missed by the
(then dip-only) raw-data gate, not by this rule. AIC alone, combined with
the two rules above and the degenerate-fit exploits already bounded out
(`t_lag`'s data-dependent lower bound, per-wave `MIN_TMAX`, `ka`=`kel`
diagonal seeding), is trusted to judge whether the extra complexity earns
its keep.

**Tried and reverted: AICc.** With only 8-9 points per curve and `k`
jumping from 3 to 5 between candidates, plain AIC's complexity penalty
(`2k`) looked too weak to lean on once `two_wave` started being attempted
much more broadly (both gates above, no R² pre-filter). AICc (the standard
small-sample correction, `+2k(k+1)/(n-k-1)`) was tried as a fix - but at
n=8 it raises the complexity gap between k=3 and k=5 from 4 points to 28,
which is *so* aggressive it re-excluded all three of the confirmed-genuine
cases above (`ER02`, `ER06` FCT2 in the k_release-fixed sense, `ER30` - all
lost to `single_wave` under AICc despite large, independently-verified R²
improvements). Overcorrecting for small n isn't the right lever here; the
raw-data gates and the bounded-out exploits are what's actually supposed
to be doing the overfitting-protection work.

### 13C6's own wave choice

13C6 does not inherit 12C's `f_delayed`. It gets its own free parameter,
`f_delayed_13C6`, fit in a second stage once 12C's `ka`/`kel`/`t_lag` are
already fixed from stage 1 - so 13C6 independently decides whether it rode
only the first wave (`f_delayed_13C6` near 0), only the second
(near 1), or split across both, sharing only `t_lag` (the timing of gastric
emptying's second wave, a systemic event) with 12C. An earlier version
forced 13C6 to inherit 12C's `f_delayed` directly, which failed concretely
on `ER03` FCT1: 12C genuinely has two waves (R² = 0.99), but 13C6's own
points show one early peak decaying monotonically with nothing at the time
the inherited second wave would place one - forcing 12C's fraction onto it
gave `r2_13C6` = 0.075. The capsule's contents don't have to split across
both waves in the same proportion as the much larger liquid 12C dose; it's
biologically plausible the capsule emptied entirely in one wave or the
other.

**Fixed:** `MIN_TMAX` now checks each wave's OWN peak (`tmax_wave1`,
`tmax_wave2`, computed directly by simulating each wave alone rather than
only the combined curve `simulate_lagged_dose()` returns), not just the
combined curve's single global maximum. The old combined-only check
(`tmax12 <- FINE_T[which.max(fine12)]`) only ever constrained whichever
wave happened to be taller, leaving the other one's timing completely
unconstrained - confirmed on real fitted results before this fix: `ER09`
FCT1/FCT2 (wave 1 taller - combined peak at t=38/33 - but wave 2 ALONE
peaked at t=151/147, never checked) and `ER16` FCT2 (wave 2 alone peaked at
t=162). This addresses the *lower* floor (a wave's peak landing
implausibly early); it does not add any *upper* constraint on how late a
wave's peak can land beyond `t_lag`'s own 150min cap plus however long
`ka`/`kel` take to rise from there - whether that's also needed is a
separate, still-open question (see below).

**Known limitation, not yet addressed (TODO):**
- **The multistart search can miss a substantially better nearby optimum
  when `ka` and `kel` are close to each other** (a known hard region for
  one-compartment models - "flip-flop" kinetics, where `bateman_conc`'s
  `ka/(ka-kel)` term makes the objective surface narrow/awkward for
  gradient-based search, and generic seeds don't reliably land near it).
  Confirmed concretely on `ER06` FCT2 (`single_wave`): the pipeline reports
  `ka`=0.016, `kel`=0.051, R2=0.945 (neither at a bound), but a direct,
  independently-verified check (multiple `optim()` starts, all converging
  cleanly to the same point) found `ka`=0.023, `kel`=0.025, R2=0.986 -
  clearly better, not noise. This isn't specific to `two_wave` or the
  gates discussed above - it's in the base `fit_curve_independent()`/joint
  multistart machinery every fit in the cohort depends on, so other
  subjects' reported "best" fits may also be suboptimal local optima, not
  true optima, which could ripple into R2-dependent decisions (the
  R2_SINGLE_WAVE_SKIP_TWO_WAVE gate, AIC comparisons). Needs investigating
  before trusting cohort-wide R2 values at face value - possible fixes
  include a denser seed grid specifically along the ka=kel diagonal, or
  reparameterizing to avoid the `ka/(ka-kel)` singularity directly.
  Deliberately not fixed yet - the post-peak-dip/two_wave gate work in
  progress should be finished first, then return to this as its own,
  separate investigation.
- **No upper constraint on how late a wave's own peak can land**, beyond
  `t_lag`'s 150min hard cap - a slow `ka` can still push a wave's actual
  peak well past `t_lag` itself (e.g. `ER16` FCT2's wave 2 above, `t_lag` =
  111min but its own peak at t=162min). t=162 itself still falls inside the
  densely-sampled (30min-spaced) 0-180min window, so it's not clearly
  unsupported the way a peak inside the 180-240/240-360min gaps would be -
  worth revisiting if a fitted wave's peak is ever found landing inside one
  of those wider gaps specifically, rather than adding a blanket upper
  bound pre-emptively.
**Fixed:** 13C6 now gets its own onset-delay option (`t_lag1_13C6`),
independent of 12C's wave choice. `simulate_delayed_release()`'s smooth
first-order dissolution (fastest right at t=0) couldn't represent a curve
whose 13C6 signal is genuinely near-zero for the first ~30min and then
rises sharply - the optimizer was forced to compromise, letting some
release happen early (overshooting the real early timepoints) to still
reach the observed later peak. Implemented via `simulate_two_lag_dose()`
(already in `pk_curves.R`) with `f_delayed` fixed at 0, which reduces it to
a pure onset-delayed single wave - the entire dose simulated as normal,
just starting `t_lag1_13C6` minutes late, with no risk of the "near-total
delay disguised as a second wave" failure mode `t_lag`/`f_delayed` had for
12C (there's no second, undelayed portion here to hide behind). Gated by
`has_onset_lag_evidence()` (a different phenomenon from 12C's
`has_peak_dip_rise()`: the FIRST observed point being a small fraction of
the curve's own peak, not a post-peak deviation), and only kept if it beats
the no-lag fit on AIC (K=2 vs K=3, computed on 13C6's own residuals, with
`ka`/`kel` fixed at the joint fit's values - never touches 12C's own
reported parameters). Validated on the two originally-confirmed subjects
plus a third found along the way: `ER25` FCT1 (R² 0.71->0.94, `t_lag1_13C6`
= 25min), `ER06` FCT2 (R² 0.64->0.94, 28min), `ER25` FCT2 (R² ->0.75,
27min, no prior baseline - not part of the original 2-subject validation).

**Fixed:** the onset-lag candidate no longer carries `k_release` at all.
It was pinned exactly at its upper bound in all 3 originally-validated
cases (`ER06` FCT2, `ER25` FCT1, `ER25` FCT2 - 3 for 3, not a coincidence)
once `t_lag1_13C6` was in the model - 30min sampling can't distinguish
"fast dissolution" from "instant" once the delay itself already explains
the flat start, so it was just an unidentifiable, boundary-pinned nuisance
parameter with no fit-quality benefit. Now modeled as an instant bolus at
`t_lag1_13C6` (plain `bateman_conc`, no capsule-release step at all),
K=2 instead of 3, compared fairly against the no-lag candidate (also K=2).

**Fixed:** 13C6 now also gets its own separate-second-hump option
(`f_delayed2_13C6`/`t_lag2_13C6`), gated on 13C6's own raw data by the same
`has_peak_dip_rise()`/`has_near_peak_neighbor()` detectors that gate 12C's
`two_wave` attempt - a different phenomenon from the onset lag directly
above (a delay before absorption starts at all): a genuinely separate wave
arriving after the first is already underway, same shape as 12C's own
`two_wave` model but decided from 13C6's own evidence, not inherited from
12C's structure.

13C6 is markedly noisier than 12C, and this needed extra caution. A
cohort-wide scan of the same two detectors against 13C6's own raw data
flagged ~69% of curves - at least as high a rate as 12C's own - which is a
warning sign given 13C6's much lower absolute signal and higher relative
noise. Standalone validation before implementing confirmed the concern
directly: `ER21` FCT1 (the strongest `dip_excess` in the cohort, 3.416)
"improved" R² 0.508->0.741, but only by pinning `f_delayed2_13C6`=0.999 -
the same near-total-delay degenerate pattern already distrusted for 12C's
`t_lag`/`f_delayed` (see the `ER06` FCT2 case in `BOUNDS_LAGGED`'s
comment) - driven by a late, near-zero crash-then-rebound (0.0117 mg/L at
t=150, 0.0354 at t=180) that reads as measurement noise, not a real second
dose. Two other candidates looked genuine by contrast: `ER11` FCT1 (R²
0.787->0.846, `f_delayed2_13C6`=0.751) and `ER18` FCT2 (R² 0.860->0.879,
`f_delayed2_13C6`=0.095) - both comfortably away from the extremes.

So two safeguards beyond what 12C's own `two_wave` needed are applied
here:
- `f_delayed2_13C6` is bounded to `[0.05, 0.95]`, not 12C's
  `[0.001, 0.999]` - directly excludes the near-total-delay trick that made
  `ER21` look good on paper.
- `t_lag2_13C6` is bounded to `[30, 150]`, not down to 5 - the onset-lag
  candidate above already owns "delay before the wave starts at all"; this
  candidate only needs to (and should only be allowed to) model a
  genuinely LATER, separate hump, not also compete for the same unsampled
  0-30min gap.

Modeled with a plain `bateman_conc` per wave (instant bolus), same as the
onset-lag candidate and for the same reason. Compared via AIC (K=3) against
whichever candidate is currently winning (no-lag or onset-lag, K=2 either
way) - mutually exclusive with the onset-lag candidate, since both explain
the same flat/slow start in different ways.

**Fixed: these two 13C6 candidates (onset-lag and independent second-wave)
are now also available when 12C itself is `two_wave`**, not just
`single_wave` - previously, whenever 12C was `two_wave`, 13C6 was stuck
with only `fit_subject_visit_two_wave()`'s stage-2 mechanism
(`simulate_delayed_release`, sharing 12C's own `t_lag`), which structurally
cannot represent an onset-delayed or independently-timed 13C6 curve.
Confirmed broken on `ER23` FCT2: 13C6's own data goes from 0.009 mg/L at
t=30 to 0.165 at t=60 (5.6% of peak, a textbook onset-lag shape), but with
12C selected as `two_wave` (correctly, on its own merits), the baseline
mechanism collapsed to `f_delayed_13C6`=0.001 and r2_13C6=0.35 - not a
search failure, a structural mismatch no amount of optimization could fix.
Now all three candidates (baseline K=3, onset-lag K=2, second-wave K=3)
are compared via AIC on 13C6's own RSS regardless of 12C's own model;
fixed, `ER23` FCT2 correctly picks up its own onset lag (~27min) and
reaches r2_13C6=0.78.

**Fixed: a boundary-pinned false positive in the second-wave candidate.**
`f_delayed2_13C6` landing at (or within 1e-3 of) its own lower bound (0.05)
means the optimizer wants an even smaller delayed fraction than the bound
allows - i.e. no real second wave at all, the exact degenerate regime that
bound exists to exclude - so this candidate is now rejected outright
(not even AIC-compared) whenever that happens, in both branches. Found on
`ER18` FCT2: won its AIC comparison at `f_delayed2_13C6`=0.05 **and**
`t_lag2_13C6`=150 simultaneously (both own bounds, at once) despite the
curve's genuine dip evidence sitting at t=120, not t=150 - a spurious
corner-of-the-box "improvement" unrelated to the actual evidence. With the
guard, `ER18` FCT2 correctly falls back to the baseline mechanism
(r2_13C6=0.66, down from the spurious 0.71, but honest).

**Tested but deliberately not adopted (yet):**
- **A two-lag onset-time extension** (`simulate_two_lag_dose()`, also in
  `pk_curves.R`) that gives the *first* wave its own fittable onset delay
  (`t_lag1`), not just the second wave - motivated by curves whose very
  first real sample is already low relative to what follows, suggesting
  absorption itself hadn't started by then (a standard, separate PK concept
  from the second-wave delay above). On `ER03` FCT1 it matched
  `former_models/MixedModel`'s historical R² almost exactly (0.71 vs 0.71),
  with `f_delayed` landing near 0 (mostly using the onset lag, not genuine
  two-wave behavior) - but this is n=1, the fit didn't converge, and its
  12C R² (0.70) was worse than the simpler unweighted baseline (0.85) for
  the same subject. Not enough validation to commit a full-cohort run to it;
  worth revisiting with more subjects and search budget.
- **Loosening `MIN_TMAX` further** (this version keeps it at 30, unchanged).
  A sweep on `ER03` FCT1 (30/20/10, keeping the *old* weighted objective)
  showed only modest improvement (R² -0.86 -> -0.67, plateauing at 20 and
  10 - the floor was not the binding constraint once weighting is the real
  problem). Whether loosening it adds anything *on top of* the unweighted
  objective and/or the lagged-dose model is untested and worth checking
  before changing it.

## Volume of distribution (Vd)

`Vd` is estimated per subject x visit as total blood volume, from that
subject's own weight, height, and sex, via Nadler's equation:

```
Men:   Vd (L) = 0.3669 x height(m)^3 + 0.03219 x weight(kg) + 0.6041
Women: Vd (L) = 0.3561 x height(m)^3 + 0.03308 x weight(kg) + 0.1833
```

> Nadler DA, Hidalgo JU, Bloch T. *Prediction of blood volume in normal
> human adults.* Surgery. 1962;51(2):224-232.

Implemented as `nadler_blood_volume()` in `scripts/assets/pk_curves.R` -
generic and project-agnostic, like the rest of that file, since a
weight+height+sex-based Vd estimate is useful well beyond this study.

**This is a deliberate change from an earlier version of this model**,
which used a flat `Vd = 0.15 L/kg` body weight for every subject - an
extracellular fluid volume (ECFV) estimate borrowed from glucose
literature, not blood volume, specifically because a small,
freely water-soluble, non-protein-bound molecule like fructose is expected
to equilibrate into interstitial fluid, not just the vascular compartment:

> van der Crabben SN et al. *Relationship between glucose volume of
> distribution and the extracellular space: a multiple tracer study.*
> Metabolism. 2011. Glucose Vd measured at 191-206 mL/kg across three
> tracers, shown to equal the extracellular fluid space rather than blood or
> plasma volume.

That argument for ECFV over blood volume was not superseded by any new
evidence - this project has not independently established that blood
volume is the more defensible choice for fructose's actual distribution
volume. The switch was made anyway, as a deliberate project decision to use
individual weight+height+sex over a flat per-kg ratio, accepting the
resulting mismatch between the Vd used and the model's own prior
physiological reasoning above.

**What changes under this switch:** blood volume (~65-75 mL/kg for a
typical adult, per Nadler) is substantially smaller than the ~150-200 mL/kg
ECFV estimate used previously. Since bioavailable fraction only enters the
model as `F * dose / Vd` (see the identifiability caveat below), a smaller
Vd at the same observed concentration implies a smaller `F` - so
`F_12C`/`F_13C6` values from this version are **not directly comparable**
to results generated before this switch. `ka`/`kel` are unaffected (they
don't depend on Vd). Within-subject, paired comparisons (e.g. baseline vs.
intervention `F`) remain valid on the same logic as before: each subject's
own Nadler-estimated Vd is used consistently within a visit, and any
systematic Vd bias still applies to both visits of the same subject (Vd
does now differ *between* baseline and intervention if weight changed,
exactly as the old flat per-kg version also would have).

## Fitting bounds (why these numbers, not others)

| Parameter | Bounds | Rationale |
|---|---|---|
| `kel` | [0.005, 0.1] /min (t½ ≈ 7-140 min) | Anchored to Hannou et al. 2018's fructose metabolism review. Early unconstrained fits gave physiologically implausible, wildly inconsistent `kel` across subjects - this bound was added for that reason, not chosen a priori. |
| `ka` | (0, 1) /min | Bounded below only, in spirit - absorption-rate differences (e.g. pre- vs. post-diet, or liquid vs. capsule route) are part of the research question, so the upper bound (1/min) is a generous ceiling, not a real constraint. |
| `F_12C`, `F_13C6` | (0, 1) | Physical: a fraction of a dose. |
| `k_release` | [0.001, 1] /min (t½ ≈ 0.7-700 min) | Wide enough to cover anything from near-instant to very slow capsule dissolution. |

A soft penalty (predicted Tmax pulled above 30 min, weight `TMAX_LAMBDA`)
rules out a specific degenerate failure mode: a spurious "early spike"
solution where the model absorbs and clears almost instantly, producing a
sharp early peak invisible between the sparse observed timepoints. A second
soft penalty (predicted Cmax pulled toward observed Cmax +/-10%) discourages
systematic over/undershoot of the real peak without rigidly constraining
the rest of the curve. Both are evaluated on a dense time grid, not just
the observed sampling times, specifically because degenerate solutions are
designed (by the optimizer, inadvertently) to look fine only at the sparse
observed points. Both are also deliberately smooth (a squared
shortfall/excess that is zero at and below the threshold, not a
discontinuous jump) - `optim()`'s L-BFGS-B relies on finite-difference
gradients, which a hard penalty boundary makes needlessly rough to search
near.

## Fitting procedure

PK objective surfaces like this one are routinely multimodal (a
fast-absorption/fast-clearance solution can fit nearly as well as a
slow/slow one) - a single optimizer run from one starting point is not a
defensible fit. Each subject x visit is fit in two stages:

1. **Independent pre-fit per curve** (closed-form for 12C, still numerical
   for 13C6 since a lone tracer curve alone has to satisfy the same ODE),
   used only to generate informed starting points for stage 2.
2. **Joint fit**, seeded from stage 1's estimates plus a systematic grid and
   random (log-uniform, for rate constants) starting points, minimizing a
   *normalized, weighted* combined objective:
   - **Across curves:** each curve contributes a weighted analog of
     `sse / ss_tot` (roughly `1 - R²` for that curve alone) rather than raw
     squared error. This normalization is essential - the 12C curve's
     concentrations are roughly 1000x the 13C6 tracer's, so combining raw
     SSE lets 12C dominate the objective and the optimizer effectively
     ignores 13C6 entirely.
   - **Within a curve:** `build_joint_objective()` (`scripts/assets/pk_fit.R`)
     supports an optional `1 / max(pred, floor)²` proportional/constant-CV
     weighting (floored at 1% of that curve's own Cmax to avoid blow-up near
     zero), added after checking residuals from an unweighted version of
     this fit against predicted concentration and finding clear
     heteroscedasticity - squared-residual scale differed ~24x between the
     top and bottom quartile of predicted concentration, for both curves -
     meaning unweighted SSE was letting each curve's own peak region
     dominate its fit even after the cross-curve normalization above, at the
     expense of the tail (which carries most of the information about
     `kel`).
     **As of this version, `02_fit_erie_model.R` calls it with
     `proportional_weighting = FALSE` (unweighted) instead.** The weighting
     fixes a real, measured problem *on average*, but it has a specific,
     serious failure mode: it discounts exactly the region a genuinely
     informative early feature lives in, by construction (a large predicted
     value gets a *small* weight). Found via `ER03` FCT1's 13C6 curve, a
     known-catastrophic fit (R² = -2.18, see "R² can be negative" below) -
     it spikes at t=30 then crashes by t=90, and no amount of loosening
     `MIN_TMAX` fixed it (a sweep from 30 to 10, keeping the weighted
     objective, only moved R² from -0.86 to -0.67 and plateaued there -
     confirming the floor was not the binding constraint). Switching that
     same fit to unweighted - same `MIN_TMAX`, same bounds, ordinary
     unseeded multistart - recovered R² = 0.37. `former_models/MixedModel`'s
     own historical fit for the *exact same subject*, independently
     re-derived and validated (see "Continuity with prior work" below),
     reached R² = 0.71 - and its objective was unweighted all along (see its
     `objective_joint()` in
     `former_models/MixedModel/Scripts/fructose_joint_model_final.R`,
     which normalizes each curve by its own `ssTot` only, with no
     `1/pred²` term). That is strong indirect evidence unweighted is sound
     at full-cohort scale too, not just for this one subject - MixedModel's
     own cohort-wide numbers (median R² ~0.9+ for 12C, ~0.7-0.9 for 13C6,
     see "Fit quality and what to trust") were already produced this way.
     Still, this is a real tradeoff, not a strict improvement - the
     heteroscedasticity that motivated weighting in the first place doesn't
     go away just because it's not the default; `weight_floor_frac` and
     `proportional_weighting = TRUE` remain available in
     `build_joint_objective()` for a sensitivity check, and `r2_12C`/`r2_13C6`
     should be checked broadly (not just on previously-known-hard subjects)
     before fully trusting this as the final word.

Both stages reuse the same generic multi-start optimization engine
(`scripts/assets/pk_fit.R`) - there is no separate, independently-maintained
"independent model" implementation to keep in sync with the joint one.

Fitting is parallelized across subject x visit pairs via
`parallel::mclapply()`. Each pair seeds its own random starting points
deterministically from its own subject/visit ID (`string_seed()` in
`pk_fit.R`), rather than from one `set.seed()` call before the parallel
loop - the latter does not actually make results reproducible under
`mclapply()` (verified: identical code, same seed, different draws on
separate runs), and even the standard workaround
(`RNGkind("L'Ecuyer-CMRG")`) only fixes reproducibility for a fixed task
order and core count, not against changes to either. Per-ID seeding makes
each subject's fit reproducible regardless of `N_CORES` or the row order of
`subject_visits`.

**Adaptive retry for flagged fits.** Any subject×visit flagged after the
standard pass above gets a second, denser retry: the prior estimate itself
as a seed (so the retry can never end up worse), a systematic grid spanning
the `kel`/`k_release` bound ranges (if either is the reason it's flagged),
and 80 more random starts, at a higher `maxit`. This exists because the
prior version of this model
(`former_models/MixedModel/Scripts/fructose_joint_model_final.R`) found a
subject stuck at the `kel` bound purely because the standard search kept
landing on a worse local optimum that happened to sit exactly on the
boundary - a substantially better solution existed comfortably inside the
bound the whole time.

The retry trigger is deliberately **not** "boundary flag alone": a fit
stuck in the same kind of worse local optimum that happens to land
comfortably *inside* the bounds would show no boundary flag at all, and
look unremarkable - but is exactly what a poor R² or a non-converged fit
would actually look like. So a subject×visit is retried if `kel_at_bound`,
`k_release_at_bound`, `r2_12C_low`, or `converged = FALSE` - any of them,
not boundary flags only. `r2_13C6_low` alone is **not** a retry trigger
(for either `single_wave` or `two_wave`), unlike the other flags: a poor
13C6 fit is often a structural mismatch - the shared `ka`/`kel`/wave-timing
13C6 inherits from 12C just isn't right for its own curve, and now that
13C6 has its own independent `f_delayed_13C6` in `two_wave` (see "13C6's
own wave choice" above) specifically to address that, a persistently poor
`r2_13C6` after that freedom means the freedom itself didn't find a good
shape - not that the search under-explored. A denser retry of the same
structure is unlikely to fix that, so retrying is reserved for cases where
it can actually help. The retry only replaces the original
result if its objective value (`objective_value` in the results table - not
comparable across subjects, only against that same subject's own
prior/retry pair) is strictly better, so broadening the trigger this way
can only improve results, at the cost of retrying more subjects (and
therefore more runtime) than a boundary-only trigger would.

**A flagged fit is not automatically a search failure**, and the results
table distinguishes this: `retried` (TRUE if this subject's fit was flagged
and a denser retry was attempted) and `retry_improved` (TRUE if that retry
found something strictly better, FALSE if the denser search was attempted
but couldn't beat the original, NA if never flagged). A low-R² 13C6 curve
in particular can be a genuinely weak, low-information fit rather than an
under-searched one - e.g. a subject with especially slow/delayed capsule
opening produces a 13C6 curve whose shape is dominated by `k_release`
rather than by `ka`/`kel`, which can leave a broad, nearly-flat region of
the objective surface where no parameter combination fits meaningfully
better than any other. `retried = TRUE, retry_improved = FALSE` is exactly
the signature of that: the denser search had every opportunity to find a
better answer and didn't, which is evidence the flag reflects real curve
weakness rather than an optimizer failure. `retried = TRUE, retry_improved
= TRUE` is the opposite signature - the original result really was
under-searched.

## Fit quality and what to trust

The 12C curve fits well for most subjects (median R² ~0.85), but the 13C6
curve does not: as of the current model (Nadler Vd, weighted objective,
adaptive retry), **45/68 (66%) of 13C6 fits still fall below R² 0.70 after
retry**, vs. 19/68 (28%) for 12C. This is a materially different picture
from earlier versions of this model, which described only "a handful of
known-hard subjects" - that language is no longer accurate and the
`r2_13C6_low` rate should not be read as rare. Per-subject R² is in
`results/fit_results_joint.csv` - always check it before trusting an
individual subject's parameters, and inspect that subject's plot in
`results/plots_individual/` if R² is low.

The likely explanation is not that these fits are under-searched: the
adaptive retry pass explicitly targets low-R² fits with a much denser
search (see "Adaptive retry" above), and this run retried 58/68 subjects
for exactly that reason. A weak, low-information 13C6 curve - e.g. from an
especially slow/delayed capsule opening, where the curve's shape is
dominated by `k_release` rather than `ka`/`kel` - can have a genuinely low
achievable R² that no amount of extra search improves. The
`retried`/`retry_improved` columns (added after this particular run;
rerun to get them) make this directly checkable per subject: `retried =
TRUE, retry_improved = FALSE` on a low-R² 13C6 fit is evidence for genuine
curve weakness rather than a search failure. Until that's been checked
across the cohort, treat the 66% figure as "13C6 fits are frequently hard
to pin down precisely," not as "45 fits are wrong."

`results/fit_results_joint.csv` also carries `r2_12C_low` / `r2_13C6_low`,
TRUE when that curve's own R² falls below `R2_RELIABLE_MIN` (0.70, set in
`02_fit_erie_model.R`). These are per-curve, not a single combined verdict
on the subject×visit: a low flag on one curve does not by itself mean the
other curve's R², or the `ka`/`kel` shared across both, are also
unreliable - though because `ka`/`kel` are fit jointly, a poor fit on one
curve can still bias them, so a low flag is a prompt to inspect that
subject's plot, not just to drop the flagged curve's own parameter.

**Trustworthy as (approximately) absolute numbers:**
- `kel`, `ka` - reasonably well-identified given the bounds and multi-start
  search.
- Within-subject, paired comparisons (e.g. baseline vs. intervention `F`) -
  a systematic Vd bias applies equally to both visits of the same subject
  and is expected to cancel in a paired comparison.
- The capsule dissolution half-life (`log(2) / k_release`) - a genuinely new
  quantity this model provides that a single-curve model cannot estimate at
  all.

**Not trustworthy as an absolute number:**
- `F` (bioavailable fraction) on its own. **This is a structural
  identifiability limit of oral-only concentration data, not a fitting
  defect**: `F` and `Vd` enter the model only as a product (`F * dose /
  Vd`), so the data cannot distinguish "small F, small Vd" from "large F,
  large Vd" - only their ratio is identified. Resolving this would require
  an independent reference (e.g. an IV tracer dose in the same subjects, or
  a directly measured individual Vd), which this protocol does not include.
  Report `F` as conditional on the Nadler blood-volume Vd assumption, not
  as a precise absolute bioavailability.

  **Sanity-checking the magnitude under the current Vd:** with Nadler
  blood-volume Vd, median `F_12C` is ~3.2% and median `F_13C6` is ~3.0%
  (range roughly 1-11% and 1-59% respectively; the 13C6 upper end is driven
  by a handful of the same weak/hard-to-fit curves discussed above). These
  are low - single-digit percent - and, unlike the earlier flat-ECFV
  version of this model, haven't yet been checked against a literature
  fructose/glucose bioavailability figure at a comparable dose. That
  comparison is worth doing before these numbers go into anything
  manuscript-facing: if literature values are substantially higher, that's
  a signal Nadler blood volume specifically (as opposed to ECFV, or some
  other Vd estimate) may be too small for fructose's actual distribution
  volume, on top of the already-acknowledged mismatch with this model's own
  prior physiological reasoning (see "Volume of distribution" above).
- Any boundary-flagged parameter (`kel_at_bound` or `k_release_at_bound` =
  TRUE in the results table) - even after the adaptive retry above, the data
  may genuinely not constrain that parameter away from the bound (e.g.
  `k_release` pinning at its ceiling simply because the first post-dose
  sample already shows near-peak tracer concentration, and nothing in the
  data argues for a slower dissolution rate - a sampling-resolution limit,
  not an error). A boundary flag surviving the retry (`retried = TRUE,
  retry_improved = FALSE`) is more trustworthy than one from a single pass,
  but still doesn't distinguish "genuinely unconstrained by the data" from
  "search still didn't find it" - inspect the subject's plot either way.
- Any fit with `converged = FALSE` - the winning multi-start result did not
  actually satisfy `optim()`'s own convergence criterion (it just had the
  lowest objective value among the seeds tried), typically because it hit
  the `maxit` cap rather than reaching a true local optimum.
- Extrapolated quantities (e.g. AUC beyond 360 min, full clearance time) -
  these are projections of the fitted curve past the last real observation,
  not confirmed by data past that point.

## Continuity with prior work

This model deliberately reconstructs the most-validated version from
`former_models/MixedModel/Scripts/fructose_joint_model_final.R` (Melany's
final joint model), rather than either of the earlier independent
single-curve models or Bas's original 3-compartment Julia model - her own
analysis showed the joint, shared-`ka`/`kel` version out-performed fitting
either curve alone (13C6 median R² improved from ~0.70 independently to
~0.78 jointly, with far fewer subjects landing below R² = 0.3), which is
itself evidence in favor of the shared-clearance biological hypothesis, not
just a modeling convenience.

What's different in this version:
- Split into a generic, project-agnostic PK/optimization engine
  (`scripts/assets/pk_curves.R`, `pk_fit.R`, `pk_diagnostics.R`) plus a
  short ERIE-specific wiring script (`02_fit_erie_model.R`), instead of one
  long self-contained script - `scripts/assets/` is meant to be copied
  as-is into other PK projects.
- The 13C6 tracer dose defaults to the confirmed administered 120 mg /
  644.78 µmol (recalculated from the confirmed mass - the xlsx's stated
  molar amount was itself inconsistent with it) rather than the ~100 mg
  originally implied by `ERIE_constants.csv`/`.xlsx` (see
  `docs/data-cleaning-notes.md`) - this shifts `F_13C6` estimates slightly
  relative to that earlier ~100 mg assumption (a larger assumed dose implies
  a somewhat smaller `F_13C6` for the same observed concentration) but does
  not change `ka`, `kel`, or any of the qualitative conclusions above.
- Fitting is parallelized across subject x visit pairs (each is
  independent) rather than run sequentially, for practical runtime.
- `Vd` is now estimated per subject x visit from weight/height/sex (Nadler
  blood volume) rather than a flat 0.15 L/kg ECFV assumption applied to
  every subject - see "Volume of distribution" above for what this changes
  and why it isn't a straightforward improvement.

### R² is lower than the former model's, and that's expected, not a regression

Running `scripts/04_compare_to_former_model.R` (`pixi run compare-former`)
against the former model's own saved results
(`former_models/MixedModel/Results/fit_results_joint.csv`) shows this
version's classical R² (`1 - SSE/SS_tot`) is **lower for most subjects**,
not higher: median 12C R² 0.93 -> 0.85 (worse for 64/68 subjects), median
13C6 R² 0.79 -> 0.54 (worse for 67/68 subjects). Full per-subject numbers
in `results/r2_comparison_vs_former_model.csv`, a scatter plot in
`results/r2_comparison_plot.pdf`.

This is not the adaptive retry pass failing, or a fitting bug - it's a
predictable consequence of the proportional weighting change (see
"Fitting procedure" above). Classical R² is computed from *raw* squared
error, which is exactly the quantity the former model's objective
minimized - its fitting target and its evaluation metric were the same
thing, so a high R² was close to guaranteed by construction. This model's
objective deliberately minimizes *proportionally weighted* error instead,
specifically to stop the peak from dominating each curve's fit at the
tail's expense (see the ~24x squared-residual-scale finding that motivated
it). That is a real change in what "good fit" means, not just a
coefficient tweak - so it should be no surprise that raw-SSE R² looks worse
under an objective that was never trying to maximize it. Whether the
proportionally-weighted fit is actually *better* for this data (more
accurate where it matters, e.g. the tail that drives `kel`) is a real
question this comparison doesn't answer - it would need a fit-quality
metric computed on the same weighted basis the model actually optimizes,
which hasn't been built yet. Until then, don't read the lower R² alone as
"the new model fits worse" - it's evaluating the new model by the old
model's yardstick.

## References

- Hannou SA, Haslam DE, McKeown NM, Herman MA. *Fructose metabolism and
  metabolic disease.* J Clin Invest. 2018;128(2):545-555.
- Jang C, Hui S, Litchfield B, et al. *The small intestine converts dietary
  fructose into glucose and organic acids.* Cell Metabolism.
  2018;27(2):351-361.
- van der Crabben SN, et al. *Relationship between glucose volume of
  distribution and the extracellular space: a multiple tracer study.*
  Metabolism. 2011.
- Nadler DA, Hidalgo JU, Bloch T. *Prediction of blood volume in normal
  human adults.* Surgery. 1962;51(2):224-232.
