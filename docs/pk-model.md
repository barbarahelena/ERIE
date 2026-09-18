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

Every subject has this pair of curves at both FCT1 (baseline) and FCT2
(post 4-week diet intervention).

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

A meaningful minority of curves (~19% by a simple threshold-free check: any
post-peak local trough followed by a rise of >=15% of that curve's own peak)
show a dip after an initial peak, then a second, often higher, rise before
the final decline - most clearly on `ER01` and `ER09` (every one of their
four curves), and at lower relative amplitude (a flattening rather than a
full second rise) on several more, including `ER03`. A smaller-amplitude
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
pool active from t=0 at a different rate. A shared fraction `f_delayed` of
each dose (12C and 13C6 both - they're ingested by the same subject at the
same time, so a plausible shared upstream trigger, e.g. transient
duodenal-brake feedback inhibition from the large osmotic/caloric load,
should affect both) contributes nothing until a shared `t_lag`, then behaves
exactly like a fresh dose given at that later time - mathematically
equivalent to two separately-timed doses. This breaks the shared-tail
constraint above because the delayed portion has a genuine gap, not just a
different decay rate, and was confirmed both analytically (a lag creates an
actual discontinuity in the delivery rate) and numerically (produces a real
local minimum followed by a second, higher peak) before being wired into
the fitting pipeline. Piloted on `ER01`/`ER03`/`ER09` (`results/model-post-peak-dip-pilot/`)
with a clear, visually-confirmed win on the genuine dip cases (`ER01`,
`ER09`) - both curves transformed from a single smooth hump to a real
double-humped fit tracking the dip - and degrades gracefully (small
`f_delayed`) on subjects that don't need it. `02_fit_erie_model.R` now fits
this model *alongside* the baseline for every subject x visit
(`fit_subject_visit_lagged()`), rather than replacing the baseline outright,
so both are available for comparison in `results/fit_results_baseline.csv`
and `results/fit_results_lagged.csv`.

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
don't depend on Vd). Within-subject, paired comparisons (e.g. FCT1 vs. FCT2
`F`) remain valid on the same logic as before: each subject's own
Nadler-estimated Vd is used consistently within a visit, and any systematic
Vd bias still applies to both visits of the same subject (Vd does now
differ *between* FCT1 and FCT2 if weight changed, exactly as the old flat
per-kg version also would have).

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
`k_release_at_bound`, `r2_12C_low`, `r2_13C6_low`, or `converged = FALSE` -
any of them, not boundary flags only. The retry only replaces the original
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

Reproducing this model against `former_models/`'s own validated results
(see "Continuity with prior work") gives fit quality in the same range
previously reported: most subjects land around R² ~0.9+ for the 12C curve
and R² ~0.7-0.9 for the 13C6 curve, with a handful of known-hard subjects
(noisy curves, e.g. an isolated early spike inconsistent with the rest of
the curve) landing much lower on either curve. Per-subject R² is in
`results/fit_results_joint.csv` - always check it before trusting an
individual subject's parameters, and inspect that subject's plot in
`results/plots_individual/` if R² is low.

`results/fit_results_joint.csv` also carries `r2_12C_low` / `r2_13C6_low`,
TRUE when that curve's own R² falls below `R2_RELIABLE_MIN` (0.70, set in
`02_fit_erie_model.R`). These are per-curve, not a single combined verdict
on the subject×visit: a low flag on one curve does not by itself mean the
other curve's R², or the `ka`/`kel` shared across both, are also
unreliable - though because `ka`/`kel` are fit jointly, a poor fit on one
curve can still bias them, so a low flag is a prompt to inspect that
subject's plot, not just to drop the flagged curve's own parameter. With
the current cohort this flags 9/68 fits on `r2_12C` and 22/68 on
`r2_13C6`.

**Trustworthy as (approximately) absolute numbers:**
- `kel`, `ka` - reasonably well-identified given the bounds and multi-start
  search.
- Within-subject, paired comparisons (e.g. FCT1 vs. FCT2 `F`) - a
  systematic Vd bias applies equally to both visits of the same subject and
  is expected to cancel in a paired comparison.
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
