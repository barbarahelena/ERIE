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

A hard rule (predicted Tmax >= 30 min) rules out a specific degenerate
failure mode: a spurious "early spike" solution where the model absorbs and
clears almost instantly, producing a sharp early peak invisible between the
sparse observed timepoints. A soft penalty (predicted Cmax pulled toward
observed Cmax +/-10%) discourages systematic over/undershoot of the real
peak without rigidly constraining the rest of the curve. Both are evaluated
on a dense time grid, not just the observed sampling times, specifically
because degenerate solutions are designed (by the optimizer, inadvertently)
to look fine only at the sparse observed points.

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
   - **Within a curve:** points are weighted `1 / max(pred, floor)²`
     (proportional/constant-CV weighting, floored at 1% of that curve's own
     Cmax to avoid blow-up near zero). Checking residuals from an earlier,
     unweighted version of this fit against predicted concentration showed
     clear heteroscedasticity - squared-residual scale differed ~24x
     between the top and bottom quartile of predicted concentration, for
     both curves - meaning unweighted SSE was letting each curve's own peak
     region dominate its fit even after the cross-curve normalization
     above, at the expense of the tail (which carries most of the
     information about `kel`). Relative residual variance was not fully
     constant across concentration either (higher at low concentration than
     pure proportional weighting assumes, consistent with a "combined"
     additive+proportional error structure) - proportional weighting
     corrects the dominant bias without fitting a full combined-error
     model, which would be a larger undertaking for a modest additional
     gain.

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
  TRUE in the results table) - the search may not have found the true
  optimum, or the data may genuinely not constrain that parameter away from
  the bound (e.g. `k_release` pinning at its ceiling simply because the
  first post-dose sample already shows near-peak tracer concentration, and
  nothing in the data argues for a slower dissolution rate - a
  sampling-resolution limit, not an error). A boundary flag does not by
  itself distinguish these two cases - that requires denser search near the
  bound (not yet automated; see below) or inspecting the subject's plot.
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
