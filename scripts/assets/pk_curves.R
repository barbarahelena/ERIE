# =============================================================================
# pk_curves.R
#
# Generic, project-agnostic building blocks for first-order oral
# pharmacokinetic (PK) compartment models. Nothing here knows about fructose,
# ERIE, or any specific study - it only knows about doses, rate constants,
# and volumes of distribution. Reuse this file as-is in other PK projects;
# put study-specific config (doses, bounds, which Vd estimate to use) in the
# analysis script instead. One reusable Vd estimator - nadler_blood_volume()
# - is provided below; whether it's the right physiological quantity for a
# given analyte is a per-study call, made in the analysis script.
#
# Two model forms are provided:
#   - A one-compartment model with first-order absorption and elimination
#     (the standard oral Bateman equation), solved in closed form.
#   - A delayed-release extension (e.g. a dissolving capsule) that adds an
#     upstream "release" compartment feeding the absorption compartment,
#     solved numerically since a closed-form solution becomes numerically
#     unstable near parameter degeneracies (ka ~= kel ~= k_release).
#
# All concentration functions share the same argument convention:
#   ka    - absorption rate constant (1/time)
#   kel   - elimination rate constant from the central compartment (1/time)
#   Fbio  - bioavailable fraction of the dose that reaches the central
#           compartment (unitless, 0-1)
#   dose  - administered dose, in mass units consistent with Vd
#   Vd    - apparent volume of distribution, in volume units such that
#           dose / Vd has the desired concentration units
# =============================================================================

#' Concentration from a one-compartment oral model (Bateman equation)
#'
#' Closed-form solution for:
#'   dA1/dt = -ka * A1                    A1(0) = dose
#'   dA2/dt =  Fbio * ka * A1 - kel * A2   A2(0) = 0
#'   C(t)   =  A2(t) / Vd
#'
#' @param t Numeric vector of times.
#' @param ka,kel,Fbio,dose,Vd See file header for conventions.
#' @return Numeric vector of concentrations at `t`.
bateman_conc <- function(t, ka, kel, Fbio, dose, Vd) {
  # ka == kel is a removable singularity in the standard Bateman formula;
  # use the limiting form instead of dividing by (ka - kel) ~= 0.
  if (abs(ka - kel) < 1e-8) {
    amount <- Fbio * ka * dose * t * exp(-ka * t)
  } else {
    amount <- (Fbio * ka * dose / (ka - kel)) * (exp(-kel * t) - exp(-ka * t))
  }
  amount / Vd
}

#' Time of peak concentration (Tmax) for the one-compartment oral model
#'
#' @param ka,kel Absorption and elimination rate constants.
#' @return Tmax, in the same time units as ka/kel.
bateman_tmax <- function(ka, kel) {
  if (abs(ka - kel) < 1e-8) return(1 / ka)
  log(ka / kel) / (ka - kel)
}

#' ODE system for a delayed-release oral model
#'
#' Three sequential first-order compartments: an undissolved "release" pool
#' (e.g. a capsule) feeds an absorption pool, which feeds the central
#' (sampled) compartment:
#'   dRelease/dt <- -k_release * Release                  Release(0) = dose
#'   dAbsorb/dt  <-  k_release * Release - ka * Absorb     Absorb(0)  = 0
#'   dCentral/dt <-  Fbio * ka * Absorb  - kel * Central   Central(0) = 0
#'
#' Not exported for direct use; called via [simulate_delayed_release()].
#'
#' @keywords internal
delayed_release_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dRelease <- -k_release * Release
    dAbsorb  <-  k_release * Release - ka * Absorb
    dCentral <-  Fbio * ka * Absorb  - kel * Central
    list(c(dRelease, dAbsorb, dCentral))
  })
}

#' Simulate a delayed-release oral model at a set of times
#'
#' Numerical integration is used deliberately instead of a closed-form
#' sum-of-exponentials solution: the closed form has removable singularities
#' whenever any two of (ka, kel, k_release) are close to equal, and a
#' hand-derived version was found to have real, sometimes large, errors near
#' those degeneracies when checked against numerical integration. Numerical
#' ODE solving avoids this at a modest speed cost.
#'
#' @param times Numeric vector of times at which to return concentration.
#' @param k_release Release-pool dissolution rate constant (1/time).
#' @param ka,kel,Fbio,dose,Vd See file header for conventions.
#' @return A tibble with columns `time` and `conc`, aligned to the
#'   requested `times` (t = 0 is added internally if not already present,
#'   but only requested times are returned).
simulate_delayed_release <- function(times, k_release, ka, kel, Fbio, dose, Vd) {
  solve_times <- sort(unique(c(0, times)))
  state0 <- c(Release = dose, Absorb = 0, Central = 0)
  parms  <- c(k_release = k_release, ka = ka, kel = kel, Fbio = Fbio)

  out <- deSolve::ode(y = state0, times = solve_times, func = delayed_release_ode,
                       parms = parms, method = "lsoda")
  out <- tibble::as_tibble(out)

  tibble::tibble(
    time = times,
    conc = out$Central[match(times, out$time)] / Vd
  )
}

#' ODE system for two simultaneously-dosed oral curves sharing a biphasic
#' gastric-emptying process
#'
#' CAVEAT (kept for the record, superseded by [simulate_lagged_dose()] for
#' reproducing an actual dip-then-rise): pilot fitting showed this structure
#' cannot represent a genuine two-humped curve for any parameter values, and
#' not for lack of searching - a targeted parameter search (thousands of
#' draws, several hand-engineered attempts, including giving the two gastric
#' pools independent absorption compartments and independent `ka`) never
#' produced one. The reason is structural: any number of pathways that all
#' still start delivering dose from t=0 and drain into one shared elimination
#' compartment converge to the same single-exponential terminal decay
#' (`kel`), so summing them can reshape/delay a single hump but not create a
#' genuine local minimum followed by a second rise. What actually works is a
#' genuine start-time offset (see [simulate_lagged_dose()]) - equivalent to a
#' second, later dose - which breaks that shared-tail constraint.
#'
#' Motivated by a post-peak dip-then-rise (or, at lower relative amplitude, a
#' mere flattening of the decline) seen consistently around t=60-90min in a
#' meaningful minority of curves, essentially always co-occurring on both
#' curves of the same subject x visit at once when it appears strongly. A
#' one-compartment absorption model is mathematically monotonic after its
#' single peak for any parameter values, so it cannot represent this shape at
#' all; a plausible mechanism is biphasic gastric emptying (e.g. transient
#' duodenal-brake feedback inhibition from a large osmotic/caloric load,
#' followed by a second emptying wave) - a process upstream of, and shared
#' by, both simultaneously-ingested doses, rather than something in either
#' curve's own absorption/clearance.
#'
#' Generalizes the "instant full dose available at t=0" assumption (curve A,
#' e.g. a liquid dose) and the single-pool delayed-release model (curve B,
#' e.g. an enteric capsule, which must itself first clear the stomach before
#' its coating can dissolve) into a shared two-pool gastric-emptying split: a
#' fraction `f_fast` of each dose empties at rate `k_ge_fast`, the remainder
#' at the slower `k_ge_slow` - same fraction and rates for both curves, since
#' both are subject to the same stomach at the same time:
#'
#'   Stomach_fast_A(0) = f_fast * doseA        Stomach_slow_A(0) = (1-f_fast) * doseA
#'   dStomach_fast_A/dt = -k_ge_fast * Stomach_fast_A
#'   dStomach_slow_A/dt = -k_ge_slow * Stomach_slow_A
#'   dAbsorbA/dt = k_ge_fast*Stomach_fast_A + k_ge_slow*Stomach_slow_A - ka*AbsorbA   AbsorbA(0) = 0
#'   dCentralA/dt = FbioA * ka * AbsorbA - kel * CentralA                            CentralA(0) = 0
#'
#'   Stomach_fast_B(0) = f_fast * doseB        Stomach_slow_B(0) = (1-f_fast) * doseB
#'   dStomach_fast_B/dt = -k_ge_fast * Stomach_fast_B
#'   dStomach_slow_B/dt = -k_ge_slow * Stomach_slow_B
#'   dReleaseB/dt = k_ge_fast*Stomach_fast_B + k_ge_slow*Stomach_slow_B - k_release*ReleaseB   ReleaseB(0) = 0
#'   dAbsorbB/dt = k_release * ReleaseB - ka * AbsorbB                                          AbsorbB(0) = 0
#'   dCentralB/dt = FbioB * ka * AbsorbB - kel * CentralB                                        CentralB(0) = 0
#'
#' `k_release` here is curve B's own post-gastric-emptying step (e.g. capsule
#' coating dissolution once past the stomach), now decoupled from the timing
#' of leaving the stomach itself, which `f_fast`/`k_ge_fast`/`k_ge_slow`
#' capture instead. `ka`, `kel` are shared across both curves exactly as in
#' [simulate_delayed_release()]; this only replaces the assumed *input*
#' shape, not the absorption/clearance structure.
#'
#' Not exported for direct use; called via [simulate_biphasic_emptying()].
#'
#' @keywords internal
biphasic_emptying_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dStomach_fast_A <- -k_ge_fast * Stomach_fast_A
    dStomach_slow_A <- -k_ge_slow * Stomach_slow_A
    dAbsorbA <- k_ge_fast * Stomach_fast_A + k_ge_slow * Stomach_slow_A - ka * AbsorbA
    dCentralA <- FbioA * ka * AbsorbA - kel * CentralA

    dStomach_fast_B <- -k_ge_fast * Stomach_fast_B
    dStomach_slow_B <- -k_ge_slow * Stomach_slow_B
    dReleaseB <- k_ge_fast * Stomach_fast_B + k_ge_slow * Stomach_slow_B - k_release * ReleaseB
    dAbsorbB <- k_release * ReleaseB - ka * AbsorbB
    dCentralB <- FbioB * ka * AbsorbB - kel * CentralB

    list(c(dStomach_fast_A, dStomach_slow_A, dAbsorbA, dCentralA,
           dStomach_fast_B, dStomach_slow_B, dReleaseB, dAbsorbB, dCentralB))
  })
}

#' Simulate two simultaneously-dosed oral curves sharing biphasic gastric
#' emptying, at a set of times
#'
#' See [biphasic_emptying_ode()] for the model. Curve A gets the instant/
#' single-compartment absorption structure (e.g. a liquid dose); curve B gets
#' the extra post-emptying `k_release` step (e.g. an enteric capsule).
#'
#' @param times Numeric vector of times at which to return concentration
#'   (same times used for both curves).
#' @param f_fast Fraction of each dose in the fast-emptying gastric pool (0-1).
#' @param k_ge_fast,k_ge_slow Fast/slow gastric-emptying rate constants (1/time).
#' @param ka,kel Shared absorption/elimination rate constants (1/time).
#' @param k_release Curve B's own post-emptying release rate constant (1/time).
#' @param FbioA,FbioB Bioavailable fraction of each dose (0-1).
#' @param doseA,doseB Administered dose of each curve, in mass units consistent with Vd.
#' @param Vd Volume of distribution (shared - same subject, same visit).
#' @return A list with `conc_A`, `conc_B` - numeric vectors of concentration
#'   at `times`, aligned to it (t = 0 is added internally if not already
#'   present, but only requested times are returned).
simulate_biphasic_emptying <- function(times, f_fast, k_ge_fast, k_ge_slow, ka, kel,
                                        k_release, FbioA, FbioB, doseA, doseB, Vd) {
  solve_times <- sort(unique(c(0, times)))
  state0 <- c(
    Stomach_fast_A = f_fast * doseA, Stomach_slow_A = (1 - f_fast) * doseA,
    AbsorbA = 0, CentralA = 0,
    Stomach_fast_B = f_fast * doseB, Stomach_slow_B = (1 - f_fast) * doseB,
    ReleaseB = 0, AbsorbB = 0, CentralB = 0
  )
  parms <- c(f_fast = f_fast, k_ge_fast = k_ge_fast, k_ge_slow = k_ge_slow,
             ka = ka, kel = kel, k_release = k_release, FbioA = FbioA, FbioB = FbioB)

  out <- deSolve::ode(y = state0, times = solve_times, func = biphasic_emptying_ode,
                       parms = parms, method = "lsoda")
  out <- tibble::as_tibble(out)

  list(
    conc_A = out$CentralA[match(times, out$time)] / Vd,
    conc_B = out$CentralB[match(times, out$time)] / Vd
  )
}

#' Split a dose into an immediate and a genuinely time-lagged second dose,
#' both absorbed/cleared through the same one-compartment kinetics
#'
#' Confirmed (pilot testing, see the caveat on [biphasic_emptying_ode()]) to
#' be able to reproduce an actual dip-then-rise: a fraction `f_delayed` of
#' the dose contributes nothing until `t_lag`, then behaves exactly like a
#' fresh dose given at that later time - mathematically equivalent to two
#' separately-timed doses, not two differently-paced continuous release
#' processes both starting at t=0 (which [biphasic_emptying_ode()] showed
#' cannot produce a real second hump). A real second gastric-emptying wave
#' triggered only after a delay (e.g. a transient duodenal-brake pause that
#' fully resolves before releasing the rest of the dose) is the motivating
#' mechanism, as opposed to two pools emptying at different constant rates
#' from the start.
#'
#' @param simulate_fn A single-dose simulator with signature
#'   `function(times, dose, ...)` returning a numeric concentration vector at
#'   `times` - e.g. a closure over [bateman_conc()] or
#'   [simulate_delayed_release()] with its other arguments (`ka`, `kel`,
#'   etc.) already fixed, leaving only `times` and `dose` free.
#' @param times Numeric vector of times at which to return concentration.
#' @param dose Total administered dose (mass units consistent with
#'   `simulate_fn`'s own Vd).
#' @param f_delayed Fraction of `dose` in the lagged portion (0-1); the
#'   remaining `1 - f_delayed` is absorbed starting at t=0 as usual.
#' @param t_lag Delay, in the same time units as `times`, before the lagged
#'   portion starts contributing at all.
#' @return Numeric vector of concentrations at `times` (immediate + lagged
#'   contributions summed). Before `t_lag`, the lagged portion is evaluated
#'   at time 0, which both `bateman_conc()` and `simulate_delayed_release()`
#'   correctly return as zero, so no explicit indicator/branch is needed.
simulate_lagged_dose <- function(simulate_fn, times, dose, f_delayed, t_lag) {
  immediate <- simulate_fn(times, dose = (1 - f_delayed) * dose)
  delayed   <- simulate_fn(pmax(times - t_lag, 0), dose = f_delayed * dose)
  immediate + delayed
}

#' Generalizes [simulate_lagged_dose()]: an onset lag for the FIRST wave too,
#' not just the second
#'
#' [simulate_lagged_dose()] assumes absorption of the non-delayed portion
#' starts exactly at t=0. That can't represent a curve whose very first real
#' sample is already low relative to what follows - evidence absorption
#' itself hadn't really started yet by then - a standard, separate PK concept
#' (an absorption lag time, sometimes called `ALAG` or `Tlag` in other PK
#' software) from the second-wave delay `simulate_lagged_dose()` targets.
#'
#' Parameterized as `t_lag1` (the first wave's own onset delay) plus `gap`
#' (an ADDITIONAL delay before the second wave, on top of `t_lag1`, so the
#' second wave's absolute start time is `t_lag1 + gap` and is always later
#' than the first wave's by construction - avoids the label-switching
#' ambiguity of fitting two independent, unordered lag times). Setting
#' `f_delayed` near 0 recovers a single-wave curve that still has its own
#' fittable onset lag (`t_lag1`) - deliberately so: an onset lag can be real
#' even for a subject with no second wave at all, and shouldn't require one
#' to be present to be fit.
#'
#' @param simulate_fn A single-dose simulator with signature
#'   `function(times, dose, ...)` - see [simulate_lagged_dose()].
#' @param times Numeric vector of times at which to return concentration.
#' @param dose Total administered dose.
#' @param f_delayed Fraction of `dose` in the second wave (0-1).
#' @param t_lag1 Onset delay for the first (majority, if `f_delayed` < 0.5)
#'   wave, in the same time units as `times`.
#' @param gap Additional delay, on top of `t_lag1`, before the second wave
#'   starts contributing (>= 0) - so the second wave's absolute start is
#'   `t_lag1 + gap`.
#' @return Numeric vector of concentrations at `times` (both waves summed).
simulate_two_lag_dose <- function(simulate_fn, times, dose, f_delayed, t_lag1, gap) {
  wave1 <- simulate_fn(pmax(times - t_lag1, 0), dose = (1 - f_delayed) * dose)
  wave2 <- simulate_fn(pmax(times - (t_lag1 + gap), 0), dose = f_delayed * dose)
  wave1 + wave2
}

#' Estimate total blood volume from weight, height, and sex (Nadler 1962)
#'
#' Nadler DA, Hidalgo JU, Bloch T. Prediction of blood volume in normal
#' human adults. Surgery. 1962;51(2):224-232. Coefficients differ by sex;
#' valid for adults and children over ~35 kg.
#'
#' A body-composition estimate like this is one candidate for a study's Vd
#' assumption (volume of distribution), but not the only one - whether
#' total blood volume is the right physiological quantity for a given
#' analyte depends on the study; see the calling script's own Vd rationale.
#'
#' @param weight_kg Body weight in kg.
#' @param height_cm Height in cm.
#' @param sex Character vector, "male"/"female" (case-insensitive) per
#'   observation - which coefficient set to use.
#' @return Estimated blood volume, in liters.
nadler_blood_volume <- function(weight_kg, height_cm, sex) {
  height_m <- height_cm / 100
  is_male  <- tolower(sex) %in% c("male", "m")
  k1 <- ifelse(is_male, 0.3669, 0.3561)
  k2 <- ifelse(is_male, 0.03219, 0.03308)
  k3 <- ifelse(is_male, 0.6041, 0.1833)
  k1 * height_m^3 + k2 * weight_kg + k3
}
