# =============================================================================
# pk_equations.R - the PK equations. Used by model_12C.R, model_13C6.R and the simulations.
#
# Building blocks for first-order oral pharmacokinetic (PK) compartment
# models: doses, rate constants and volumes of distribution. Doses, bounds and
# Vd are set in the analysis scripts.
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
  # (exp(-kel t) - exp(-ka t)) / (ka - kel) written as
  # t * exp(-kel t) * (-expm1(-d t)) / (d t), d = ka - kel. This avoids the
  # cancellation of two nearly equal exponentials as ka -> kel and is smooth
  # through d = 0, where it equals t * exp(-ka t) (the limiting form).
  d <- ka - kel
  dt <- d * t
  shape <- ifelse(abs(dt) < 1e-10, 1, -expm1(-dt) / ifelse(dt == 0, 1, dt))
  Fbio * ka * dose * t * exp(-kel * t) * shape / Vd
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
#' Numerical integration is used instead of a closed-form
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

#' Split a dose into an immediate and a time-lagged second dose,
#' both absorbed/cleared through the same one-compartment kinetics
#'
#' Reproduces an actual dip-then-rise: a fraction `f_delayed` of the dose
#' contributes nothing until `t_lag`, then behaves exactly like a fresh dose
#' given at that later time - mathematically equivalent to two separately-timed
#' doses. Two differently-paced release processes that both start at t=0
#' cannot produce a real second hump (see "The post-peak dip and the
#' lagged-dose model" in docs/pk-model.md). A real second gastric-emptying wave
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
