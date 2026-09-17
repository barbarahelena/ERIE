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
