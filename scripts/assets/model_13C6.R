# =============================================================================
# model_13C6.R - the 13C6 model: ka and kel fixed at the 12C shared-kel fit, three candidate mechanisms for the
# tracer's own absorption, chosen by AIC. Used by 04_fit_13C6.R (fit) and 05_diet_summary_13C6.R (prediction).
# Mechanisms, bounds and selection: docs/pk-model.md, "The 13C6 tracer with fixed rates".
# =============================================================================
suppressPackageStartupMessages(library(purrr))
source("scripts/assets/pk_equations.R")
source("scripts/assets/optimizer.R")
source("scripts/assets/curve_checks.R")

BOUNDS_RELEASE <- list(F_13C6 = c(1e-4, 1), k_release = c(0.001, 1), f_delayed_13C6 = c(0.001, 0.999))
BOUNDS_ONSET   <- list(F_13C6 = c(1e-4, 1), t_lag1_13C6 = c(0, 90))
BOUNDS_WAVE    <- list(F_13C6 = c(1e-4, 1), f_delayed2_13C6 = c(0.05, 0.95), t_lag2_13C6 = c(30, 150))
N_RANDOM_13C6 <- 20; MAXIT_13C6 <- 60

#' 13C6 prediction with ka and kel fixed
#'
#' @param mech `release`, `onset_lag` or `second_wave`.
#' @param t Times (min).
#' @param th Named 13C6 parameters: F_13C6 plus the ones the mechanism uses.
#' @param r List or one-row data frame with `ka`, `kel`, `model_12C` and `t_lag_12C`.
#' @param dose 13C6 dose (mg). @param Vd Volume of distribution (L).
pred_13C6 <- function(mech, t, th, r, dose, Vd) {
  bolus <- function(tt, dose) bateman_conc(tt, r$ka, r$kel, th[["F_13C6"]], dose, Vd)
  switch(mech,
    release = {
      simB <- function(tt, dose) simulate_delayed_release(tt, th[["k_release"]], r$ka, r$kel, th[["F_13C6"]], dose, Vd)$conc
      if (r$model_12C == "two_wave") simulate_lagged_dose(simB, t, dose, th[["f_delayed_13C6"]], r$t_lag_12C) else simB(t, dose)
    },
    onset_lag   = bolus(pmax(t - th[["t_lag1_13C6"]], 0), dose),
    second_wave = simulate_lagged_dose(bolus, t, dose, th[["f_delayed2_13C6"]], th[["t_lag2_13C6"]]))
}

#' Fit one 13C6 curve: `release` always, `onset_lag` and `second_wave` when the raw data show evidence for them
#'
#' @param key Stable ID of the curve, for the random seeds.
#' @param times,conc The observed 13C6 curve (t > 0).
#' @param r The 12C fit of this visit, as in [pred_13C6()], plus `F_12C`.
#' @return NULL when nothing could be fitted, otherwise a list with `best` and `cands` (all candidates, by
#'   mechanism). Each candidate has `mech`, `par`, `rss`, `r2`, `aic` and `converged`.
fit_13C6 <- function(key, times, conc, r, dose, Vd) {
  set.seed(string_seed(key))
  n <- length(conc)
  if (n < 3) return(NULL)

  run <- function(mech, bounds, seeds) {
    curve <- list(times = times, conc = conc,
                  simulate = function(th) pred_13C6(mech, times, th, r, dose, Vd),
                  fine_simulate = function(th) list(time = FINE_T, conc = pred_13C6(mech, FINE_T, th, r, dose, Vd)))
    objective <- build_objective(list(curve), min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                 tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
    fit <- fit_multistart(objective, vapply(bounds, `[`, numeric(1), 1), vapply(bounds, `[`, numeric(1), 2),
                          seeds = seeds, control = list(maxit = MAXIT_13C6))
    if (is.null(fit)) return(NULL)
    pred <- curve$simulate(fit$par)
    rss <- sum((conc - pred)^2)
    list(mech = mech, par = fit$par, rss = rss, r2 = r_squared(conc, pred), aic = aic_rss(rss, n, length(bounds)),
         converged = fit$convergence == 0)
  }

  F0 <- r$F_12C
  two_wave <- r$model_12C == "two_wave"
  b_rel <- if (two_wave) BOUNDS_RELEASE else BOUNDS_RELEASE[c("F_13C6", "k_release")]
  seeds_rel <- c(
    unlist(lapply(c(0.01, 0.03, 0.1, 0.3), function(k) lapply(c(F0, 0.02), function(f)
      c(F_13C6 = f, k_release = k, f_delayed_13C6 = 0.3)[names(b_rel)])), recursive = FALSE),
    if (two_wave) list(c(F_13C6 = F0, k_release = 0.05, f_delayed_13C6 = 0.01), c(F_13C6 = F0, k_release = 0.05, f_delayed_13C6 = 0.99)),
    random_seeds(N_RANDOM_13C6, b_rel, log_scale = "k_release"))
  cands <- list(release = run("release", b_rel, seeds_rel))
  if (has_onset_lag_evidence(times, conc))
    cands$onset_lag <- run("onset_lag", BOUNDS_ONSET, lapply(c(10, 20, 30, 40, 50), function(tl) c(F_13C6 = F0, t_lag1_13C6 = tl)))
  if (has_peak_dip_rise(times, conc) || has_near_peak_neighbor(times, conc)) {
    sw <- run("second_wave", BOUNDS_WAVE, lapply(c(40, 60, 90, 120), function(tl) c(F_13C6 = F0, f_delayed2_13C6 = 0.5, t_lag2_13C6 = tl)))
    # pinned at the lower bound of f_delayed2_13C6 means no real second wave
    if (!is.null(sw) && sw$par[["f_delayed2_13C6"]] > BOUNDS_WAVE$f_delayed2_13C6[1] + 1e-3) cands$second_wave <- sw
  }
  cands <- compact(cands)
  if (length(cands) == 0) return(NULL)
  list(best = cands[[which.min(map_dbl(cands, "aic"))]], cands = cands)
}
