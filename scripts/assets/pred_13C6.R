# =============================================================================
# pred_13C6.R - 13C6 prediction with ka and kel fixed, shared by 07 and 08. Needs pk_curves.R.
# =============================================================================

#' @param mech `release`, `onset_lag` or `second_wave` (see 07_fit_13C6.R).
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
