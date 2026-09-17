# =============================================================================
# pk_diagnostics.R
#
# Generic, project-agnostic diagnostics for evaluating a fitted PK curve.
# Nothing here is study-specific.
# =============================================================================

#' Coefficient of determination (R^2) between observed and predicted values
#'
#' @param observed,predicted Numeric vectors of equal length.
#' @return Scalar R^2. Can be negative if `predicted` fits worse than the
#'   mean of `observed`.
r_squared <- function(observed, predicted) {
  1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
}

#' Time at which a simulated curve has cleared to a fraction of its own peak
#'
#' Searches only after the peak (so an early, sub-peak dip is not mistaken
#' for clearance). Intended for choosing a sensible x-axis upper limit when
#' plotting a fitted curve well past the last real observation, not for use
#' inside any fitting objective.
#'
#' @param fine_time,fine_conc Numeric vectors describing a simulated curve on
#'   a dense time grid (dense enough that the peak and post-peak decay are
#'   well resolved).
#' @param frac Fraction of the peak concentration considered "cleared"
#'   (default 1%).
#' @return The time of clearance, or `max(fine_time)` if the curve never
#'   drops below `frac` of its peak within the supplied grid.
time_to_clearance <- function(fine_time, fine_conc, frac = 0.01) {
  peak_idx <- which.max(fine_conc)
  post_time <- fine_time[peak_idx:length(fine_conc)]
  post_conc <- fine_conc[peak_idx:length(fine_conc)]
  hit <- which(post_conc <= frac * fine_conc[peak_idx])[1]
  if (is.na(hit)) return(max(fine_time))
  post_time[hit]
}
