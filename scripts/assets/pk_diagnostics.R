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

#' Whether a curve's own raw observations show a post-peak dip followed by a
#' meaningful rise
#'
#' A simple, threshold-free check (no curve fitting involved) for whether
#' the data itself gives any direct evidence of a genuine second wave,
#' meant to gate a more complex model BEFORE it's even attempted - fitting
#' first and relying on AIC/R^2 to notice afterward that there was nothing
#' to explain lets a flexible model "explain" ordinary sampling noise around
#' a ordinary single peak, not just a genuine second one.
#'
#' @param times,conc Numeric vectors of equal length, one curve's own
#'   observed sampling times and concentrations (not necessarily sorted).
#' @param min_rise_frac After the post-peak trough, require a subsequent
#'   rise of at least this fraction of the curve's own overall peak
#'   concentration (default 0.15, i.e. 15%) before counting it as a genuine
#'   second wave rather than noise around a flat decline.
#' @param min_first_peak_frac The candidate "first wave" peak itself must be
#'   at least this fraction of the curve's own overall peak (default 0.5) -
#'   otherwise a tiny early wobble while concentration is still low and
#'   genuinely on its way up to the real peak (ordinary sampling noise)
#'   would get mistaken for a first wave.
#' @return TRUE if a first peak, a later local trough, and a later rise of
#'   at least `min_rise_frac * overall peak` all exist in that order; FALSE
#'   otherwise.
has_peak_dip_rise <- function(times, conc, min_rise_frac = 0.15, min_first_peak_frac = 0.5) {
  ord <- order(times)
  t <- times[ord]; c <- conc[ord]
  n <- length(c)
  if (n < 4) return(FALSE)
  overall_peak <- max(c)

  # The first "wave 1" candidate: the first point after which concentration
  # turns down, that's still a substantial fraction of the curve's overall
  # peak. Deliberately NOT anchored on the curve's global max (unlike an
  # earlier version of this function) - a genuine second wave is often
  # TALLER than the first (see docs/pk-model.md, e.g. ER01/ER09), so
  # requiring the trough to come after the tallest point structurally
  # misses exactly that case (there's nothing left to rise back up to).
  candidates <- which(diff(c) < 0 & c[-n] >= min_first_peak_frac * overall_peak)
  if (length(candidates) == 0) return(FALSE)
  peak_idx <- candidates[1]
  if (peak_idx >= n - 1) return(FALSE)   # need a trough AND a later rise point after the peak

  # The first LOCAL trough after that peak - the first point after which
  # concentration starts rising again - not the post-peak segment's global
  # minimum, which on any normal declining curve is almost always its very
  # last (most-decayed) point, not the dip actually being looked for.
  post_peak <- c[(peak_idx + 1):n]
  m <- length(post_peak)
  falls <- diff(post_peak) < 0
  trough_rel <- which(!falls)[1]
  if (is.na(trough_rel) || trough_rel >= m) return(FALSE)
  trough_val <- post_peak[trough_rel]
  after_trough <- post_peak[(trough_rel + 1):m]
  (max(after_trough) - trough_val) >= min_rise_frac * overall_peak
}
