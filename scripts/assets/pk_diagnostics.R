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

#' Whether a curve's own raw observations show a post-peak deviation from a
#' smooth decline - a plateau, wobble, or outright second rise
#'
#' A simple, threshold-free check (no curve fitting involved) for whether
#' the data itself gives any direct evidence of a genuine second wave,
#' meant to gate a more complex model BEFORE it's even attempted - fitting
#' first and relying on AIC/R^2 to notice afterward that there was nothing
#' to explain lets a flexible model "explain" ordinary sampling noise around
#' an ordinary single peak, not just a genuine second one.
#'
#' The test is NOT "does some later point exceed some earlier local
#' minimum" - that misses a plateau or slowed decline that never actually
#' rises above what came immediately before it, but is still real evidence
#' of extra absorption countering ongoing elimination at that moment. A
#' genuine single-exponential elimination phase is convex (decays ever more
#' slowly, never abruptly), so any real point on it always sits AT OR BELOW
#' the straight line connecting its two neighbors - that's what convexity
#' means. A point sitting MEASURABLY ABOVE that interpolated line is
#' therefore inherent evidence of non-monotonic-decay behavior, whether or
#' not it becomes an outright new local maximum - so every point after the
#' first peak is checked against its own neighbors' interpolation, not just
#' compared to the single deepest trough.
#'
#' @param times,conc Numeric vectors of equal length, one curve's own
#'   observed sampling times and concentrations (not necessarily sorted).
#' @param min_rise_frac A point counts as evidence of a second wave if it
#'   sits at least this fraction ABOVE the value linearly interpolated
#'   (in time) between its immediate neighbors (default 0.15, i.e. 15%).
#' @param min_first_peak_frac The candidate "first wave" peak itself must be
#'   at least this fraction of the curve's own overall peak (default 0.5) -
#'   otherwise a tiny early wobble while concentration is still low and
#'   genuinely on its way up to the real peak (ordinary sampling noise)
#'   would get mistaken for a first wave.
#' @return A list with `detected` (TRUE if, after a first peak, any later
#'   point sits at least `min_rise_frac` above its neighbors' linear
#'   interpolation) and `trigger_time` (the time of the FIRST point that
#'   triggered it, or NA if `detected` is FALSE) - when multiple points
#'   deviate from a smooth decline (multiple humps/wobbles), this is the
#'   one with the LARGEST deviation, not just the first one found, since
#'   that's the strongest single piece of evidence and the most useful
#'   point to anchor a second-wave model's search on. The trigger time is
#'   useful for seeding that model's search near the real evidence rather
#'   than leaving it to find an unrelated local optimum, and `excess` (the
#'   deviation's own magnitude) is useful for judging how strong that
#'   evidence is - e.g. gating a "definitely two_wave" fast path for
#'   unambiguous cases (see fit_one() in 02_fit_erie_model.R).
peak_dip_rise_info <- function(times, conc, min_rise_frac = 0.15, min_first_peak_frac = 0.5) {
  none <- list(detected = FALSE, trigger_time = NA_real_, excess = NA_real_)
  ord <- order(times)
  t <- times[ord]; c <- conc[ord]
  n <- length(c)
  if (n < 4) return(none)
  overall_peak <- max(c)

  # The first "wave 1" candidate: the first point after which concentration
  # turns down, that's still a substantial fraction of the curve's overall
  # peak. Deliberately NOT anchored on the curve's global max - a genuine
  # second wave is often TALLER than the first (see docs/pk-model.md, e.g.
  # ER01/ER09), so requiring evidence to come after the tallest point would
  # structurally miss exactly that case.
  candidates <- which(diff(c) < 0 & c[-n] >= min_first_peak_frac * overall_peak)
  if (length(candidates) == 0) return(none)
  peak_idx <- candidates[1]
  # Require at least one point of separation between the first peak and any
  # candidate second peak - the point immediately after the peak (i =
  # peak_idx+1) is still describing the first wave's own decline shape
  # (how sharply it turns over), not a genuinely separate second wave.
  if (peak_idx >= n - 2) return(none)

  best_excess <- -Inf
  best_time <- NA_real_
  for (i in (peak_idx + 2):(n - 1)) {
    interp <- c[i - 1] + (c[i + 1] - c[i - 1]) * (t[i] - t[i - 1]) / (t[i + 1] - t[i - 1])
    if (!is.finite(interp) || interp <= 0) next
    excess <- (c[i] - interp) / interp
    if (excess > best_excess) { best_excess <- excess; best_time <- t[i] }
  }
  if (best_excess >= min_rise_frac) return(list(detected = TRUE, trigger_time = best_time, excess = best_excess))
  none
}

#' Whether a curve's own first real observation suggests a genuine onset
#' delay (absorption hadn't really started yet), not just an ordinarily
#' slow rise
#'
#' A different phenomenon from [peak_dip_rise_info()] - a delay before a
#' SINGLE wave starts at all (a standard PK "absorption lag time" concept),
#' not a second, separate wave. Simple and threshold-free: if the first
#' observed point is only a small fraction of the curve's own eventual
#' peak, that's consistent with "nothing had happened yet" rather than an
#' ordinary gradual rise.
#'
#' @param times,conc Numeric vectors of equal length, one curve's own
#'   observed sampling times and concentrations.
#' @param min_first_frac The first point must be below this fraction of the
#'   curve's own peak to count as onset-lag evidence (default 0.25) -
#'   comfortably above both confirmed cases this was validated against
#'   (ER25 FCT1: 18.4%, ER06 FCT2: 2.2%).
#' @return TRUE if the first observed point is below `min_first_frac` of
#'   the curve's own peak.
has_onset_lag_evidence <- function(times, conc, min_first_frac = 0.25) {
  ord <- order(times)
  c <- conc[ord]
  if (length(c) < 3) return(FALSE)
  peak <- max(c)
  if (peak <= 0) return(FALSE)
  c[1] / peak < min_first_frac
}

#' Boolean-only wrapper around [peak_dip_rise_info()] - see there for details.
has_peak_dip_rise <- function(times, conc, min_rise_frac = 0.15, min_first_peak_frac = 0.5) {
  peak_dip_rise_info(times, conc, min_rise_frac, min_first_peak_frac)$detected
}
