# =============================================================================
# curve_checks.R - checks on a curve's raw observations (no fitting) that decide which structures a model
# tries: a second wave, a near-equal neighbor of the peak, an onset lag. Used by model_12C.R and model_13C6.R.
# =============================================================================

#' Whether a curve's own raw observations show a post-peak deviation from a
#' smooth decline - a plateau, wobble, or outright second rise
#'
#' A simple check (no curve fitting involved) meant to gate a more complex model
#' BEFORE it's attempted: fitting first and relying on AIC/R^2 to notice there
#' was nothing to explain lets a flexible model fit ordinary sampling noise.
#'
#' The test is whether a post-peak point sits measurably ABOVE the straight line
#' joining its two neighbors. A genuine single-exponential elimination phase is
#' convex, so any real point on it lies at or below that chord; a point above it
#' is evidence of non-monotonic decay even if it never becomes a new local
#' maximum (a plateau or slowed decline). See "How the dip is detected" in
#' docs/pk-model.md.
#'
#' @param times,conc Numeric vectors of equal length, one curve's own
#'   observed sampling times and concentrations (not necessarily sorted).
#' @param min_rise_frac A point counts as evidence of a second wave if it
#'   sits at least this fraction ABOVE the value linearly interpolated
#'   (in time) between its immediate neighbors (default 0.15, i.e. 15%).
#' @param min_first_peak_frac The candidate "first wave" peak itself must be
#'   at least this fraction of the curve's own overall peak (default 0.5) -
#'   otherwise a tiny early wobble while concentration is still low and
#'   on its way up to the real peak (ordinary sampling noise)
#'   would get mistaken for a first wave.
#' @param max_time Only points up to this time count (e.g. the latest time a
#'   second wave allowed by the t_lag bound can become visible).
#' @param min_peak_frac Only rises of at least this fraction of the curve's
#'   peak count, so noise near zero in the tail cannot qualify.
#' @return A list with `detected` (TRUE if, after a first peak, any later point
#'   sits at least `min_rise_frac` above its neighbors' linear interpolation),
#'   `trigger_time` (the time of the point with the LARGEST such deviation, or NA)
#'   and `excess` (that deviation's magnitude, or NA). The trigger time is useful
#'   for seeding a second-wave model's search near the evidence, `excess` for
#'   judging how strong the evidence is.
peak_dip_rise_info <- function(times, conc, min_rise_frac = 0.15, min_first_peak_frac = 0.5,
                               max_time = Inf, min_peak_frac = -Inf) {
  none <- list(detected = FALSE, trigger_time = NA_real_, excess = NA_real_)
  ord <- order(times)
  t <- times[ord]; c <- conc[ord]
  n <- length(c)
  if (n < 4) return(none)
  overall_peak <- max(c)

  # The first "wave 1" candidate: the first point after which concentration
  # turns down, that's still a substantial fraction of the curve's overall
  # peak. Not anchored on the global max: a genuine second wave is
  # often taller than the first, so that would miss exactly that case.
  candidates <- which(diff(c) < 0 & c[-n] >= min_first_peak_frac * overall_peak)
  if (length(candidates) == 0) return(none)
  peak_idx <- candidates[1]
  # Require at least one point of separation between the first peak and any
  # candidate second peak - the point immediately after the peak (i =
  # peak_idx+1) is still describing the first wave's own decline shape
  # (how sharply it turns over), not a separate second wave.
  if (peak_idx >= n - 2) return(none)

  best_excess <- -Inf
  best_time <- NA_real_
  for (i in (peak_idx + 2):(n - 1)) {
    if (t[i] > max_time) next
    interp <- c[i - 1] + (c[i + 1] - c[i - 1]) * (t[i] - t[i - 1]) / (t[i + 1] - t[i - 1])
    if (!is.finite(interp) || interp <= 0 || (c[i] - interp) / overall_peak < min_peak_frac) next
    excess <- (c[i] - interp) / interp
    if (excess > best_excess) { best_excess <- excess; best_time <- t[i] }
  }
  if (best_excess >= min_rise_frac) return(list(detected = TRUE, trigger_time = best_time, excess = best_excess))
  none
}

#' Whether a curve's peak has a near-equal-height neighbor - a plateau or
#' two closely-overlapping waves, rather than one clean single maximum
#'
#' A different, complementary signature from [peak_dip_rise_info()]: two waves
#' close enough together in time don't necessarily produce a visible dip - they
#' can blend into a flat top or an irregular, non-smoothly-decelerating rise. A
#' genuine single-compartment absorption curve has one clean maximum whose
#' adjacent points are not usually THIS close to it. Catches cases
#' `peak_dip_rise_info()` structurally cannot (no trough exists to detect); see
#' "Choosing between single_wave and two_wave" in docs/pk-model.md.
#'
#' @param times,conc Numeric vectors of equal length, one curve's own
#'   observed sampling times and concentrations.
#' @param near_peak_frac A neighbor immediately before or after the peak
#'   counts as evidence if it's at least this fraction of the peak's own
#'   height (default 0.85, i.e. 85%).
#' @return TRUE if the point immediately before or immediately after the
#'   curve's own maximum (EXCLUDING t=30, see below) is at least
#'   `near_peak_frac` of the peak's height, and the maximum itself is not at
#'   t=30 (see below).
has_near_peak_neighbor <- function(times, conc, near_peak_frac = 0.85) {
  ord <- order(times)
  t <- times[ord]; c <- conc[ord]
  n <- length(c)
  if (n < 3) return(FALSE)
  peak_idx <- which.max(c)
  peak <- c[peak_idx]
  if (peak <= 0) return(FALSE)
  # A peak AT t=30 is excluded outright, not just as a neighbor (see below):
  # the true peak may sit anywhere at or before the first sample, so a nearly
  # as tall point at t=60 reflects an early single peak's slow decline over
  # that window, not evidence of a second wave arriving around t=60. Mirrors
  # the neighbor-at-30 exclusion below, for the mirror-image case.
  if (t[peak_idx] == 30) return(FALSE)
  neighbor_idx <- c(peak_idx - 1, peak_idx + 1)
  neighbor_idx <- neighbor_idx[neighbor_idx >= 1 & neighbor_idx <= n]
  # t=30 is excluded as a valid neighbor: it's the first real sample for
  # virtually every curve in this study, so being close to the peak there just
  # reflects fast absorption, not a second wave. (Study-specific: change it if
  # the first sample is at a different time.)
  neighbor_idx <- neighbor_idx[t[neighbor_idx] != 30]
  if (length(neighbor_idx) == 0) return(FALSE)
  any(c[neighbor_idx] / peak >= near_peak_frac)
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
#'   curve's own peak to count as onset-lag evidence (default 0.25; see "How the
#'   dip is detected" in docs/pk-model.md for the cases it was validated on).
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
