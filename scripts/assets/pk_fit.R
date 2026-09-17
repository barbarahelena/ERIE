# =============================================================================
# pk_fit.R
#
# Generic optimization engine for fitting PK models to observed concentration
# curves - single-curve or jointly across several curves that share some
# parameters and not others. Nothing here is study-specific: it operates on
# plain "curve" lists (see `build_joint_objective()`) and named parameter
# vectors. Reuse this file as-is in other projects; the analysis script
# supplies the model-specific `simulate` closures, bounds, and doses.
#
# WHY MULTI-START: first-order PK objective surfaces are routinely
# multimodal (e.g. a fast-absorption/fast-clearance solution can fit nearly
# as well as a slow/slow one). A single optim() call from one starting point
# is not a defensible fit; always use fit_multistart() with a diverse set of
# seeds (informed + systematic grid + random) rather than optim() directly.
#
# WHY PER-TASK SEEDING, NOT ONE set.seed() BEFORE A PARALLEL LOOP: a single
# set.seed() call before parallel::mclapply() does NOT make the random
# starting points reproducible. mclapply reseeds each forked worker using
# something not deterministically derived from that seed (verified: the
# same script, same set.seed(), produces different random draws on two
# separate runs). Even the documented fix (RNGkind("L'Ecuyer-CMRG")) only
# fixes *that* symptom - which stream a given task draws from is still tied
# to its fork order/position, so changing core count or task order can
# silently change which random seeds a given item's fit uses. string_seed()
# sidesteps all of this: call set.seed(string_seed(<task's own stable ID>))
# inside the worker function, once per task, and its random draws depend
# only on that ID - invariant to core count, task order, and how many other
# tasks exist.
# =============================================================================

#' Deterministic integer seed from a string, for reproducible per-task
#' set.seed() calls inside a parallel loop (see note above)
#'
#' Not cryptographic - just needs to decorrelate a modest number of short
#' keys (e.g. "<subject>_<visit>") well enough for seeding `set.seed()`.
#' A standard base-31 polynomial rolling hash, computed in double precision
#' and reduced mod a prime at every step to avoid integer overflow.
#'
#' @param key A single string (e.g. paste(subject_id, visit)).
#' @return An integer suitable for set.seed().
string_seed <- function(key) {
  codes <- utf8ToInt(key)
  h <- Reduce(function(acc, code) (acc * 31 + code) %% 2147483647, codes, init = 7)
  as.integer(h)
}

#' Build a normalized multi-curve objective function
#'
#' Each curve contributes `sse / ss_tot` (i.e. 1 - that curve's own R^2) to
#' the total, rather than raw squared error. This is essential whenever
#' curves differ in concentration scale (e.g. a tracer dose vs. a much
#' larger unlabeled dose): without normalizing, the larger-scale curve's
#' error dominates the objective and the optimizer effectively ignores the
#' smaller curve. Do not combine raw SSE across curves of different scale.
#'
#' Optionally adds two penalty terms per curve, evaluated on a fine time
#' grid (not just the observed sampling times) so they can catch degenerate
#' solutions that look fine only at the sparse observed points:
#'   - a soft floor on predicted Tmax (discourages spurious early-spike fits
#'     where the peak occurs before absorption could plausibly complete)
#'   - a soft band around the curve's own observed Cmax (discourages
#'     systematic over/undershoot of the peak without rigidly constraining
#'     the rest of the curve)
#' Both are smooth (a squared shortfall/excess, zero at and below the
#' threshold) rather than a discontinuous jump - `optim()`'s L-BFGS-B relies
#' on finite-difference gradients, which a hard penalty boundary makes
#' needlessly rough to search near.
#'
#' @param curves A list of curve specifications. Each element must have:
#'   - `times`: numeric vector of observed times
#'   - `conc`: numeric vector of observed concentrations, same length
#'   - `simulate`: function(theta) -> predicted concentration at `times`
#'   and may optionally have:
#'   - `fine_simulate`: function(theta) -> list(time = ..., conc = ...) on a
#'     dense time grid, used only for the Tmax/Cmax penalties below
#' @param min_tmax Soft floor on predicted Tmax (time units matching the
#'   data), or NULL to disable. Requires `fine_simulate` on each curve.
#' @param cmax_tol Fractional tolerance band around each curve's own
#'   observed Cmax (e.g. 0.10 = +/-10%), or NULL to disable.
#' @param cmax_lambda Penalty weight applied to Cmax band violations.
#' @param tmax_lambda Penalty weight applied to Tmax shortfalls.
#' @return A function(theta) -> scalar objective value to minimize.
build_joint_objective <- function(curves, min_tmax = NULL, cmax_tol = NULL, cmax_lambda = 20, tmax_lambda = 50) {
  curves <- lapply(curves, function(cv) {
    if (is.null(cv$ss_tot)) cv$ss_tot <- sum((cv$conc - mean(cv$conc))^2)
    cv
  })

  function(theta) {
    total <- 0
    for (cv in curves) {
      pred <- tryCatch(cv$simulate(theta), error = function(e) NA_real_)
      if (length(pred) != length(cv$conc) || any(!is.finite(pred))) return(1e10)
      total <- total + sum((pred - cv$conc)^2) / cv$ss_tot

      if (!is.null(cv$fine_simulate) && (!is.null(min_tmax) || !is.null(cmax_tol))) {
        fine <- tryCatch(cv$fine_simulate(theta), error = function(e) NULL)
        if (is.null(fine) || any(!is.finite(fine$conc))) return(1e10)

        if (!is.null(min_tmax)) {
          tmax <- fine$time[which.max(fine$conc)]
          shortfall <- max(0, min_tmax - tmax)
          total <- total + tmax_lambda * shortfall^2
        }
        if (!is.null(cmax_tol)) {
          obs_cmax <- max(cv$conc)
          pred_cmax <- max(fine$conc)
          excess <- max(0, abs(pred_cmax - obs_cmax) / obs_cmax - cmax_tol)
          total <- total + cmax_lambda * excess^2
        }
      }
    }
    total
  }
}

#' Multi-start bounded optimization
#'
#' Runs `optim(method = "L-BFGS-B")` from every seed in `seeds` and returns
#' the best (lowest-objective) result. Failed starts (errors, or objective
#' values at the escape-value ceiling used by [build_joint_objective()]) are
#' silently skipped rather than aborting the whole fit.
#'
#' @param objective_fn function(theta) -> scalar, e.g. from
#'   [build_joint_objective()].
#' @param lower,upper Named numeric bound vectors.
#' @param seeds List of named numeric starting vectors (same names/order as
#'   `lower`/`upper`).
#' @param control List passed through to `optim()`.
#' @return The best `optim()` result (a list with `par`, `value`, ...), or
#'   NULL if every start failed.
fit_multistart <- function(objective_fn, lower, upper, seeds, control = list(maxit = 100)) {
  best <- NULL
  for (start in seeds) {
    start <- pmax(pmin(start, upper * 0.99), lower * 1.01)
    fit <- tryCatch(
      stats::optim(start, objective_fn, method = "L-BFGS-B",
                   lower = lower, upper = upper, control = control),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < 1e6 &&
        (is.null(best) || fit$value < best$value)) {
      best <- fit
    }
  }
  best
}

#' Generate random starting seeds, optionally log-uniform per parameter
#'
#' Rate constants are usually better sampled log-uniformly (equal weight per
#' order of magnitude) than uniformly; fractions/proportions are usually
#' fine sampled uniformly.
#'
#' @param n Number of seeds to generate.
#' @param bounds Named list, one `c(lower, upper)` per parameter.
#' @param log_scale Character vector of parameter names to sample
#'   log-uniformly (must have lower > 0 for those parameters).
#' @return A list of `n` named numeric vectors.
random_seeds <- function(n, bounds, log_scale = character(0)) {
  param_names <- names(bounds)
  replicate(n, {
    stats::setNames(vapply(param_names, function(p) {
      lo <- bounds[[p]][1]; hi <- bounds[[p]][2]
      if (p %in% log_scale) exp(stats::runif(1, log(lo), log(hi))) else stats::runif(1, lo, hi)
    }, numeric(1)), param_names)
  }, simplify = FALSE)
}

#' Generate systematic grid seeds, crossed over a subset of parameters
#'
#' The remaining parameters (not in `grid`) are held fixed at the values
#' given in `fixed`. Useful for densely covering a region of concern (e.g.
#' near a bound that a previous fit got stuck at) without the combinatorial
#' blow-up of gridding every parameter.
#'
#' @param grid Named list, one numeric vector of grid values per parameter.
#' @param fixed Named list of fixed values for the other parameters.
#' @return A list of named numeric vectors, one per grid combination.
grid_seeds <- function(grid, fixed = list()) {
  combos <- expand.grid(grid, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  lapply(seq_len(nrow(combos)), function(i) {
    unlist(c(as.list(combos[i, , drop = FALSE]), fixed))
  })
}
