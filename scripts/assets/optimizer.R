# =============================================================================
# optimizer.R
#
# Optimization engine for fitting PK models to observed concentration curves -
# single-curve or jointly across several curves that share some parameters and
# not others. It operates on plain "curve" lists (see `build_objective()`)
# and named parameter vectors; model_12C.R and model_13C6.R supply the
# `simulate` closures, bounds, and doses. Also holds the penalty settings every
# fit uses and the fit-quality helpers (R², AIC).
#
# WHY MULTI-START: first-order PK objective surfaces are routinely
# multimodal (e.g. a fast-absorption/fast-clearance solution can fit nearly
# as well as a slow/slow one). A single optim() call from one starting point
# is not a defensible fit; always use fit_multistart() with a diverse set of
# seeds (informed + systematic grid + random) rather than optim() directly.
#
# WHY PER-TASK SEEDING, NOT ONE set.seed() BEFORE A PARALLEL LOOP: a single
# set.seed() before parallel::mclapply() does NOT make the random starting
# points reproducible (verified: same script, same seed, different draws on two
# runs), and even RNGkind("L'Ecuyer-CMRG") only fixes that for a fixed core
# count and task order. string_seed() sidesteps this: call
# set.seed(string_seed(<task's own stable ID>)) inside the worker, so a task's
# draws depend only on its ID, not on core count, task order or other tasks.
# =============================================================================

# Penalty settings of every fit (see build_objective()): the Tmax floor (min), the Cmax band (fraction), their
# weights, and the fine time grid the penalties are evaluated on (min). See "Fitting procedure" in docs/pk-model.md.
MIN_TMAX <- 10; CMAX_TOL <- 0.10; CMAX_LAMBDA <- 20; TMAX_LAMBDA <- 50
FINE_T <- seq(0, 400, length.out = 50)

#' Coefficient of determination (R^2) between observed and predicted values
#'
#' @param observed,predicted Numeric vectors of equal length.
#' @return Scalar R^2. Can be negative if `predicted` fits worse than the
#'   mean of `observed`.
r_squared <- function(observed, predicted) {
  1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
}

#' AIC from residual sums of squares, with a residual variance per curve
#'
#' @param rss,n Residual sum of squares and number of points, one per curve.
#' @param k Number of fitted parameters.
aic_rss <- function(rss, n, k) sum(n * log(rss / n)) + 2 * k

#' A fitted parameter, or NA when the fit is NULL or has no such parameter
#'
#' @param fit A list with `par` (a named numeric vector), or NULL.
par_or_na <- function(fit, name) if (is.null(fit) || !(name %in% names(fit$par))) NA_real_ else fit$par[[name]]

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

#' Build a normalized, proportionally-weighted multi-curve objective function
#'
#' Each curve contributes a weighted analog of `sse / ss_tot` (i.e. 1 - that
#' curve's own weighted R^2) to the total, rather than raw squared error.
#' Two separate normalizations are combined here:
#'
#'   1. ACROSS curves: dividing by each curve's own (weighted) total
#'      variance is essential whenever curves differ in concentration scale
#'      (e.g. a tracer dose vs. a much larger unlabeled dose) - without it,
#'      the larger-scale curve's error dominates the objective and the
#'      optimizer effectively ignores the smaller curve.
#'   2. WITHIN a curve: points are weighted `1 / max(pred, floor)^2`
#'      (proportional/constant-CV weighting, floored to avoid blow-up near
#'      zero). Unweighted SSE lets the peak region dominate a curve's own fit
#'      when residual scale grows with concentration, under-weighting the tail.
#'      `floor` is `weight_floor_frac` times that curve's own observed Cmax.
#'
#' Because weights depend on `pred`, which changes every evaluation, the
#' within-curve normalizer (`weighted_ss_tot`) is recomputed at every call
#' rather than precomputed once from the observed data alone.
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
#' @param weight_floor_frac Floor on the within-curve weighting, as a
#'   fraction of that curve's own observed Cmax.
#' @param proportional_weighting If TRUE (default), apply the within-curve
#'   `1/pred^2` weighting described above. If FALSE, use plain per-curve SSE
#'   normalized only by that curve's own total variance - every observed point
#'   counts equally, including the peak/early region the weighted version
#'   discounts. Use FALSE when a curve's most informative feature is its
#'   peak/early shape (e.g. a sharp early spike a sparse tail can't pin down).
#'   See "Fitting procedure" in docs/pk-model.md for why the ERIE fits use FALSE.
#' @return A function(theta) -> scalar objective value to minimize.
build_objective <- function(curves, min_tmax = NULL, cmax_tol = NULL, cmax_lambda = 20,
                                   tmax_lambda = 50, weight_floor_frac = 0.01,
                                   proportional_weighting = TRUE) {
  function(theta) {
    total <- 0
    for (cv in curves) {
      pred <- tryCatch(cv$simulate(theta), error = function(e) NA_real_)
      if (length(pred) != length(cv$conc) || any(!is.finite(pred))) return(1e10)

      if (proportional_weighting) {
        floor_val <- weight_floor_frac * max(cv$conc)
        w <- 1 / pmax(pred, floor_val)^2
        w_mean_obs <- sum(w * cv$conc) / sum(w)
        weighted_ss_tot <- sum(w * (cv$conc - w_mean_obs)^2)
        total <- total + sum(w * (pred - cv$conc)^2) / weighted_ss_tot
      } else {
        ss_tot <- sum((cv$conc - mean(cv$conc))^2)
        total <- total + sum((pred - cv$conc)^2) / ss_tot
      }

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
#' values at the escape-value ceiling used by [build_objective()]) are
#' silently skipped rather than aborting the whole fit.
#'
#' @param objective_fn function(theta) -> scalar, e.g. from
#'   [build_objective()].
#' @param lower,upper Named numeric bound vectors.
#' @param seeds List of named numeric starting vectors (same names/order as
#'   `lower`/`upper`).
#' @param control List passed through to `optim()`. If it doesn't already
#'   set `parscale`, one is computed PER START from that start's own
#'   magnitude (see below) - L-BFGS-B's internal step sizing assumes
#'   roughly unit-scaled parameters, and fitted PK parameters routinely
#'   span very different magnitudes (e.g. rate constants ~1e-4-1 alongside
#'   fractions 0-1) without it.
#' @return The best `optim()` result (a list with `par`, `value`, ...), or
#'   NULL if every start failed.
fit_multistart <- function(objective_fn, lower, upper, seeds, control = list(maxit = 100)) {
  bound_range <- upper - lower

  best <- NULL
  for (start in seeds) {
    start <- pmax(pmin(start, upper * 0.99), lower * 1.01)
    start_control <- control
    if (is.null(start_control$parscale)) {
      # parscale from each start's own magnitude (floored at 1% of the bound
      # range), not from upper - lower: when bounds are wide
      # relative to typical values, bound-width scaling is far too large and
      # L-BFGS-B's step sizing gets stuck near the start. A multistart seed is
      # meant to be in the right neighborhood, so its own magnitude is the better
      # local-scale estimate. See "Parameter scaling" under "Fitting
      # procedure" in docs/pk-model.md.
      start_control$parscale <- pmax(abs(start), 0.01 * bound_range)
    }
    fit <- tryCatch(
      stats::optim(start, objective_fn, method = "L-BFGS-B",
                   lower = lower, upper = upper, control = start_control),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$value) && fit$value < 1e6 &&
        (is.null(best) || fit$value < best$value)) {
      best <- fit
    }
  }
  best
}

#' Seeds for a targeted retry near a fit's flagged bound(s)
#'
#' A boundary-flagged fit (a parameter landing essentially at its lower or
#' upper bound) is ambiguous: either that bound is a real constraint, or the
#' standard search just missed a better solution sitting inside it. This
#' builds a seed set concentrated on exactly that possibility: the prior
#' estimate itself (so a retry can never end up worse than not retrying),
#' a systematic grid spanning each flagged parameter's own bound range
#' (holding every other parameter at its prior value), and `n_random` more
#' random starts across the full bound space.
#'
#' @param prior_par Named numeric vector, the fit being retried.
#' @param bounds Named list, one `c(lower, upper)` per parameter (the same
#'   set `fit_multistart()` was called with).
#' @param flagged_params Character vector of parameter names currently at a
#'   bound - only these get a dedicated grid.
#' @param grid_fracs Fractions of each flagged parameter's bound range to
#'   grid over.
#' @param n_random Number of additional random starts.
#' @param log_scale Passed through to [random_seeds()].
#' @return A list of named numeric vectors, ready for `fit_multistart()`.
retry_seeds_near_bounds <- function(prior_par, bounds, flagged_params,
                                     grid_fracs = c(0.05, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 0.99),
                                     n_random = 80, log_scale = character(0)) {
  grid <- unlist(lapply(flagged_params, function(p) {
    vals <- bounds[[p]][1] + diff(bounds[[p]]) * grid_fracs
    lapply(vals, function(v) { par <- prior_par; par[[p]] <- v; par })
  }), recursive = FALSE)
  c(list(prior_par), grid, random_seeds(n_random, bounds, log_scale = log_scale))
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
