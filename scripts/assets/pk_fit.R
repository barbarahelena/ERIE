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
#'      zero). Unweighted SSE lets the peak region dominate a single curve's
#'      own fit - checking residuals from an unweighted fit against
#'      predicted concentration showed squared-residual scale differing
#'      ~24x between the top and bottom quartile of predicted concentration
#'      (see docs/pk-model.md), which under-weights the tail even after the
#'      cross-curve normalization above. `floor` is `weight_floor_frac`
#'      times that curve's own observed Cmax.
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
#' @return A function(theta) -> scalar objective value to minimize.
build_joint_objective <- function(curves, min_tmax = NULL, cmax_tol = NULL, cmax_lambda = 20,
                                   tmax_lambda = 50, weight_floor_frac = 0.01) {
  function(theta) {
    total <- 0
    for (cv in curves) {
      pred <- tryCatch(cv$simulate(theta), error = function(e) NA_real_)
      if (length(pred) != length(cv$conc) || any(!is.finite(pred))) return(1e10)

      floor_val <- weight_floor_frac * max(cv$conc)
      w <- 1 / pmax(pred, floor_val)^2
      w_mean_obs <- sum(w * cv$conc) / sum(w)
      weighted_ss_tot <- sum(w * (cv$conc - w_mean_obs)^2)
      total <- total + sum(w * (pred - cv$conc)^2) / weighted_ss_tot

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
#' @param control List passed through to `optim()`. If it doesn't already
#'   set `parscale`, defaults to `upper - lower` per parameter - L-BFGS-B's
#'   internal step sizing assumes roughly unit-scaled parameters, and
#'   fitted PK parameters routinely span very different magnitudes (e.g.
#'   rate constants ~1e-4-1 alongside fractions 0-1) without it.
#' @return The best `optim()` result (a list with `par`, `value`, ...), or
#'   NULL if every start failed.
fit_multistart <- function(objective_fn, lower, upper, seeds, control = list(maxit = 100)) {
  if (is.null(control$parscale)) control$parscale <- upper - lower

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

#' Adaptive retry pass: densely re-search any still-flagged fit, keep only
#' if strictly better
#'
#' Generic driver for the failure mode `retry_seeds_near_bounds()` targets -
#' a fit stuck in a worse local optimum than one the standard search should
#' have found (see a calling project's own docs for the specific incident
#' that motivated this). A boundary flag is only ONE symptom of that
#' failure, and an unreliable one: it only shows up when the worse local
#' optimum happens to sit exactly on a bound. A fit stuck in a worse local
#' optimum that lands comfortably *inside* the bounds looks unremarkable -
#' no boundary flag - but is just as wrong, and typically shows up instead
#' as a poor fit (low R^2) or a search that didn't converge - hence
#' `extra_flag_cols`/`convergence_col` below, not boundary flags only.
#'
#' Retries every row flagged by any of `bound_cols`, `extra_flag_cols`, or
#' `convergence_col`, and replaces that row only if the retry's
#' `objective_value` is strictly lower - so this can only improve `results`,
#' never make it worse, however many rows are retried. A caller supplying
#' only `bound_cols` still works (boundary-only retry), but should also
#' pass its poor-fit/non-convergence columns where it has them, for the
#' reason above.
#'
#' @param results A data frame with one row per fit, including an
#'   `objective_value` column, the parameter columns named in `par_cols`,
#'   and every column named in `bound_cols`/`extra_flag_cols`/`convergence_col`.
#' @param bound_cols Named character vector: each name is one of `results`'
#'   boundary-flag columns, each value the parameter name (matching
#'   `bounds`/`par_cols`) that column is about - e.g.
#'   `c(kel_at_bound = "kel", k_release_at_bound = "k_release")`. Triggers a
#'   retry AND decides which parameter(s) get a targeted grid for it.
#' @param extra_flag_cols Character vector of other `results` columns that
#'   should also trigger a retry when TRUE (e.g. `"r2_12C_low"`), without
#'   implying any particular parameter to grid over - a row retried only for
#'   one of these gets the prior estimate plus random seeds, no grid.
#' @param convergence_col A single `results` column name that should also
#'   trigger a retry when FALSE (e.g. `"converged"`), or NULL to skip this.
#' @param bounds Named list, one `c(lower, upper)` per fitted parameter.
#' @param par_cols Character vector naming which `results` columns hold the
#'   fitted parameter values (same names as `bounds`).
#' @param refit_fn `function(i, seeds, maxit)` - refit just row `i` of
#'   `results` with the given retry `seeds` and `maxit`, and return a
#'   one-row data frame shaped like a row of `results`. This is the
#'   caller's own per-item fitting function (e.g. a project's
#'   `fit_subject_visit()`), so this stays a driver, not a second,
#'   separately-maintained fitting implementation.
#' @param maxit_retry `maxit` passed to `refit_fn` for the retry.
#' @param log_scale Passed through to `retry_seeds_near_bounds()`.
#' @param mc.cores Cores to retry flagged rows across
#'   (`parallel::mclapply()`); default 1 (sequential, works everywhere).
#' @param verbose Print a one-line progress/outcome summary.
#' @return `results`, with any strictly-improved rows replaced, plus two new
#'   columns: `retried` (TRUE for every row this pass attempted) and
#'   `retry_improved` (TRUE if the retry replaced the row, FALSE if it was
#'   attempted but didn't beat the original, NA if never retried). A row
#'   that was retried but not improved is informative on its own - it's
#'   evidence the original result wasn't an under-searched local optimum,
#'   e.g. a genuinely weak/low-information curve (rather than a search
#'   failure) is expected to stay flagged even after a denser retry.
adaptive_retry <- function(results, bound_cols, bounds, par_cols, refit_fn,
                            extra_flag_cols = character(0), convergence_col = NULL,
                            maxit_retry = 60, log_scale = character(0),
                            mc.cores = 1, verbose = TRUE) {
  results$retried <- FALSE
  results$retry_improved <- NA

  flag_cols <- c(names(bound_cols), extra_flag_cols)
  is_flagged <- if (length(flag_cols)) Reduce(`|`, lapply(flag_cols, function(col) results[[col]] %in% TRUE)) else FALSE
  not_converged <- if (!is.null(convergence_col)) results[[convergence_col]] %in% FALSE else FALSE
  retry_mask <- is_flagged | not_converged

  flagged <- which(retry_mask)
  if (length(flagged) == 0) return(results)
  if (verbose) cat("\nRetrying", length(flagged), "flagged fit(s) (boundary and/or poor fit) with a denser search...\n")

  retry_one <- function(i) {
    prior_par <- stats::setNames(as.numeric(results[i, par_cols]), par_cols)
    flagged_params <- unname(bound_cols[vapply(names(bound_cols), function(col) isTRUE(results[[col]][i]), logical(1))])
    seeds <- retry_seeds_near_bounds(prior_par, bounds, flagged_params, log_scale = log_scale)
    refit_fn(i, seeds, maxit_retry)
  }
  retry_results <- parallel::mclapply(flagged, retry_one, mc.cores = mc.cores)

  n_improved <- 0
  for (j in seq_along(flagged)) {
    i <- flagged[j]; new_row <- retry_results[[j]]
    results$retried[i] <- TRUE
    improved <- !is.na(new_row$objective_value) && new_row$objective_value < results$objective_value[i]
    results$retry_improved[i] <- improved
    # new_row only carries the fit columns (ka, kel, ..., objective_value),
    # not retried/retry_improved, so this can't clobber what was just set above.
    if (improved) {
      results[i, names(new_row)] <- new_row
      n_improved <- n_improved + 1
    }
  }
  if (verbose) cat("Retry improved", n_improved, "of", length(flagged), "flagged fit(s).\n")
  results
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
