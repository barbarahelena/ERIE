# =============================================================================
# tests.R
#
# Lightweight regression tests for the generic PK engine (pk_curves.R,
# pk_fit.R, pk_diagnostics.R). Base R (stopifnot()) rather than a testing
# framework, to keep scripts/assets/ dependency-light and copyable as-is
# into another project. Run directly: Rscript scripts/assets/tests.R
# (or `pixi run test-engine`).
#
# These exist because several of this engine's design choices were
# validated once, by hand, during development (see docs/pk-model.md) but not
# previously captured anywhere re-runnable:
#   - the closed-form Bateman equation's ka==kel degeneracy
#   - simulate_delayed_release() reducing to bateman_conc() for very fast
#     release, which is what "no dissolution delay" should mean
#   - string_seed()'s determinism and (for this project's actual keys)
#     distinctness, which the whole reproducibility argument in
#     docs/pk-model.md rests on
# =============================================================================

suppressMessages({
  library(tibble)
})

source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_fit.R")
source("scripts/assets/pk_diagnostics.R")

set.seed(1)   # deterministic test run - single-threaded here, so a plain set.seed() (not string_seed()) is fine

n_pass <- 0
check <- function(desc, expr) {
  ok <- isTRUE(tryCatch(expr, error = function(e) { message("ERROR in '", desc, "': ", conditionMessage(e)); FALSE }))
  if (!ok) stop("FAILED: ", desc)
  n_pass <<- n_pass + 1
  cat("PASS:", desc, "\n")
}

# ---- bateman_conc(): ka == kel removable singularity -----------------------

check("bateman_conc() is continuous across the ka==kel singularity", {
  t <- seq(1, 300, by = 1)
  kel <- 0.03
  c_exact  <- bateman_conc(t, ka = kel,            kel = kel, Fbio = 0.2, dose = 100, Vd = 10)
  c_close  <- bateman_conc(t, ka = kel * (1 + 1e-6), kel = kel, Fbio = 0.2, dose = 100, Vd = 10)
  max(abs(c_exact - c_close)) / max(c_exact) < 1e-3
})

# ---- bateman_tmax() matches the actual numerical peak of bateman_conc() ----

check("bateman_tmax() matches bateman_conc()'s numerical peak", {
  ka <- 0.05; kel <- 0.02
  analytic_tmax <- bateman_tmax(ka, kel)
  fine_t <- seq(0, 400, by = 0.1)
  numeric_tmax <- fine_t[which.max(bateman_conc(fine_t, ka, kel, 0.2, 100, 10))]
  abs(analytic_tmax - numeric_tmax) < 0.2
})

# ---- simulate_delayed_release() reduces to bateman_conc() for fast release -

check("simulate_delayed_release() ~ bateman_conc() as k_release -> fast (near-instant dissolution)", {
  t <- seq(10, 300, by = 10)
  ka <- 0.05; kel <- 0.02; Fbio <- 0.2; dose <- 100; Vd <- 10
  direct  <- bateman_conc(t, ka, kel, Fbio, dose, Vd)
  delayed <- simulate_delayed_release(t, k_release = 50, ka = ka, kel = kel, Fbio = Fbio, dose = dose, Vd = Vd)$conc
  max(abs(direct - delayed)) / max(direct) < 0.02
})

# ---- r_squared() ------------------------------------------------------------

check("r_squared() is 1 for a perfect prediction", {
  obs <- c(1, 5, 3, 8, 2)
  r_squared(obs, obs) > 1 - 1e-10
})

check("r_squared() is ~0 for predicting the observed mean", {
  obs <- c(1, 5, 3, 8, 2)
  abs(r_squared(obs, rep(mean(obs), length(obs)))) < 1e-10
})

# ---- time_to_clearance() ----------------------------------------------------

check("time_to_clearance() finds the point a decaying curve drops below frac*peak", {
  t <- seq(0, 100, by = 1)
  conc <- exp(-0.1 * t) * 10   # peak at t=0, decays
  tc <- time_to_clearance(t, conc, frac = 0.01)
  cleared <- conc[t == tc] <= 0.01 * max(conc) + 1e-9
  cleared && tc > 0
})

# ---- string_seed(): deterministic and, for this project's real keys, distinct

check("string_seed() is deterministic (same key -> same seed)", {
  string_seed("ER07 intervention") == string_seed("ER07 intervention")
})

check("string_seed() gives distinct seeds for all 70 ERIE subject x visit keys", {
  keys <- as.vector(outer(sprintf("ER%02d", 1:35), c("baseline", "intervention"), paste))
  seeds <- vapply(keys, string_seed, integer(1))
  length(unique(seeds)) == length(keys)
})

# ---- fit_multistart(): recovers the true CURVE from noiseless synthetic data
#
# Checks curve fit, not exact parameter recovery: the
# one-compartment oral model has an exact "flip-flop" identifiability twin
# - (ka, kel, F) and (kel, ka, F*ka/kel) produce IDENTICAL concentration
# curves (verified directly: with wide, symmetric ka/kel bounds, this test
# originally asserted exact parameter recovery and failed deterministically
# under set.seed(1), converging instead to the exact flip-flop twin of the
# true parameters). This is real model structure, not an engine bug - it's
# exactly why 02_fit_erie_model.R's kel bound is much tighter than ka's,
# to break this degeneracy in production. A generic engine test shouldn't
# assume a specific parameterization is uniquely recoverable when the model
# itself doesn't guarantee that.

check("fit_multistart() recovers the true curve (not necessarily the true parameterization) from noiseless synthetic data", {
  true_par <- c(ka = 0.05, kel = 0.02, F = 0.2)
  t <- c(30, 60, 90, 120, 180, 240, 360)
  dose <- 100; Vd <- 10
  synthetic_conc <- bateman_conc(t, true_par[["ka"]], true_par[["kel"]], true_par[["F"]], dose, Vd)
  curve <- list(times = t, conc = synthetic_conc,
                simulate = function(theta) bateman_conc(t, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd))
  objective <- build_joint_objective(list(curve))
  lower <- c(ka = 1e-4, kel = 1e-4, F = 1e-4)
  upper <- c(ka = 1,    kel = 1,    F = 1)
  seeds <- c(list(true_par * 1.5), random_seeds(20, list(ka = c(1e-4, 1), kel = c(1e-4, 1), F = c(1e-4, 1)), log_scale = c("ka", "kel")))
  fit <- fit_multistart(objective, lower, upper, seeds = seeds)
  refit_conc <- bateman_conc(t, fit$par[["ka"]], fit$par[["kel"]], fit$par[["F"]], dose, Vd)
  !is.null(fit) && max(abs(refit_conc - synthetic_conc) / synthetic_conc) < 0.01
})

# ---- adaptive_retry(): only replaces a flagged row, and only if strictly better

check("adaptive_retry() replaces a flagged row only when the retry is strictly better", {
  bounds <- list(a = c(0, 10), b = c(0, 1))
  results <- tibble(a = c(1, 9.999, 3), b = c(0.1, 0.5, 0.3),
                     a_at_bound = c(FALSE, TRUE, FALSE),
                     objective_value = c(0.5, 0.9, 0.3))
  refit_better <- function(i, seeds, maxit) tibble(a = 5, b = 0.5, a_at_bound = FALSE, objective_value = 0.1)
  out_better <- adaptive_retry(results, bound_cols = c(a_at_bound = "a"), bounds = bounds,
                                par_cols = c("a", "b"), refit_fn = refit_better, mc.cores = 1, verbose = FALSE)
  refit_worse <- function(i, seeds, maxit) tibble(a = 5, b = 0.5, a_at_bound = FALSE, objective_value = 5.0)
  out_worse <- adaptive_retry(results, bound_cols = c(a_at_bound = "a"), bounds = bounds,
                               par_cols = c("a", "b"), refit_fn = refit_worse, mc.cores = 1, verbose = FALSE)
  out_better$a[2] == 5 && out_better$objective_value[2] == 0.1 &&    # improved row replaced
    out_worse$a[2] == 9.999 && out_worse$objective_value[2] == 0.9 && # worse retry discarded
    out_better$a[1] == 1 && out_better$a[3] == 3                      # unflagged rows untouched
})

check("adaptive_retry() also retries a row with no bound flag via extra_flag_cols/convergence_col (e.g. poor fit)", {
  bounds <- list(a = c(0, 10), b = c(0, 1))
  results <- tibble(a = c(1, 4, 3, 2), b = c(0.1, 0.5, 0.3, 0.2),
                     a_at_bound = c(FALSE, FALSE, FALSE, FALSE),  # no bound flags at all
                     r2_low = c(FALSE, TRUE, FALSE, FALSE),       # row 2: flagged via extra_flag_cols
                     converged = c(TRUE, TRUE, TRUE, FALSE),      # row 4: flagged via convergence_col
                     objective_value = c(0.5, 0.9, 0.3, 0.4))
  refit <- function(i, seeds, maxit) tibble(a = 9, b = 0.5, a_at_bound = FALSE, r2_low = FALSE, converged = TRUE, objective_value = 0.05)
  out <- adaptive_retry(results, bound_cols = c(a_at_bound = "a"), bounds = bounds, par_cols = c("a", "b"),
                         extra_flag_cols = "r2_low", convergence_col = "converged",
                         refit_fn = refit, mc.cores = 1, verbose = FALSE)
  out$a[2] == 9 && out$a[4] == 9 &&      # both the poor-fit row and the non-converged row got retried
    out$a[1] == 1 && out$a[3] == 3       # unflagged rows untouched
})

check("adaptive_retry()'s retried/retry_improved columns distinguish never-flagged, improved, and retried-but-not-improved", {
  bounds <- list(a = c(0, 10), b = c(0, 1))
  results <- tibble(a = c(1, 9.999, 3, 4), b = c(0.1, 0.5, 0.3, 0.5),
                     a_at_bound = c(FALSE, TRUE, FALSE, FALSE),
                     r2_low = c(FALSE, FALSE, FALSE, TRUE),
                     objective_value = c(0.5, 0.9, 0.3, 0.9))
  refit <- function(i, seeds, maxit) {
    if (i == 2) return(tibble(a = 5, b = 0.5, a_at_bound = FALSE, r2_low = FALSE, objective_value = 0.1))
    tibble(a = 99, b = 0.9, a_at_bound = FALSE, r2_low = TRUE, objective_value = 99)  # never beats the original
  }
  out <- adaptive_retry(results, bound_cols = c(a_at_bound = "a"), bounds = bounds, par_cols = c("a", "b"),
                         extra_flag_cols = "r2_low", refit_fn = refit, mc.cores = 1, verbose = FALSE)
  out$retried[1] == FALSE && is.na(out$retry_improved[1]) &&              # never flagged
    out$retried[2] == TRUE && isTRUE(out$retry_improved[2]) && out$a[2] == 5 &&   # flagged + improved
    out$retried[3] == FALSE && is.na(out$retry_improved[3]) &&            # never flagged
    out$retried[4] == TRUE && isFALSE(out$retry_improved[4]) && out$a[4] == 4     # retried, not improved, original kept
})

# ---- watson_ecf_volume(): hand-computed TBW / 3 ------------------------------

check("watson_ecf_volume() matches hand-computed values, and only males use age", {
  male   <- watson_ecf_volume(80, 180, 40, "Male")     # (2.447 - 3.8064 + 19.332 + 26.896) / 3
  female <- watson_ecf_volume(60, 165, 40, "female")   # (-2.097 + 17.6385 + 14.796) / 3
  abs(male - 14.9562) < 1e-4 && abs(female - 10.1125) < 1e-4 &&
    watson_ecf_volume(60, 165, NA, "female") == female &&
    is.na(watson_ecf_volume(80, 180, NA, "male"))
})

cat("\nAll", n_pass, "tests passed.\n")
