# =============================================================================
# tests.R
#
# Lightweight regression tests for the PK engine (pk_equations.R, optimizer.R,
# curve_checks.R, stats.R), in base R without a testing framework.
# Run directly: Rscript scripts/assets/tests.R (or `pixi run test-engine`).
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

source("scripts/assets/pk_equations.R")
source("scripts/assets/optimizer.R")
source("scripts/assets/curve_checks.R")

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

check("bateman_conc() matches the textbook Bateman form away from ka==kel, for ka < kel and ka > kel", {
  t <- c(0, 15, 30, 60, 120, 240, 360)
  textbook <- function(ka, kel) 0.2 * ka * 100 / (ka - kel) * (exp(-kel * t) - exp(-ka * t)) / 10
  all(abs(bateman_conc(t, 0.05, 0.02, 0.2, 100, 10) - textbook(0.05, 0.02)) < 1e-12,
      abs(bateman_conc(t, 0.02, 0.05, 0.2, 100, 10) - textbook(0.02, 0.05)) < 1e-12)
})

check("bateman_conc() is smooth for ka-kel from 1e-12 to 1e-6 (no cancellation noise)", {
  t <- c(30, 60, 120, 240, 360)
  ref <- bateman_conc(t, 0.03, 0.03, 0.2, 100, 10)
  worst <- max(vapply(10^(-12:-6), function(d)
    max(abs(bateman_conc(t, 0.03 + d, 0.03, 0.2, 100, 10) - ref) / ref), numeric(1)))
  worst < 1e-3 && all(is.finite(bateman_conc(t, 0.03 + 1e-14, 0.03, 0.2, 100, 10)))
})

# ---- mixed_visit_tests(): recovers known within-group and between-group effects ----

check("mixed_visit_tests() recovers a known visit effect per group and the interaction, and returns NA when nothing varies within subject", {
  source("scripts/assets/stats.R")
  set.seed(1)
  n <- 30
  d <- expand.grid(subject_id = paste0("S", 1:(2 * n)), visit = c("baseline", "intervention"), stringsAsFactors = FALSE)
  d$diet <- ifelse(as.integer(sub("S", "", d$subject_id)) <= n, "low_fructose", "high_fructose")
  subj <- rnorm(2 * n, 0, 0.3)[as.integer(sub("S", "", d$subject_id))]
  true_effect <- ifelse(d$diet == "high_fructose", 0.4, 0) * (d$visit == "intervention")
  d$value <- exp(2 + subj + true_effect + rnorm(nrow(d), 0, 0.1))
  r <- mixed_visit_tests(d)
  ok_effects <- abs(r$lmm_estimate[r$diet == "high_fructose"] - 0.4) < 0.1 && abs(r$lmm_estimate[r$diet == "low_fructose"]) < 0.1
  ok_interaction <- r$p_interaction[1] < 0.001 && abs(r$interaction_estimate[1] - 0.4) < 0.1
  ok_ci <- r$lmm_lo[r$diet == "high_fructose"] < 0.4 && r$lmm_hi[r$diet == "high_fructose"] > 0.4
  const <- d; const$value <- exp(subj[match(const$subject_id, d$subject_id)])   # identical at both visits
  ok_na <- all(is.na(mixed_visit_tests(const)$p_lmm))
  ok_effects && ok_interaction && ok_ci && ok_na
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
# why the ka and kel labels are only resolved with one kel per subject (see
# "Why one kel per subject" in docs/pk-model.md). An engine test shouldn't
# assume a specific parameterization is uniquely recoverable when the model
# itself doesn't guarantee that.

check("fit_multistart() recovers the true curve (not necessarily the true parameterization) from noiseless synthetic data", {
  true_par <- c(ka = 0.05, kel = 0.02, F = 0.2)
  t <- c(30, 60, 90, 120, 180, 240, 360)
  dose <- 100; Vd <- 10
  synthetic_conc <- bateman_conc(t, true_par[["ka"]], true_par[["kel"]], true_par[["F"]], dose, Vd)
  curve <- list(times = t, conc = synthetic_conc,
                simulate = function(theta) bateman_conc(t, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd))
  objective <- build_objective(list(curve))
  lower <- c(ka = 1e-4, kel = 1e-4, F = 1e-4)
  upper <- c(ka = 1,    kel = 1,    F = 1)
  seeds <- c(list(true_par * 1.5), random_seeds(20, list(ka = c(1e-4, 1), kel = c(1e-4, 1), F = c(1e-4, 1)), log_scale = c("ka", "kel")))
  fit <- fit_multistart(objective, lower, upper, seeds = seeds)
  refit_conc <- bateman_conc(t, fit$par[["ka"]], fit$par[["kel"]], fit$par[["F"]], dose, Vd)
  !is.null(fit) && max(abs(refit_conc - synthetic_conc) / synthetic_conc) < 0.01
})

# ---- retry_seeds_near_bounds(): prior first, grid only over the flagged parameter ----

check("retry_seeds_near_bounds() starts from the prior and grids only the flagged parameter within its bounds", {
  bounds <- list(a = c(0, 10), b = c(0, 1))
  prior <- c(a = 9.999, b = 0.3)
  seeds <- retry_seeds_near_bounds(prior, bounds, "a", n_random = 5)
  grid <- seeds[2:9]
  identical(seeds[[1]], prior) && length(seeds) == 1 + 8 + 5 &&
    all(vapply(grid, function(p) p[["b"]] == 0.3 && p[["a"]] >= 0 && p[["a"]] <= 10, logical(1)))
})

# ---- has_near_peak_neighbor(): peak-at-t=30 is excluded, not just neighbor-at-t=30 ----

check("has_near_peak_neighbor() excludes a peak AT t=30, mirroring the neighbor-at-30 exclusion", {
  # Peak at t=30, near-tall neighbor at t=60: an early single peak (at or
  # before t=30) declining slowly can produce exactly this shape, so it must
  # NOT count as second-wave evidence (this is the case the fix targets).
  peak_at_30 <- !has_near_peak_neighbor(c(0, 30, 60, 90), c(1, 100, 90, 40))
  # Same shape shifted one sample later (peak at t=60, tall neighbor at
  # t=90) is NOT excluded - a real signature the detector should still catch.
  peak_at_60 <- has_near_peak_neighbor(c(0, 30, 60, 90), c(1, 50, 100, 90))
  peak_at_30 && peak_at_60
})

cat("\nAll", n_pass, "tests passed.\n")
