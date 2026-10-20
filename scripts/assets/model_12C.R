# =============================================================================
# model_12C.R - the 12C model: settings, stage 1 (each visit alone, structure chosen per visit) and stage 2
# (one kel per subject). Used by 02_fit_shared_kel.R and simulations/04_structure_selection.R, so the simulation
# tests the same code. Model: docs/pk-model.md.
# =============================================================================
suppressPackageStartupMessages(library(purrr))
source("scripts/assets/pk_equations.R")
source("scripts/assets/optimizer.R")
source("scripts/assets/curve_checks.R")

KEL <- log(2) / c(70, 7); KA <- c(1e-4, 1); FBIO <- c(1e-4, 1); F_DELAYED <- c(0.001, 0.999); T_LAG <- c(5, 150)
EARLY_LAG_OK_MGL <- 5             # t=30 at or below this: t_lag may go down to 5 min, otherwise 30
# How stage 1 chooses the structure of a visit (simulations/04_structure_selection.R compares the rules):
#   "gated"       two_wave is fitted only when the dip or plateau check fires or single_wave R2 is below 0.50,
#                 the lower AIC decides, and a large dip selects it without the AIC comparison
#   "tail_limit"  as gated, but the dip check only counts points up to 180 min (a wave starting by t_lag's
#                 150-min bound is visible by then) that rise at least 5% of the peak, and single_wave R2 below
#                 0.95 also triggers two_wave
#   "always_try"  two_wave is fitted whenever t=30 exists and the lower AIC decides
STRUCTURE_RULE <- "tail_limit"
R2_ALWAYS_TRY_TWO_WAVE <- if (STRUCTURE_RULE == "tail_limit") 0.95 else 0.50   # single_wave R2 below this: two_wave is tried
DIP_MAX_TIME <- if (STRUCTURE_RULE == "tail_limit") 180 else Inf                # dip check: latest point that counts (min)
DIP_MIN_PEAK_FRAC <- if (STRUCTURE_RULE == "tail_limit") 0.05 else -Inf         # dip check: smallest rise, as a share of the peak
EXCESS_DEFINITE_TWO_WAVE <- 0.90  # gated, tail_limit: peak_dip_rise_info()$excess at or above this: two_wave without the AIC comparison
N_RANDOM <- 30; MAXIT <- 80
R2_RETRY <- 0.70; MAXIT_RETRY <- 100   # retry a fit with kel at a bound, a visit below this R2, or no convergence
RATE_SEEDS <- c(list(c(ka = 0.03, kel = 0.02, F = 0.3), c(ka = 0.08, kel = 0.05, F = 0.15),
                     c(ka = 0.01, kel = 0.01, F = 0.5), c(ka = 0.05, kel = 0.08, F = 0.2)),
                lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) c(ka = v, kel = v, F = 0.03)))
WAVE_SEEDS <- list(c(0.01, 60), c(0.01, 90), c(0.3, 60), c(0.4, 45), c(0.5, 60), c(0.35, 75), c(0.45, 90))   # (f_delayed, t_lag)

# One visit: data, dose, Vd and the raw-data evidence for a second wave
make_visit <- function(times, conc, dose, Vd) {
  dip <- peak_dip_rise_info(times, conc, max_time = DIP_MAX_TIME, min_peak_frac = DIP_MIN_PEAK_FRAC)
  t30 <- conc[times == 30]
  list(times = times, conc = conc, dose = dose, Vd = Vd, t30_present = length(t30) == 1,
       evidence = dip$detected || has_near_peak_neighbor(times, conc),
       definite = isTRUE(dip$excess >= EXCESS_DEFINITE_TWO_WAVE), trigger_time = dip$trigger_time,
       t_lag_lower = if (length(t30) == 1 && t30 > EARLY_LAG_OK_MGL) 30 else T_LAG[1])
}

# Parameters are named by visit suffix (ka_b, F_i, t_lag_b, ...); kel has no suffix and is shared by every visit fitted together
pred_visit <- function(v, s, two_wave, th, t) {
  simA <- function(tt, dose) bateman_conc(tt, th[[paste0("ka_", s)]], th[["kel"]], th[[paste0("F_", s)]], dose, v$Vd)
  if (two_wave) simulate_lagged_dose(simA, t, v$dose, th[[paste0("f_delayed_", s)]], th[[paste0("t_lag_", s)]]) else simA(t, v$dose)
}

param_bounds <- function(vs, waves) {
  b <- list(kel = KEL)
  for (s in names(vs)) {
    b[[paste0("ka_", s)]] <- KA; b[[paste0("F_", s)]] <- FBIO
    if (waves[[s]]) { b[[paste0("f_delayed_", s)]] <- F_DELAYED; b[[paste0("t_lag_", s)]] <- c(vs[[s]]$t_lag_lower, T_LAG[2]) }
  }
  b
}

# Starting points: the fixed and ka = kel rate seeds, each with one wave seed (cycled), every `nested` fit with
# each wave seed, and random seeds. The wave seeds include the trigger time of the dip and 20 and 40 min before it.
make_seeds <- function(vs, waves, nested = list()) {
  wave_opts <- setNames(lapply(names(vs), function(s) {
    tt <- vs[[s]]$trigger_time
    c(WAVE_SEEDS, if (!is.na(tt)) list(c(0.4, tt), c(0.3, tt - 20), c(0.35, tt - 40)))
  }), names(vs))
  add_waves <- function(p, j) {
    for (s in names(vs)) if (waves[[s]] && !(paste0("t_lag_", s) %in% names(p))) {
      w <- wave_opts[[s]]; p[paste0(c("f_delayed_", "t_lag_"), s)] <- w[[(j - 1) %% length(w) + 1]]
    }
    p
  }
  rates <- lapply(RATE_SEEDS, function(r) {
    p <- c(kel = r[["kel"]])
    for (s in names(vs)) p[paste0(c("ka_", "F_"), s)] <- r[c("ka", "F")]
    p
  })
  n_wave <- max(lengths(wave_opts))
  bounds <- param_bounds(vs, waves)
  seeds <- c(imap(rates, add_waves),
             unlist(lapply(nested, function(p) lapply(seq_len(n_wave), function(j) add_waves(p, j))), recursive = FALSE),
             random_seeds(N_RANDOM, bounds, log_scale = c("kel", paste0("ka_", names(vs)))))
  unique(lapply(seeds, function(p) p[names(bounds)]))
}

# Fit the visits in `vs` (a list named by suffix) together; `waves` says per visit whether it is two_wave
fit_visits <- function(vs, waves, seeds) {
  curves <- imap(vs, function(v, s) list(times = v$times, conc = v$conc,
    simulate = function(th) pred_visit(v, s, waves[[s]], th, v$times),
    fine_simulate = function(th) list(time = FINE_T, conc = pred_visit(v, s, waves[[s]], th, FINE_T))))
  joint <- build_objective(curves, min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                 tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
  # The Tmax floor of the combined curve only reaches the taller wave. Both waves have the shape of the
  # plain Bateman curve, so the floor on that shape applies it to each wave.
  two <- names(waves)[unlist(waves)]
  objective <- function(th) {
    total <- joint(th)
    for (s in two) {
      shape <- bateman_conc(FINE_T, th[[paste0("ka_", s)]], th[["kel"]], 1, 1, 1)
      total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - FINE_T[which.max(shape)])^2
    }
    total
  }
  bounds <- param_bounds(vs, waves); lower <- map_dbl(bounds, 1); upper <- map_dbl(bounds, 2)
  fit <- fit_multistart(objective, lower, upper, seeds = seeds, control = list(maxit = MAXIT))
  if (is.null(fit)) return(NULL)
  rss_of <- function(par) map_dbl(curves, ~ sum((.x$conc - .x$simulate(par))^2))
  r2_of <- function(par) map_dbl(curves, ~ r_squared(.x$conc, .x$simulate(par)))
  # A suspect fit (kel at a bound, a visit with a poor R2, or stopped at maxit) may sit in a worse local optimum.
  # It is searched again from the estimate, a grid over kel's range when kel is at a bound and random seeds,
  # and replaced only by a lower objective.
  kel_at_bound <- fit$par[["kel"]] < KEL[1] + 1e-4 || fit$par[["kel"]] > KEL[2] - 1e-4
  retried <- kel_at_bound || any(r2_of(fit$par) < R2_RETRY) || fit$convergence != 0
  retry_improved <- NA
  if (retried) {
    retry <- fit_multistart(objective, lower, upper, control = list(maxit = MAXIT_RETRY),
                            seeds = retry_seeds_near_bounds(fit$par, bounds, if (kel_at_bound) "kel" else character(0),
                                                            log_scale = c("kel", paste0("ka_", names(vs)))))
    retry_improved <- !is.null(retry) && retry$value < fit$value
    if (retry_improved) fit <- retry
  }
  list(par = fit$par, waves = waves, converged = fit$convergence == 0, rss = rss_of(fit$par), r2 = r2_of(fit$par),
       n = map_int(curves, ~ length(.x$conc)), retried = retried, retry_improved = retry_improved)
}

# AIC with a residual variance per visit, so each visit is judged on its own residuals
aic_of <- function(f) aic_rss(f$rss, f$n, length(f$par))
VISITS <- c(b = "baseline", i = "intervention")

# Stage 1: one visit alone, with its own kel. single_wave always, and two_wave as STRUCTURE_RULE says; the lower
# AIC on this curve's residuals decides (with the gated rules, definite raw-data evidence keeps two_wave).
fit_alone <- function(sid, vs, s) {
  v <- vs[s]
  run <- function(two, nested = list()) {
    set.seed(string_seed(paste(sid, s, two, "alone")))
    w <- setNames(list(two), s)
    fit_visits(v, w, make_seeds(v, w, nested))
  }
  sw <- run(FALSE)
  if (is.null(sw)) return(NULL)
  gated <- STRUCTURE_RULE != "always_try"
  tried <- v[[s]]$t30_present && (!gated || v[[s]]$evidence || sw$r2[[s]] < R2_ALWAYS_TRY_TWO_WAVE)
  tw <- if (tried) run(TRUE, list(sw$par)) else NULL
  aic <- c(single_wave = aic_of(sw), two_wave = if (is.null(tw)) NA_real_ else aic_of(tw))
  best <- if (!is.null(tw) && ((gated && v[[s]]$definite) || aic[["two_wave"]] < aic[["single_wave"]])) tw else sw
  c(best, list(tried = tried, aic = aic))
}

# Stage 2: both visits with one kel, the structure of stage 1 fixed and f_delayed and t_lag re-estimated.
# Seeds from stage 1: kel averaged, and each visit's slower rate taken as kel with the faster rates as ka.
fit_shared <- function(key, vs, alone) {
  w <- map(alone, ~ .x$waves[[1]])
  rate <- function(s, f) f(alone[[s]]$par[["kel"]], alone[[s]]$par[[paste0("ka_", s)]])
  rest <- c(alone$b$par[names(alone$b$par) != "kel"], alone$i$par[names(alone$i$par) != "kel"])
  as_ka <- function(p) { p[c("ka_b", "ka_i")] <- c(rate("b", max), rate("i", max)); p }
  informed <- list(c(kel = sqrt(alone$b$par[["kel"]] * alone$i$par[["kel"]]), rest),
                   as_ka(c(kel = rate("b", min), rest)), as_ka(c(kel = rate("i", min), rest)))
  set.seed(string_seed(paste(key, "shared")))
  fit_visits(vs, w, make_seeds(vs, w, informed))
}
