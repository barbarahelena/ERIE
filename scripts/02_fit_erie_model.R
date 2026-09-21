# Fitting the joint 12C/13C6 fructose PK model
# Barbara Verhaar
#
# Fits two candidate models per subject x visit and picks one by AIC on 12C's
# own residuals:
#   - single_wave: the plain joint delayed-release model.
#   - two_wave: a lagged-second-dose extension (simulate_lagged_dose, in
#     pk_curves.R) that can represent the post-peak dip, which single_wave
#     structurally cannot.
# two_wave is only attempted when 12C's t=30 sample exists and its raw data
# show evidence of a second wave (or single_wave fits very poorly). 13C6 then
# makes its own independent choice of mechanism. Both models use an unweighted
# objective (proportional_weighting = FALSE).
#
# The reasoning, validation cases and known limitations behind all of this are
# in docs/pk-model.md, mainly "The post-peak dip and the lagged-dose model" and
# "Fitting procedure".

# Libraries
suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(grid)
  library(ggthemes)
  library(stringr)
})

# Source function scripts
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_fit.R")
source("scripts/assets/pk_diagnostics.R")

# House plotting theme
theme_Publication <- function(base_size=14, base_family="sans") {
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0),
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                strip.text.y = element_text(size = rel(1.05)),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))

}

# Plot labels
ISOTOPE_LABELS <- c("12C" = "Fructose 12C", "13C6" = "Fructose 13C6")
ISOTOPE_COLORS <- c("12C" = "steelblue", "13C6" = "firebrick")
VISIT_LABELS   <- c(baseline = "Baseline (FCT1)", intervention = "Intervention (FCT2)")

# Config
MIN_TMAX    <- 30     # min - soft floor on predicted Tmax for both curves - see file header
CMAX_TOL    <- 0.10   # +/-10% soft band around each curve's own observed Cmax
CMAX_LAMBDA <- 20     # penalty weight for the Cmax band
TMAX_LAMBDA <- 50     # penalty weight for the Tmax floor
# t_lag's lower bound is a step function of the subject's own 12C t=30 sample:
# 5 at or below this threshold, else 30 (see "The post-peak dip and the
# lagged-dose model" in docs/pk-model.md).
EARLY_LAG_OK_MGL  <- 5    # t=30 at or below this: early t_lag (as low as 5) is fine
# kel: Hannou et al. 2018 (t1/2 ~ 7-140 min). ka: bounded only below in
# spirit (absorption-rate differences are part of the research question).
BOUNDS_INDEP <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F = c(1e-4, 1))
BOUNDS_JOINT <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1),
                     F_12C = c(1e-4, 1), F_13C6 = c(1e-4, 1),
                     k_release = c(0.001, 1))   # capsule dissolution t1/2 ~ 0.7-700 min
BOUNDS_LAGGED <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1),
                      F_12C = c(1e-4, 1), F_13C6 = c(1e-4, 1),
                      k_release = c(0.001, 1),
                      # t_lag [5, 150]: why these ends (sampling gaps, per-subject lower bound)
                      # is in "The post-peak dip and the lagged-dose model" in docs/pk-model.md.
                      f_delayed = c(0.001, 0.999), t_lag = c(5, 150),
                      # 13C6's own delayed fraction, separate from 12C's f_delayed (see
                      # "13C6's own wave choice" in docs/pk-model.md).
                      f_delayed_13C6 = c(0.001, 0.999))
# BOUNDS_LAGGED above is the combined 8-parameter view adaptive_retry() works
# over; the two-stage fit uses these two subsets. t_lag's lower bound here (5)
# is the permissive floor - fit_subject_visit_two_wave() raises it to 30 for
# subjects whose t=30 is above EARLY_LAG_OK_MGL.
BOUNDS_LAGGED_12C <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F_12C = c(1e-4, 1),
                          f_delayed = c(0.001, 0.999), t_lag = c(5, 150))
BOUNDS_13C6_GIVEN_LAG <- list(F_13C6 = c(1e-4, 1), k_release = c(0.001, 1),
                              f_delayed_13C6 = c(0.001, 0.999))
N_RANDOM_INDEP <- 40
N_RANDOM_JOINT <- 16
N_RANDOM_LAGGED <- 60   # pilot-validated budget (ER01/ER03/ER09) - see results/model-post-peak-dip-pilot/
N_RANDOM_13C6_GIVEN_LAG <- 60  # stage 2 is 3-dim and cheap (no 12C simulation), so a generous budget
FINE_T <- seq(0, 400, length.out = 50)    # grid for Tmax/Cmax penalty checks
CLEARANCE_FRAC <- 0.01                    # "fully cleared", for plot x-axis limits only
R2_RELIABLE_MIN <- 0.70                   # below this, that curve's fit is flagged unreliable
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = max(1, parallel::detectCores() - 2))) # To parallelize loop

# NOTE: no top-level set.seed() here - it wouldn't make the multi-start seeds
# reproducible under parallel::mclapply(). Each subject x visit fit seeds itself
# from its own ID in fit_one() (see "Fitting procedure" in docs/pk-model.md).

# Load cleaned data
concentrations <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE)
covariates     <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants      <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)

# Prepare data
MW_12C  <- constants$value[constants$constant == "MW_12C"]
MW_13C6 <- constants$value[constants$constant == "MW_13C6"]
dose_13C6_mg <- constants$value[constants$constant == "tracer_13C6_dose_mg"]

data <- concentrations %>%
  left_join(covariates, by = c("subject_id", "visit")) %>%
  mutate(conc_mgL = conc_umol_L * if_else(isotope == "12C", MW_12C, MW_13C6) / 1000) %>%
  filter(!is.na(bw_kg), !is.na(conc_mgL))

# Baseline correction for 12C - deducting the baseline fasted fructose level
baseline_12C <- data %>% filter(isotope == "12C", time_min == 0) %>%
  select(subject_id, visit, baseline_mgL = conc_mgL)

data_corrected <- data %>%
  left_join(baseline_12C, by = c("subject_id", "visit")) %>%
  mutate(conc_mgL = if_else(isotope == "12C", conc_mgL - baseline_mgL, conc_mgL)) %>%
  select(-baseline_mgL)

# t=0 carries no fitting information (both models predict exactly 0 there), so
# it's excluded from the fit; data_corrected keeps it so the plots still show
# the observed point.
data_fit <- data_corrected %>% filter(time_min > 0)

subject_visits <- data_fit %>% distinct(subject_id, visit) %>% arrange(subject_id, visit)

# Fast-testing mode: restrict to a handful of subjects for quick local
# iteration without waiting on the full cohort run. Set via
# ERIE_TEST_SUBJECTS="ER01,ER02,ER03" Rscript scripts/02_fit_erie_model.R -
# unset (the default) runs every subject.
test_subjects <- Sys.getenv("ERIE_TEST_SUBJECTS", unset = "")
if (nzchar(test_subjects)) {
  keep <- trimws(strsplit(test_subjects, ",")[[1]])
  subject_visits <- subject_visits %>% filter(subject_id %in% keep)
  cat("ERIE_TEST_SUBJECTS set - restricting to:", paste(keep, collapse = ", "), "\n")
}

# ---------------------------------------------------------------------------
# Independent pre-fit (unweighted) - starting values for both joint models
# ---------------------------------------------------------------------------
fit_curve_independent <- function(obs_time, obs_conc, dose, Vd) {
  curve <- list(
    times = obs_time, conc = obs_conc,
    simulate      = function(theta) bateman_conc(obs_time, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd),
    fine_simulate = function(theta) list(time = FINE_T, conc = bateman_conc(FINE_T, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd))
  )
  objective <- build_joint_objective(list(curve), min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL,
                                      cmax_lambda = CMAX_LAMBDA, tmax_lambda = TMAX_LAMBDA,
                                      proportional_weighting = FALSE)

  # ka=kel diagonal seeds: the objective is awkward near ka ~ kel and generic
  # seeds miss the optimum there (see "Known limitations and open questions" in docs/pk-model.md).
  diagonal_seeds <- lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07),
                            function(v) c(ka = v, kel = v, F = 0.03))
  seeds <- c(
    list(c(ka = 0.03, kel = 0.02, F = 0.3), c(ka = 0.08, kel = 0.05, F = 0.15),
         c(ka = 0.01, kel = 0.01, F = 0.5), c(ka = 0.05, kel = 0.08, F = 0.2)),
    diagonal_seeds,
    random_seeds(N_RANDOM_INDEP, BOUNDS_INDEP, log_scale = c("ka", "kel"))
  )
  fit_multistart(objective,
                  lower = c(ka = BOUNDS_INDEP$ka[1], kel = BOUNDS_INDEP$kel[1], F = BOUNDS_INDEP$F[1]),
                  upper = c(ka = BOUNDS_INDEP$ka[2], kel = BOUNDS_INDEP$kel[2], F = BOUNDS_INDEP$F[2]),
                  seeds = seeds)
}

# ---------------------------------------------------------------------------
# Baseline joint fit: shared ka/kel, separate F per curve, 13C6 gets a
# delayed-release step. extra_seeds/maxit let the retry pass reuse this.
# ---------------------------------------------------------------------------
fit_subject_visit_single_wave <- function(sid, vis, extra_seeds = list(), maxit = 40) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA, t_lag1_13C6 = NA,
                   f_delayed2_13C6 = NA, t_lag2_13C6 = NA,
                   r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                   converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA,
                   rss_12C = NA, n_12C = NA, t30_present = NA)

  obs12 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  obs13 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "13C6")
  if (nrow(obs12) < 3 || nrow(obs13) < 3) return(empty)

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- watson_ecf_volume(cov$bw_kg, cov$height_cm, cov$age_years, cov$sex)
  dose_12C <- cov$dose_12C_mg

  pre12 <- fit_curve_independent(obs12$time_min, obs12$conc_mgL, dose_12C, Vd)
  pre13 <- fit_curve_independent(obs13$time_min, obs13$conc_mgL, dose_13C6_mg, Vd)
  if (is.null(pre12) || is.null(pre13)) return(empty)

  curve_12C <- list(
    times = obs12$time_min, conc = obs12$conc_mgL,
    simulate      = function(theta) bateman_conc(obs12$time_min, theta[["ka"]], theta[["kel"]], theta[["F_12C"]], dose_12C, Vd),
    fine_simulate = function(theta) list(time = FINE_T, conc = bateman_conc(FINE_T, theta[["ka"]], theta[["kel"]], theta[["F_12C"]], dose_12C, Vd))
  )
  curve_13C6 <- list(
    times = obs13$time_min, conc = obs13$conc_mgL,
    simulate      = function(theta) simulate_delayed_release(obs13$time_min, theta[["k_release"]], theta[["ka"]], theta[["kel"]], theta[["F_13C6"]], dose_13C6_mg, Vd)$conc,
    fine_simulate = function(theta) simulate_delayed_release(FINE_T, theta[["k_release"]], theta[["ka"]], theta[["kel"]], theta[["F_13C6"]], dose_13C6_mg, Vd)
  )

  objective <- build_joint_objective(list(curve_12C, curve_13C6),
                                      min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL,
                                      cmax_lambda = CMAX_LAMBDA, tmax_lambda = TMAX_LAMBDA,
                                      proportional_weighting = FALSE)

  mean_ka  <- mean(c(pre12$par[["ka"]],  pre13$par[["ka"]]))
  mean_kel <- mean(c(pre12$par[["kel"]], pre13$par[["kel"]]))
  informed_seeds <- list(
    c(ka = pre12$par[["ka"]], kel = mean_kel, F_12C = pre12$par[["F"]], F_13C6 = pre13$par[["F"]], k_release = 0.05),
    c(ka = pre13$par[["ka"]], kel = mean_kel, F_12C = pre12$par[["F"]], F_13C6 = pre13$par[["F"]], k_release = 0.02),
    c(ka = mean_ka, kel = pre12$par[["kel"]], F_12C = pre12$par[["F"]], F_13C6 = pre13$par[["F"]], k_release = 0.10),
    c(ka = mean_ka, kel = pre13$par[["kel"]], F_12C = pre12$par[["F"]], F_13C6 = pre13$par[["F"]], k_release = 0.03)
  )
  grid_seeds_joint <- grid_seeds(
    grid  = list(ka = c(0.02, 0.08), kel = c(0.01, 0.03, 0.06, 0.09)),
    fixed = list(F_12C = 0.1, F_13C6 = 0.1, k_release = 0.05)
  )
  # ka=kel diagonal seeds - see fit_curve_independent().
  diagonal_seeds_joint <- lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07),
                                  function(v) c(ka = v, kel = v, F_12C = 0.03, F_13C6 = 0.03, k_release = 0.05))
  jitter_seeds <- random_seeds(N_RANDOM_JOINT, BOUNDS_JOINT, log_scale = c("ka", "kel", "k_release"))

  lower <- vapply(BOUNDS_JOINT, `[`, numeric(1), 1)
  upper <- vapply(BOUNDS_JOINT, `[`, numeric(1), 2)
  fit <- fit_multistart(objective, lower, upper,
                         seeds = c(informed_seeds, grid_seeds_joint, diagonal_seeds_joint, jitter_seeds, extra_seeds),
                         control = list(maxit = maxit))
  if (is.null(fit)) return(empty)

  par <- fit$par
  pred12 <- curve_12C$simulate(par)
  pred13 <- curve_13C6$simulate(par)

  r2_12C  <- r_squared(curve_12C$conc, pred12)
  r2_13C6 <- r_squared(curve_13C6$conc, pred13)

  # 13C6's own onset lag: if its raw data show absorption hadn't started yet
  # (has_onset_lag_evidence()), refit F_13C6 plus t_lag1_13C6 with ka/kel fixed
  # at the joint fit's values, as an instant bolus (no k_release). Kept only if
  # it beats the no-lag fit on AIC (K=2 either way; computed inline since aic()
  # below isn't defined yet). Never touches 12C's parameters. See "13C6's
  # onset-lag and second-wave candidates" in docs/pk-model.md.
  n13 <- length(curve_13C6$conc)
  rss_no_lag <- sum((curve_13C6$conc - pred13)^2)
  best_aic <- n13 * log(rss_no_lag / n13) + 2 * 2   # no-lag baseline: F_13C6/k_release, K=2
  t_lag1_13C6 <- NA_real_
  t_lag2_13C6 <- NA_real_
  f_delayed2_13C6 <- NA_real_

  if (has_onset_lag_evidence(curve_13C6$times, curve_13C6$conc)) {
    obj_13C6_lag <- function(theta) {
      simB <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], theta[["F_13C6"]], dose, Vd)
      pred <- tryCatch(simulate_two_lag_dose(simB, curve_13C6$times, dose_13C6_mg, f_delayed = 0, t_lag1 = theta[["t_lag1_13C6"]], gap = 0), error = function(e) NULL)
      if (is.null(pred) || any(!is.finite(pred))) return(1e10)
      sum((pred - curve_13C6$conc)^2)
    }
    seeds_lag <- lapply(c(10, 20, 30, 40, 50), function(tl) c(F_13C6 = par[["F_13C6"]], t_lag1_13C6 = tl))
    lower_lag <- c(F_13C6 = BOUNDS_JOINT$F_13C6[1], t_lag1_13C6 = 0)
    upper_lag <- c(F_13C6 = BOUNDS_JOINT$F_13C6[2], t_lag1_13C6 = 90)
    fit_lag <- fit_multistart(obj_13C6_lag, lower_lag, upper_lag, seeds = seeds_lag, control = list(maxit = 60))
    if (!is.null(fit_lag)) {
      aic_lag <- n13 * log(fit_lag$value / n13) + 2 * 2   # F_13C6/t_lag1_13C6, K=2
      if (aic_lag < best_aic) {
        best_aic <- aic_lag
        par[["F_13C6"]] <- fit_lag$par[["F_13C6"]]
        par[["k_release"]] <- NA_real_   # instant bolus - no meaningful capsule dissolution rate for this candidate
        t_lag1_13C6 <- fit_lag$par[["t_lag1_13C6"]]
        simB_final <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], par[["F_13C6"]], dose, Vd)
        pred13 <- simulate_two_lag_dose(simB_final, curve_13C6$times, dose_13C6_mg, f_delayed = 0, t_lag1 = t_lag1_13C6, gap = 0)
        r2_13C6 <- r_squared(curve_13C6$conc, pred13)
      }
    }
  }

  # 13C6's own second, later hump: a separate wave after the first is underway,
  # decided from 13C6's own raw data (the detectors that gate 12C's two_wave),
  # not inherited from 12C. Bounded more tightly than 12C's (f_delayed2_13C6 in
  # [0.05, 0.95], t_lag2_13C6 in [30, 150]) because 13C6 is much noisier. Instant
  # bolus, compared by AIC (K=3) against the current winner, and mutually
  # exclusive with the onset lag. See "13C6's onset-lag and second-wave
  # candidates" in docs/pk-model.md.
  if (has_peak_dip_rise(curve_13C6$times, curve_13C6$conc) ||
      has_near_peak_neighbor(curve_13C6$times, curve_13C6$conc)) {
    obj_13C6_tw <- function(theta) {
      simB <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], theta[["F_13C6"]], dose, Vd)
      pred <- tryCatch(simulate_lagged_dose(simB, curve_13C6$times, dose_13C6_mg, theta[["f_delayed2_13C6"]], theta[["t_lag2_13C6"]]), error = function(e) NULL)
      if (is.null(pred) || any(!is.finite(pred))) return(1e10)
      sum((pred - curve_13C6$conc)^2)
    }
    seeds_tw <- lapply(c(40, 60, 90, 120), function(tl)
      c(F_13C6 = par[["F_13C6"]], f_delayed2_13C6 = 0.5, t_lag2_13C6 = tl))
    lower_tw <- c(F_13C6 = BOUNDS_JOINT$F_13C6[1], f_delayed2_13C6 = 0.05, t_lag2_13C6 = 30)
    upper_tw <- c(F_13C6 = BOUNDS_JOINT$F_13C6[2], f_delayed2_13C6 = 0.95, t_lag2_13C6 = 150)
    fit_tw <- fit_multistart(obj_13C6_tw, lower_tw, upper_tw, seeds = seeds_tw, control = list(maxit = 60))
    # A result pinned at f_delayed2_13C6's lower bound means "no real second
    # wave" - the degenerate regime the bound excludes - so it is rejected
    # outright rather than compared on AIC.
    at_lower_bound <- !is.null(fit_tw) && fit_tw$par[["f_delayed2_13C6"]] <= lower_tw[["f_delayed2_13C6"]] + 1e-3
    if (!is.null(fit_tw) && !at_lower_bound) {
      aic_tw <- n13 * log(fit_tw$value / n13) + 2 * 3
      if (aic_tw < best_aic) {
        best_aic <- aic_tw
        par[["F_13C6"]] <- fit_tw$par[["F_13C6"]]
        par[["k_release"]] <- NA_real_
        t_lag1_13C6 <- NA_real_   # mutually exclusive with the onset-lag candidate above
        f_delayed2_13C6 <- fit_tw$par[["f_delayed2_13C6"]]
        t_lag2_13C6 <- fit_tw$par[["t_lag2_13C6"]]
        simB_final <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], par[["F_13C6"]], dose, Vd)
        pred13 <- simulate_lagged_dose(simB_final, curve_13C6$times, dose_13C6_mg, f_delayed2_13C6, t_lag2_13C6)
        r2_13C6 <- r_squared(curve_13C6$conc, pred13)
      }
    }
  }

  tibble(
    ka = par[["ka"]], kel = par[["kel"]], F_12C = par[["F_12C"]], F_13C6 = par[["F_13C6"]],
    k_release = par[["k_release"]], t_lag1_13C6 = t_lag1_13C6,
    f_delayed2_13C6 = f_delayed2_13C6, t_lag2_13C6 = t_lag2_13C6,
    r2_12C  = r2_12C,
    r2_13C6 = r2_13C6,
    kel_at_bound = kel > (BOUNDS_JOINT$kel[2] - 1e-4),
    k_release_at_bound = isTRUE(par[["k_release"]] > (BOUNDS_JOINT$k_release[2] - 1e-4)),
    converged = fit$convergence == 0,   # see "Fit quality and what to trust" in docs/pk-model.md
    r2_12C_low  = r2_12C  < R2_RELIABLE_MIN,
    r2_13C6_low = r2_13C6 < R2_RELIABLE_MIN,
    objective_value = fit$value,   # for comparing against a retry; not meaningful across subjects
    rss_12C = sum((curve_12C$conc - pred12)^2),   # raw (unweighted, native mg/L^2) - for the single_wave-vs-two_wave AIC comparison below
    n_12C = length(curve_12C$conc),
    t30_present = any(obs12$time_min == 30)   # ka is poorly anchored without this sample - see fit_subject_visit_two_wave()
  )
}

# ---------------------------------------------------------------------------
# Two-wave (lagged-second-dose) fit, in two stages rather than one joint optimization:
#   Stage 1: fit ka, kel, f_delayed, t_lag, F_12C to 12C alone.
#   Stage 2: with those fixed, fit F_13C6/k_release/f_delayed_13C6 to 13C6 alone.
# See "Two-stage fit: 12C first, then 13C6" and "13C6's own wave choice" in
# docs/pk-model.md.
# ---------------------------------------------------------------------------
fit_subject_visit_two_wave <- function(sid, vis, extra_seeds = list(), maxit = 80, single_wave_r2_12C = NA_real_) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                   f_delayed = NA, t_lag = NA, f_delayed_13C6 = NA,
                   t_lag1_13C6 = NA, f_delayed2_13C6 = NA, t_lag2_13C6 = NA,
                   r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                   converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA,
                   rss_12C = NA, n_12C = NA, t30_present = NA)

  obs12 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  obs13 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "13C6")
  if (nrow(obs12) < 3 || nrow(obs13) < 3) return(empty)
  # No two_wave without a t=30 sample: the 0-60min window would be
  # under-identified and AIC would rubber-stamp whatever 5 sparse points fit.
  # Model selection then falls back to single_wave.
  if (!any(obs12$time_min == 30)) return(mutate(empty, t30_present = FALSE))
  # Require evidence of a second wave in 12C's own raw data: a post-peak dip
  # (peak_dip_rise_info()) or a near-equal-height neighbor of the peak
  # (has_near_peak_neighbor()). dip_evidence$trigger_time is kept to seed the
  # search below. See "How the dip is detected" and "Choosing between
  # single_wave and two_wave" in docs/pk-model.md.
  dip_evidence <- peak_dip_rise_info(obs12$time_min, obs12$conc_mgL)
  plateau_evidence <- has_near_peak_neighbor(obs12$time_min, obs12$conc_mgL)
  # A third path in: single_wave fitting 12C very poorly (R2 below
  # R2_ALWAYS_TRY_TWO_WAVE) also triggers an attempt, to catch smoothly
  # accelerating rises neither detector sees (ER21 baseline). Set well below
  # R2_RELIABLE_MIN so it only fires when single_wave doesn't fit at all.
  R2_ALWAYS_TRY_TWO_WAVE <- 0.50
  poor_single_wave <- !is.na(single_wave_r2_12C) && single_wave_r2_12C < R2_ALWAYS_TRY_TWO_WAVE
  if (!dip_evidence$detected && !plateau_evidence && !poor_single_wave) return(mutate(empty, t30_present = TRUE))

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- watson_ecf_volume(cov$bw_kg, cov$height_cm, cov$age_years, cov$sex)
  dose_12C <- cov$dose_12C_mg

  # ---- Stage 1: ka, kel, f_delayed, t_lag, F_12C from 12C alone ----
  objective_12C <- function(theta) {
    simA <- function(t, dose) bateman_conc(t, theta[["ka"]], theta[["kel"]], theta[["F_12C"]], dose, Vd)
    pred12 <- tryCatch(simulate_lagged_dose(simA, obs12$time_min, dose_12C, theta[["f_delayed"]], theta[["t_lag"]]), error = function(e) NULL)
    if (is.null(pred12) || any(!is.finite(pred12))) return(1e10)
    ss_tot12 <- sum((obs12$conc_mgL - mean(obs12$conc_mgL))^2)
    total <- sum((pred12 - obs12$conc_mgL)^2) / ss_tot12

    fine12 <- simulate_lagged_dose(simA, FINE_T, dose_12C, theta[["f_delayed"]], theta[["t_lag"]])
    if (any(!is.finite(fine12))) return(1e10)
    tmax12 <- FINE_T[which.max(fine12)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax12)^2
    obs_cmax12 <- max(obs12$conc_mgL)
    excess12 <- max(0, abs(max(fine12) - obs_cmax12) / obs_cmax12 - CMAX_TOL)
    total <- total + CMAX_LAMBDA * excess12^2

    # Each wave's own peak must also clear MIN_TMAX, not just the combined
    # curve's, or the non-tallest wave's timing is unconstrained (see "13C6's
    # own wave choice" in docs/pk-model.md). Computed per wave because
    # simulate_lagged_dose() only returns the summed curve.
    fine_wave1 <- simA(FINE_T, (1 - theta[["f_delayed"]]) * dose_12C)
    fine_wave2 <- simA(pmax(FINE_T - theta[["t_lag"]], 0), theta[["f_delayed"]] * dose_12C)
    if (any(!is.finite(fine_wave1)) || any(!is.finite(fine_wave2))) return(1e10)
    tmax_wave1 <- FINE_T[which.max(fine_wave1)]
    tmax_wave2 <- FINE_T[which.max(fine_wave2)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax_wave1)^2
    total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax_wave2)^2
  }

  pre12 <- fit_curve_independent(obs12$time_min, obs12$conc_mgL, dose_12C, Vd)
  informed_seed_12C <- if (!is.null(pre12)) {
    list(c(ka = pre12$par[["ka"]], kel = pre12$par[["kel"]], F_12C = pre12$par[["F"]], f_delayed = 0.01, t_lag = 60))
  } else list()
  no_lag_seeds <- list(  # f_delayed near 0 should recover ~a plain Bateman fit
    c(ka = 0.02, kel = 0.02, F_12C = 0.02, f_delayed = 0.01, t_lag = 60),
    c(ka = 0.05, kel = 0.03, F_12C = 0.01, f_delayed = 0.01, t_lag = 90)
  )
  dip_seeds <- list(  # pilot-informed: real dip cases (ER01/ER09) converged near t_lag~90
    c(ka = 0.05, kel = 0.02, F_12C = 0.02, f_delayed = 0.5, t_lag = 60),
    c(ka = 0.03, kel = 0.02, F_12C = 0.02, f_delayed = 0.4, t_lag = 45),
    c(ka = 0.06, kel = 0.03, F_12C = 0.015, f_delayed = 0.35, t_lag = 75),
    c(ka = 0.04, kel = 0.025, F_12C = 0.02, f_delayed = 0.45, t_lag = 90)
  )
  # Seeds anchored on the evidence that justified trying two_wave
  # (dip_evidence$trigger_time), at and before it: generic seeds can land on an
  # unrelated near-total-delay optimum instead (ER04 baseline; see "Two-stage fit"
  # in docs/pk-model.md).
  evidence_seeds <- if (!is.na(dip_evidence$trigger_time)) {
    tt <- dip_evidence$trigger_time
    list(
      c(ka = 0.03, kel = 0.02, F_12C = 0.02, f_delayed = 0.4, t_lag = tt),
      c(ka = 0.04, kel = 0.025, F_12C = 0.02, f_delayed = 0.3, t_lag = max(5, tt - 20)),
      c(ka = 0.05, kel = 0.03, F_12C = 0.015, f_delayed = 0.35, t_lag = max(5, tt - 40))
    )
  } else list()
  # ka=kel diagonal seeds - see the comment in fit_curve_independent() for
  # why this region needs explicit seeding.
  diagonal_seeds_12C <- lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) {
    c(ka = v, kel = v, F_12C = 0.02, f_delayed = 0.3, t_lag = 60)
  })
  jitter_seeds_12C <- random_seeds(N_RANDOM_LAGGED, BOUNDS_LAGGED_12C, log_scale = c("ka", "kel"))
  # extra_seeds (from the retry driver) carry all 7 combined-model parameter
  # names - only the stage-1-relevant subset is used here.
  extra_seeds_12C <- lapply(extra_seeds, function(s) s[c("ka", "kel", "F_12C", "f_delayed", "t_lag")])
  seeds_12C <- c(informed_seed_12C, no_lag_seeds, dip_seeds, evidence_seeds, diagonal_seeds_12C, jitter_seeds_12C, extra_seeds_12C)

  # t_lag's lower bound for this subject: a step on their own t=30 sample (5 if
  # at or below EARLY_LAG_OK_MGL, else 30). Nothing is sampled between t=0 and
  # t=30, so a partial floor is as unsupported as none (see "The post-peak dip
  # and the lagged-dose model" in docs/pk-model.md).
  obs30_12C <- obs12$conc_mgL[obs12$time_min == 30]
  t_lag_lower <- if (obs30_12C <= EARLY_LAG_OK_MGL) BOUNDS_LAGGED_12C$t_lag[1] else 30

  lower12 <- vapply(BOUNDS_LAGGED_12C, `[`, numeric(1), 1)
  upper12 <- vapply(BOUNDS_LAGGED_12C, `[`, numeric(1), 2)
  lower12[["t_lag"]] <- t_lag_lower
  fit12 <- fit_multistart(objective_12C, lower12, upper12, seeds = seeds_12C, control = list(maxit = maxit))
  if (is.null(fit12)) return(empty)
  par12 <- fit12$par

  simA_fixed <- function(t, dose) bateman_conc(t, par12[["ka"]], par12[["kel"]], par12[["F_12C"]], dose, Vd)
  pred12 <- simulate_lagged_dose(simA_fixed, obs12$time_min, dose_12C, par12[["f_delayed"]], par12[["t_lag"]])
  r2_12C <- r_squared(obs12$conc_mgL, pred12)

  # ---- Stage 2: F_13C6, k_release and 13C6's own delayed fraction, from 13C6
  # alone, with ka/kel/t_lag fixed at stage 1's estimate (see "13C6's own wave
  # choice" in docs/pk-model.md).
  objective_13C6 <- function(theta) {
    simB <- function(t, dose) simulate_delayed_release(t, theta[["k_release"]], par12[["ka"]], par12[["kel"]], theta[["F_13C6"]], dose, Vd)$conc
    pred13 <- tryCatch(simulate_lagged_dose(simB, obs13$time_min, dose_13C6_mg, theta[["f_delayed_13C6"]], par12[["t_lag"]]), error = function(e) NULL)
    if (is.null(pred13) || any(!is.finite(pred13))) return(1e10)
    ss_tot13 <- sum((obs13$conc_mgL - mean(obs13$conc_mgL))^2)
    total <- sum((pred13 - obs13$conc_mgL)^2) / ss_tot13

    fine13 <- simulate_lagged_dose(simB, FINE_T, dose_13C6_mg, theta[["f_delayed_13C6"]], par12[["t_lag"]])
    if (any(!is.finite(fine13))) return(1e10)
    tmax13 <- FINE_T[which.max(fine13)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax13)^2
    obs_cmax13 <- max(obs13$conc_mgL)
    excess13 <- max(0, abs(max(fine13) - obs_cmax13) / obs_cmax13 - CMAX_TOL)
    total + CMAX_LAMBDA * excess13^2
  }

  fixed_seeds_13C6 <- list(  # explicitly span wave-1-only, wave-2-only, and split
    c(F_13C6 = 0.05, k_release = 0.05, f_delayed_13C6 = 0.01),
    c(F_13C6 = 0.05, k_release = 0.05, f_delayed_13C6 = 0.99),
    c(F_13C6 = 0.1,  k_release = 0.02, f_delayed_13C6 = 0.5),
    c(F_13C6 = 0.02, k_release = 0.3,  f_delayed_13C6 = 0.3)
  )
  extra_seeds_13C6 <- lapply(extra_seeds, function(s) s[c("F_13C6", "k_release", "f_delayed_13C6")])
  jitter_seeds_13C6 <- random_seeds(N_RANDOM_13C6_GIVEN_LAG, BOUNDS_13C6_GIVEN_LAG, log_scale = c("k_release"))
  seeds_13C6 <- c(fixed_seeds_13C6, jitter_seeds_13C6, extra_seeds_13C6)

  lower13 <- vapply(BOUNDS_13C6_GIVEN_LAG, `[`, numeric(1), 1)
  upper13 <- vapply(BOUNDS_13C6_GIVEN_LAG, `[`, numeric(1), 2)
  fit13 <- fit_multistart(objective_13C6, lower13, upper13, seeds = seeds_13C6, control = list(maxit = maxit))
  if (is.null(fit13)) return(empty)
  par13 <- fit13$par

  simB_fixed <- function(t, dose) simulate_delayed_release(t, par13[["k_release"]], par12[["ka"]], par12[["kel"]], par13[["F_13C6"]], dose, Vd)$conc
  pred13 <- simulate_lagged_dose(simB_fixed, obs13$time_min, dose_13C6_mg, par13[["f_delayed_13C6"]], par12[["t_lag"]])
  r2_13C6 <- r_squared(obs13$conc_mgL, pred13)

  # 13C6's own onset-lag and second-wave candidates, as in
  # fit_subject_visit_single_wave(), so 13C6's mechanism doesn't depend on
  # 12C's model (ER23 intervention needed this). All three candidates (baseline K=3,
  # onset-lag K=2, second-wave K=3) are compared by AIC on 13C6's own RSS.
  n13 <- length(obs13$conc_mgL)
  best_aic13 <- n13 * log(fit13$value / n13) + 2 * 3   # baseline: F_13C6/k_release/f_delayed_13C6, K=3
  t_lag1_13C6 <- NA_real_
  t_lag2_13C6 <- NA_real_
  f_delayed2_13C6 <- NA_real_

  if (has_onset_lag_evidence(obs13$time_min, obs13$conc_mgL)) {
    obj_13C6_lag <- function(theta) {
      simB <- function(t, dose) bateman_conc(t, par12[["ka"]], par12[["kel"]], theta[["F_13C6"]], dose, Vd)
      pred <- tryCatch(simulate_two_lag_dose(simB, obs13$time_min, dose_13C6_mg, f_delayed = 0, t_lag1 = theta[["t_lag1_13C6"]], gap = 0), error = function(e) NULL)
      if (is.null(pred) || any(!is.finite(pred))) return(1e10)
      sum((pred - obs13$conc_mgL)^2)
    }
    seeds_lag <- lapply(c(10, 20, 30, 40, 50), function(tl) c(F_13C6 = par13[["F_13C6"]], t_lag1_13C6 = tl))
    lower_lag <- c(F_13C6 = BOUNDS_JOINT$F_13C6[1], t_lag1_13C6 = 0)
    upper_lag <- c(F_13C6 = BOUNDS_JOINT$F_13C6[2], t_lag1_13C6 = 90)
    fit_lag <- fit_multistart(obj_13C6_lag, lower_lag, upper_lag, seeds = seeds_lag, control = list(maxit = 60))
    if (!is.null(fit_lag)) {
      aic_lag <- n13 * log(fit_lag$value / n13) + 2 * 2   # F_13C6/t_lag1_13C6, K=2
      if (aic_lag < best_aic13) {
        best_aic13 <- aic_lag
        par13[["F_13C6"]] <- fit_lag$par[["F_13C6"]]
        par13[["k_release"]] <- NA_real_
        par13[["f_delayed_13C6"]] <- NA_real_
        t_lag1_13C6 <- fit_lag$par[["t_lag1_13C6"]]
        simB_final <- function(t, dose) bateman_conc(t, par12[["ka"]], par12[["kel"]], par13[["F_13C6"]], dose, Vd)
        pred13 <- simulate_two_lag_dose(simB_final, obs13$time_min, dose_13C6_mg, f_delayed = 0, t_lag1 = t_lag1_13C6, gap = 0)
        r2_13C6 <- r_squared(obs13$conc_mgL, pred13)
      }
    }
  }

  if (has_peak_dip_rise(obs13$time_min, obs13$conc_mgL) ||
      has_near_peak_neighbor(obs13$time_min, obs13$conc_mgL)) {
    obj_13C6_tw <- function(theta) {
      simB <- function(t, dose) bateman_conc(t, par12[["ka"]], par12[["kel"]], theta[["F_13C6"]], dose, Vd)
      pred <- tryCatch(simulate_lagged_dose(simB, obs13$time_min, dose_13C6_mg, theta[["f_delayed2_13C6"]], theta[["t_lag2_13C6"]]), error = function(e) NULL)
      if (is.null(pred) || any(!is.finite(pred))) return(1e10)
      sum((pred - obs13$conc_mgL)^2)
    }
    seeds_tw <- lapply(c(40, 60, 90, 120), function(tl)
      c(F_13C6 = par13[["F_13C6"]], f_delayed2_13C6 = 0.5, t_lag2_13C6 = tl))
    lower_tw <- c(F_13C6 = BOUNDS_JOINT$F_13C6[1], f_delayed2_13C6 = 0.05, t_lag2_13C6 = 30)
    upper_tw <- c(F_13C6 = BOUNDS_JOINT$F_13C6[2], f_delayed2_13C6 = 0.95, t_lag2_13C6 = 150)
    fit_tw <- fit_multistart(obj_13C6_tw, lower_tw, upper_tw, seeds = seeds_tw, control = list(maxit = 60))
    # As in fit_subject_visit_single_wave(): pinned at the lower bound means no
    # real second wave, so reject rather than compare on AIC.
    at_lower_bound <- !is.null(fit_tw) && fit_tw$par[["f_delayed2_13C6"]] <= lower_tw[["f_delayed2_13C6"]] + 1e-3
    if (!is.null(fit_tw) && !at_lower_bound) {
      aic_tw <- n13 * log(fit_tw$value / n13) + 2 * 3   # F_13C6/f_delayed2_13C6/t_lag2_13C6, K=3
      if (aic_tw < best_aic13) {
        best_aic13 <- aic_tw
        par13[["F_13C6"]] <- fit_tw$par[["F_13C6"]]
        par13[["k_release"]] <- NA_real_
        par13[["f_delayed_13C6"]] <- NA_real_
        t_lag1_13C6 <- NA_real_   # mutually exclusive with the onset-lag candidate above
        f_delayed2_13C6 <- fit_tw$par[["f_delayed2_13C6"]]
        t_lag2_13C6 <- fit_tw$par[["t_lag2_13C6"]]
        simB_final <- function(t, dose) bateman_conc(t, par12[["ka"]], par12[["kel"]], par13[["F_13C6"]], dose, Vd)
        pred13 <- simulate_lagged_dose(simB_final, obs13$time_min, dose_13C6_mg, f_delayed2_13C6, t_lag2_13C6)
        r2_13C6 <- r_squared(obs13$conc_mgL, pred13)
      }
    }
  }

  tibble(
    ka = par12[["ka"]], kel = par12[["kel"]], F_12C = par12[["F_12C"]], F_13C6 = par13[["F_13C6"]],
    k_release = par13[["k_release"]], f_delayed = par12[["f_delayed"]], t_lag = par12[["t_lag"]],
    f_delayed_13C6 = par13[["f_delayed_13C6"]],
    t_lag1_13C6 = t_lag1_13C6, f_delayed2_13C6 = f_delayed2_13C6, t_lag2_13C6 = t_lag2_13C6,
    r2_12C = r2_12C, r2_13C6 = r2_13C6,
    kel_at_bound = par12[["kel"]] > (BOUNDS_LAGGED_12C$kel[2] - 1e-4),
    k_release_at_bound = isTRUE(par13[["k_release"]] > (BOUNDS_13C6_GIVEN_LAG$k_release[2] - 1e-4)),
    converged = fit12$convergence == 0 && fit13$convergence == 0,
    r2_12C_low = r2_12C < R2_RELIABLE_MIN,
    r2_13C6_low = r2_13C6 < R2_RELIABLE_MIN,
    objective_value = fit12$value + fit13$value,   # not comparable across subjects, only against this subject's own prior/retry pair
    rss_12C = sum((obs12$conc_mgL - pred12)^2),   # raw (unweighted, native mg/L^2) - for the single_wave-vs-two_wave AIC comparison below
    n_12C = length(obs12$conc_mgL),
    t30_present = TRUE   # the t30-missing check above already returned early otherwise
  )
}

# ---------------------------------------------------------------------------
# Run both candidate models for every subject x visit
# ---------------------------------------------------------------------------
# peak_dip_rise_info()'s $excess at or above which the raw data count as
# unambiguous evidence of a second wave, overriding the AIC comparison (see
# "Definite two_wave override" in docs/pk-model.md).
EXCESS_DEFINITE_TWO_WAVE <- 0.90

fit_one <- function(i) {
  sid <- subject_visits$subject_id[i]; vis <- subject_visits$visit[i]
  set.seed(string_seed(paste(sid, vis)))   # reproducible regardless of N_CORES or row order - see pk_fit.R
  cat(sid, vis, "(single_wave)\n")
  sw <- fit_subject_visit_single_wave(sid, vis)

  obs12_raw <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  dip_evidence <- if (nrow(obs12_raw) >= 4) peak_dip_rise_info(obs12_raw$time_min, obs12_raw$conc_mgL) else list(detected = FALSE, excess = NA_real_)
  definite_two_wave <- isTRUE(dip_evidence$detected) && !is.na(dip_evidence$excess) && dip_evidence$excess >= EXCESS_DEFINITE_TWO_WAVE

  # two_wave is not skipped on single_wave's own R2: once the raw-data gates in
  # fit_subject_visit_two_wave() pass, AIC decides (see "Choosing between
  # single_wave and two_wave" in docs/pk-model.md).
  set.seed(string_seed(paste(sid, vis, "two_wave")))
  cat(sid, vis, if (definite_two_wave) "(two_wave - definite evidence)" else "(two_wave)", "\n")
  tw <- fit_subject_visit_two_wave(sid, vis, single_wave_r2_12C = sw$r2_12C)
  list(single_wave = sw, two_wave = tw, definite_two_wave = definite_two_wave)
}

fit_list <- parallel::mclapply(seq_len(nrow(subject_visits)), fit_one, mc.cores = N_CORES)
single_wave_results <- bind_cols(subject_visits, bind_rows(lapply(fit_list, `[[`, "single_wave")))
two_wave_results    <- bind_cols(subject_visits, bind_rows(lapply(fit_list, `[[`, "two_wave")))
definite_two_wave_flags <- bind_cols(subject_visits, tibble(definite_two_wave = vapply(fit_list, `[[`, logical(1), "definite_two_wave")))

# Adaptive retry for flagged fits (boundary/poor-fit/non-convergence) - see
# "Adaptive retry" in docs/pk-model.md for why. Run separately per model.
refit_flagged_single_wave <- function(i, seeds, maxit) {
  sid <- single_wave_results$subject_id[i]; vis <- single_wave_results$visit[i]
  set.seed(string_seed(paste(sid, vis, "retry")))
  cat(sid, vis, "(single_wave retry)\n")
  fit_subject_visit_single_wave(sid, vis, extra_seeds = seeds, maxit = maxit)
}
single_wave_results <- adaptive_retry(
  single_wave_results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  # r2_13C6_low is excluded: a poor 13C6 fit is usually a structural
  # mismatch, not an under-searched optimum (see "Adaptive retry" in
  # docs/pk-model.md).
  extra_flag_cols = c("r2_12C_low"),
  convergence_col = "converged",
  bounds = BOUNDS_JOINT,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release"),
  refit_fn = refit_flagged_single_wave,
  log_scale = c("ka", "kel", "k_release"),
  mc.cores = N_CORES
)

refit_flagged_two_wave <- function(i, seeds, maxit) {
  sid <- two_wave_results$subject_id[i]; vis <- two_wave_results$visit[i]
  set.seed(string_seed(paste(sid, vis, "two_wave retry")))
  cat(sid, vis, "(two_wave retry)\n")
  # single_wave_results is already computed here (its retry runs first). Looked
  # up rather than threaded through adaptive_retry()'s generic signature, so a
  # fit that passed the raw-data gate via the R2_ALWAYS_TRY_TWO_WAVE override
  # isn't reset to "not attempted" on retry.
  sw_r2 <- single_wave_results$r2_12C[single_wave_results$subject_id == sid & single_wave_results$visit == vis]
  fit_subject_visit_two_wave(sid, vis, extra_seeds = seeds, maxit = maxit,
                              single_wave_r2_12C = if (length(sw_r2) == 1) sw_r2 else NA_real_)
}
two_wave_results <- adaptive_retry(
  two_wave_results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  # r2_13C6_low excluded here too, for the same reason as in single_wave's
  # retry above (13C6 has its own f_delayed_13C6, so a poor r2_13C6 isn't a
  # search failure).
  extra_flag_cols = c("r2_12C_low"),
  convergence_col = "converged",
  bounds = BOUNDS_LAGGED,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release", "f_delayed", "t_lag", "f_delayed_13C6"),
  refit_fn = refit_flagged_two_wave,
  log_scale = c("ka", "kel", "k_release"),
  maxit_retry = 100,
  mc.cores = N_CORES
)

# ---------------------------------------------------------------------------
# Model selection: single_wave vs two_wave per subject x visit, by AIC on 12C's
# own raw residuals only, so 13C6's fit quality can't leak into whether 12C has
# a second wave. 13C6 follows whichever structure wins and makes its own
# mechanism choice inside the fit functions. See "Choosing between single_wave
# and two_wave" in docs/pk-model.md.
# ---------------------------------------------------------------------------
K_12C_SINGLE_WAVE <- 3   # ka, kel, F_12C
K_12C_TWO_WAVE    <- 5   # + f_delayed, t_lag

# Plain AIC (no AICc): at n=8, AICc excludes every confirmed plateau case (see
# "Choosing between single_wave and two_wave" in docs/pk-model.md).
aic <- function(rss, n, k) n * log(rss / n) + 2 * k

sw <- single_wave_results %>% select(subject_id, visit, rss_12C, n_12C) %>%
  rename(rss_sw = rss_12C, n_sw = n_12C)
tw <- two_wave_results %>% select(subject_id, visit, rss_12C, n_12C) %>%
  rename(rss_tw = rss_12C, n_tw = n_12C)
selection <- sw %>% left_join(tw, by = c("subject_id", "visit")) %>%
  left_join(definite_two_wave_flags, by = c("subject_id", "visit")) %>%
  mutate(
    aic_single_wave = aic(rss_sw, n_sw, K_12C_SINGLE_WAVE),
    aic_two_wave    = aic(rss_tw, n_tw, K_12C_TWO_WAVE),
    model = case_when(
      is.na(aic_single_wave) & is.na(aic_two_wave) ~ NA_character_,
      is.na(aic_two_wave)                          ~ "single_wave",
      # definite_two_wave (see EXCESS_DEFINITE_TWO_WAVE) overrides the AIC
      # comparison; only a failed two_wave fit (NA, caught above) falls back.
      definite_two_wave %in% TRUE                  ~ "two_wave",
      is.na(aic_single_wave)                       ~ "two_wave",
      aic_two_wave < aic_single_wave                ~ "two_wave",
      TRUE                                          ~ "single_wave"
    )
  ) %>%
  select(subject_id, visit, model, aic_single_wave, aic_two_wave)

results <- selection %>%
  left_join(single_wave_results, by = c("subject_id", "visit"), suffix = c("", ".sw")) %>%
  left_join(two_wave_results, by = c("subject_id", "visit"), suffix = c("", ".tw"))

# Pick each output column from whichever model was selected for that row.
pick <- function(col) {
  sw_col <- results[[col]]; tw_col <- results[[paste0(col, ".tw")]]
  if_else(results$model == "two_wave", tw_col, sw_col)
}
for (col in c("ka", "kel", "F_12C", "F_13C6", "k_release", "r2_12C", "r2_13C6",
              "kel_at_bound", "k_release_at_bound", "converged", "r2_12C_low", "r2_13C6_low",
              # 13C6's own mechanism columns exist under both models, so they are
              # picked like any other shared column, not nulled by model.
              "t_lag1_13C6", "f_delayed2_13C6", "t_lag2_13C6")) {
  results[[col]] <- pick(col)
}
# f_delayed/t_lag/f_delayed_13C6 only exist in two_wave_results, so the join
# brings them in unsuffixed; null them where two_wave lost the AIC comparison so
# the losing fit's values don't leak through.
results$f_delayed <- if_else(results$model == "two_wave", results$f_delayed, NA_real_)
results$t_lag      <- if_else(results$model == "two_wave", results$t_lag, NA_real_)
results$f_delayed_13C6 <- if_else(results$model == "two_wave", results$f_delayed_13C6, NA_real_)
# t30_present comes from single_wave's copy via the join (always attempted). It
# describes the data, not the winning model - see "Output files" in
# docs/pk-model.md.
results <- results %>% select(subject_id, visit, model, ka, kel, F_12C, F_13C6, k_release,
                               f_delayed, t_lag, f_delayed_13C6, t_lag1_13C6,
                               f_delayed2_13C6, t_lag2_13C6, t30_present,
                               r2_12C, r2_13C6, kel_at_bound, k_release_at_bound, converged,
                               r2_12C_low, r2_13C6_low, aic_single_wave, aic_two_wave)

cat("\n=== Model selection ===\n")
print(table(results$model, useNA = "ifany"))

dir.create("results", showWarnings = FALSE)
write_csv(single_wave_results, "results/fit_results_single_wave.csv")
write_csv(two_wave_results, "results/fit_results_two_wave.csv")
write_csv(results, "results/fit_results.csv")

# Plots: one figure per subject (both isotopes x both visits) using each subject
# x visit's selected model, against the observed points (including t=0, unlike
# data_fit). 12C's panel also shows both candidates' AIC. See "Output files" in
# docs/pk-model.md.
simulate_fit <- function(sid, vis, r) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- watson_ecf_volume(cov$bw_kg, cov$height_cm, cov$age_years, cov$sex)
  dose_12C <- cov$dose_12C_mg
  fine <- seq(0, 400, length.out = 400)

  if (r$model == "two_wave") {
    simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_12C, dose, Vd)
    sim12 <- simulate_lagged_dose(simA, fine, dose_12C, r$f_delayed, r$t_lag)
  } else {
    sim12 <- bateman_conc(fine, r$ka, r$kel, r$F_12C, dose_12C, Vd)
  }

  # 13C6's mechanism is chosen independently of 12C's model and can be any of
  # the four below; this must dispatch exactly as the fit did, or the curve won't
  # match the r2_13C6 it is annotated with.
  simB_instant <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_13C6, dose, Vd)
  sim13 <- if (!is.na(r$t_lag1_13C6)) {
    simulate_two_lag_dose(simB_instant, fine, dose_13C6_mg, f_delayed = 0, t_lag1 = r$t_lag1_13C6, gap = 0)
  } else if (!is.na(r$t_lag2_13C6)) {
    simulate_lagged_dose(simB_instant, fine, dose_13C6_mg, r$f_delayed2_13C6, r$t_lag2_13C6)
  } else if (r$model == "two_wave") {
    simB <- function(t, dose) simulate_delayed_release(t, r$k_release, r$ka, r$kel, r$F_13C6, dose, Vd)$conc
    simulate_lagged_dose(simB, fine, dose_13C6_mg, r$f_delayed_13C6, r$t_lag)
  } else {
    simulate_delayed_release(fine, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, Vd)$conc
  }

  t_end <- max(time_to_clearance(fine, sim12, CLEARANCE_FRAC),
               time_to_clearance(fine, sim13, CLEARANCE_FRAC))
  keep <- fine <= t_end

  bind_rows(
    tibble(time_min = fine[keep], conc_mgL = sim12[keep], isotope = "12C"),
    tibble(time_min = fine[keep], conc_mgL = sim13[keep], isotope = "13C6")
  )
}

plot_subject_fit <- function(sid) {
  obs <- data_corrected %>% filter(subject_id == sid) %>% select(visit, isotope, time_min, conc_mgL)
  fits <- results %>% filter(subject_id == sid, !is.na(ka))
  if (nrow(fits) == 0) return(NULL)

  sim <- fits %>% pmap_dfr(function(...) {
    row <- tibble(...)
    bind_cols(visit = row$visit, simulate_fit(sid, row$visit, row))
  })

  fmt_aic <- function(x) if_else(is.na(x), "NA", sprintf("%.0f", x))
  ann <- bind_rows(
    fits %>% transmute(visit, isotope = "12C",
                        label = sprintf("R2=%.3f\nAIC sw/tw=%s/%s", r2_12C,
                                        fmt_aic(aic_single_wave), fmt_aic(aic_two_wave))),
    fits %>% transmute(visit, isotope = "13C6", label = sprintf("R2=%.3f", r2_13C6))
  )

  ggplot(obs, aes(time_min, conc_mgL)) +
    geom_point(color = "black", size = 1.4) +
    geom_line(data = sim, aes(color = isotope), linewidth = 0.8, show.legend = FALSE) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.2, size = 2.5, lineheight = 0.9) +
    facet_grid(rows = vars(isotope), cols = vars(visit), scales = "free",
               labeller = labeller(visit = VISIT_LABELS, isotope = ISOTOPE_LABELS)) +
    scale_color_manual(values = ISOTOPE_COLORS) +
    labs(title = sid, x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
}

dir.create("results/plots_individual", showWarnings = FALSE, recursive = TRUE)
for (sid in unique(subject_visits$subject_id)) {
  p <- plot_subject_fit(sid)
  if (!is.null(p)) ggsave(file.path("results/plots_individual", sprintf("%s_joint.pdf", sid)),
                           p, width = 8, height = 5, dpi = 120)
}

cat("\nDone. results/fit_results.csv (selected model per subject x visit),",
    "results/fit_results_single_wave.csv, results/fit_results_two_wave.csv (both candidates),",
    "results/plots_individual/*.pdf\n")
