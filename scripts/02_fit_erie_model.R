# Fitting the joint 12C/13C6 fructose PK model
# Barbara Verhaar
#
# Fits TWO models per subject x visit: the baseline joint delayed-release
# model, and a lagged-second-dose extension (simulate_lagged_dose, added to
# pk_curves.R) that can represent the post-peak dip seen in a meaningful
# minority of curves - the baseline model is structurally incapable of it (a
# one-compartment absorption curve is monotonic after its single peak for
# any parameter values). The lagged model degrades gracefully when a subject
# doesn't need it (f_delayed fits near 0), so both are reported side by side
# rather than picking one model for the whole cohort.
#
# Both now use an UNWEIGHTED objective (proportional_weighting = FALSE),
# changed from the previous default. Confirmed by direct comparison on ER03
# FCT1's 13C6 curve - catastrophic under the old proportionally-weighted
# objective (R2 = -2.18) despite MIN_TMAX and bounds unchanged - that
# unweighted alone (ordinary multistart, no special seeding) recovers R2 =
# 0.37, and that former_models/MixedModel's own historical fit for the exact
# same subject, built with an unweighted objective all along, reached R2 =
# 0.71. See the proportional_weighting doc in scripts/assets/pk_fit.R for
# the full comparison and the tradeoff being made (the weighted version was
# separately added to fix a measured ~24x peak-vs-tail heteroscedasticity;
# removing it cohort-wide is a real tradeoff, not a strict improvement -
# check r2_12C/r2_13C6 broadly, not just on known-hard subjects, before
# trusting this as the final word).
#
# NOT included here (deliberately, pending further validation): the
# MIN_TMAX floor is left at 30 (unchanged) rather than loosened further,
# since once the weighting is fixed there's no confirmed evidence it's still
# the binding constraint; and a two-lag onset-time extension
# (simulate_two_lag_dose, also in pk_curves.R) that matched
# former_models/MixedModel's R2 on ER03 FCT1 almost exactly but on n=1,
# didn't converge, and fit the 12C curve worse than the simpler unweighted
# baseline - not enough validation yet to commit a full cohort run to it.

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
VISIT_LABELS   <- c(FCT1 = "FCT1 (before diet)", FCT2 = "FCT2 (after diet)")
MODEL_COLORS   <- c(baseline = "grey45", lagged = "firebrick")
MODEL_LABELS   <- c(baseline = "Baseline", lagged = "Lagged second dose")

# Config
MIN_TMAX    <- 30     # min - soft floor on predicted Tmax for both curves - see file header
CMAX_TOL    <- 0.10   # +/-10% soft band around each curve's own observed Cmax
CMAX_LAMBDA <- 20     # penalty weight for the Cmax band
TMAX_LAMBDA <- 50     # penalty weight for the Tmax floor
# kel: Hannou et al. 2018 (t1/2 ~ 7-140 min). ka: bounded only below in
# spirit (absorption-rate differences are part of the research question).
BOUNDS_INDEP <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F = c(1e-4, 1))
BOUNDS_JOINT <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1),
                     F_12C = c(1e-4, 1), F_13C6 = c(1e-4, 1),
                     k_release = c(0.001, 1))   # capsule dissolution t1/2 ~ 0.7-700 min
BOUNDS_LAGGED <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1),
                      F_12C = c(1e-4, 1), F_13C6 = c(1e-4, 1),
                      k_release = c(0.001, 1),
                      # t_lag upper bound: sampling is 30-min-spaced through
                      # t=180, then widens to 60 (180->240) and 120 (240->360)
                      # min gaps. A t_lag landing past 180 lets the optimizer
                      # place an entire invented second peak inside one of
                      # those wide, unsampled gaps - fits nothing, is
                      # penalized by nothing, and can silently make the fit
                      # worse than baseline (found on ER12 FCT2: t_lag=179,
                      # r2_13C6 dropped from baseline's 0.95 to 0.75, and the
                      # resulting curve shows a large peak with zero
                      # supporting data around t=200-220). Both validated
                      # genuine-dip subjects (ER01, ER09) converged to
                      # t_lag in 60-114 min, well inside this bound, so this
                      # doesn't constrain any real case found so far.
                      f_delayed = c(0.001, 0.999), t_lag = c(5, 150))
N_RANDOM_INDEP <- 40
N_RANDOM_JOINT <- 16
N_RANDOM_LAGGED <- 60   # pilot-validated budget (ER01/ER03/ER09) - see results/model-post-peak-dip-pilot/
FINE_T <- seq(0, 400, length.out = 50)    # grid for Tmax/Cmax penalty checks
CLEARANCE_FRAC <- 0.01                    # "fully cleared", for plot x-axis limits only
R2_RELIABLE_MIN <- 0.70                   # below this, that curve's fit is flagged unreliable
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = max(1, parallel::detectCores() - 2))) # To parallelize loop

# NOTE: no top-level set.seed() here - it would not actually make the
# multi-start random seeds reproducible under parallel::mclapply() (see the
# note in scripts/assets/pk_fit.R). Each subject x visit fit instead seeds
# itself deterministically from its own ID in fit_one() below.

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

# t=0 carries no fitting information (12C is 0 by construction after baseline
# correction; 13C6 hasn't been dosed yet - both models predict exactly 0
# there regardless of parameters), so it's excluded here - but
# data_corrected itself keeps t=0, so the plots below (which read from
# data_corrected, not data_fit) still show the observed point.
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

  seeds <- c(
    list(c(ka = 0.03, kel = 0.02, F = 0.3), c(ka = 0.08, kel = 0.05, F = 0.15),
         c(ka = 0.01, kel = 0.01, F = 0.5), c(ka = 0.05, kel = 0.08, F = 0.2)),
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
fit_subject_visit_baseline <- function(sid, vis, extra_seeds = list(), maxit = 40) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                   r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                   converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA)

  obs12 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  obs13 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "13C6")
  if (nrow(obs12) < 3 || nrow(obs13) < 3) return(empty)

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
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
  jitter_seeds <- random_seeds(N_RANDOM_JOINT, BOUNDS_JOINT, log_scale = c("ka", "kel", "k_release"))

  lower <- vapply(BOUNDS_JOINT, `[`, numeric(1), 1)
  upper <- vapply(BOUNDS_JOINT, `[`, numeric(1), 2)
  fit <- fit_multistart(objective, lower, upper,
                         seeds = c(informed_seeds, grid_seeds_joint, jitter_seeds, extra_seeds),
                         control = list(maxit = maxit))
  if (is.null(fit)) return(empty)

  par <- fit$par
  pred12 <- curve_12C$simulate(par)
  pred13 <- curve_13C6$simulate(par)

  r2_12C  <- r_squared(curve_12C$conc, pred12)
  r2_13C6 <- r_squared(curve_13C6$conc, pred13)

  tibble(
    ka = par[["ka"]], kel = par[["kel"]], F_12C = par[["F_12C"]], F_13C6 = par[["F_13C6"]],
    k_release = par[["k_release"]],
    r2_12C  = r2_12C,
    r2_13C6 = r2_13C6,
    kel_at_bound = kel > (BOUNDS_JOINT$kel[2] - 1e-4),
    k_release_at_bound = par[["k_release"]] > (BOUNDS_JOINT$k_release[2] - 1e-4),
    converged = fit$convergence == 0,   # see "Fit quality and what to trust" in docs/pk-model.md
    r2_12C_low  = r2_12C  < R2_RELIABLE_MIN,
    r2_13C6_low = r2_13C6 < R2_RELIABLE_MIN,
    objective_value = fit$value   # for comparing against a retry; not meaningful across subjects
  )
}

# ---------------------------------------------------------------------------
# Lagged-second-dose joint fit: same shared ka/kel/F/k_release structure as
# baseline, plus a shared f_delayed/t_lag (see simulate_lagged_dose() in
# pk_curves.R) - a fraction of each dose contributes nothing until t_lag,
# then behaves like a fresh dose given at that later time. Shared across
# curves because both doses are ingested by the same subject at the same
# time and (per the pilot investigation) the trigger is plausibly an
# upstream, shared gastric-emptying event, not something curve-specific.
# ---------------------------------------------------------------------------
fit_subject_visit_lagged <- function(sid, vis, extra_seeds = list(), maxit = 80) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                   f_delayed = NA, t_lag = NA,
                   r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                   converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA)

  obs12 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  obs13 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "13C6")
  if (nrow(obs12) < 3 || nrow(obs13) < 3) return(empty)

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
  dose_12C <- cov$dose_12C_mg

  objective <- function(theta) {
    simA <- function(t, dose) bateman_conc(t, theta[["ka"]], theta[["kel"]], theta[["F_12C"]], dose, Vd)
    simB <- function(t, dose) simulate_delayed_release(t, theta[["k_release"]], theta[["ka"]], theta[["kel"]], theta[["F_13C6"]], dose, Vd)$conc
    pred12 <- tryCatch(simulate_lagged_dose(simA, obs12$time_min, dose_12C, theta[["f_delayed"]], theta[["t_lag"]]), error = function(e) NULL)
    pred13 <- tryCatch(simulate_lagged_dose(simB, obs13$time_min, dose_13C6_mg, theta[["f_delayed"]], theta[["t_lag"]]), error = function(e) NULL)
    if (is.null(pred12) || is.null(pred13) || any(!is.finite(pred12)) || any(!is.finite(pred13))) return(1e10)
    ss_tot12 <- sum((obs12$conc_mgL - mean(obs12$conc_mgL))^2)
    ss_tot13 <- sum((obs13$conc_mgL - mean(obs13$conc_mgL))^2)
    total <- sum((pred12 - obs12$conc_mgL)^2) / ss_tot12 + sum((pred13 - obs13$conc_mgL)^2) / ss_tot13

    fine12 <- simulate_lagged_dose(simA, FINE_T, dose_12C, theta[["f_delayed"]], theta[["t_lag"]])
    fine13 <- simulate_lagged_dose(simB, FINE_T, dose_13C6_mg, theta[["f_delayed"]], theta[["t_lag"]])
    if (any(!is.finite(fine12)) || any(!is.finite(fine13))) return(1e10)
    tmax12 <- FINE_T[which.max(fine12)]; tmax13 <- FINE_T[which.max(fine13)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax12)^2 + TMAX_LAMBDA * max(0, MIN_TMAX - tmax13)^2
    obs_cmax12 <- max(obs12$conc_mgL); obs_cmax13 <- max(obs13$conc_mgL)
    excess12 <- max(0, abs(max(fine12) - obs_cmax12) / obs_cmax12 - CMAX_TOL)
    excess13 <- max(0, abs(max(fine13) - obs_cmax13) / obs_cmax13 - CMAX_TOL)
    total + CMAX_LAMBDA * (excess12^2 + excess13^2)
  }

  no_lag_seeds <- list(  # f_delayed near 0 should recover ~baseline
    c(ka = 0.02, kel = 0.02, F_12C = 0.02, F_13C6 = 0.05, k_release = 0.05, f_delayed = 0.01, t_lag = 60),
    c(ka = 0.05, kel = 0.03, F_12C = 0.01, F_13C6 = 0.10, k_release = 0.02, f_delayed = 0.01, t_lag = 90)
  )
  dip_seeds <- list(  # pilot-informed: real dip cases (ER01/ER09) converged near t_lag~90
    c(ka = 0.05, kel = 0.02, F_12C = 0.02, F_13C6 = 0.05, k_release = 0.05, f_delayed = 0.5, t_lag = 60),
    c(ka = 0.03, kel = 0.02, F_12C = 0.02, F_13C6 = 0.05, k_release = 0.05, f_delayed = 0.4, t_lag = 45),
    c(ka = 0.06, kel = 0.03, F_12C = 0.015, F_13C6 = 0.08, k_release = 0.03, f_delayed = 0.35, t_lag = 75),
    c(ka = 0.04, kel = 0.025, F_12C = 0.02, F_13C6 = 0.06, k_release = 0.04, f_delayed = 0.45, t_lag = 90)
  )
  jitter_seeds <- random_seeds(N_RANDOM_LAGGED, BOUNDS_LAGGED, log_scale = c("ka", "kel", "k_release"))
  seeds <- c(no_lag_seeds, dip_seeds, jitter_seeds, extra_seeds)

  lower <- vapply(BOUNDS_LAGGED, `[`, numeric(1), 1)
  upper <- vapply(BOUNDS_LAGGED, `[`, numeric(1), 2)
  fit <- fit_multistart(objective, lower, upper, seeds = seeds, control = list(maxit = maxit))
  if (is.null(fit)) return(empty)

  par <- fit$par
  simA <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], par[["F_12C"]], dose, Vd)
  simB <- function(t, dose) simulate_delayed_release(t, par[["k_release"]], par[["ka"]], par[["kel"]], par[["F_13C6"]], dose, Vd)$conc
  pred12 <- simulate_lagged_dose(simA, obs12$time_min, dose_12C, par[["f_delayed"]], par[["t_lag"]])
  pred13 <- simulate_lagged_dose(simB, obs13$time_min, dose_13C6_mg, par[["f_delayed"]], par[["t_lag"]])

  r2_12C  <- r_squared(obs12$conc_mgL, pred12)
  r2_13C6 <- r_squared(obs13$conc_mgL, pred13)

  tibble(
    ka = par[["ka"]], kel = par[["kel"]], F_12C = par[["F_12C"]], F_13C6 = par[["F_13C6"]],
    k_release = par[["k_release"]], f_delayed = par[["f_delayed"]], t_lag = par[["t_lag"]],
    r2_12C = r2_12C, r2_13C6 = r2_13C6,
    kel_at_bound = par[["kel"]] > (BOUNDS_LAGGED$kel[2] - 1e-4),
    k_release_at_bound = par[["k_release"]] > (BOUNDS_LAGGED$k_release[2] - 1e-4),
    converged = fit$convergence == 0,
    r2_12C_low = r2_12C < R2_RELIABLE_MIN,
    r2_13C6_low = r2_13C6 < R2_RELIABLE_MIN,
    objective_value = fit$value
  )
}

# ---------------------------------------------------------------------------
# Run both models for every subject x visit
# ---------------------------------------------------------------------------
fit_one <- function(i) {
  sid <- subject_visits$subject_id[i]; vis <- subject_visits$visit[i]
  set.seed(string_seed(paste(sid, vis)))   # reproducible regardless of N_CORES or row order - see pk_fit.R
  cat(sid, vis, "(baseline)\n")
  base <- fit_subject_visit_baseline(sid, vis)
  set.seed(string_seed(paste(sid, vis, "lagged")))
  cat(sid, vis, "(lagged)\n")
  lag <- fit_subject_visit_lagged(sid, vis)
  list(baseline = base, lagged = lag)
}

fit_list <- parallel::mclapply(seq_len(nrow(subject_visits)), fit_one, mc.cores = N_CORES)
baseline_results <- bind_cols(subject_visits, bind_rows(lapply(fit_list, `[[`, "baseline")))
lagged_results   <- bind_cols(subject_visits, bind_rows(lapply(fit_list, `[[`, "lagged")))

# Adaptive retry for flagged fits (boundary/poor-fit/non-convergence) - see
# "Adaptive retry" in docs/pk-model.md for why. Run separately per model.
refit_flagged_baseline <- function(i, seeds, maxit) {
  sid <- baseline_results$subject_id[i]; vis <- baseline_results$visit[i]
  set.seed(string_seed(paste(sid, vis, "retry")))
  cat(sid, vis, "(baseline retry)\n")
  fit_subject_visit_baseline(sid, vis, extra_seeds = seeds, maxit = maxit)
}
baseline_results <- adaptive_retry(
  baseline_results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  extra_flag_cols = c("r2_12C_low", "r2_13C6_low"),
  convergence_col = "converged",
  bounds = BOUNDS_JOINT,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release"),
  refit_fn = refit_flagged_baseline,
  log_scale = c("ka", "kel", "k_release"),
  mc.cores = N_CORES
)

refit_flagged_lagged <- function(i, seeds, maxit) {
  sid <- lagged_results$subject_id[i]; vis <- lagged_results$visit[i]
  set.seed(string_seed(paste(sid, vis, "lagged retry")))
  cat(sid, vis, "(lagged retry)\n")
  fit_subject_visit_lagged(sid, vis, extra_seeds = seeds, maxit = maxit)
}
lagged_results <- adaptive_retry(
  lagged_results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  extra_flag_cols = c("r2_12C_low", "r2_13C6_low"),
  convergence_col = "converged",
  bounds = BOUNDS_LAGGED,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release", "f_delayed", "t_lag"),
  refit_fn = refit_flagged_lagged,
  log_scale = c("ka", "kel", "k_release"),
  maxit_retry = 100,
  mc.cores = N_CORES
)

baseline_results$capsule_dissolution_halflife_min <- log(2) / baseline_results$k_release
lagged_results$capsule_dissolution_halflife_min <- log(2) / lagged_results$k_release

dir.create("results", showWarnings = FALSE)
write_csv(baseline_results, "results/fit_results_baseline.csv")
write_csv(lagged_results, "results/fit_results_lagged.csv")

# Plots: one figure per subject, both isotopes x both visits, both models
# overlaid against the observed points (which include t=0, unlike the
# t>0-only data_fit used for fitting).
simulate_fit <- function(sid, vis, r, model) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
  dose_12C <- cov$dose_12C_mg
  fine <- seq(0, 400, length.out = 400)

  if (model == "baseline") {
    sim12 <- bateman_conc(fine, r$ka, r$kel, r$F_12C, dose_12C, Vd)
    sim13 <- simulate_delayed_release(fine, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, Vd)$conc
  } else {
    simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_12C, dose, Vd)
    simB <- function(t, dose) simulate_delayed_release(t, r$k_release, r$ka, r$kel, r$F_13C6, dose, Vd)$conc
    sim12 <- simulate_lagged_dose(simA, fine, dose_12C, r$f_delayed, r$t_lag)
    sim13 <- simulate_lagged_dose(simB, fine, dose_13C6_mg, r$f_delayed, r$t_lag)
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
  fits_base <- baseline_results %>% filter(subject_id == sid, !is.na(ka))
  fits_lag  <- lagged_results   %>% filter(subject_id == sid, !is.na(ka))
  if (nrow(fits_base) == 0 && nrow(fits_lag) == 0) return(NULL)

  sim <- bind_rows(
    fits_base %>% pmap_dfr(function(...) {
      row <- tibble(...)
      bind_cols(visit = row$visit, model = "baseline", simulate_fit(sid, row$visit, row, "baseline"))
    }),
    fits_lag %>% pmap_dfr(function(...) {
      row <- tibble(...)
      bind_cols(visit = row$visit, model = "lagged", simulate_fit(sid, row$visit, row, "lagged"))
    })
  )

  ggplot(obs, aes(time_min, conc_mgL)) +
    geom_point(color = "black", size = 1.4) +
    geom_line(data = sim, aes(color = model, linetype = model), linewidth = 0.8) +
    facet_grid(rows = vars(isotope), cols = vars(visit), scales = "free",
               labeller = labeller(visit = VISIT_LABELS, isotope = ISOTOPE_LABELS)) +
    scale_color_manual(values = MODEL_COLORS, labels = MODEL_LABELS, name = NULL) +
    scale_linetype_manual(values = c(baseline = "dashed", lagged = "solid"), labels = MODEL_LABELS, name = NULL) +
    labs(title = sid, x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
}

dir.create("results/plots_individual", showWarnings = FALSE, recursive = TRUE)
for (sid in unique(subject_visits$subject_id)) {
  p <- plot_subject_fit(sid)
  if (!is.null(p)) ggsave(file.path("results/plots_individual", sprintf("%s_joint.pdf", sid)),
                           p, width = 8, height = 5, dpi = 120)
}

cat("\nDone. results/fit_results_baseline.csv, results/fit_results_lagged.csv,",
    "results/plots_individual/*.pdf\n")
