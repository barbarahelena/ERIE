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
DIET_LABELS    <- c(low_fructose = "Diet A: low fructose", high_fructose = "Diet B: high fructose")
DIET_COLORS    <- c(low_fructose = "#1b9e77", high_fructose = "#d95f02")

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
# BOUNDS_LAGGED above stays as the "combined" 7-parameter view the generic
# adaptive_retry() driver operates over (bound_cols/par_cols); the actual
# two-stage fit below (fit_subject_visit_lagged()) uses these two subsets.
BOUNDS_LAGGED_12C <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F_12C = c(1e-4, 1),
                          f_delayed = c(0.001, 0.999), t_lag = c(5, 150))
BOUNDS_13C6_GIVEN_LAG <- list(F_13C6 = c(1e-4, 1), k_release = c(0.001, 1))
N_RANDOM_INDEP <- 40
N_RANDOM_JOINT <- 16
N_RANDOM_LAGGED <- 60   # pilot-validated budget (ER01/ER03/ER09) - see results/model-post-peak-dip-pilot/
N_RANDOM_13C6_GIVEN_LAG <- 40  # stage 2 is only 2-dim and cheap (no 12C simulation), so a generous budget
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
fit_subject_visit_single_wave <- function(sid, vis, extra_seeds = list(), maxit = 40) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                   r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                   converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA,
                   rss_12C = NA, n_12C = NA, t30_present = NA)

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
    objective_value = fit$value,   # for comparing against a retry; not meaningful across subjects
    rss_12C = sum((curve_12C$conc - pred12)^2),   # raw (unweighted, native mg/L^2) - for the single_wave-vs-two_wave AIC comparison below
    n_12C = length(curve_12C$conc),
    t30_present = any(obs12$time_min == 30)   # ka is poorly anchored without this sample - see fit_subject_visit_two_wave()
  )
}

# ---------------------------------------------------------------------------
# Two-wave (lagged-second-dose) fit, in TWO STAGES rather than one joint optimization:
#   Stage 1: fit ka, kel, f_delayed, t_lag, F_12C to 12C ALONE.
#   Stage 2: with those FIXED, fit only F_13C6/k_release to 13C6.
#
# An earlier version shared f_delayed/t_lag symmetrically, informed by both
# curves' own residuals in one joint objective - this turned out to be
# unreliable. 13C6 already has its own dedicated onset-delay parameter
# (k_release, the enteric capsule's dissolution rate); letting the SHARED
# lag also be jointly informed by 13C6's residuals let 13C6's genuinely
# different, k_release-explained slow onset get misattributed to a "second
# wave" and imposed onto 12C even when 12C's own data gave zero support for
# one. Confirmed concretely on ER32 FCT1 (12C's t=30 sample is already near
# its eventual peak - no plausible onset delay - yet the shared fit still
# pinned f_delayed~1, t_lag~19min, and made BOTH curves' R2 worse than the
# plain baseline model) and ER35 (same pattern). 12C is the right curve to
# derive the lag from on its own merits, not just to route around this
# failure: it has a far larger, cleaner signal, and it's also the
# physiologically primary trigger for a real biphasic-emptying event (the
# large 1g/kg osmotic/caloric liquid load, not the tiny fixed 13C6 tracer),
# so a genuine dip should show up in 12C's own data if it's real at all.
# ---------------------------------------------------------------------------
fit_subject_visit_two_wave <- function(sid, vis, extra_seeds = list(), maxit = 80) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                   f_delayed = NA, t_lag = NA,
                   r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                   converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA,
                   rss_12C = NA, n_12C = NA, t30_present = NA)

  obs12 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  obs13 <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "13C6")
  if (nrow(obs12) < 3 || nrow(obs13) < 3) return(empty)
  # A genuine second wave in the 0-60min window can only be confirmed (or
  # ruled out) against an actual t=30 sample - without it, that whole window
  # has only its t=0 and t=60 endpoints to anchor a 5-parameter fit, which
  # both under-identifies the model and removes the only data point that
  # could tell a real early dip apart from an unconstrained one. Rather than
  # let AIC quietly rubber-stamp whatever shape 5 sparse points can fit
  # almost exactly, don't attempt two_wave at all when t=30 is missing -
  # model selection below then has no two_wave candidate to prefer, so
  # single_wave is used automatically. Found via ER05 FCT1 (t=30, 90, 150,
  # 240 all missing - only 5 real points for 12C).
  if (!any(obs12$time_min == 30)) return(mutate(empty, t30_present = FALSE))

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
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
    total + CMAX_LAMBDA * excess12^2
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
  jitter_seeds_12C <- random_seeds(N_RANDOM_LAGGED, BOUNDS_LAGGED_12C, log_scale = c("ka", "kel"))
  # extra_seeds (from the retry driver) carry all 7 combined-model parameter
  # names - only the stage-1-relevant subset is used here.
  extra_seeds_12C <- lapply(extra_seeds, function(s) s[c("ka", "kel", "F_12C", "f_delayed", "t_lag")])
  seeds_12C <- c(informed_seed_12C, no_lag_seeds, dip_seeds, jitter_seeds_12C, extra_seeds_12C)

  lower12 <- vapply(BOUNDS_LAGGED_12C, `[`, numeric(1), 1)
  upper12 <- vapply(BOUNDS_LAGGED_12C, `[`, numeric(1), 2)
  fit12 <- fit_multistart(objective_12C, lower12, upper12, seeds = seeds_12C, control = list(maxit = maxit))
  if (is.null(fit12)) return(empty)
  par12 <- fit12$par

  simA_fixed <- function(t, dose) bateman_conc(t, par12[["ka"]], par12[["kel"]], par12[["F_12C"]], dose, Vd)
  pred12 <- simulate_lagged_dose(simA_fixed, obs12$time_min, dose_12C, par12[["f_delayed"]], par12[["t_lag"]])
  r2_12C <- r_squared(obs12$conc_mgL, pred12)

  # ---- Stage 2: F_13C6, k_release from 13C6, with ka/kel/f_delayed/t_lag fixed at stage 1's estimate ----
  objective_13C6 <- function(theta) {
    simB <- function(t, dose) simulate_delayed_release(t, theta[["k_release"]], par12[["ka"]], par12[["kel"]], theta[["F_13C6"]], dose, Vd)$conc
    pred13 <- tryCatch(simulate_lagged_dose(simB, obs13$time_min, dose_13C6_mg, par12[["f_delayed"]], par12[["t_lag"]]), error = function(e) NULL)
    if (is.null(pred13) || any(!is.finite(pred13))) return(1e10)
    ss_tot13 <- sum((obs13$conc_mgL - mean(obs13$conc_mgL))^2)
    total <- sum((pred13 - obs13$conc_mgL)^2) / ss_tot13

    fine13 <- simulate_lagged_dose(simB, FINE_T, dose_13C6_mg, par12[["f_delayed"]], par12[["t_lag"]])
    if (any(!is.finite(fine13))) return(1e10)
    tmax13 <- FINE_T[which.max(fine13)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax13)^2
    obs_cmax13 <- max(obs13$conc_mgL)
    excess13 <- max(0, abs(max(fine13) - obs_cmax13) / obs_cmax13 - CMAX_TOL)
    total + CMAX_LAMBDA * excess13^2
  }

  fixed_seeds_13C6 <- list(c(F_13C6 = 0.05, k_release = 0.05), c(F_13C6 = 0.1, k_release = 0.02),
                            c(F_13C6 = 0.02, k_release = 0.3))
  extra_seeds_13C6 <- lapply(extra_seeds, function(s) s[c("F_13C6", "k_release")])
  jitter_seeds_13C6 <- random_seeds(N_RANDOM_13C6_GIVEN_LAG, BOUNDS_13C6_GIVEN_LAG, log_scale = c("k_release"))
  seeds_13C6 <- c(fixed_seeds_13C6, jitter_seeds_13C6, extra_seeds_13C6)

  lower13 <- vapply(BOUNDS_13C6_GIVEN_LAG, `[`, numeric(1), 1)
  upper13 <- vapply(BOUNDS_13C6_GIVEN_LAG, `[`, numeric(1), 2)
  fit13 <- fit_multistart(objective_13C6, lower13, upper13, seeds = seeds_13C6, control = list(maxit = maxit))
  if (is.null(fit13)) return(empty)
  par13 <- fit13$par

  simB_fixed <- function(t, dose) simulate_delayed_release(t, par13[["k_release"]], par12[["ka"]], par12[["kel"]], par13[["F_13C6"]], dose, Vd)$conc
  pred13 <- simulate_lagged_dose(simB_fixed, obs13$time_min, dose_13C6_mg, par12[["f_delayed"]], par12[["t_lag"]])
  r2_13C6 <- r_squared(obs13$conc_mgL, pred13)

  tibble(
    ka = par12[["ka"]], kel = par12[["kel"]], F_12C = par12[["F_12C"]], F_13C6 = par13[["F_13C6"]],
    k_release = par13[["k_release"]], f_delayed = par12[["f_delayed"]], t_lag = par12[["t_lag"]],
    r2_12C = r2_12C, r2_13C6 = r2_13C6,
    kel_at_bound = par12[["kel"]] > (BOUNDS_LAGGED_12C$kel[2] - 1e-4),
    k_release_at_bound = par13[["k_release"]] > (BOUNDS_13C6_GIVEN_LAG$k_release[2] - 1e-4),
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
R2_SINGLE_WAVE_SKIP_TWO_WAVE <- 0.95   # see the comment on skip_two_wave below

fit_one <- function(i) {
  sid <- subject_visits$subject_id[i]; vis <- subject_visits$visit[i]
  set.seed(string_seed(paste(sid, vis)))   # reproducible regardless of N_CORES or row order - see pk_fit.R
  cat(sid, vis, "(single_wave)\n")
  sw <- fit_subject_visit_single_wave(sid, vis)

  # If single_wave already fits 12C very well, don't try two_wave at all -
  # not a computational shortcut (AIC on ER06 FCT2 showed a >0.95 R2 single
  # fit CAN still be decisively beaten by two_wave, a real ~10x RSS
  # reduction, not noise), but a deliberate choice to not trust a purely
  # statistical fit-improvement criterion over biological plausibility once
  # the data is already well explained. AIC has no concept of whether a
  # shape is physiologically real; with only 8-9 sparse points, 2 extra
  # degrees of freedom can find a "better" fit that's just exploiting
  # flexibility rather than a genuine second absorption wave, and that
  # risk is judged not worth taking once the simple model already works.
  skip_two_wave <- !is.na(sw$r2_12C) && sw$r2_12C > R2_SINGLE_WAVE_SKIP_TWO_WAVE
  if (skip_two_wave) {
    cat(sid, vis, "(two_wave skipped - single_wave R2 =", round(sw$r2_12C, 3), "already >",
        R2_SINGLE_WAVE_SKIP_TWO_WAVE, ")\n")
    tw <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                 f_delayed = NA, t_lag = NA,
                 r2_12C = NA, r2_13C6 = NA, kel_at_bound = NA, k_release_at_bound = NA,
                 converged = NA, r2_12C_low = NA, r2_13C6_low = NA, objective_value = NA,
                 rss_12C = NA, n_12C = NA, t30_present = sw$t30_present)
  } else {
    set.seed(string_seed(paste(sid, vis, "two_wave")))
    cat(sid, vis, "(two_wave)\n")
    tw <- fit_subject_visit_two_wave(sid, vis)
  }
  list(single_wave = sw, two_wave = tw)
}

fit_list <- parallel::mclapply(seq_len(nrow(subject_visits)), fit_one, mc.cores = N_CORES)
single_wave_results <- bind_cols(subject_visits, bind_rows(lapply(fit_list, `[[`, "single_wave")))
two_wave_results    <- bind_cols(subject_visits, bind_rows(lapply(fit_list, `[[`, "two_wave")))

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
  extra_flag_cols = c("r2_12C_low", "r2_13C6_low"),
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
  fit_subject_visit_two_wave(sid, vis, extra_seeds = seeds, maxit = maxit)
}
two_wave_results <- adaptive_retry(
  two_wave_results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  extra_flag_cols = c("r2_12C_low", "r2_13C6_low"),
  convergence_col = "converged",
  bounds = BOUNDS_LAGGED,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release", "f_delayed", "t_lag"),
  refit_fn = refit_flagged_two_wave,
  log_scale = c("ka", "kel", "k_release"),
  maxit_retry = 100,
  mc.cores = N_CORES
)

single_wave_results$capsule_dissolution_halflife_min <- log(2) / single_wave_results$k_release
two_wave_results$capsule_dissolution_halflife_min <- log(2) / two_wave_results$k_release

# ---------------------------------------------------------------------------
# Model selection: single_wave vs two_wave, per subject x visit, by AIC on
# 12C's own raw residuals only. The two-stage design means 13C6 never gains
# extra free parameters between the two candidates - it's always fit with
# exactly F_13C6/k_release, just conditioned on different inherited
# ka/kel/timing from 12C. All the actual complexity difference (f_delayed,
# t_lag: 2 extra parameters) lives in 12C's own stage-1 fit, so that's where
# "does the added complexity earn its keep" should be judged - comparing
# AIC on a single curve's own residuals avoids the cross-curve concentration-
# scale mixing problem a combined-curve AIC would have. 13C6 (and everything
# else) simply follows whichever structure wins for 12C.
# ---------------------------------------------------------------------------
K_12C_SINGLE_WAVE <- 3   # ka, kel, F_12C
K_12C_TWO_WAVE    <- 5   # + f_delayed, t_lag

aic <- function(rss, n, k) n * log(rss / n) + 2 * k

sw <- single_wave_results %>% select(subject_id, visit, rss_12C, n_12C) %>%
  rename(rss_sw = rss_12C, n_sw = n_12C)
tw <- two_wave_results %>% select(subject_id, visit, rss_12C, n_12C) %>%
  rename(rss_tw = rss_12C, n_tw = n_12C)
selection <- sw %>% left_join(tw, by = c("subject_id", "visit")) %>%
  mutate(
    aic_single_wave = aic(rss_sw, n_sw, K_12C_SINGLE_WAVE),
    aic_two_wave    = aic(rss_tw, n_tw, K_12C_TWO_WAVE),
    model = case_when(
      is.na(aic_single_wave) & is.na(aic_two_wave) ~ NA_character_,
      is.na(aic_two_wave)                          ~ "single_wave",
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
              "capsule_dissolution_halflife_min")) {
  results[[col]] <- pick(col)
}
# f_delayed/t_lag only exist in two_wave_results (single_wave has no such
# columns), so the join above brings them in unsuffixed, not as ".tw" -
# still need to null them out when two_wave was attempted but lost the AIC
# comparison (their values would otherwise leak through from that losing fit).
results$f_delayed <- if_else(results$model == "two_wave", results$f_delayed, NA_real_)
results$t_lag      <- if_else(results$model == "two_wave", results$t_lag, NA_real_)
# t30_present (kept from single_wave's own copy via the join above, since
# that model is always attempted regardless of whether two_wave was skipped)
# is a property of the data, not of which model won - ka is poorly anchored
# without it even under single_wave; downstream reliability filtering should
# take this into account alongside r2/converged/bound flags.
results <- results %>% select(subject_id, visit, model, ka, kel, F_12C, F_13C6, k_release,
                               f_delayed, t_lag, capsule_dissolution_halflife_min, t30_present,
                               r2_12C, r2_13C6, kel_at_bound, k_release_at_bound, converged,
                               r2_12C_low, r2_13C6_low, aic_single_wave, aic_two_wave)

cat("\n=== Model selection ===\n")
print(table(results$model, useNA = "ifany"))

dir.create("results", showWarnings = FALSE)
write_csv(single_wave_results, "results/fit_results_single_wave.csv")
write_csv(two_wave_results, "results/fit_results_two_wave.csv")
write_csv(results, "results/fit_results.csv")

# Plots: one figure per subject, both isotopes x both visits, using each
# subject x visit's SELECTED model (single_wave or two_wave, per
# results$model above), colored by that subject's diet arm, against the
# observed points (which include t=0, unlike the t>0-only data_fit used for
# fitting).
simulate_fit <- function(sid, vis, r) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
  dose_12C <- cov$dose_12C_mg
  fine <- seq(0, 400, length.out = 400)

  if (r$model == "two_wave") {
    simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_12C, dose, Vd)
    simB <- function(t, dose) simulate_delayed_release(t, r$k_release, r$ka, r$kel, r$F_13C6, dose, Vd)$conc
    sim12 <- simulate_lagged_dose(simA, fine, dose_12C, r$f_delayed, r$t_lag)
    sim13 <- simulate_lagged_dose(simB, fine, dose_13C6_mg, r$f_delayed, r$t_lag)
  } else {
    sim12 <- bateman_conc(fine, r$ka, r$kel, r$F_12C, dose_12C, Vd)
    sim13 <- simulate_delayed_release(fine, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, Vd)$conc
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

  diet_val <- covariates %>% filter(subject_id == sid) %>% pull(diet) %>% first()

  sim <- fits %>% pmap_dfr(function(...) {
    row <- tibble(...)
    bind_cols(visit = row$visit, diet = diet_val, simulate_fit(sid, row$visit, row))
  })

  ggplot(obs, aes(time_min, conc_mgL)) +
    geom_point(color = "black", size = 1.4) +
    geom_line(data = sim, aes(color = diet), linewidth = 0.8) +
    facet_grid(rows = vars(isotope), cols = vars(visit), scales = "free",
               labeller = labeller(visit = VISIT_LABELS, isotope = ISOTOPE_LABELS)) +
    scale_color_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
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
