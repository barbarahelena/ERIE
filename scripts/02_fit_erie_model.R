# Fitting the joint 12C/13C6 fructose PK model
# Barbara Verhaar
#
# Fits TWO candidate models per subject x visit and picks one:
#   - single_wave: the plain joint delayed-release model.
#   - two_wave: a lagged-second-dose extension (simulate_lagged_dose, in
#     pk_curves.R) that can represent the post-peak dip seen in a meaningful
#     minority of curves - single_wave is structurally incapable of it (a
#     one-compartment absorption curve is monotonic after its single peak
#     for any parameter values).
# Selection is by AIC on 12C's own residuals, gated by two rules that
# override a pure fit-quality comparison - two_wave is never attempted when:
#   1. 12C's t=30 sample is missing (the 0-60min window is then
#      under-identified);
#   2. 12C's own raw observations show no evidence of a second wave at all
#      - either a post-peak dip-then-rise (has_peak_dip_rise() in
#      pk_diagnostics.R) or a near-equal-height neighbor right next to the
#      peak (has_near_peak_neighbor() - two overlapping waves close enough
#      together in time can blend into a flat top or an irregular rise
#      without ever producing a visible dip). There has to be direct
#      evidence of the phenomenon in the data, not just AIC noticing after
#      the fact that a flexible model found SOME improvement.
# A THIRD rule, single_wave already fitting 12C with R2 > 0.95, was tried
# and then removed: it wasn't actually the thing blocking genuine two_wave
# improvements (several, e.g. ER30 FCT1 at R2=0.846, were already well
# below it) and was redundant with the AIC comparison itself once the
# raw-data gate above was broadened to catch plateaus and the degenerate-
# fit exploits were bounded out (the data-dependent t_lag floor, per-wave
# MIN_TMAX, ka=kel diagonal seeding) - AIC alone is now trusted to judge
# whether the extra complexity earns its keep, once evidence justifies
# trying at all. See the comments on fit_subject_visit_two_wave(), fit_one(),
# and the "Model selection" block below for the evidence behind each.
# 13C6 gets its own independent choice of how much of ITS dose rode the
# second wave (f_delayed_13C6, separate from 12C's f_delayed) once 12C's
# structure is decided - see the comment on fit_subject_visit_two_wave().
#
# Both models use an UNWEIGHTED objective (proportional_weighting = FALSE),
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

# Config
MIN_TMAX    <- 30     # min - soft floor on predicted Tmax for both curves - see file header
CMAX_TOL    <- 0.10   # +/-10% soft band around each curve's own observed Cmax
CMAX_LAMBDA <- 20     # penalty weight for the Cmax band
TMAX_LAMBDA <- 50     # penalty weight for the Tmax floor
# t_lag's effective lower bound (see BOUNDS_LAGGED comment below) depends on
# that subject's OWN observed 12C concentration at t=30: fully permissive
# (5) when it's at or below this threshold (a genuine early delay is
# plausible - and exactly where within the unsampled 0-30min window can't
# be pinned down from data anyway, so no finer-grained floor is used), else
# fully restricted (30) - a STEP, not a linear interpolation: every value
# strictly between 5 and 30 sits in that same unsampled gap regardless of
# how far above this threshold t=30 is, so a "partial" floor is just as
# unsupported as no floor at all (found on ER04 FCT1: t=30=51.4mg/L, only
# moderately above this threshold, still produced an unjustified
# flat-then-rise artifact under the old interpolated version).
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
                      # t_lag bounds: sampling is 30-min-spaced through t=180,
                      # then widens to 60 (180->240) and 120 (240->360) min
                      # gaps. Both ends of this range let the optimizer place
                      # an entire invented second wave inside a gap with no
                      # sample to confirm or rule it out:
                      #  - UPPER (150): past 180 lets it hide in the wide
                      #    180-240/240-360 gaps - fits nothing, penalized by
                      #    nothing, can silently make the fit worse than
                      #    single_wave (found on ER12 FCT2: t_lag=179,
                      #    r2_13C6 dropped from baseline's 0.95 to 0.75, a
                      #    large invented peak with zero supporting data
                      #    around t=200-220).
                      #  - LOWER (5 here, but see EARLY_LAG_OK/BAD_MGL above):
                      #    below 30 lets it hide in the first, always-
                      #    unsampled 0-30min gap - there is no data point
                      #    between t=0 and t=30 to tell a genuine early
                      #    second wave apart from the optimizer just faking
                      #    an onset delay (small t_lag + f_delayed near 1,
                      #    i.e. "almost the whole dose is late") that has no
                      #    support either. Found on ER06 FCT2: t_lag=19.7min,
                      #    f_delayed=0.999, even though 12C's own t=30 sample
                      #    there was already near its peak (~80 mg/L) -
                      #    nothing in the data suggested any delay that
                      #    early. But a flat higher lower bound would equally
                      #    wrongly forbid a GENUINE early delay when t=30
                      #    really is near 0 (a real, separately-discussed
                      #    phenomenon) - so the actual bound used in
                      #    fit_subject_visit_two_wave() is computed per
                      #    subject from that subject's own t=30 value, not a
                      #    fixed constant; 5 here is just the floor of that
                      #    computed range, used as-is only for the retry
                      #    driver's generic seed generation (bounds =
                      #    BOUNDS_LAGGED below), where seeds get clipped to
                      #    the real, subject-specific bound by
                      #    fit_multistart() regardless.
                      # Both validated genuine-dip subjects (ER01, ER09)
                      # converged to t_lag in 60-114 min, well inside this
                      # range either way, so this doesn't constrain any real
                      # case found so far.
                      f_delayed = c(0.001, 0.999), t_lag = c(5, 150),
                      # 13C6's OWN delayed fraction - see BOUNDS_13C6_GIVEN_LAG
                      # and the comment on fit_subject_visit_two_wave() below
                      # for why this is separate from 12C's f_delayed.
                      f_delayed_13C6 = c(0.001, 0.999))
# BOUNDS_LAGGED above stays as the "combined" 8-parameter view the generic
# adaptive_retry() driver operates over (bound_cols/par_cols); the actual
# two-stage fit below (fit_subject_visit_lagged()) uses these two subsets.
# t_lag's lower bound here (5) is the permissive floor - the bound actually
# passed to fit_multistart() in fit_subject_visit_two_wave() jumps to 30 for
# any subject whose own t=30 concentration is above EARLY_LAG_OK_MGL above.
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

  # ka=kel diagonal seeds: bateman_conc's ka/(ka-kel) term makes the
  # objective surface narrow/awkward right where ka and kel are close
  # ("flip-flop" kinetics, a known hard region for one-compartment models) -
  # confirmed concretely on ER06 FCT2, where the production search reported
  # ka=0.016/kel=0.051 (R2=0.945) but the true nearby optimum is ka=0.023/
  # kel=0.025 (R2=0.986, verified from multiple starts). Generic random/grid
  # seeds don't reliably land close enough to this region for L-BFGS-B to
  # find it, so it's seeded explicitly here.
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
  # ka=kel diagonal seeds - see the comment in fit_curve_independent() for
  # why this region needs explicit seeding (verified on ER06 FCT2: the true
  # nearby optimum, ka=0.023/kel=0.025, was missed by the seeds below alone).
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

  # 13C6's own onset lag: if its raw data shows evidence absorption hadn't
  # started yet (has_onset_lag_evidence() - a different phenomenon from
  # 12C's two_wave dip, a delay before 13C6's SINGLE wave starts at all),
  # refit F_13C6 plus a new t_lag1_13C6, with ka/kel FIXED at the joint
  # fit's values above - never touches 12C's own reported parameters,
  # mirroring how fit_subject_visit_two_wave()'s stage 2 fixes ka/kel
  # before fitting 13C6's own extras.
  #
  # Modeled as an INSTANT BOLUS at t_lag1_13C6 (plain bateman_conc, not
  # simulate_delayed_release/k_release) rather than a gradual post-delay
  # release: k_release consistently pinned to its own upper bound in every
  # validated case once an onset lag was already in the model (ER25 FCT1,
  # ER06 FCT2, ER25 FCT2 - 3 for 3, not a coincidence) - 30min sampling
  # can't distinguish "fast dissolution" from "instant" once the delay
  # itself already explains the flat start, so carrying k_release as a
  # third free parameter was just an unidentifiable, boundary-pinned
  # nuisance dimension with no fit-quality benefit. Compared via AIC on
  # 13C6's own RSS (K=2 either way - F_13C6/k_release without lag vs
  # F_13C6/t_lag1_13C6 with - computed inline since the shared aic() helper
  # below isn't defined yet at the point this function is actually called)
  # so it's only used when it earns its keep. See "13C6's capsule release
  # has no genuine onset lag" in docs/pk-model.md - validated standalone on
  # ER25 FCT1 (R2 0.71->0.94) and ER06 FCT2 (R2 0.64->0.95).
  t_lag1_13C6 <- NA_real_
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
      n13 <- length(curve_13C6$conc)
      rss_no_lag <- sum((curve_13C6$conc - pred13)^2)
      rss_lag <- fit_lag$value
      aic_no_lag <- n13 * log(rss_no_lag / n13) + 2 * 2
      aic_lag    <- n13 * log(rss_lag    / n13) + 2 * 2
      if (aic_lag < aic_no_lag) {
        par[["F_13C6"]] <- fit_lag$par[["F_13C6"]]
        par[["k_release"]] <- NA_real_   # instant bolus - no meaningful capsule dissolution rate for this candidate
        t_lag1_13C6 <- fit_lag$par[["t_lag1_13C6"]]
        simB_final <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], par[["F_13C6"]], dose, Vd)
        pred13 <- simulate_two_lag_dose(simB_final, curve_13C6$times, dose_13C6_mg, f_delayed = 0, t_lag1 = t_lag1_13C6, gap = 0)
        r2_13C6 <- r_squared(curve_13C6$conc, pred13)
      }
    }
  }

  tibble(
    ka = par[["ka"]], kel = par[["kel"]], F_12C = par[["F_12C"]], F_13C6 = par[["F_13C6"]],
    k_release = par[["k_release"]], t_lag1_13C6 = t_lag1_13C6,
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
#
# 13C6's OWN delayed fraction (f_delayed_13C6, stage 2) is fit independently
# of 12C's f_delayed, sharing only t_lag (the timing of gastric emptying's
# second wave, a systemic event) and ka/kel. The capsule's contents don't
# have to split across both waves in the same proportion as the much larger
# liquid 12C dose - it's plausible the capsule emptied entirely in the first
# wave (f_delayed_13C6 ~ 0) or entirely in the second (~1), and forcing it
# to inherit 12C's own fraction produced exactly that failure on ER03 FCT1:
# 12C genuinely has two waves (R2 = 0.99), but 13C6's own points show a
# single early peak decaying monotonically with nothing at the time the
# inherited second wave would place one - forcing 12C's fraction onto it
# gave r2_13C6 = 0.075. Letting it fit its own fraction lets the data say
# which wave (if not both) 13C6 actually rode.
# ---------------------------------------------------------------------------
fit_subject_visit_two_wave <- function(sid, vis, extra_seeds = list(), maxit = 80) {
  empty <- tibble(ka = NA, kel = NA, F_12C = NA, F_13C6 = NA, k_release = NA,
                   f_delayed = NA, t_lag = NA, f_delayed_13C6 = NA,
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
  # Require actual evidence of a second wave in 12C's own raw observations -
  # EITHER a post-peak deviation from a smooth decline (peak_dip_rise_info(),
  # >=15% above what linear interpolation between neighbors predicts) OR a
  # near-equal-height neighbor right next to the peak (has_near_peak_neighbor()
  # - two waves close enough together in time don't always produce a visible
  # dip; they can blend into a flat top or an irregular rise instead, which
  # peak_dip_rise_info() structurally cannot see since there's no trough to
  # find). Without this, two_wave would be tried on every curve regardless
  # of shape, leaving AIC to notice only after the fact whether there was
  # anything to explain - checked directly against the data one step
  # earlier instead. Found via direct verification (grid/multistart search
  # outside the normal bounds-checked pipeline) that this was blocking real,
  # well-bounded, non-degenerate improvements: ER02 FCT1 (a flat top - both
  # neighbors within 91-93% of the peak - R2 0.955->0.997), ER06 FCT2 (peak
  # neighbor at 86% - R2 0.981->0.9985), ER30 FCT1 (irregular pre-peak rise,
  # neighbor at 93% - R2 0.859->0.957) - none of these show a post-peak dip
  # at all, so peak_dip_rise_info() alone can never catch them.
  # dip_evidence$trigger_time (when the dip check specifically fired) is
  # kept and used to seed the search below - found on ER04 FCT1 that
  # passing this gate does NOT guarantee the optimizer's generic seeds find
  # a fit anywhere near the evidence: it landed on an unrelated, spuriously-
  # slightly-better-AIC local optimum (t_lag=21.6, an early near-total-delay
  # trick) instead of the genuine t=150 plateau that triggered the gate in
  # the first place.
  dip_evidence <- peak_dip_rise_info(obs12$time_min, obs12$conc_mgL)
  plateau_evidence <- has_near_peak_neighbor(obs12$time_min, obs12$conc_mgL)
  if (!dip_evidence$detected && !plateau_evidence) return(mutate(empty, t30_present = TRUE))

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
    total <- total + CMAX_LAMBDA * excess12^2

    # Each wave's OWN peak must also clear MIN_TMAX, not just the combined
    # curve's tallest point above - the combined check alone only ever
    # constrains whichever wave happens to be taller, leaving the other
    # one's timing completely free. Confirmed on real fitted results: ER09
    # FCT1/FCT2 (wave 1 taller - combined peak at t=38/33, but wave 2 ALONE
    # peaks at t=151/147, never checked) and ER16 FCT2 (wave 2 alone peaks
    # at t=162). Computed directly (not via simulate_lagged_dose, which
    # only returns the summed curve) so each wave's own peak can be
    # checked independently.
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
  # Seeds anchored on the ACTUAL evidence that justified trying two_wave at
  # all (dip_evidence$trigger_time, from peak_dip_rise_info() above) - the
  # generic seeds above have no reason to explore near where the real
  # evidence is, and the optimizer can land on an unrelated, spuriously-
  # better-AIC local optimum instead. Found on ER04 FCT1: without these,
  # the search converged to t_lag=21.6 (an early near-total-delay trick),
  # completely missing the genuine t=150 plateau that triggered this fit
  # being attempted in the first place. Seeded at and somewhat before the
  # trigger time, since the trigger point is where the second wave's
  # contribution becomes visibly evident, not necessarily exactly when it
  # started.
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

  # t_lag's actual lower bound for THIS subject: a STEP function of their
  # own t=30 sample, not a linear interpolation (an earlier version of this
  # bound interpolated between 5 and 30 min as t=30 rose from
  # EARLY_LAG_OK_MGL to EARLY_LAG_BAD_MGL - wrong, because there is no
  # sample ANYWHERE between t=0 and t=30, so every value strictly between 5
  # and 30 sits in exactly the same unsampled gap regardless of how close
  # t=30 is to either threshold. A "partial credit" floor is just as
  # unsupported as the original flat 5min one - confirmed on ER04 FCT1:
  # t=30 = 51.4 mg/L, only "moderately" elevated (not near
  # EARLY_LAG_BAD_MGL), so the interpolated floor came out to 21.6min - and
  # the fit landed EXACTLY there, reproducing the same unjustified
  # flat-then-rise artifact the whole mechanism exists to prevent, just
  # shifted from 5 to 21.6). Either t=30 is genuinely near zero (a delay
  # somewhere in the unsampled window is plausible, and 5 is used as the
  # permissive floor since exactly where within it can't be pinned down
  # from data anyway), or it isn't - in which case NO point in that window
  # is defensible, so the floor jumps straight to 30, not partway there.
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

  # ---- Stage 2: F_13C6, k_release, AND 13C6's own delayed fraction, from
  # 13C6 alone - ka/kel/t_lag fixed at stage 1's estimate, but f_delayed_13C6
  # is free (see the comment above this function for why it isn't inherited
  # from 12C's f_delayed).
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

  tibble(
    ka = par12[["ka"]], kel = par12[["kel"]], F_12C = par12[["F_12C"]], F_13C6 = par13[["F_13C6"]],
    k_release = par13[["k_release"]], f_delayed = par12[["f_delayed"]], t_lag = par12[["t_lag"]],
    f_delayed_13C6 = par13[["f_delayed_13C6"]],
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
# Deviation magnitude (peak_dip_rise_info()'s $excess) above which the raw
# data is treated as UNAMBIGUOUS evidence of a genuine second wave, strong
# enough to override the AIC comparison itself (see the model-selection
# block below). Set well above the 15% detection floor and above ER01's
# own 0.40 (a clearly genuine, but unremarkable-by-comparison, case) -
# chosen from the actual cohort distribution: only ER08 FCT1 (1.57) and
# ER11 FCT1 (1.01) clear this bar, both independently confirmed genuine
# (ER08: corroborated by 13C6 peaking at the identical timepoint; ER11: a
# clean, sustained dip-then-peak, not just a single-point spike).
# Deliberately conservative - most real dip cases (including ER01, ER09)
# still go through the standard AIC comparison, which already handles them
# correctly.
EXCESS_DEFINITE_TWO_WAVE <- 0.90

fit_one <- function(i) {
  sid <- subject_visits$subject_id[i]; vis <- subject_visits$visit[i]
  set.seed(string_seed(paste(sid, vis)))   # reproducible regardless of N_CORES or row order - see pk_fit.R
  cat(sid, vis, "(single_wave)\n")
  sw <- fit_subject_visit_single_wave(sid, vis)

  obs12_raw <- data_fit %>% filter(subject_id == sid, visit == vis, isotope == "12C")
  dip_evidence <- if (nrow(obs12_raw) >= 4) peak_dip_rise_info(obs12_raw$time_min, obs12_raw$conc_mgL) else list(detected = FALSE, excess = NA_real_)
  definite_two_wave <- isTRUE(dip_evidence$detected) && !is.na(dip_evidence$excess) && dip_evidence$excess >= EXCESS_DEFINITE_TWO_WAVE

  # two_wave is no longer skipped based on single_wave's own R2 - that
  # skip existed to guard against AIC "making stuff up" on curves already
  # well explained, but single_wave's R2 turned out to be the wrong signal
  # for that: several curves with single_wave R2 well below 0.95 (ER30
  # FCT1: 0.846) were STILL only blocked by the raw-data gate, not by this
  # skip, while genuinely good two_wave improvements (ER02 FCT1, ER06 FCT2)
  # were being missed entirely by the OLD, dip-only raw-data gate rather
  # than correctly caught and then fairly judged by AIC. Now that the
  # raw-data gate also catches plateaus (has_near_peak_neighbor(), see
  # fit_subject_visit_two_wave()) and the degenerate-fit exploits are
  # bounded out (the data-dependent t_lag floor, per-wave MIN_TMAX, ka=kel
  # diagonal seeding), the AIC comparison itself - not a pre-emptive R2
  # veto - is trusted to decide whether the extra complexity earns its keep.
  set.seed(string_seed(paste(sid, vis, "two_wave")))
  cat(sid, vis, if (definite_two_wave) "(two_wave - definite evidence)" else "(two_wave)", "\n")
  tw <- fit_subject_visit_two_wave(sid, vis)
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
  # r2_13C6_low deliberately excluded: a poor 13C6 fit alone is often a
  # structural mismatch (the shared ka/kel/wave-timing 13C6 inherited from
  # 12C just isn't right for its own curve), not a search that missed a
  # better local optimum - a denser search with the SAME structure won't
  # fix that. Retrying only on r2_12C_low keeps retries targeted at fits a
  # denser search can actually help.
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
  fit_subject_visit_two_wave(sid, vis, extra_seeds = seeds, maxit = maxit)
}
two_wave_results <- adaptive_retry(
  two_wave_results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  # r2_13C6_low excluded here too - same reasoning as single_wave's retry
  # above, and even more directly applicable: 13C6 now has its own
  # f_delayed_13C6 (see fit_subject_visit_two_wave()) specifically so it
  # isn't stuck with 12C's wave split, so a poor r2_13C6 here means that
  # freedom still didn't find a good shape, not that the search under-
  # explored - a denser retry of the same stage-2 optimization is unlikely
  # to help.
  extra_flag_cols = c("r2_12C_low"),
  convergence_col = "converged",
  bounds = BOUNDS_LAGGED,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release", "f_delayed", "t_lag", "f_delayed_13C6"),
  refit_fn = refit_flagged_two_wave,
  log_scale = c("ka", "kel", "k_release"),
  maxit_retry = 100,
  mc.cores = N_CORES
)

single_wave_results$capsule_dissolution_halflife_min <- log(2) / single_wave_results$k_release
two_wave_results$capsule_dissolution_halflife_min <- log(2) / two_wave_results$k_release

# ---------------------------------------------------------------------------
# Model selection: single_wave vs two_wave, per subject x visit, by AIC on
# 12C's own raw residuals only - deciding purely whether 12C's OWN curve
# earns the extra complexity (f_delayed, t_lag: 2 extra parameters in its
# stage-1 fit), independent of whatever 13C6 ends up doing. This is
# necessary, not just convenient: 13C6 now has its own independent
# f_delayed_13C6 (see fit_subject_visit_two_wave()), so it's no longer true
# that two_wave adds the same fixed complexity to 13C6 every time - basing
# selection on a combined-curve AIC would let 13C6's own fit quality leak
# into whether 12C is judged to have a second wave, which isn't a question
# about 13C6 at all. Comparing AIC on 12C's own residuals only avoids that,
# as well as the cross-curve concentration-scale mixing problem a combined
# AIC would have. 13C6 (and everything else) simply follows whichever
# structure wins for 12C; its own wave choice is decided independently,
# inside fit_subject_visit_two_wave(), once 12C's structure is already fixed.
# ---------------------------------------------------------------------------
K_12C_SINGLE_WAVE <- 3   # ka, kel, F_12C
K_12C_TWO_WAVE    <- 5   # + f_delayed, t_lag

# Plain AIC, not AICc: AICc (the standard small-sample correction,
# +2k(k+1)/(n-k-1)) was tried and reverted - at n=8 with k jumping from 3
# to 5, its correction is so aggressive (raising the complexity gap from 4
# points to 28) that it re-excluded ALL THREE of the confirmed-genuine
# plateau cases that motivated widening the raw-data gate in the first
# place (ER02 FCT1, ER06 FCT2, ER30 FCT1 - all lost to single_wave under
# AICc despite large, independently-verified R2 improvements of
# 0.04-0.10+). Overcorrecting for small n isn't the right lever here - the
# raw-data gate itself (requiring real dip-or-plateau evidence, not just
# AIC finding SOME improvement) and the bounded-out degenerate-fit exploits
# (the data-dependent t_lag floor, per-wave MIN_TMAX, ka=kel diagonal
# seeding) are what's actually supposed to be doing the overfitting-
# protection work, not an extra-harsh information criterion that throws
# out real signal along with noise.
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
      # definite_two_wave (see EXCESS_DEFINITE_TWO_WAVE above) overrides
      # the AIC comparison itself: if the raw data makes a second wave
      # unambiguous, a marginal AIC edge for single_wave (e.g. from an
      # unlucky local optimum) shouldn't be allowed to override direct
      # visual evidence - only a failed two_wave fit (aic_two_wave NA,
      # caught above) does.
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
              "capsule_dissolution_halflife_min")) {
  results[[col]] <- pick(col)
}
# f_delayed/t_lag/f_delayed_13C6 only exist in two_wave_results (single_wave
# has no such columns), so the join above brings them in unsuffixed, not as
# ".tw" - still need to null them out when two_wave was attempted but lost
# the AIC comparison (their values would otherwise leak through from that
# losing fit).
results$f_delayed <- if_else(results$model == "two_wave", results$f_delayed, NA_real_)
results$t_lag      <- if_else(results$model == "two_wave", results$t_lag, NA_real_)
results$f_delayed_13C6 <- if_else(results$model == "two_wave", results$f_delayed_13C6, NA_real_)
# t_lag1_13C6 is the mirror image - it only exists in single_wave_results
# (not yet wired into two_wave's own 13C6 stage), so mask it out when
# two_wave wins instead.
results$t_lag1_13C6 <- if_else(results$model == "single_wave", results$t_lag1_13C6, NA_real_)
# t30_present (kept from single_wave's own copy via the join above, since
# that model is always attempted regardless of whether two_wave was skipped)
# is a property of the data, not of which model won - ka is poorly anchored
# without it even under single_wave; downstream reliability filtering should
# take this into account alongside r2/converged/bound flags.
results <- results %>% select(subject_id, visit, model, ka, kel, F_12C, F_13C6, k_release,
                               f_delayed, t_lag, f_delayed_13C6, t_lag1_13C6, capsule_dissolution_halflife_min, t30_present,
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
# results$model above), colored by isotope (12C/13C6 - no legend needed,
# already labelled by the facet rows), against the observed points (which
# include t=0, unlike the t>0-only data_fit used for fitting). Each panel is
# annotated with that curve's own R2, and 12C's panel also gets both
# candidate models' AIC (the value the model-selection choice above was
# actually based on).
simulate_fit <- function(sid, vis, r) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
  dose_12C <- cov$dose_12C_mg
  fine <- seq(0, 400, length.out = 400)

  if (r$model == "two_wave") {
    simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_12C, dose, Vd)
    simB <- function(t, dose) simulate_delayed_release(t, r$k_release, r$ka, r$kel, r$F_13C6, dose, Vd)$conc
    sim12 <- simulate_lagged_dose(simA, fine, dose_12C, r$f_delayed, r$t_lag)
    sim13 <- simulate_lagged_dose(simB, fine, dose_13C6_mg, r$f_delayed_13C6, r$t_lag)
  } else {
    sim12 <- bateman_conc(fine, r$ka, r$kel, r$F_12C, dose_12C, Vd)
    # 13C6's own onset lag (t_lag1_13C6, see fit_subject_visit_single_wave())
    # is an instant bolus at that delayed start, not simulate_delayed_release
    # - must be applied here too, or the plotted curve silently doesn't match
    # what r2_13C6 was actually computed against.
    sim13 <- if (!is.na(r$t_lag1_13C6)) {
      simB <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_13C6, dose, Vd)
      simulate_two_lag_dose(simB, fine, dose_13C6_mg, f_delayed = 0, t_lag1 = r$t_lag1_13C6, gap = 0)
    } else {
      simulate_delayed_release(fine, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, Vd)$conc
    }
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
