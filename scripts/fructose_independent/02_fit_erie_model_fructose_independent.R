# Fitting the fructose PK model independently of the 13C6 tracer
# Barbara Verhaar
#
# A standalone counterpart to 02_fit_erie_model.R's joint fructose/13C6 model:
# fits fructose alone, with no shared ka/kel and no coupling to 13C6 at all.
# This reproduces the "fit either curve alone" alternative that
# docs/pk-model.md ("Origin of the model") reports as giving poor cross-curve
# kel agreement and motivated the joint model in the first place - kept here
# as a separate script so that comparison stays possible without touching the
# joint pipeline.
#
# Fits two candidate models per subject x visit and picks one by AIC on
# fructose's own residuals, exactly as the joint script does for it:
#   - single_wave: plain Bateman (one-compartment oral) model.
#   - two_wave: a lagged-second-dose extension (simulate_lagged_dose, in
#     pk_curves.R) that can represent the post-peak dip, which single_wave
#     structurally cannot.
# two_wave is only attempted when the t=30 sample exists and the raw data
# show evidence of a second wave (or single_wave fits very poorly). Both
# models use an unweighted objective (proportional_weighting = FALSE).
#
# See docs/pk-model.md for the shared reasoning (post-peak dip detection,
# fitting procedure, AIC choice) - only the 13C6-specific sections don't
# apply here.

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
VISIT_LABELS <- c(baseline = "Baseline (FCT1)", intervention = "Intervention (FCT2)")

# Config
# MIN_TMAX: soft floor on predicted Tmax (see file header). 30 min (the first
# sampling time) was tried and found too strict here: some subjects' real
# Tmax falls before the first sample (e.g. ER12 intervention spikes to 180
# mg/L at t=30 then crashes to 57 by t=60), and a 30-min floor structurally
# prevents the optimizer from ever matching that. 10 min was tried once
# before, under the joint script's now-abandoned proportionally-weighted
# objective, and found inconclusive there (see the loosen-tmax-and-test-mode
# branch and docs/pk-model.md's "Parameter weighting"). Re-tested here under
# the current unweighted objective on the full cohort: median R2 0.976 ->
# 0.980, convergence 49/68 -> 64/68, zero regressions (worst delta -2e-11,
# floating-point noise), and the improved fits' own Tmax lands well clear of
# the new floor (e.g. ER12 intervention: 17.9 min) rather than collapsing to
# a degenerate near-instant spike - see docs/pk-model.md for the full check.
MIN_TMAX    <- 10     # min - soft floor on predicted Tmax - see file header
CMAX_TOL    <- 0.10   # +/-10% soft band around the curve's own observed Cmax
CMAX_LAMBDA <- 20     # penalty weight for the Cmax band
TMAX_LAMBDA <- 50     # penalty weight for the Tmax floor
# t_lag's lower bound is a step function of the subject's own t=30 sample:
# 5 at or below this threshold, else 30 (see "The post-peak dip and the
# lagged-dose model" in docs/pk-model.md).
EARLY_LAG_OK_MGL  <- 5    # t=30 at or below this: early t_lag (as low as 5) is fine
# kel: Hannou et al. 2018 (t1/2 ~ 7-140 min). ka: bounded only below in
# spirit (absorption-rate differences are part of the research question).
BOUNDS_INDEP <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F = c(1e-4, 1))
BOUNDS_LAGGED <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F = c(1e-4, 1),
                          # t_lag [5, 150]: why these ends (sampling gaps, per-subject lower bound)
                          # is in "The post-peak dip and the lagged-dose model" in docs/pk-model.md.
                          f_delayed = c(0.001, 0.999), t_lag = c(5, 150))
N_RANDOM_INDEP <- 40
N_RANDOM_LAGGED <- 60   # pilot-validated budget in the joint script - see results/model-post-peak-dip-pilot/
FIT_SEED <- 2   # salt mixed into every per-fit string_seed() key; change it to redraw all random multi-start points
FINE_T <- seq(0, 400, length.out = 50)    # grid for Tmax/Cmax penalty checks
CLEARANCE_FRAC <- 0.01                    # "fully cleared", for plot x-axis limits only
R2_RELIABLE_MIN <- 0.70                   # below this, the fit is flagged unreliable
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = max(1, parallel::detectCores() - 2))) # To parallelize loop

# NOTE: no top-level set.seed() here - it wouldn't make the multi-start seeds
# reproducible under parallel::mclapply(). Each subject x visit fit seeds itself
# from its own ID in fit_one() (see "Fitting procedure" in docs/pk-model.md).

# Load cleaned data
concentrations <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE)
covariates     <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants      <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)

# Prepare data - fructose (12C isotope) only
MW_fructose <- constants$value[constants$constant == "MW_12C"]

data <- concentrations %>%
  filter(isotope == "12C") %>%
  left_join(covariates, by = c("subject_id", "visit")) %>%
  mutate(conc_mgL = conc_umol_L * MW_fructose / 1000) %>%
  filter(!is.na(bw_kg), !is.na(conc_mgL))

# Baseline correction - deducting the baseline fasted fructose level
baseline_fructose <- data %>% filter(time_min == 0) %>%
  select(subject_id, visit, baseline_mgL = conc_mgL)

data_corrected <- data %>%
  left_join(baseline_fructose, by = c("subject_id", "visit")) %>%
  mutate(conc_mgL = conc_mgL - baseline_mgL) %>%
  select(-baseline_mgL)

# t=0 carries no fitting information (both models predict exactly 0 there), so
# it's excluded from the fit; data_corrected keeps it so the plots still show
# the observed point.
data_fit <- data_corrected %>% filter(time_min > 0)

subject_visits <- data_fit %>% distinct(subject_id, visit) %>% arrange(subject_id, visit)

# Fast-testing mode: restrict to a handful of subjects for quick local
# iteration without waiting on the full cohort run. Set via
# ERIE_TEST_SUBJECTS="ER01,ER02,ER03" Rscript scripts/fructose_independent/02_fit_erie_model_fructose_independent.R -
# unset (the default) runs every subject.
test_subjects <- Sys.getenv("ERIE_TEST_SUBJECTS", unset = "")
if (nzchar(test_subjects)) {
  keep <- trimws(strsplit(test_subjects, ",")[[1]])
  subject_visits <- subject_visits %>% filter(subject_id %in% keep)
  cat("ERIE_TEST_SUBJECTS set - restricting to:", paste(keep, collapse = ", "), "\n")
}

# ---------------------------------------------------------------------------
# single_wave: plain Bateman fit to fructose alone
# ---------------------------------------------------------------------------
fit_subject_visit_single_wave <- function(sid, vis, extra_seeds = list(), maxit = 40) {
  empty <- tibble(ka = NA, kel = NA, F = NA, r2 = NA, kel_at_bound = NA,
                   converged = NA, r2_low = NA, objective_value = NA,
                   rss = NA, n = NA, t30_present = NA)

  obs <- data_fit %>% filter(subject_id == sid, visit == vis)
  if (nrow(obs) < 3) return(empty)

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- cov$vd_L
  dose_mg <- cov$dose_12C_mg

  curve <- list(
    times = obs$time_min, conc = obs$conc_mgL,
    simulate      = function(theta) bateman_conc(obs$time_min, theta[["ka"]], theta[["kel"]], theta[["F"]], dose_mg, Vd),
    fine_simulate = function(theta) list(time = FINE_T, conc = bateman_conc(FINE_T, theta[["ka"]], theta[["kel"]], theta[["F"]], dose_mg, Vd))
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
    random_seeds(N_RANDOM_INDEP, BOUNDS_INDEP, log_scale = c("ka", "kel")),
    extra_seeds
  )

  lower <- vapply(BOUNDS_INDEP, `[`, numeric(1), 1)
  upper <- vapply(BOUNDS_INDEP, `[`, numeric(1), 2)
  fit <- fit_multistart(objective, lower, upper, seeds = seeds, control = list(maxit = maxit))
  if (is.null(fit)) return(empty)

  par <- fit$par
  pred <- curve$simulate(par)
  r2 <- r_squared(curve$conc, pred)

  tibble(
    ka = par[["ka"]], kel = par[["kel"]], F = par[["F"]],
    r2 = r2,
    kel_at_bound = kel > (BOUNDS_INDEP$kel[2] - 1e-4),
    converged = fit$convergence == 0,   # see "Fit quality and what to trust" in docs/pk-model.md
    r2_low = r2 < R2_RELIABLE_MIN,
    objective_value = fit$value,   # for comparing against a retry; not meaningful across subjects
    rss = sum((curve$conc - pred)^2),   # raw (unweighted, native mg/L^2) - for the single_wave-vs-two_wave AIC comparison below
    n = length(curve$conc),
    t30_present = any(obs$time_min == 30)   # ka is poorly anchored without this sample - see fit_subject_visit_two_wave()
  )
}

# ---------------------------------------------------------------------------
# two_wave (lagged-second-dose) fit to fructose alone - identical to "Stage 1" of
# the joint script's fit_subject_visit_two_wave(), minus the 13C6 stage that
# follows it there.
# ---------------------------------------------------------------------------
fit_subject_visit_two_wave <- function(sid, vis, extra_seeds = list(), maxit = 80, single_wave_r2 = NA_real_) {
  empty <- tibble(ka = NA, kel = NA, F = NA, f_delayed = NA, t_lag = NA,
                   r2 = NA, kel_at_bound = NA, converged = NA, r2_low = NA,
                   objective_value = NA, rss = NA, n = NA, t30_present = NA)

  obs <- data_fit %>% filter(subject_id == sid, visit == vis)
  if (nrow(obs) < 3) return(empty)
  # No two_wave without a t=30 sample: the 0-60min window would be
  # under-identified and AIC would rubber-stamp whatever 5 sparse points fit.
  # Model selection then falls back to single_wave.
  if (!any(obs$time_min == 30)) return(mutate(empty, t30_present = FALSE))
  # Require evidence of a second wave in fructose's own raw data: a post-peak dip
  # (peak_dip_rise_info()) or a near-equal-height neighbor of the peak
  # (has_near_peak_neighbor()). dip_evidence$trigger_time is kept to seed the
  # search below. See "How the dip is detected" and "Choosing between
  # single_wave and two_wave" in docs/pk-model.md.
  dip_evidence <- peak_dip_rise_info(obs$time_min, obs$conc_mgL)
  plateau_evidence <- has_near_peak_neighbor(obs$time_min, obs$conc_mgL)
  # A third path in: single_wave fitting fructose very poorly (R2 below
  # R2_ALWAYS_TRY_TWO_WAVE) also triggers an attempt, to catch smoothly
  # accelerating rises neither detector sees (ER21 baseline, in the joint
  # script). Set well below R2_RELIABLE_MIN so it only fires when single_wave
  # doesn't fit at all.
  R2_ALWAYS_TRY_TWO_WAVE <- 0.50
  poor_single_wave <- !is.na(single_wave_r2) && single_wave_r2 < R2_ALWAYS_TRY_TWO_WAVE
  if (!dip_evidence$detected && !plateau_evidence && !poor_single_wave) return(mutate(empty, t30_present = TRUE))

  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- cov$vd_L
  dose_mg <- cov$dose_12C_mg

  objective_fn <- function(theta) {
    simA <- function(t, dose) bateman_conc(t, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd)
    pred <- tryCatch(simulate_lagged_dose(simA, obs$time_min, dose_mg, theta[["f_delayed"]], theta[["t_lag"]]), error = function(e) NULL)
    if (is.null(pred) || any(!is.finite(pred))) return(1e10)
    ss_tot12 <- sum((obs$conc_mgL - mean(obs$conc_mgL))^2)
    total <- sum((pred - obs$conc_mgL)^2) / ss_tot12

    fine <- simulate_lagged_dose(simA, FINE_T, dose_mg, theta[["f_delayed"]], theta[["t_lag"]])
    if (any(!is.finite(fine))) return(1e10)
    tmax12 <- FINE_T[which.max(fine)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax12)^2
    obs_cmax12 <- max(obs$conc_mgL)
    excess12 <- max(0, abs(max(fine) - obs_cmax12) / obs_cmax12 - CMAX_TOL)
    total <- total + CMAX_LAMBDA * excess12^2

    # Each wave's own peak must also clear MIN_TMAX, not just the combined
    # curve's, or the non-tallest wave's timing is unconstrained (see "13C6's
    # own wave choice" in docs/pk-model.md - the same reasoning applies to
    # fructose's two waves here). Computed per wave because simulate_lagged_dose()
    # only returns the summed curve.
    fine_wave1 <- simA(FINE_T, (1 - theta[["f_delayed"]]) * dose_mg)
    fine_wave2 <- simA(pmax(FINE_T - theta[["t_lag"]], 0), theta[["f_delayed"]] * dose_mg)
    if (any(!is.finite(fine_wave1)) || any(!is.finite(fine_wave2))) return(1e10)
    tmax_wave1 <- FINE_T[which.max(fine_wave1)]
    tmax_wave2 <- FINE_T[which.max(fine_wave2)]
    total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax_wave1)^2
    total + TMAX_LAMBDA * max(0, MIN_TMAX - tmax_wave2)^2
  }

  no_lag_seeds <- list(  # f_delayed near 0 should recover ~a plain Bateman fit
    c(ka = 0.02, kel = 0.02, F = 0.02, f_delayed = 0.01, t_lag = 60),
    c(ka = 0.05, kel = 0.03, F = 0.01, f_delayed = 0.01, t_lag = 90)
  )
  dip_seeds <- list(  # pilot-informed: real dip cases (ER01/ER09) converged near t_lag~90
    c(ka = 0.05, kel = 0.02, F = 0.02, f_delayed = 0.5, t_lag = 60),
    c(ka = 0.03, kel = 0.02, F = 0.02, f_delayed = 0.4, t_lag = 45),
    c(ka = 0.06, kel = 0.03, F = 0.015, f_delayed = 0.35, t_lag = 75),
    c(ka = 0.04, kel = 0.025, F = 0.02, f_delayed = 0.45, t_lag = 90)
  )
  # Seeds anchored on the evidence that justified trying two_wave
  # (dip_evidence$trigger_time), at and before it: generic seeds can land on an
  # unrelated near-total-delay optimum instead (ER04 baseline, in the joint
  # script; see "Two-stage fit" in docs/pk-model.md).
  evidence_seeds <- if (!is.na(dip_evidence$trigger_time)) {
    tt <- dip_evidence$trigger_time
    list(
      c(ka = 0.03, kel = 0.02, F = 0.02, f_delayed = 0.4, t_lag = tt),
      c(ka = 0.04, kel = 0.025, F = 0.02, f_delayed = 0.3, t_lag = max(5, tt - 20)),
      c(ka = 0.05, kel = 0.03, F = 0.015, f_delayed = 0.35, t_lag = max(5, tt - 40))
    )
  } else list()
  # ka=kel diagonal seeds - see the comment in fit_subject_visit_single_wave()
  # for why this region needs explicit seeding.
  diagonal_seeds <- lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) {
    c(ka = v, kel = v, F = 0.02, f_delayed = 0.3, t_lag = 60)
  })
  jitter_seeds <- random_seeds(N_RANDOM_LAGGED, BOUNDS_LAGGED, log_scale = c("ka", "kel"))
  seeds <- c(no_lag_seeds, dip_seeds, evidence_seeds, diagonal_seeds, jitter_seeds, extra_seeds)

  # t_lag's lower bound for this subject: a step on their own t=30 sample (5 if
  # at or below EARLY_LAG_OK_MGL, else 30). Nothing is sampled between t=0 and
  # t=30, so a partial floor is as unsupported as none (see "The post-peak dip
  # and the lagged-dose model" in docs/pk-model.md).
  obs30 <- obs$conc_mgL[obs$time_min == 30]
  t_lag_lower <- if (obs30 <= EARLY_LAG_OK_MGL) BOUNDS_LAGGED$t_lag[1] else 30

  lower <- vapply(BOUNDS_LAGGED, `[`, numeric(1), 1)
  upper <- vapply(BOUNDS_LAGGED, `[`, numeric(1), 2)
  lower[["t_lag"]] <- t_lag_lower
  fit <- fit_multistart(objective_fn, lower, upper, seeds = seeds, control = list(maxit = maxit))
  if (is.null(fit)) return(empty)
  par <- fit$par

  simA_fixed <- function(t, dose) bateman_conc(t, par[["ka"]], par[["kel"]], par[["F"]], dose, Vd)
  pred <- simulate_lagged_dose(simA_fixed, obs$time_min, dose_mg, par[["f_delayed"]], par[["t_lag"]])
  r2 <- r_squared(obs$conc_mgL, pred)

  tibble(
    ka = par[["ka"]], kel = par[["kel"]], F = par[["F"]],
    f_delayed = par[["f_delayed"]], t_lag = par[["t_lag"]],
    r2 = r2,
    kel_at_bound = par[["kel"]] > (BOUNDS_LAGGED$kel[2] - 1e-4),
    converged = fit$convergence == 0,
    r2_low = r2 < R2_RELIABLE_MIN,
    objective_value = fit$value,   # not comparable across subjects, only against this subject's own prior/retry pair
    rss = sum((obs$conc_mgL - pred)^2),   # raw (unweighted, native mg/L^2) - for the single_wave-vs-two_wave AIC comparison below
    n = length(obs$conc_mgL),
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
  set.seed(string_seed(paste(sid, vis, FIT_SEED, "fructose_indep")))   # reproducible regardless of N_CORES or row order - see pk_fit.R
  cat(sid, vis, "(single_wave)\n")
  sw <- fit_subject_visit_single_wave(sid, vis)

  obs_raw <- data_fit %>% filter(subject_id == sid, visit == vis)
  dip_evidence <- if (nrow(obs_raw) >= 4) peak_dip_rise_info(obs_raw$time_min, obs_raw$conc_mgL) else list(detected = FALSE, excess = NA_real_)
  definite_two_wave <- isTRUE(dip_evidence$detected) && !is.na(dip_evidence$excess) && dip_evidence$excess >= EXCESS_DEFINITE_TWO_WAVE

  # two_wave is not skipped on single_wave's own R2: once the raw-data gates in
  # fit_subject_visit_two_wave() pass, AIC decides (see "Choosing between
  # single_wave and two_wave" in docs/pk-model.md).
  set.seed(string_seed(paste(sid, vis, FIT_SEED, "fructose_indep two_wave")))
  cat(sid, vis, if (definite_two_wave) "(two_wave - definite evidence)" else "(two_wave)", "\n")
  tw <- fit_subject_visit_two_wave(sid, vis, single_wave_r2 = sw$r2)
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
  set.seed(string_seed(paste(sid, vis, FIT_SEED, "fructose_indep retry")))
  cat(sid, vis, "(single_wave retry)\n")
  fit_subject_visit_single_wave(sid, vis, extra_seeds = seeds, maxit = maxit)
}
single_wave_results <- adaptive_retry(
  single_wave_results,
  bound_cols = c(kel_at_bound = "kel"),
  extra_flag_cols = c("r2_low"),
  convergence_col = "converged",
  bounds = BOUNDS_INDEP,
  par_cols = c("ka", "kel", "F"),
  refit_fn = refit_flagged_single_wave,
  log_scale = c("ka", "kel"),
  mc.cores = N_CORES
)

refit_flagged_two_wave <- function(i, seeds, maxit) {
  sid <- two_wave_results$subject_id[i]; vis <- two_wave_results$visit[i]
  set.seed(string_seed(paste(sid, vis, FIT_SEED, "fructose_indep two_wave retry")))
  cat(sid, vis, "(two_wave retry)\n")
  # single_wave_results is already computed here (its retry runs first). Looked
  # up rather than threaded through adaptive_retry()'s generic signature, so a
  # fit that passed the raw-data gate via the R2_ALWAYS_TRY_TWO_WAVE override
  # isn't reset to "not attempted" on retry.
  sw_r2 <- single_wave_results$r2[single_wave_results$subject_id == sid & single_wave_results$visit == vis]
  fit_subject_visit_two_wave(sid, vis, extra_seeds = seeds, maxit = maxit,
                              single_wave_r2 = if (length(sw_r2) == 1) sw_r2 else NA_real_)
}
two_wave_results <- adaptive_retry(
  two_wave_results,
  bound_cols = c(kel_at_bound = "kel"),
  extra_flag_cols = c("r2_low"),
  convergence_col = "converged",
  bounds = BOUNDS_LAGGED,
  par_cols = c("ka", "kel", "F", "f_delayed", "t_lag"),
  refit_fn = refit_flagged_two_wave,
  log_scale = c("ka", "kel"),
  maxit_retry = 100,
  mc.cores = N_CORES
)

# ---------------------------------------------------------------------------
# Model selection: single_wave vs two_wave per subject x visit, by AIC on
# fructose's own raw residuals. See "Choosing between single_wave and two_wave" in
# docs/pk-model.md.
# ---------------------------------------------------------------------------
K_SINGLE_WAVE <- 3   # ka, kel, F
K_TWO_WAVE    <- 5   # + f_delayed, t_lag

# Plain AIC (no AICc): at n=8, AICc excludes every confirmed plateau case (see
# "Choosing between single_wave and two_wave" in docs/pk-model.md).
aic <- function(rss, n, k) n * log(rss / n) + 2 * k

sw <- single_wave_results %>% select(subject_id, visit, rss, n) %>%
  rename(rss_sw = rss, n_sw = n)
tw <- two_wave_results %>% select(subject_id, visit, rss, n) %>%
  rename(rss_tw = rss, n_tw = n)
selection <- sw %>% left_join(tw, by = c("subject_id", "visit")) %>%
  left_join(definite_two_wave_flags, by = c("subject_id", "visit")) %>%
  mutate(
    aic_single_wave = aic(rss_sw, n_sw, K_SINGLE_WAVE),
    aic_two_wave    = aic(rss_tw, n_tw, K_TWO_WAVE),
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
for (col in c("ka", "kel", "F", "r2", "kel_at_bound", "converged", "r2_low")) {
  results[[col]] <- pick(col)
}
# f_delayed/t_lag only exist in two_wave_results, so the join brings them in
# unsuffixed; null them where two_wave lost the AIC comparison so the losing
# fit's values don't leak through.
results$f_delayed <- if_else(results$model == "two_wave", results$f_delayed, NA_real_)
results$t_lag      <- if_else(results$model == "two_wave", results$t_lag, NA_real_)
# t30_present comes from single_wave's copy via the join (always attempted). It
# describes the data, not the winning model - see "Output files" in
# docs/pk-model.md.
results <- results %>% select(subject_id, visit, model, ka, kel, F,
                               f_delayed, t_lag, t30_present,
                               r2, kel_at_bound, converged,
                               r2_low, aic_single_wave, aic_two_wave)

cat("\n=== Model selection ===\n")
print(table(results$model, useNA = "ifany"))

dir.create("results/fructose_independent", showWarnings = FALSE, recursive = TRUE)
write_csv(single_wave_results, "results/fructose_independent/fit_results_single_wave.csv")
write_csv(two_wave_results, "results/fructose_independent/fit_results_two_wave.csv")
write_csv(results, "results/fructose_independent/fit_results.csv")

# Plots: one figure per subject (both visits) using each subject x visit's
# selected model, against the observed points (including t=0, unlike
# data_fit). Annotates both candidates' AIC. See "Output files" in
# docs/pk-model.md.
simulate_fit <- function(sid, vis, r) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- cov$vd_L
  dose_mg <- cov$dose_12C_mg
  fine <- seq(0, 400, length.out = 400)

  if (r$model == "two_wave") {
    simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F, dose, Vd)
    sim <- simulate_lagged_dose(simA, fine, dose_mg, r$f_delayed, r$t_lag)
  } else {
    sim <- bateman_conc(fine, r$ka, r$kel, r$F, dose_mg, Vd)
  }

  t_end <- time_to_clearance(fine, sim, CLEARANCE_FRAC)
  keep <- fine <= t_end

  tibble(time_min = fine[keep], conc_mgL = sim[keep])
}

plot_subject_fit <- function(sid) {
  obs <- data_corrected %>% filter(subject_id == sid) %>% select(visit, time_min, conc_mgL)
  fits <- results %>% filter(subject_id == sid, !is.na(ka))
  if (nrow(fits) == 0) return(NULL)

  sim <- fits %>% pmap_dfr(function(...) {
    row <- tibble(...)
    bind_cols(visit = row$visit, simulate_fit(sid, row$visit, row))
  })

  fmt_aic <- function(x) if_else(is.na(x), "NA", sprintf("%.0f", x))
  ann <- fits %>% transmute(visit,
                             label = sprintf("R2=%.3f\nAIC sw/tw=%s/%s", r2,
                                             fmt_aic(aic_single_wave), fmt_aic(aic_two_wave)))

  ggplot(obs, aes(time_min, conc_mgL)) +
    geom_point(color = "black", size = 1.4) +
    geom_line(data = sim, color = "steelblue", linewidth = 0.8) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.2, size = 2.5, lineheight = 0.9) +
    facet_grid(cols = vars(visit), scales = "free",
               labeller = labeller(visit = VISIT_LABELS)) +
    labs(title = sid, x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
}

dir.create("results/fructose_independent/plots_individual", showWarnings = FALSE, recursive = TRUE)
for (sid in unique(subject_visits$subject_id)) {
  p <- plot_subject_fit(sid)
  if (!is.null(p)) ggsave(file.path("results/fructose_independent/plots_individual", sprintf("%s_fructose_indep.pdf", sid)),
                           p, width = 6, height = 4, dpi = 120)
}

cat("\nDone. results/fructose_independent/fit_results.csv (selected model per subject x visit),",
    "results/fructose_independent/fit_results_single_wave.csv, fit_results_two_wave.csv (both candidates),",
    "results/fructose_independent/plots_individual/*.pdf\n")
