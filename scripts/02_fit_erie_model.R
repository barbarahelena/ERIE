# Fitting the joint 12C/13C6 fructose PK model
# Barbara Verhaar

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
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))

}

# Plot labels
ISOTOPE_LABELS <- c("12C" = "Fructose 12C", "13C6" = "Fructose 13C6")
ISOTOPE_COLORS <- c("12C" = "steelblue", "13C6" = "firebrick")
VISIT_LABELS   <- c(FCT1 = "FCT1 (before diet)", FCT2 = "FCT2 (after diet)")

# Config
# MIN_TMAX loosened from 30 to 10 min: some subjects' real Tmax looks like
# it's at or before the very first sample (t=30) - e.g. ER03 baseline 13C6
# spikes at t=30 then crashes by t=90 - and a floor of 30 was structurally
# preventing the optimizer from ever matching that, forcing a slower
# compromise fit that was wrong almost everywhere (see the R2<0 discussion
# in "Fit quality and what to trust" in docs/pk-model.md). 10 min still
# rules out the genuinely degenerate near-instant-spike solutions the floor
# was originally added for, without penalizing a real fast peak the sparse
# 30-min-resolution sampling can't itself distinguish from an even earlier one.
MIN_TMAX    <- 10     # min - soft floor on predicted Tmax for both curves
CMAX_TOL    <- 0.10   # +/-10% soft band around each curve's own observed Cmax
CMAX_LAMBDA <- 20     # penalty weight for the Cmax band
TMAX_LAMBDA <- 50     # penalty weight for the Tmax floor
# kel: Hannou et al. 2018 (t1/2 ~ 7-140 min). ka: bounded only below in
# spirit (absorption-rate differences are part of the research question).
BOUNDS_INDEP <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F = c(1e-4, 1))
BOUNDS_JOINT <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1),
                     F_12C = c(1e-4, 1), F_13C6 = c(1e-4, 1),
                     k_release = c(0.001, 1))   # capsule dissolution t1/2 ~ 0.7-700 min
N_RANDOM_INDEP <- 40
N_RANDOM_JOINT <- 16
FINE_T <- seq(0, 400, length.out = 50)    # grid for Tmax/Cmax penalty checks
CLEARANCE_FRAC <- 0.01                    # "fully cleared", for plot x-axis limits only
R2_RELIABLE_MIN <- 0.70                   # below this, that curve's fit is flagged unreliable
N_CORES <- max(1, parallel::detectCores() - 2) # To parallelize loop

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
# correction; 13C6 hasn't been dosed yet), so it's excluded here - but
# data_corrected itself keeps t=0, so plot_subject_fit() below (which reads
# from data_corrected, not data_fit) still shows the observed point.
data_fit <- data_corrected %>% filter(time_min > 0)

subject_visits <- data_fit %>% distinct(subject_id, visit) %>% arrange(subject_id, visit)

# Fast-testing mode: restrict to a handful of subjects for quick local
# iteration (e.g. checking a config change) without waiting on the full
# ~68-fit run. Set via `pixi run fit-model-test` or
# ERIE_TEST_SUBJECTS="ER01,ER02,ER03" Rscript scripts/02_fit_erie_model.R -
# unset (the default) runs every subject.
test_subjects <- Sys.getenv("ERIE_TEST_SUBJECTS", unset = "")
if (nzchar(test_subjects)) {
  keep <- trimws(strsplit(test_subjects, ",")[[1]])
  subject_visits <- subject_visits %>% filter(subject_id %in% keep)
  cat("ERIE_TEST_SUBJECTS set - restricting to:", paste(keep, collapse = ", "), "\n")
}

# Fit each curve on its own (ignoring the other curve) to get a quick
# ka/kel/F estimate. These are used as starting values for the slower joint
# optimization below, instead of starting it from arbitrary guesses.
# "On its own" just means calling build_joint_objective() with a list
# containing only that one curve - the same function used for the joint
# fit, so there's one fitting engine rather than two versions to keep in sync.
fit_curve_independent <- function(obs_time, obs_conc, dose, Vd) {
  curve <- list(
    times = obs_time, conc = obs_conc,
    simulate      = function(theta) bateman_conc(obs_time, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd),
    fine_simulate = function(theta) list(time = FINE_T, conc = bateman_conc(FINE_T, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd))
  )
  objective <- build_joint_objective(list(curve), min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL,
                                      cmax_lambda = CMAX_LAMBDA, tmax_lambda = TMAX_LAMBDA)

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

# Joint fit: shared ka/kel, separate F per curve, 13C6 gets a delayed-release
# step. extra_seeds/maxit let the retry pass below reuse this function.
fit_subject_visit <- function(sid, vis, extra_seeds = list(), maxit = 40) {
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
                                      cmax_lambda = CMAX_LAMBDA, tmax_lambda = TMAX_LAMBDA)

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

# Run for all subjects x visits
fit_one <- function(i) {
  sid <- subject_visits$subject_id[i]; vis <- subject_visits$visit[i]
  set.seed(string_seed(paste(sid, vis)))   # reproducible regardless of N_CORES or row order - see pk_fit.R
  cat(sid, vis, "\n")
  fit_subject_visit(sid, vis)
}

fit_list <- parallel::mclapply(seq_len(nrow(subject_visits)), fit_one, mc.cores = N_CORES)
results <- bind_cols(subject_visits, bind_rows(fit_list))

# Adaptive retry for flagged fits (boundary/poor-fit/non-convergence) - see
# "Adaptive retry" in docs/pk-model.md for why. Driver is generic
# (scripts/assets/pk_fit.R); refit_flagged is the only ERIE-specific glue.
refit_flagged <- function(i, seeds, maxit) {
  sid <- results$subject_id[i]; vis <- results$visit[i]
  set.seed(string_seed(paste(sid, vis, "retry")))
  cat(sid, vis, "(retry)\n")   # same per-item progress visibility as fit_one() above
  fit_subject_visit(sid, vis, extra_seeds = seeds, maxit = maxit)
}

results <- adaptive_retry(
  results,
  bound_cols = c(kel_at_bound = "kel", k_release_at_bound = "k_release"),
  extra_flag_cols = c("r2_12C_low", "r2_13C6_low"),
  convergence_col = "converged",
  bounds = BOUNDS_JOINT,
  par_cols = c("ka", "kel", "F_12C", "F_13C6", "k_release"),
  refit_fn = refit_flagged,
  log_scale = c("ka", "kel", "k_release"),
  mc.cores = N_CORES
)

results$capsule_dissolution_halflife_min <- log(2) / results$k_release

dir.create("results", showWarnings = FALSE)
write_csv(results, "results/fit_results_joint.csv")

# Plots: one figure per subject, both isotopes x both visits, each curve
# extended until it drops below CLEARANCE_FRAC of its own peak
simulate_fit <- function(sid, vis, r) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
  dose_12C <- cov$dose_12C_mg

  fine <- seq(0, 400, length.out = 400)
  sim12 <- bateman_conc(fine, r$ka, r$kel, r$F_12C, dose_12C, Vd)
  sim13 <- simulate_delayed_release(fine, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, Vd)$conc

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

  # facet_grid, not facet_wrap: isotope as rows / visit as columns puts each
  # variable's label on its own margin once (facet_wrap with two variables
  # stacks both labels on every individual panel - redundant and cluttered).
  # scales = "free" still gives each ROW (isotope) its own y-scale shared
  # only across its two visit columns - the same pattern 03_diet_summary.R's
  # diet plot already uses - so the ~1000x 12C/13C6 scale difference still
  # can't squash one isotope flat against the other; only the two visits of
  # the *same* isotope (a comparable scale) ever share an axis.
  # No color legend: isotope is already fully conveyed by the facet rows
  # (each row is exactly one isotope), so a color legend for the same
  # variable would just be a second, redundant label - the color itself
  # (via ISOTOPE_COLORS) stays, only its guide/legend is suppressed.
  ggplot(obs, aes(time_min, conc_mgL, color = isotope)) +
    geom_point() +
    geom_line(data = sim, linewidth = 0.8) +
    facet_grid(rows = vars(isotope), cols = vars(visit), scales = "free",
               labeller = labeller(visit = VISIT_LABELS, isotope = ISOTOPE_LABELS)) +
    scale_color_manual(values = ISOTOPE_COLORS, guide = "none") +
    labs(title = sid, x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication() +
    # theme_Publication()'s strip.text isn't size-reduced like axis text is;
    # at full base_size the isotope labels ("Fructose 13C6") clip against
    # the narrow right-margin row strip facet_grid uses. Only override
    # here, not in the shared theme, since facet_wrap panels elsewhere
    # don't have this narrow-strip problem.
    theme(strip.text.y = element_text(size = rel(0.75)))
}

dir.create("results/plots_individual", showWarnings = FALSE, recursive = TRUE)
for (sid in unique(results$subject_id[!is.na(results$ka)])) {
  p <- plot_subject_fit(sid)
  if (!is.null(p)) ggsave(file.path("results/plots_individual", sprintf("%s_joint.pdf", sid)),
                           p, width = 9, height = 5, dpi = 120)
}
