# =============================================================================
# 06_simulate_ka_kel_recovery.R - can the pipeline tell a ka change from a kel change?
#
# Simulates the design of the study (same subjects' sampling times, Vd and
# doses, both visits) from single_wave truth with known ka and kel, adds noise
# at the level of the real fits, refits every curve exactly as
# 02_fit_erie_model_fructose_independent.R does for single_wave (same
# objective, penalties and bounds; fewer random seeds), and runs the paired
# Wilcoxon on ka, kel, max(ka, kel) and min(ka, kel). Per scenario it reports:
#   - rejection rate at p < 0.05 for each quantity (power if the quantity
#     truly changed, false-positive rate if not)
#   - bias of the fitted ka and kel, and how often the labels come out swapped
#   - how often a fit lands on the ka = kel diagonal
# Scenarios (true ka / kel at baseline, in 1/min; "up" = x1.5 at intervention):
#   null_sep, ka_up_sep, kel_up_sep       0.040 / 0.025  well separated
#   null_diag, ka_up_diag, kel_up_diag    0.025 / 0.025  on the diagonal
#   ka_up_flip                            0.020 / 0.050  kel > ka
# Limitation: truth has no second wave; real curves with two waves are harder.
# Env: ERIE_SIM_REPS (default 50), ERIE_N_CORES. Output: results/diagnostics_ka_kel/simulation/
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_fit.R")

OUT_DIR <- "results/diagnostics_ka_kel/simulation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
N_REPS  <- as.integer(Sys.getenv("ERIE_SIM_REPS", unset = "50"))
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = max(1, parallel::detectCores() - 2)))
SIM_SEED <- 11

# Same settings as 02_fit_erie_model_fructose_independent.R
MIN_TMAX <- 10; CMAX_TOL <- 0.10; CMAX_LAMBDA <- 20; TMAX_LAMBDA <- 50
BOUNDS_INDEP <- list(ka = c(1e-4, 1), kel = c(0.005, 0.1), F = c(1e-4, 1))
N_RANDOM_INDEP <- 15   # the independent fit uses 40; fewer here to keep the simulation affordable
FINE_T <- seq(0, 400, length.out = 50)

# ---- design: the real completers' sampling times, Vd and dose ---------------
cov <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE) %>% filter(!is.na(bw_kg))
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", time_min > 0, !is.na(conc_umol_L))
completers <- cov %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)
design <- conc %>% filter(subject_id %in% completers) %>%
  group_by(subject_id, visit) %>% summarise(times = list(time_min), .groups = "drop") %>%
  left_join(cov %>% select(subject_id, visit, vd_L, dose_12C_mg), by = c("subject_id", "visit"))

# Noise: additive, at the level of the real single-fit residuals (the objective is unweighted)
# Taken from the better of each curve's single_wave and two_wave fit, so unmodelled second waves don't inflate it.
rss_of <- function(f) read_csv(file.path("results/fructose_independent", f), show_col_types = FALSE) %>% select(subject_id, visit, rss, n)
real <- bind_rows(rss_of("fit_results_single_wave.csv"), rss_of("fit_results_two_wave.csv")) %>%
  filter(!is.na(rss)) %>% group_by(subject_id, visit) %>% slice_min(rss, n = 1, with_ties = FALSE) %>% ungroup()
NOISE_SD <- median(sqrt(real$rss / real$n))
cat("Residual SD of the real fits (median):", round(NOISE_SD, 2), "mg/L; completers:", length(completers), "\n")

scenarios <- tribble(
  ~scenario,    ~ka,   ~kel,  ~ka_x, ~kel_x,
  "null_sep",   0.040, 0.025, 1,     1,
  "ka_up_sep",  0.040, 0.025, 1.5,   1,
  "kel_up_sep", 0.040, 0.025, 1,     1.5,
  "null_diag",  0.025, 0.025, 1,     1,
  "ka_up_diag", 0.025, 0.025, 1.5,   1,
  "kel_up_diag",0.025, 0.025, 1,     1.5,
  "ka_up_flip", 0.020, 0.050, 1.5,   1
)

fit_curve <- function(times, y, dose, Vd, seed) {
  set.seed(seed)
  curve <- list(times = times, conc = y,
    simulate      = function(theta) bateman_conc(times, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd),
    fine_simulate = function(theta) list(time = FINE_T, conc = bateman_conc(FINE_T, theta[["ka"]], theta[["kel"]], theta[["F"]], dose, Vd)))
  objective <- build_joint_objective(list(curve), min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                     tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
  seeds <- c(list(c(ka = 0.03, kel = 0.02, F = 0.3), c(ka = 0.08, kel = 0.05, F = 0.15),
                  c(ka = 0.01, kel = 0.01, F = 0.5), c(ka = 0.05, kel = 0.08, F = 0.2)),
             lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) c(ka = v, kel = v, F = 0.03)),
             random_seeds(N_RANDOM_INDEP, BOUNDS_INDEP, log_scale = c("ka", "kel")))
  fit <- fit_multistart(objective, vapply(BOUNDS_INDEP, `[`, numeric(1), 1), vapply(BOUNDS_INDEP, `[`, numeric(1), 2),
                        seeds = seeds, control = list(maxit = 40))
  if (is.null(fit)) return(c(ka = NA_real_, kel = NA_real_, F = NA_real_))
  fit$par[c("ka", "kel", "F")]
}

# One replicate of one scenario: simulate both visits for every subject, refit, summarise.
run_rep <- function(sc, rep) {
  set.seed(SIM_SEED * 100000 + rep * 100 + match(sc$scenario, scenarios$scenario))
  subj <- tibble(subject_id = completers,
                 ka_s  = sc$ka  * exp(rnorm(length(completers), 0, 0.2)),
                 kel_s = sc$kel * exp(rnorm(length(completers), 0, 0.2)),
                 cmax_s = 80 * exp(rnorm(length(completers), 0, 0.2)))
  d <- design %>% left_join(subj, by = "subject_id") %>%
    mutate(ka_true  = ka_s  * if_else(visit == "intervention", sc$ka_x, 1),
           kel_true = kel_s * if_else(visit == "intervention", sc$kel_x, 1))
  fits <- pmap_dfr(d, function(subject_id, visit, times, vd_L, dose_12C_mg, ka_s, kel_s, cmax_s, ka_true, kel_true) {
    shape <- bateman_conc(FINE_T, ka_true, kel_true, 1, dose_12C_mg, vd_L)
    F_true <- cmax_s / max(shape)
    y <- bateman_conc(times, ka_true, kel_true, F_true, dose_12C_mg, vd_L) + rnorm(length(times), 0, NOISE_SD)
    f <- fit_curve(times, y, dose_12C_mg, vd_L, seed = sample.int(1e6, 1))
    tibble(subject_id, visit, ka_true, kel_true, ka = f[["ka"]], kel = f[["kel"]])
  })
  fits <- fits %>% mutate(k_fast = pmax(ka, kel), k_slow = pmin(ka, kel),
                          k_fast_true = pmax(ka_true, kel_true), k_slow_true = pmin(ka_true, kel_true))
  pw <- function(col) {
    w <- fits %>% select(subject_id, visit, v = all_of(col)) %>% pivot_wider(names_from = visit, values_from = v) %>%
      filter(is.finite(baseline), is.finite(intervention))
    tryCatch(wilcox.test(log(w$intervention / w$baseline), mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
  }
  tibble(scenario = sc$scenario, rep = rep,
         p_ka = pw("ka"), p_kel = pw("kel"), p_k_fast = pw("k_fast"), p_k_slow = pw("k_slow"),
         bias_log_ka = median(log(fits$ka / fits$ka_true), na.rm = TRUE),
         bias_log_kel = median(log(fits$kel / fits$kel_true), na.rm = TRUE),
         swapped = mean((fits$ka < fits$kel) != (fits$ka_true < fits$kel_true) & abs(fits$ka_true - fits$kel_true) > 1e-9, na.rm = TRUE),
         on_diagonal = mean(abs(fits$ka - fits$kel) / pmax(fits$ka, fits$kel) < 1e-3, na.rm = TRUE),
         ka_above_0.1 = mean(fits$ka > 0.1, na.rm = TRUE))
}

# options(erie.sim.lib = TRUE) before source() loads the design and helpers without running the simulation
if (!isTRUE(getOption("erie.sim.lib"))) {
jobs <- expand_grid(s = seq_len(nrow(scenarios)), rep = seq_len(N_REPS))
cat("Running", nrow(jobs), "replicates on", N_CORES, "cores (", nrow(design), "curves each )\n")
t0 <- Sys.time()
out <- bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(i) run_rep(scenarios[jobs$s[i], ], jobs$rep[i]),
                                    mc.cores = N_CORES))
cat("Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
write_csv(out, file.path(OUT_DIR, "simulation_replicates.csv"))

summ <- out %>% group_by(scenario) %>%
  summarise(reps = n(),
            rej_ka = mean(p_ka < 0.05, na.rm = TRUE), rej_kel = mean(p_kel < 0.05, na.rm = TRUE),
            rej_k_fast = mean(p_k_fast < 0.05, na.rm = TRUE), rej_k_slow = mean(p_k_slow < 0.05, na.rm = TRUE),
            bias_ka = median(bias_log_ka), bias_kel = median(bias_log_kel),
            swapped = mean(swapped), on_diagonal = mean(on_diagonal), ka_above_0.1 = mean(ka_above_0.1),
            .groups = "drop") %>% arrange(match(scenario, scenarios$scenario))
write_csv(summ, file.path(OUT_DIR, "simulation_summary.csv"))
cat("\n=== Rejection rates (p < 0.05), bias (median log fitted/true), label swaps ===\n")
print(as.data.frame(summ), digits = 2)
}
