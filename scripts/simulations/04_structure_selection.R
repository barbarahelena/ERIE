# =============================================================================
# simulations/04_structure_selection.R - simulated study with single_wave and two_wave curves, fitted as in
# 02_fit_shared_kel.R (structure chosen per visit, then one kel per subject), comparing rules for choosing
# the structure. Each visit is fitted alone with single_wave and, when t=30 exists, two_wave; the rules
# differ only in which of the two they keep:
#   current           02 (gated): two_wave tried when the dip check (rise relative to the neighbours' line
#                     >= 15%), the plateau check or single_wave R2 < 0.50 fires; the lower AIC decides, and a
#                     rise >= 90% selects two_wave without the AIC comparison
#   tail_limit        the same, with the dip check limited to points up to 180 min (t_lag is capped at 150) and
#                     to rises of at least 5% of the curve's peak
#   tail_limit_r2_90  tail_limit, and two_wave also tried when single_wave R2 < 0.90
#   tail_limit_r2_95  tail_limit, and two_wave also tried when single_wave R2 < 0.95
#   always_try        two_wave tried whenever t=30 exists, the lower AIC decides, no override
#   true_structure    the true structure (the best any rule can do)
# Truth: the 33 completers' sampling times, Vd and doses, 20% between-subject variation in ka, kel and Cmax,
# two_wave with probability p_two_wave, (f_delayed, t_lag) drawn from 02's two_wave fits, additive noise at
# the median residual SD of 02, and the same kel at both visits.
# Env: ERIE_SIM_REPS (default 20), ERIE_N_CORES. Output: results/02_diagnostics_ka_kel/structure_selection/
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
source("scripts/assets/model_12C.R")   # the model code of 02

OUT_DIR <- "results/02_diagnostics_ka_kel/structure_selection"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
N_REPS  <- as.integer(Sys.getenv("ERIE_SIM_REPS", unset = "20"))
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = max(1, parallel::detectCores() - 2)))
SIM_SEED <- 17
METHODS <- c("current", "tail_limit", "tail_limit_r2_90", "tail_limit_r2_95", "always_try", "true_structure")

# ---- design and truth from the data and from 02 -----------------------------
cov <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE) %>% filter(!is.na(bw_kg))
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", time_min > 0, !is.na(conc_umol_L))
completers <- cov %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)
design <- conc %>% filter(subject_id %in% completers) %>% arrange(time_min) %>%
  group_by(subject_id, visit) %>% summarise(times = list(time_min), .groups = "drop") %>%
  left_join(cov %>% select(subject_id, visit, vd_L, dose_12C_mg), by = c("subject_id", "visit"))
visit_fits <- read_csv("results/01_shared_kel/visit_fits.csv", show_col_types = FALSE)
shared <- read_csv("results/01_shared_kel/shared_kel_fits.csv", show_col_types = FALSE)
NOISE_SD <- median(sqrt(c(shared$rss_b / shared$n_b, shared$rss_i / shared$n_i)), na.rm = TRUE)
P_TWO_WAVE <- mean(visit_fits$model == "two_wave")
WAVES <- visit_fits %>% filter(model == "two_wave") %>% select(f_delayed, t_lag)
scenarios <- tribble(
  ~scenario,         ~ka,   ~kel,  ~ka_x, ~p_two_wave,
  "null",            0.040, 0.025, 1,     P_TWO_WAVE,
  "ka_up",           0.040, 0.025, 1.5,   P_TWO_WAVE,
  "null_fewer_waves", 0.040, 0.025, 1,    0.35
)
cat("Noise SD", round(NOISE_SD, 2), "mg/L | two_wave share in 02", round(P_TWO_WAVE, 2), "|", nrow(WAVES), "(f_delayed, t_lag) pairs\n")

# The dip check of peak_dip_rise_info(): the largest rise of a point above the line between its neighbours,
# relative to that line, optionally only for points up to max_time and rising at least min_peak_frac of the peak
dip_rise <- function(t, c, max_time = Inf, min_peak_frac = -Inf) {
  n <- length(c); peak <- max(c)
  first <- which(diff(c) < 0 & c[-n] >= 0.5 * peak)
  if (n < 4 || length(first) == 0 || first[1] >= n - 2) return(-Inf)
  rise <- -Inf
  for (i in (first[1] + 2):(n - 1)) {
    if (t[i] > max_time) next
    line <- c[i - 1] + (c[i + 1] - c[i - 1]) * (t[i] - t[i - 1]) / (t[i + 1] - t[i - 1])
    if (!is.finite(line) || line <= 0 || (c[i] - line) / peak < min_peak_frac) next
    rise <- max(rise, (c[i] - line) / line)
  }
  rise
}

# Which structure each rule keeps for one visit, given both fits
choose <- function(v, sw, tw) {
  if (is.null(tw)) return(setNames(rep(FALSE, length(METHODS) - 1), METHODS[-length(METHODS)]))
  aic_two_better <- aic_of(tw) < aic_of(sw)
  plateau <- has_near_peak_neighbor(v$times, v$conc); r2_one <- sw$r2[[1]]
  rule <- function(rise, r2_try) (rise >= EXCESS_DEFINITE_TWO_WAVE) || ((rise >= 0.15 || plateau || r2_one < r2_try) && aic_two_better)
  tail <- dip_rise(v$times, v$conc, max_time = 180, min_peak_frac = 0.05)
  c(current          = rule(dip_rise(v$times, v$conc), 0.50),
    tail_limit       = rule(tail, 0.50),
    tail_limit_r2_90 = rule(tail, 0.90),
    tail_limit_r2_95 = rule(tail, 0.95),
    always_try       = aic_two_better)
}

sim_subject <- function(sc, rep, sid) {
  key <- paste(sc$scenario, rep, sid)
  set.seed(string_seed(paste("truth", SIM_SEED, key)))
  ka_s <- sc$ka * exp(rnorm(1, 0, 0.2)); kel_s <- sc$kel * exp(rnorm(1, 0, 0.2)); cmax_s <- 80 * exp(rnorm(1, 0, 0.2))
  truth <- map(c(b = "baseline", i = "intervention"), function(vis) {
    dd <- design %>% filter(subject_id == sid, visit == vis)
    two <- runif(1) < sc$p_two_wave; w <- WAVES[sample.int(nrow(WAVES), 1), ]
    th <- c(ka_x = ka_s * if (vis == "intervention") sc$ka_x else 1, kel = kel_s, F_x = 1,
            f_delayed_x = w$f_delayed, t_lag_x = w$t_lag)
    v <- list(Vd = dd$vd_L, dose = dd$dose_12C_mg); times <- dd$times[[1]]
    th[["F_x"]] <- cmax_s / max(pred_visit(v, "x", two, th, seq(0, 400, by = 1)))
    list(two = two, ka = th[["ka_x"]], f_delayed = w$f_delayed, t_lag = w$t_lag,
         v = make_visit(times, pred_visit(v, "x", two, th, times) + rnorm(length(times), 0, NOISE_SD), dd$dose_12C_mg, dd$vd_L))
  })
  vs <- map(truth, "v")
  # Stage 1 for both structures, with the seeds of fit_alone()
  alone <- imap(vs, function(v, s) {
    run <- function(two, nested = list()) {
      set.seed(string_seed(paste(key, s, two, "alone")))
      fit_visits(vs[s], setNames(list(two), s), make_seeds(vs[s], setNames(list(two), s), nested))
    }
    sw <- run(FALSE)
    tw <- if (v$t30_present && !is.null(sw)) run(TRUE, list(sw$par)) else NULL
    list(sw = sw, tw = tw)
  })
  if (any(map_lgl(alone, ~ is.null(.x$sw)))) return(NULL)
  chosen <- imap(alone, function(a, s) c(choose(vs[[s]], a$sw, a$tw), true_structure = truth[[s]]$two && !is.null(a$tw)))
  # Stage 2 once per distinct structure
  stage2 <- list()
  bind_rows(lapply(METHODS, function(m) {
    w <- map(chosen, m); combo <- paste(unlist(w), collapse = "_")
    if (is.null(stage2[[combo]])) stage2[[combo]] <<- fit_shared(key, vs, imap(alone, ~ if (w[[.y]]) .x$tw else .x$sw))
    f <- stage2[[combo]]
    if (is.null(f)) return(NULL)
    tibble(scenario = sc$scenario, rep = rep, subject_id = sid, method = m,
           true_two_b = truth$b$two, true_two_i = truth$i$two, two_b = w$b, two_i = w$i,
           f_delayed_true_b = truth$b$f_delayed, t_lag_true_b = truth$b$t_lag,
           f_delayed_true_i = truth$i$f_delayed, t_lag_true_i = truth$i$t_lag,
           ka_b_true = truth$b$ka, ka_i_true = truth$i$ka, kel_true = kel_s,
           ka_b = f$par[["ka_b"]], ka_i = f$par[["ka_i"]], kel = f$par[["kel"]])
  }))
}

jobs <- expand_grid(s = seq_len(nrow(scenarios)), rep = seq_len(N_REPS), sid = completers)
cat("Running", nrow(jobs), "subject simulations (", N_REPS, "replicates x", nrow(scenarios), "scenarios ) on", N_CORES, "cores\n")
t0 <- Sys.time()
out <- bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(j) sim_subject(scenarios[jobs$s[j], ], jobs$rep[j], jobs$sid[j]),
                                    mc.cores = N_CORES, mc.preschedule = FALSE))
cat("Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
write_csv(out, file.path(OUT_DIR, "structure_rules_subjects.csv"))

# ---- summaries ----------------------------------------------------------------
out <- out %>% mutate(err = log(ka_i / ka_b) - log(ka_i_true / ka_b_true))
curves <- out %>% select(scenario, method, rep, subject_id, true_two_b, true_two_i, two_b, two_i) %>%
  pivot_longer(c(true_two_b, true_two_i, two_b, two_i), names_to = c(".value", "visit"), names_pattern = "(true_two|two)_(b|i)")
structure <- curves %>% group_by(scenario, method) %>%
  summarise(structure_correct = mean(two == true_two), two_wave_added = mean(two[!true_two]), two_wave_missed = mean(!two[true_two]),
            .groups = "drop")
per_rep <- out %>% group_by(scenario, method, rep) %>%
  summarise(p_ka = tryCatch(wilcox.test(log(ka_i / ka_b), mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
            bias_ka_change = median(log(ka_i / ka_b)) - median(log(ka_i_true / ka_b_true)), .groups = "drop")
summ <- out %>% group_by(scenario, method) %>%
  summarise(subjects = n(), median_abs_err = median(abs(err)), err_above_0.5 = mean(abs(err) > 0.5), err_above_1 = mean(abs(err) > 1),
            .groups = "drop") %>%
  left_join(per_rep %>% group_by(scenario, method) %>%
              summarise(reps = n(), rej_ka = mean(p_ka < 0.05, na.rm = TRUE), bias_ka_change = median(bias_ka_change), .groups = "drop"),
            by = c("scenario", "method")) %>%
  left_join(structure, by = c("scenario", "method")) %>%
  arrange(match(scenario, scenarios$scenario), match(method, METHODS))
write_csv(summ, file.path(OUT_DIR, "structure_rules_summary.csv"))
cat("\n=== Structure recovery, ka change error per subject (log scale) and the ka test across subjects ===\n")
print(as.data.frame(summ), digits = 3)
