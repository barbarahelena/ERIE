# =============================================================================
# 05_identified_quantities.R - baseline vs intervention on quantities the data
# identify, without assigning ka and kel to a role
#
# ka and kel are interchangeable in the one-compartment oral model (see "Known
# limitations" in docs/pk-model.md), so a ka-only or kel-only comparison mixes
# real change with which rate the optimizer happened to label ka. This script
# tests instead:
#   - from the fits:        k_fast = max(ka, kel), k_slow = min(ka, kel)
#   - from observed data:   Cmax, Tmax, AUC(0-360), time to 50% of Cmax on the
#                           rise, and the terminal rate from the t >= 180 samples
# The observed-data metrics use no model. A change in the absorption rate moves
# the rise (t50_rise, Tmax); a change in elimination moves the terminal rate.
#
# Paired Wilcoxon signed-rank, intervention vs baseline, within each diet arm
# (log ratio for rates, Cmax, AUC and t50_rise; minutes for Tmax). No
# multiplicity correction. Outputs: results/03_diagnostics_ka_kel/identified_quantities/
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })

RESULTS_DIR <- "results/01_independent_fit"
OUT_DIR     <- "results/03_diagnostics_ka_kel/identified_quantities"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
R2_INCLUDE_MIN <- 0.9   # same gate as scripts/archive/fructose_independent/03_diet_summary_fructose_independent.R

covariates <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants  <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
MW <- constants$value[constants$constant == "MW_12C"]
fits <- read_csv(file.path(RESULTS_DIR, "fit_results.csv"), show_col_types = FALSE) %>%
  left_join(covariates %>% select(subject_id, visit, bw_kg, diet), by = c("subject_id", "visit")) %>%
  filter(!is.na(ka), !is.na(diet))

completers <- covariates %>% filter(!is.na(bw_kg)) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)

# ---- observed-data metrics --------------------------------------------------
obs_all <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", !is.na(conc_umol_L)) %>%
  mutate(conc = conc_umol_L * MW / 1000) %>%
  group_by(subject_id, visit) %>%
  mutate(conc = conc - conc[time_min == 0][1]) %>%      # fasted level subtracted, as in the fits
  ungroup() %>% filter(!is.na(conc))

observed_metrics <- function(d) {
  d <- d %>% arrange(time_min)
  t <- d$time_min; y <- d$conc
  i <- which.max(y)
  half <- 0.5 * y[i]
  # first crossing of 50% of Cmax on the rise, linear interpolation
  j <- which(y[seq_len(i)] >= half)[1]
  t50 <- if (is.na(j) || j == 1) NA_real_ else t[j - 1] + (half - y[j - 1]) * (t[j] - t[j - 1]) / (y[j] - y[j - 1])
  late <- d %>% filter(time_min >= 180, conc > 0)
  k_term <- if (nrow(late) >= 3) -unname(coef(lm(log(conc) ~ time_min, late))[2]) else NA_real_
  tibble(Cmax_obs = y[i], Tmax_obs = t[i], AUC_obs = sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2),
         t50_rise = t50, k_term_obs = k_term)
}
obs_metrics <- obs_all %>% filter(time_min >= 0) %>%
  group_by(subject_id, visit) %>% group_modify(~ observed_metrics(.x)) %>% ungroup()

# ---- assemble ---------------------------------------------------------------
fit_part <- fits %>%
  mutate(reliable = !kel_at_bound & r2 >= R2_INCLUDE_MIN, k_fast = pmax(ka, kel), k_slow = pmin(ka, kel)) %>%
  select(subject_id, visit, diet, reliable, ka, kel, k_fast, k_slow)
d <- fit_part %>% inner_join(obs_metrics, by = c("subject_id", "visit")) %>% filter(subject_id %in% completers)
write_csv(d, file.path(OUT_DIR, "identified_quantities_per_visit.csv"))

# Fit-derived rows need a reliable fit at both visits (the diet-summary rule);
# observed-data metrics need only the data, so they use every completer.
reliable_pairs <- d %>% filter(reliable) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)

metrics <- tribble(
  ~metric,       ~source,   ~scale,
  "ka",          "fit",     "log",
  "kel",         "fit",     "log",
  "k_fast",      "fit",     "log",
  "k_slow",      "fit",     "log",
  "Tmax_obs",    "observed", "diff",
  "t50_rise",    "observed", "log",
  "Cmax_obs",    "observed", "log",
  "AUC_obs",     "observed", "log",
  "k_term_obs",  "observed", "log"
)

res <- pmap_dfr(metrics, function(metric, source, scale) {
  dd <- if (source == "fit") d %>% filter(subject_id %in% reliable_pairs) else d
  dd %>% select(subject_id, diet, visit, value = all_of(metric)) %>%
    pivot_wider(names_from = visit, values_from = value) %>%
    filter(is.finite(baseline), is.finite(intervention), baseline > 0 | scale == "diff") %>%
    mutate(delta = if (scale == "log") log(intervention / baseline) else intervention - baseline) %>%
    group_by(diet) %>%
    summarise(metric = metric, source = source, n_paired = n(), median_delta = median(delta),
              p = if (n() >= 3) tryCatch(wilcox.test(delta, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_) else NA_real_,
              .groups = "drop")
}) %>% select(metric, source, diet, n_paired, median_delta, p) %>% arrange(match(metric, metrics$metric), diet)
write_csv(res, file.path(OUT_DIR, "identified_quantities_wilcoxon.csv"))
cat("\n=== Paired Wilcoxon, intervention vs baseline (median delta: log ratio, or minutes for Tmax_obs) ===\n")
print(as.data.frame(res), digits = 3)
