# Diet-arm summary plots - fructose independent fit
# Barbara Verhaar
#
# Counterpart to 03_diet_summary.R for the standalone fructose fit
# (02_fit_erie_model_fructose_independent.R): same statistics and plot choices,
# with every 13C6-specific piece (F_13C6, isotope faceting, dose_13C6_mg,
# f_delayed_13C6) dropped, since that fit has no 13C6 curve at all.

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
  library(ggpubr)
})

# Get functions
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_stats.R")

# Theme
theme_Publication <- function(base_size=14, base_family="sans") {
    suppressWarnings(theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                plot.subtitle = element_text(size = rel(0.6), hjust = 0.5),
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
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))

}

# Plot labels: readable names and consistent colors.
VISIT_LABELS <- c(baseline = "Baseline (FCT1)", intervention = "Intervention (FCT2)")
DIET_LABELS  <- c(low_fructose = "Low fructose diet", high_fructose = "High fructose diet")
DIET_COLORS  <- c(low_fructose = "#1b9e77", high_fructose = "#d95f02")
VISIT_COLORS <- c(baseline = "#4477AA", intervention = "#CC6677")
FINE_T_DIET <- seq(0, 400, by = 2)
# Inclusion cutoff for this script's plots/statistics - separate from, and
# stricter than, 02_fit_erie_model_fructose_independent.R's R2_RELIABLE_MIN.
R2_INCLUDE_MIN <- 0.9

# Open data
results    <- read_csv("results/fructose_independent/fit_results.csv", show_col_types = FALSE)
covariates <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)

# All of this script's own outputs go in one subfolder, alongside
# 02_fit_erie_model_fructose_independent.R's plots_individual/, rather than loose
# in results/fructose_independent/.
OUT_DIR <- "results/fructose_independent/diet_summary"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Subjects who completed both visits (a body weight at each). ER33 and ER34 were
# lost to follow-up after baseline: they stay in the per-subject CSVs, flagged
# `completed = FALSE`, and are left out of every plot and statistic.
completers <- covariates %>% filter(!is.na(bw_kg)) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)

fits <- results %>%
  left_join(covariates %>% select(subject_id, visit, bw_kg, vd_L, sex, age_years, diet), by = c("subject_id", "visit")) %>%
  filter(!is.na(ka))

fits <- fits %>%
  filter(!is.na(diet)) %>%
  mutate(
    Vd = vd_L,
    completed = subject_id %in% completers,
    dose_mg = 1000 * bw_kg,
    # Reliable = the fit isn't stuck at its kel bound and its R2 clears
    # R2_INCLUDE_MIN. `converged` is not required.
    reliable = !kel_at_bound & r2 >= R2_INCLUDE_MIN
  )

reliable_fits <- fits %>% filter(reliable)
# What the plots and statistics use: completers with a reliable fit at BOTH visits,
# so every summarised subject contributes a baseline/intervention pair.
paired_ids   <- reliable_fits %>% filter(completed) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)
summary_fits <- reliable_fits %>% filter(subject_id %in% paired_ids)

# ---- Volume of distribution (Vd) by diet arm -------------------------------
# Vd is the measured ECW, not a fit result, so it isn't gated on `reliable`.
# Split by visit (dodged, matching the R2 plot below) as a covariate-balance check.

vd_data <- fits %>% filter(completed) %>% distinct(subject_id, diet, visit, sex, Vd)

p_vd <- ggplot(vd_data, aes(diet, Vd, fill = visit)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, position = position_dodge(width = 0.6)) +
  # group = visit, not the default (fill x shape): without it, adding shape=sex here
  # splits the dodge into 4 slots (visit x sex) instead of matching the boxplot's 2
  # (visit), and points land outside their own box.
  geom_jitter(aes(shape = sex, group = visit), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.6),
              size = 1.8, alpha = 0.8, color = "black") +
  scale_x_discrete(labels = DIET_LABELS) +
  scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
  labs(title = "Extracellular fluid volume",
       x = NULL, y = "Vd (L)", shape = "Sex") +
  theme_Publication()

ggsave(file.path(OUT_DIR, "diet_vd_boxplot.pdf"), p_vd, width = 6, height = 5, dpi = 150)
cat("\nSaved", file.path(OUT_DIR, "diet_vd_boxplot.pdf"), "\n")

# ---- Fit quality (R2) ------------------------------------------------------
# Not gated on `reliable` (which filters on r2 itself): a whole-cohort QC view.
R2_RELIABLE_MIN <- 0.70   # matches 02_fit_erie_model_fructose_independent.R's own threshold

r2_data <- fits %>% filter(completed, !is.na(r2)) %>% transmute(subject_id, diet, visit, r2)

p_r2 <- ggplot(r2_data, aes(visit, r2, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.6),
             size = 1.2, alpha = 0.6, color = "black") +
  geom_hline(yintercept = R2_RELIABLE_MIN, linetype = "dashed", color = "grey40") +
  scale_x_discrete(labels = VISIT_LABELS) +
  scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  labs(title = "R2",
       x = NULL, y = expression(R^2)) +
  theme_Publication()

ggsave(file.path(OUT_DIR, "diet_r2_boxplot.pdf"), p_r2, width = 7, height = 5, dpi = 150)
cat("\nSaved", file.path(OUT_DIR, "diet_r2_boxplot.pdf"), "\n")

# ---- Parameter summary table (mean/SD/SEM/n per diet x visit) -------------

summarise_param <- function(param) {
  summary_fits %>%
    group_by(diet, visit) %>%
    summarise(
      mean = mean(.data[[param]], na.rm = TRUE),
      sd   = sd(.data[[param]], na.rm = TRUE),
      n    = n(),
      .groups = "drop"
    ) %>%
    mutate(sem = sd / sqrt(n), parameter = param, .before = 1)
}

param_table <- bind_rows(
  summarise_param("ka"),
  summarise_param("kel"),
  summarise_param("F")
) %>% arrange(parameter, diet, visit)

write_csv(param_table, file.path(OUT_DIR, "diet_parameter_summary.csv"))
cat("=== Parameter summary by diet x visit (reliable fits only) ===\n")
print(as.data.frame(param_table), digits = 3)

# ---- Boxplot: diet x visit distributions for each parameter ---------------
PARAM_LABELS <- c(ka = "ka (1/min)", kel = "kel (1/min)", F = "F")

box_data <- bind_rows(
  summary_fits %>% transmute(subject_id, diet, visit, parameter = "ka", value = ka),
  summary_fits %>% transmute(subject_id, diet, visit, parameter = "kel", value = kel),
  summary_fits %>% transmute(subject_id, diet, visit, parameter = "F", value = F)
)

# Paired Wilcoxon signed-rank test, baseline vs intervention within each diet
# arm (subjects with both visits reliable; needs >= 3 pairs). Signed differences
# are on the log scale (fold-change), matching the curve-derived metrics test
# further down.
param_deltas <- box_data %>%
  pivot_wider(names_from = visit, values_from = value) %>%
  filter(!is.na(baseline), !is.na(intervention)) %>%
  mutate(delta = log(intervention) - log(baseline))

param_wilcoxon <- param_deltas %>%
  group_by(parameter, diet) %>%
  summarise(n_paired = n(), median_delta = median(delta),
            p = if (n() >= 3) tryCatch(wilcox.test(delta, mu = 0)$p.value, error = function(e) NA_real_) else NA_real_,
            .groups = "drop") %>%
  mutate(label = case_when(is.na(p) ~ sprintf("n=%d", n_paired),
                           p < 0.001 ~ sprintf("p<0.001 (n=%d)", n_paired),
                           TRUE ~ sprintf("p=%.3f (n=%d)", p, n_paired)))
write_csv(param_wilcoxon, file.path(OUT_DIR, "diet_parameter_wilcoxon.csv"))
cat("\n=== Paired Wilcoxon, intervention vs baseline within each arm (fitted parameters) ===\n")
print(as.data.frame(param_wilcoxon %>% select(-label)), digits = 3)

# Random-intercept mixed model per parameter (log scale), shown next to the Wilcoxon:
# the within-arm change and the diet x visit interaction (difference between arms).
param_mixed <- box_data %>%
  group_by(parameter) %>% group_modify(~ mixed_visit_tests(.x, log_scale = TRUE)) %>% ungroup()
param_tests <- param_wilcoxon %>% select(-label) %>% rename(p_wilcoxon = p) %>%
  left_join(param_mixed %>% select(-n_subjects), by = c("parameter", "diet")) %>%
  mutate(label = test_label(n_paired, p_wilcoxon, p_lmm))
write_csv(param_tests, file.path(OUT_DIR, "diet_parameter_tests.csv"))
cat("\n=== Wilcoxon and LMM side by side (estimates are log ratios, intervention vs baseline) ===\n")
print(as.data.frame(param_tests %>% select(-label)), digits = 3)

# One PDF per parameter, faceted by diet arm, annotated with the paired
# Wilcoxon test above. Jitter is plain black (not colored by visit) so it stays
# visible against the fill, matching the curve-derived metrics boxplots below.
plot_param_boxplot <- function(param_name) {
  d <- box_data %>% filter(parameter == param_name)
  if (nrow(d) == 0) return(NULL)
  ann <- param_tests %>% filter(parameter == param_name)
  ggplot(d, aes(visit, value, fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.6, color = "black") +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
    facet_wrap(vars(diet), nrow = 1, labeller = labeller(diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    # ka spans a wide range with a long upper tail, so it is drawn on a log scale
    (if (param_name == "ka") scale_y_log10(expand = expansion(mult = c(0.05, 0.15)))
     else scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(subtitle = interaction_subtitle(ann), title = PARAM_LABELS[[param_name]], x = NULL, y = NULL) +
    theme_Publication()
}

for (param_name in c("ka", "kel", "F")) {
  p <- plot_param_boxplot(param_name)
  if (!is.null(p)) {
    out_path <- file.path(OUT_DIR, sprintf("diet_parameter_boxplot_%s.pdf", param_name))
    ggsave(out_path, p, width = 7, height = 5, dpi = 150)
    cat("\nSaved", out_path, "\n")
  }
}

# ---- Two-wave ("second peak") characteristics by diet arm -----------------
# Exploratory only.

two_wave_rate <- summary_fits %>%
  group_by(diet, visit) %>%
  summarise(n = n(), n_two_wave = sum(model == "two_wave"), .groups = "drop") %>%
  mutate(pct_two_wave = round(100 * n_two_wave / n, 1))

write_csv(two_wave_rate, file.path(OUT_DIR, "diet_two_wave_rate.csv"))
cat("\n=== two_wave selection rate by diet x visit (reliable fits only) ===\n")
print(as.data.frame(two_wave_rate), digits = 3)

two_wave_tbl <- table(summary_fits$diet, summary_fits$model == "two_wave")
if (all(dim(two_wave_tbl) == c(2, 2))) {
  fisher_p <- fisher.test(two_wave_tbl)$p.value
  cat("\nFisher's exact test, diet x two_wave selection: p =", round(fisher_p, 3), "\n")
}

two_wave_params <- reliable_fits %>% filter(model == "two_wave") %>%
  select(subject_id, diet, visit, completed, t_lag, f_delayed)
write_csv(two_wave_params, file.path(OUT_DIR, "diet_two_wave_params.csv"))
two_wave_completed <- two_wave_params %>% filter(subject_id %in% paired_ids)

if (nrow(two_wave_completed) >= 4) {
  p_two_wave <- two_wave_completed %>%
    pivot_longer(c(t_lag, f_delayed), names_to = "parameter", values_to = "value") %>%
    filter(!is.na(value)) %>%
    ggplot(aes(diet, value, fill = diet)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.08, size = 1.4, alpha = 0.7, color = "black") +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    facet_wrap(vars(parameter), scales = "free_y",
               labeller = labeller(parameter = c(t_lag = "t_lag (min)", f_delayed = "f_delayed"))) +
    scale_x_discrete(labels = DIET_LABELS) +
    scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
    labs(title = sprintf("Two-wave characteristics by diet arm (exploratory, n=%d curves)", nrow(two_wave_completed)),
         x = NULL, y = NULL) +
    theme_Publication()
  ggsave(file.path(OUT_DIR, "diet_two_wave_boxplot.pdf"), p_two_wave, width = 7, height = 4.5, dpi = 150)
  cat("\nSaved", file.path(OUT_DIR, "diet_two_wave_boxplot.pdf"), ",", file.path(OUT_DIR, "diet_two_wave_rate.csv"), ",", file.path(OUT_DIR, "diet_two_wave_params.csv"), "\n")
} else {
  cat("\nToo few two_wave fits to plot a diet comparison (n =", nrow(two_wave_completed), ")\n")
}

# ---- Average concentration-time curves per diet x visit --------------------
# Must mirror simulate_fit() in 02_fit_erie_model_fructose_independent.R exactly:
# the mechanism that won for each subject x visit determines the curve.
curve_matrix <- function(diet_val, vis, completed_only = TRUE) {
  sub <- (if (completed_only) summary_fits else reliable_fits) %>% filter(diet == diet_val, visit == vis)
  if (nrow(sub) == 0) return(NULL)

  conc_list <- sub %>% pmap(function(...) {
    r <- list(...)
    if (r$model == "two_wave") {
      simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F, dose, r$Vd)
      simulate_lagged_dose(simA, FINE_T_DIET, r$dose_mg, r$f_delayed, r$t_lag)
    } else {
      bateman_conc(FINE_T_DIET, r$ka, r$kel, r$F, r$dose_mg, r$Vd)
    }
  })
  conc_mat <- do.call(cbind, conc_list)
  colnames(conc_mat) <- sub$subject_id
  conc_mat
}

# Pointwise Harrell-Davis median (rows = time points, columns = subjects): a
# weighted average of the sorted subject values, so it changes smoothly where
# subjects swap rank.
hd_median_rows <- function(m) {
  n <- ncol(m)
  w <- diff(pbeta(0:n / n, (n + 1) / 2, (n + 1) / 2))
  drop(crossprod(apply(m, 1, sort), w))
}

BOOT_N    <- 1000
BOOT_SEED <- 1

average_curve <- function(diet_val, vis) {
  conc_mat <- curve_matrix(diet_val, vis)
  if (is.null(conc_mat)) return(NULL)
  n <- ncol(conc_mat)

  # 95% percentile bootstrap CI of the median: resample subjects (columns).
  set.seed(BOOT_SEED)
  boot <- replicate(BOOT_N, hd_median_rows(conc_mat[, sample.int(n, replace = TRUE), drop = FALSE]))
  ci   <- apply(boot, 1, quantile, probs = c(0.025, 0.975))
  sem  <- apply(conc_mat, 1, sd) / sqrt(n)
  mean_conc <- rowMeans(conc_mat)

  tibble(
    time_min    = FINE_T_DIET,
    median_conc = hd_median_rows(conc_mat),
    median_lo   = ci[1, ],
    median_hi   = ci[2, ],
    mean_conc   = mean_conc,
    mean_lo     = mean_conc - sem,
    mean_hi     = mean_conc + sem,
    n           = n
  )
}

# Faceted by diet arm with baseline/intervention overlaid by colour, one PDF per
# statistic (median with bootstrap CI band, mean with SEM band).
CURVE_STATS <- list(
  median = list(cols = c("median_conc", "median_lo", "median_hi"), title = "Median fructose concentration (95% bootstrap CI)"),
  mean   = list(cols = c("mean_conc", "mean_lo", "mean_hi"),       title = "Mean fructose concentration (SEM band)")
)

plot_diet_curves <- function(stat) {
  curves <- expand_grid(diet = c("low_fructose", "high_fructose"), visit = c("baseline", "intervention")) %>%
    pmap_dfr(function(diet, visit) {
      ac <- average_curve(diet, visit)
      if (is.null(ac)) return(NULL)
      ac %>% mutate(diet = diet, visit = visit)
    })
  if (nrow(curves) == 0) return(NULL)

  cols <- CURVE_STATS[[stat]]$cols
  ggplot(curves, aes(time_min, .data[[cols[1]]], color = visit, fill = visit)) +
    geom_ribbon(aes(ymin = .data[[cols[2]]], ymax = .data[[cols[3]]]), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(vars(diet), nrow = 1, labeller = labeller(diet = DIET_LABELS)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(title = paste0(CURVE_STATS[[stat]]$title, ": baseline vs intervention"),
         x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
}

for (stat in names(CURVE_STATS)) {
  p <- plot_diet_curves(stat)
  if (!is.null(p)) {
    out_path <- file.path(OUT_DIR, sprintf("diet_summary_curves_%s.pdf", stat))
    ggsave(out_path, p, width = 9, height = 5, dpi = 150)
    cat("\nSaved", out_path, "\n")
  }
}

# ---- Curve-derived PK metrics (AUC, Cmax, Tmax) per subject x visit -------
# Cmax and Tmax come from each subject's simulated fitted curve on FINE_T_DIET.
curve_metrics <- expand_grid(diet = c("low_fructose", "high_fructose"), visit = c("baseline", "intervention")) %>%
  pmap_dfr(function(diet, visit) {
    m <- curve_matrix(diet, visit, completed_only = FALSE)
    if (is.null(m)) return(NULL)
    f <- reliable_fits[match(paste(colnames(m), visit), paste(reliable_fits$subject_id, reliable_fits$visit)), ]
    # AUC to infinity: every mechanism delivers F * dose to the blood and kel clears it.
    tibble(
      subject_id = colnames(m), diet = diet, visit = visit, completed = f$completed,
      AUC  = f$F * f$dose_mg / (f$Vd * f$kel),
      Cmax = apply(m, 2, max),
      Tmax = FINE_T_DIET[apply(m, 2, which.max)]
    )
  })
write_csv(curve_metrics, file.path(OUT_DIR, "diet_curve_metrics.csv"))
metrics_completed <- curve_metrics %>% filter(completed)

curve_metrics_summary <- metrics_completed %>%
  pivot_longer(c(AUC, Cmax, Tmax), names_to = "metric", values_to = "value") %>%
  group_by(diet, visit, metric) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75), .groups = "drop") %>%
  arrange(metric, diet, visit)
write_csv(curve_metrics_summary, file.path(OUT_DIR, "diet_curve_metrics_summary.csv"))
cat("\n=== Curve-derived metrics: median (IQR) per diet x visit ===\n")
print(as.data.frame(curve_metrics_summary), digits = 3)

# Paired Wilcoxon signed-rank test, baseline vs intervention within each arm
# (subjects with both visits; needs >= 3 pairs). The signed differences are on
# the log scale for AUC and Cmax and in minutes for Tmax, and no multiplicity
# correction is applied.
metric_deltas <- metrics_completed %>%
  pivot_longer(c(AUC, Cmax, Tmax), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = visit, values_from = value) %>%
  filter(!is.na(baseline), !is.na(intervention)) %>%
  mutate(delta = if_else(metric == "Tmax", intervention - baseline, log(intervention / baseline)))

metric_wilcoxon <- metric_deltas %>%
  group_by(metric, diet) %>%
  summarise(n_paired = n(), median_delta = median(delta),
            p = if (n() >= 3) tryCatch(wilcox.test(delta, mu = 0, exact = if (first(metric) == "Tmax") FALSE else NULL)$p.value,
                                       error = function(e) NA_real_) else NA_real_,
            .groups = "drop") %>%
  mutate(label = case_when(is.na(p) ~ sprintf("n=%d", n_paired),
                           p < 0.001 ~ sprintf("p<0.001 (n=%d)", n_paired),
                           TRUE ~ sprintf("p=%.3f (n=%d)", p, n_paired)))
write_csv(metric_wilcoxon, file.path(OUT_DIR, "diet_curve_metrics_wilcoxon.csv"))
cat("\n=== Paired Wilcoxon, intervention vs baseline within each arm (curve-derived metrics) ===\n")
print(as.data.frame(metric_wilcoxon %>% select(-label)), digits = 3)

# Same mixed model for the curve metrics (log scale for AUC and Cmax, minutes for Tmax).
metric_mixed <- metrics_completed %>%
  pivot_longer(c(AUC, Cmax, Tmax), names_to = "metric", values_to = "value") %>%
  group_by(metric) %>% group_modify(~ mixed_visit_tests(.x, log_scale = .y$metric != "Tmax")) %>% ungroup()
metric_tests <- metric_wilcoxon %>% select(-label) %>% rename(p_wilcoxon = p) %>%
  left_join(metric_mixed %>% select(-n_subjects), by = c("metric", "diet")) %>%
  mutate(label = test_label(n_paired, p_wilcoxon, p_lmm))
write_csv(metric_tests, file.path(OUT_DIR, "diet_curve_metrics_tests.csv"))
cat("\n=== Wilcoxon and LMM side by side (curve-derived metrics) ===\n")
print(as.data.frame(metric_tests %>% select(-label)), digits = 3)

METRIC_LABELS <- c(AUC = "AUC (mg*min/L)", Cmax = "Cmax (mg/L)", Tmax = "Tmax (min)")

for (metric in names(METRIC_LABELS)) {
  ann <- metric_tests %>% filter(metric == !!metric)
  p <- ggplot(metrics_completed, aes(visit, .data[[metric]], fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.6, color = "black") +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
    facet_wrap(vars(diet), nrow = 1, scales = "free_y", labeller = labeller(diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(subtitle = interaction_subtitle(ann), title = metric, x = NULL, y = METRIC_LABELS[[metric]]) +
    theme_Publication()
  out_path <- file.path(OUT_DIR, sprintf("diet_curve_%s_boxplot.pdf", tolower(metric)))
  ggsave(out_path, p, width = 7, height = 5, dpi = 150)
  cat("\nSaved", out_path, "\n")
}

cat("\nSaved", file.path(OUT_DIR, "diet_parameter_summary.csv"), "\n")
