# Diet-arm summary plots - 13C6 fit with fixed rates (07_fit_13C6.R)
# Barbara Verhaar
#
# ka and kel are 12C inputs and are not tested. See docs/pk-model.md, "Output files".

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

# Get functions
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_diagnostics.R")
source("scripts/assets/pk_stats.R")
source("scripts/assets/pred_13C6.R")

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

# Plot labels: readable names and consistent colors (as in 03).
VISIT_LABELS <- c(baseline = "Baseline (FCT1)", intervention = "Intervention (FCT2)")
DIET_LABELS  <- c(low_fructose = "Low fructose diet", high_fructose = "High fructose diet")
DIET_COLORS  <- c(low_fructose = "#1b9e77", high_fructose = "#d95f02")
VISIT_COLORS <- c(baseline = "#4477AA", intervention = "#CC6677")
FINE_T_DIET <- seq(0, 400, by = 2)
R2_INCLUDE_MIN <- 0.7
CAPTION <- "13C6 fitted with ka and kel of the 12C shared-kel fit: conditional on a stable kel and on the tracer following the 12C rates"

OUT_DIR <- "results/03_13C6_fixed_rates/diet_summary"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Fits from 07, plus kel_at_bound of the 12C fit
covariates <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
shared <- read_csv("results/01_shared_kel/shared_kel_fits.csv", show_col_types = FALSE) %>% select(subject_id, kel_at_bound)
fits <- read_csv("results/03_13C6_fixed_rates/fits_13C6.csv", show_col_types = FALSE) %>%
  filter(!is.na(mechanism_13C6), !is.na(diet)) %>%
  left_join(shared, by = "subject_id") %>%
  left_join(covariates %>% select(subject_id, visit, vd_L, dose_13C6_mg), by = c("subject_id", "visit")) %>%
  # Reliable: 13C6 R2 >= R2_INCLUDE_MIN and a reliable 12C fit
  mutate(reliable = r2_13C6 >= R2_INCLUDE_MIN & r2_12C >= R2_INCLUDE_MIN & !kel_at_bound)
summary_fits <- fits %>% filter(reliable)

# Stop if an R2 recomputed from the drawn curve differs from the one 07 saved
MW_13C6 <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE) %>% filter(constant == "MW_13C6") %>% pull(value)
obs13 <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "13C6", time_min > 0, !is.na(conc_umol_L)) %>% mutate(conc_mgL = conc_umol_L * MW_13C6 / 1000)
th_of <- function(r) unlist(r[c("F_13C6", "k_release", "f_delayed_13C6", "t_lag1_13C6", "f_delayed2_13C6", "t_lag2_13C6")])
r2_check <- map_dbl(seq_len(nrow(fits)), function(j) {
  r <- fits[j, ]; o <- obs13 %>% filter(subject_id == r$subject_id, visit == r$visit) %>% arrange(time_min)
  r_squared(o$conc_mgL, pred_13C6(r$mechanism_13C6, o$time_min, th_of(r), r, r$dose_13C6_mg, r$vd_L))
})
stopifnot(max(abs(r2_check - fits$r2_13C6)) < 1e-6)

cat("Curves:", nrow(fits), "| reliable:", nrow(summary_fits),
    "| excluded for 13C6 R2:", sum(fits$r2_13C6 < R2_INCLUDE_MIN), "| for the 12C fit:", sum(fits$r2_13C6 >= R2_INCLUDE_MIN & !fits$reliable), "\n")

# ---- Fit quality (R2) ------------------------------------------------------
# Not gated on `reliable` (which filters on r2 itself): a whole-cohort QC view.
p_r2 <- ggplot(fits, aes(visit, r2_13C6, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.6),
             size = 1.2, alpha = 0.6, color = "black") +
  geom_hline(yintercept = R2_INCLUDE_MIN, linetype = "dashed", color = "grey40") +
  scale_x_discrete(labels = VISIT_LABELS) +
  scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  labs(title = "R2 (13C6)", x = NULL, y = expression(R^2)) +
  theme_Publication() + labs(caption = CAPTION)
ggsave(file.path(OUT_DIR, "diet_r2_boxplot.pdf"), p_r2, width = 7, height = 5, dpi = 150)

# ---- F_13C6 and F_13C6 / F_12C: summary, tests and boxplots ----------------
# F_ratio compares the tracer with the drink in the same visit; Vd cancels in it.
PARAM_LABELS <- c(F_13C6 = "F (13C6)", F_ratio = "F (13C6) / F (12C)")
box_data <- bind_rows(
  summary_fits %>% transmute(subject_id, diet, visit, parameter = "F_13C6", value = F_13C6),
  summary_fits %>% transmute(subject_id, diet, visit, parameter = "F_ratio", value = F_13C6 / F_12C))

param_table <- box_data %>% group_by(parameter, diet, visit) %>%
  summarise(mean = mean(value), sd = sd(value), n = n(), median = median(value), .groups = "drop") %>%
  mutate(sem = sd / sqrt(n))
write_csv(param_table, file.path(OUT_DIR, "diet_parameter_summary.csv"))
cat("\n=== Parameter summary by diet x visit (reliable fits only) ===\n")
print(as.data.frame(param_table), digits = 3)

# Paired Wilcoxon on the log ratio (>= 3 pairs) and the mixed model, as in 03
param_wilcoxon <- box_data %>%
  pivot_wider(names_from = visit, values_from = value) %>%
  filter(!is.na(baseline), !is.na(intervention)) %>%
  mutate(delta = log(intervention) - log(baseline)) %>%
  group_by(parameter, diet) %>%
  summarise(n_paired = n(), median_delta = median(delta),
            p_wilcoxon = if (n() >= 3) tryCatch(wilcox.test(delta)$p.value, error = function(e) NA_real_) else NA_real_,
            .groups = "drop")
param_mixed <- box_data %>% group_by(parameter) %>% group_modify(~ mixed_visit_tests(.x, log_scale = TRUE)) %>% ungroup()
param_tests <- param_wilcoxon %>%
  left_join(param_mixed %>% select(-n_subjects), by = c("parameter", "diet")) %>%
  mutate(label = test_label(n_paired, p_wilcoxon, p_lmm))
write_csv(param_tests %>% select(-label), file.path(OUT_DIR, "diet_parameter_tests.csv"))
cat("\n=== Wilcoxon and LMM side by side (estimates are log ratios, intervention vs baseline) ===\n")
print(as.data.frame(param_tests %>% select(-label)), digits = 3)

for (param_name in names(PARAM_LABELS)) {
  ann <- param_tests %>% filter(parameter == param_name)
  p <- ggplot(box_data %>% filter(parameter == param_name), aes(visit, value, fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.6, color = "black") +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
    facet_wrap(vars(diet), nrow = 1, labeller = labeller(diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(subtitle = interaction_subtitle(ann), title = PARAM_LABELS[[param_name]], x = NULL, y = NULL) +
    theme_Publication() + labs(caption = CAPTION)
  ggsave(file.path(OUT_DIR, sprintf("diet_parameter_boxplot_%s.pdf", param_name)), p, width = 7, height = 5, dpi = 150)
}

# ---- 13C6 mechanism by diet arm (exploratory) ------------------------------
mechanism_rate <- summary_fits %>% count(diet, visit, mechanism_13C6) %>%
  group_by(diet, visit) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
write_csv(mechanism_rate, file.path(OUT_DIR, "diet_mechanism_rate.csv"))
cat("\n=== 13C6 mechanism by diet x visit (reliable fits only) ===\n")
print(as.data.frame(mechanism_rate))
mech_tbl <- table(summary_fits$diet, summary_fits$mechanism_13C6)
if (nrow(mech_tbl) == 2 && ncol(mech_tbl) >= 2)
  cat("\nFisher's exact test, diet x 13C6 mechanism: p =", round(fisher.test(mech_tbl)$p.value, 3), "\n")

# ---- Average concentration-time curves per diet x visit --------------------
curve_matrix <- function(diet_val, vis) {
  sub <- summary_fits %>% filter(diet == diet_val, visit == vis)
  if (nrow(sub) == 0) return(NULL)
  conc_mat <- sapply(seq_len(nrow(sub)), function(j) {
    r <- sub[j, ]; pred_13C6(r$mechanism_13C6, FINE_T_DIET, th_of(r), r, r$dose_13C6_mg, r$vd_L)
  })
  colnames(conc_mat) <- sub$subject_id
  conc_mat
}

# Pointwise Harrell-Davis median with a 95% percentile bootstrap CI over subjects (as in 03).
hd_median_rows <- function(m) {
  n <- ncol(m)
  w <- diff(pbeta(0:n / n, (n + 1) / 2, (n + 1) / 2))
  drop(crossprod(apply(m, 1, sort), w))
}
BOOT_N <- 1000; BOOT_SEED <- 1

average_curve <- function(conc_mat) {
  n <- ncol(conc_mat)
  set.seed(BOOT_SEED)
  boot <- replicate(BOOT_N, hd_median_rows(conc_mat[, sample.int(n, replace = TRUE), drop = FALSE]))
  ci <- apply(boot, 1, quantile, probs = c(0.025, 0.975))
  mean_conc <- rowMeans(conc_mat); sem <- apply(conc_mat, 1, sd) / sqrt(n)
  tibble(time_min = FINE_T_DIET, median_conc = hd_median_rows(conc_mat), median_lo = ci[1, ], median_hi = ci[2, ],
         mean_conc = mean_conc, mean_lo = mean_conc - sem, mean_hi = mean_conc + sem, n = n)
}

mats <- expand_grid(diet = c("low_fructose", "high_fructose"), visit = c("baseline", "intervention")) %>%
  mutate(mat = map2(diet, visit, curve_matrix))
curves <- mats %>% filter(!map_lgl(mat, is.null)) %>%
  mutate(curve = map(mat, average_curve)) %>% select(diet, visit, curve) %>% unnest(curve)

CURVE_STATS <- list(
  median = list(cols = c("median_conc", "median_lo", "median_hi"), title = "Median 13C6 concentration (95% bootstrap CI)"),
  mean   = list(cols = c("mean_conc", "mean_lo", "mean_hi"),       title = "Mean 13C6 concentration (SEM band)"))
for (stat in names(CURVE_STATS)) {
  cols <- CURVE_STATS[[stat]]$cols
  p <- ggplot(curves, aes(time_min, .data[[cols[1]]], color = visit, fill = visit)) +
    geom_ribbon(aes(ymin = .data[[cols[2]]], ymax = .data[[cols[3]]]), alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.9) +
    facet_wrap(vars(diet), nrow = 1, labeller = labeller(diet = DIET_LABELS)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(title = paste0(CURVE_STATS[[stat]]$title, ": baseline vs intervention"), x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication() + labs(caption = CAPTION)
  ggsave(file.path(OUT_DIR, sprintf("diet_summary_curves_%s.pdf", stat)), p, width = 9, height = 5, dpi = 150)
}

# ---- Curve-derived PK metrics (AUC, Cmax, Tmax) per subject x visit -------
# AUC to infinity: every mechanism delivers F_13C6 * dose to the blood and kel clears it.
curve_metrics <- mats %>% filter(!map_lgl(mat, is.null)) %>%
  pmap_dfr(function(diet, visit, mat) {
    f <- summary_fits[match(paste(colnames(mat), visit), paste(summary_fits$subject_id, summary_fits$visit)), ]
    tibble(subject_id = colnames(mat), diet = diet, visit = visit,
           AUC = f$F_13C6 * f$dose_13C6_mg / (f$vd_L * f$kel),
           Cmax = apply(mat, 2, max), Tmax = FINE_T_DIET[apply(mat, 2, which.max)])
  })
write_csv(curve_metrics, file.path(OUT_DIR, "diet_curve_metrics.csv"))

metrics_long <- curve_metrics %>% pivot_longer(c(AUC, Cmax, Tmax), names_to = "metric", values_to = "value")
curve_metrics_summary <- metrics_long %>% group_by(diet, visit, metric) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75), .groups = "drop") %>%
  arrange(metric, diet, visit)
write_csv(curve_metrics_summary, file.path(OUT_DIR, "diet_curve_metrics_summary.csv"))
cat("\n=== Curve-derived metrics: median (IQR) per diet x visit ===\n")
print(as.data.frame(curve_metrics_summary), digits = 3)

# Log scale for AUC and Cmax, minutes for Tmax (as in 03).
metric_wilcoxon <- metrics_long %>%
  pivot_wider(names_from = visit, values_from = value) %>%
  filter(!is.na(baseline), !is.na(intervention)) %>%
  mutate(delta = if_else(metric == "Tmax", intervention - baseline, log(intervention / baseline))) %>%
  group_by(metric, diet) %>%
  summarise(n_paired = n(), median_delta = median(delta),
            p_wilcoxon = if (n() >= 3) tryCatch(wilcox.test(delta, mu = 0, exact = if (first(metric) == "Tmax") FALSE else NULL)$p.value,
                                                error = function(e) NA_real_) else NA_real_,
            .groups = "drop")
metric_mixed <- metrics_long %>% group_by(metric) %>%
  group_modify(~ mixed_visit_tests(.x, log_scale = .y$metric != "Tmax")) %>% ungroup()
metric_tests <- metric_wilcoxon %>%
  left_join(metric_mixed %>% select(-n_subjects), by = c("metric", "diet")) %>%
  mutate(label = test_label(n_paired, p_wilcoxon, p_lmm))
write_csv(metric_tests %>% select(-label), file.path(OUT_DIR, "diet_curve_metrics_tests.csv"))
cat("\n=== Wilcoxon and LMM side by side (curve-derived metrics) ===\n")
print(as.data.frame(metric_tests %>% select(-label)), digits = 3)

METRIC_LABELS <- c(AUC = "AUC (mg*min/L)", Cmax = "Cmax (mg/L)", Tmax = "Tmax (min)")
for (metric in names(METRIC_LABELS)) {
  ann <- metric_tests %>% filter(metric == !!metric)
  p <- ggplot(curve_metrics, aes(visit, .data[[metric]], fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.08, size = 1.2, alpha = 0.6, color = "black") +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
    facet_wrap(vars(diet), nrow = 1, scales = "free_y", labeller = labeller(diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(subtitle = interaction_subtitle(ann), title = paste(metric, "(13C6)"), x = NULL, y = METRIC_LABELS[[metric]]) +
    theme_Publication() + labs(caption = CAPTION)
  ggsave(file.path(OUT_DIR, sprintf("diet_curve_%s_boxplot.pdf", tolower(metric))), p, width = 7, height = 5, dpi = 150)
}

cat("\nDone:", OUT_DIR, "\n")
