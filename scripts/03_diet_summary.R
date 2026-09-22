# Diet-arm summary plots
# Barbara Verhaar
#
# Which fits are included, the statistics used and the plot choices are
# explained in docs/diet-summary.md.

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
  library(lmerTest)
  library(ggpubr)
})

# Get functions
source("scripts/assets/pk_curves.R")

# Theme
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

# Plot labels: readable names and consistent colors.
ISOTOPE_LABELS <- c("12C" = "Fructose 12C", "13C6" = "Fructose 13C6")
VISIT_LABELS   <- c(baseline = "Baseline (FCT1)", intervention = "Intervention (FCT2)")
DIET_LABELS    <- c(low_fructose = "Low fructose diet", high_fructose = "High fructose diet")
DIET_COLORS    <- c(low_fructose = "#1b9e77", high_fructose = "#d95f02")
VISIT_COLORS   <- c(baseline = "#4477AA", intervention = "#CC6677")
FINE_T_DIET <- seq(0, 400, by = 2)
# Inclusion cutoff for this script's plots/statistics - separate from, and
# stricter than, 02_fit_erie_model.R's R2_RELIABLE_MIN (see "Which fits are
# included" in docs/diet-summary.md).
R2_INCLUDE_MIN <- 0.75

# Open data
results      <- read_csv("results/fit_results.csv", show_col_types = FALSE)
covariates   <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants    <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
dose_13C6_mg <- constants$value[constants$constant == "tracer_13C6_dose_mg"]

fits <- results %>%
  left_join(covariates %>% select(subject_id, visit, bw_kg, vd_L, sex, age_years, diet), by = c("subject_id", "visit")) %>%
  filter(!is.na(ka))

fits <- fits %>%
  filter(!is.na(diet)) %>%
  mutate(
    Vd = vd_L,
    dose_12C_mg = 1000 * bw_kg,
    # Reliable = 12C's fit isn't stuck at the shared kel bound and its R2 clears
    # R2_INCLUDE_MIN, applied to every parameter. `converged` is not
    # required. See "Which fits are included" in docs/diet-summary.md.
    reliable = !kel_at_bound & r2_12C >= R2_INCLUDE_MIN
  )

# ---- Volume of distribution (Vd) by diet arm -------------------------------
# Vd is 0.4 x the measured TBW, not a fit result, so it isn't gated on `reliable`
# and is shown once per subject (mean over visits) as a covariate-balance check. See "Outputs" in docs/diet-summary.md.

vd_data <- fits %>% group_by(subject_id, diet, sex) %>%
  summarise(Vd = mean(Vd), .groups = "drop")

p_vd <- ggplot(vd_data, aes(diet, Vd, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(shape = sex), width = 0.08, size = 1.8, alpha = 0.8) +
  scale_x_discrete(labels = DIET_LABELS) +
  scale_fill_manual(values = DIET_COLORS, guide = "none") +
  labs(title = "Extracellular fluid volume (Vd = 0.4 x measured TBW) by diet arm",
       x = NULL, y = "Vd (L)", shape = "Sex") +
  theme_Publication()

ggsave("results/diet_vd_boxplot.pdf", p_vd, width = 6, height = 5, dpi = 150)
cat("\nSaved results/diet_vd_boxplot.pdf\n")

# ---- Fit quality (R2) by diet arm and visit --------------------------------
# Not gated on `reliable` (which filters on r2_12C itself): a
# whole-cohort QC view. See "Outputs" in docs/diet-summary.md.
R2_RELIABLE_MIN <- 0.70   # matches 02_fit_erie_model.R's own threshold

r2_data <- bind_rows(
  fits %>% filter(!is.na(r2_12C)) %>% transmute(subject_id, diet, visit, isotope = "12C", r2 = r2_12C),
  fits %>% filter(!is.na(r2_13C6)) %>% transmute(subject_id, diet, visit, isotope = "13C6", r2 = r2_13C6)
)

p_r2 <- ggplot(r2_data, aes(visit, r2, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(aes(color = diet), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.6),
             size = 1.2, alpha = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = R2_RELIABLE_MIN, linetype = "dashed", color = "grey40") +
  facet_wrap(vars(isotope), nrow = 1, labeller = labeller(isotope = ISOTOPE_LABELS)) +
  scale_x_discrete(labels = VISIT_LABELS) +
  scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  scale_color_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  labs(title = "Fit quality (R2) by diet arm and visit, all fitted subjects",
       x = NULL, y = expression(R^2)) +
  theme_Publication()

ggsave("results/diet_r2_boxplot.pdf", p_r2, width = 9, height = 5, dpi = 150)
cat("\nSaved results/diet_r2_boxplot.pdf\n")

# ---- Parameter summary table (mean/SD/SEM/n per diet x visit) -------------

summarise_param <- function(param) {
  fits %>%
    filter(reliable) %>%
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
  summarise_param("F_12C"),
  summarise_param("F_13C6")
) %>% arrange(parameter, diet, visit)

dir.create("results", showWarnings = FALSE)
write_csv(param_table, "results/diet_parameter_summary.csv")
cat("=== Parameter summary by diet x visit (reliable fits only) ===\n")
print(as.data.frame(param_table), digits = 3)

# ---- Statistical comparison: diet x time linear mixed models --------------
# One LMM per parameter on reliable fits, with subject_id as a random
# intercept for the paired baseline/intervention observations. See "Linear mixed models" in
# docs/diet-summary.md.

# Modelled on the log scale, with sex, age and body weight as fixed-effect
# covariates (rationale in docs/diet-summary.md). Shared by the fitted
# parameters here and the curve-derived metrics (AUC, Cmax, Tmax) further down,
# once those are computed.
run_lmm <- function(d, label, log_scale = TRUE) {
  if (log_scale) d <- d %>% mutate(value = log(value))
  n_subjects_both <- d %>% count(subject_id) %>% filter(n == 2) %>% nrow()
  if (n_subjects_both < 3) {
    warning(label, ": fewer than 3 subjects with both visits reliable - skipping LMM")
    return(NULL)
  }
  model <- lmer(value ~ diet * visit + sex + age_years + bw_kg + (1 | subject_id), data = d)
  a <- anova(model)  # Type III, Satterthwaite df (lmerTest default)
  tibble(parameter = label, term = rownames(a), `F` = a$`F value`, df1 = a$NumDF, df2 = a$DenDF, p = a$`Pr(>F)`)
}

fit_lmm <- function(param) {
  d <- fits %>% filter(reliable) %>% select(subject_id, diet, visit, sex, age_years, bw_kg, value = all_of(param))
  run_lmm(d, param)
}

lmm_results_params <- bind_rows(
  fit_lmm("ka"),
  fit_lmm("kel"),
  fit_lmm("F_12C"),
  fit_lmm("F_13C6")
)

# ---- Boxplot: diet x visit distributions for each parameter ---------------
PARAM_LABELS <- c(ka = "ka (1/min)", kel = "kel (1/min)", F_12C = "F[12C]", F_13C6 = "F[13C6]")

box_data <- bind_rows(
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "ka", value = ka),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "kel", value = kel),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "F_12C", value = F_12C),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "F_13C6", value = F_13C6)
)

# One PDF per parameter, faceted by diet arm and colored by visit (VISIT_COLORS,
# as in the curve plots).
plot_param_boxplot <- function(param_name) {
  d <- box_data %>% filter(parameter == param_name)
  if (nrow(d) == 0) return(NULL)
  ggplot(d, aes(visit, value, fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(aes(color = visit), width = 0.08, size = 1.2, alpha = 0.6, show.legend = FALSE) +
    facet_wrap(vars(diet), nrow = 1, labeller = labeller(diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    scale_color_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(title = sprintf("%s by diet arm and visit", PARAM_LABELS[[param_name]]), x = NULL, y = NULL) +
    theme_Publication()
}

for (param_name in c("ka", "kel", "F_12C", "F_13C6")) {
  p <- plot_param_boxplot(param_name)
  if (!is.null(p)) {
    out_path <- sprintf("results/diet_parameter_boxplot_%s.pdf", param_name)
    ggsave(out_path, p, width = 7, height = 5, dpi = 150)
    cat("\nSaved", out_path, "\n")
  }
}

# ---- Delta plots: within-subject before/after diet change, by diet arm ----
# Per-subject log fold-change (log(intervention) - log(baseline)) for subjects with both
# visits reliable, compared between arms with a Wilcoxon rank-sum test. See
# "Within-subject change (delta)" in docs/diet-summary.md.

delta_param <- function(param) {
  fits %>% filter(reliable) %>%
    select(subject_id, diet, visit, value = all_of(param)) %>%
    pivot_wider(names_from = visit, values_from = value) %>%
    filter(!is.na(baseline), !is.na(intervention)) %>%
    transmute(subject_id, diet, parameter = param, delta = log(intervention) - log(baseline))
}

delta_data <- bind_rows(
  delta_param("ka"), delta_param("kel"), delta_param("F_12C"),
  delta_param("F_13C6")
)

write_csv(delta_data, "results/diet_parameter_deltas.csv")

p_delta <- ggplot(delta_data, aes(diet, delta, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(color = diet), width = 0.08, size = 1.2, alpha = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  facet_wrap(vars(parameter), scales = "free_y", nrow = 2,
             labeller = labeller(parameter = PARAM_LABELS)) +
  scale_x_discrete(labels = DIET_LABELS) +
  scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  scale_color_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  labs(title = "Within-subject diet change (log fold-change, intervention vs baseline) by diet arm",
       x = NULL, y = "log(intervention / baseline)") +
  theme_Publication()

ggsave("results/diet_parameter_delta_boxplot.pdf", p_delta, width = 11, height = 7, dpi = 150)
cat("\nSaved results/diet_parameter_delta_boxplot.pdf and results/diet_parameter_deltas.csv\n")

# ---- baseline vs intervention ("before/after") WITHIN each diet arm --------------------
# Paired Wilcoxon signed-rank test (a one-sample Wilcoxon of each subject's
# delta against 0), computed explicitly rather than via ggpubr's pairing
# detection. The boxplot shows every reliable value, so its n differs from the
# test's n. See "baseline vs intervention within each arm" in docs/diet-summary.md.
before_after_p <- box_data %>%
  group_by(parameter, diet) %>%
  summarise(y.position = max(value, na.rm = TRUE) * 1.08, .groups = "drop") %>%
  left_join(
    delta_data %>% group_by(parameter, diet) %>%
      summarise(n_paired = n(),
                p = if (n() >= 3) wilcox.test(delta, mu = 0)$p.value else NA_real_,
                .groups = "drop"),
    by = c("parameter", "diet")
  ) %>%
  mutate(group1 = "baseline", group2 = "intervention",
         label = if_else(is.na(p), sprintf("n=%d", n_paired), sprintf("p=%.3f (n=%d)", p, n_paired)))

write_csv(before_after_p, "results/diet_before_after_wilcoxon.csv")

# One PDF per parameter here too, faceted by diet arm, colored by visit.
plot_before_after <- function(param_name) {
  d <- box_data %>% filter(parameter == param_name)
  if (nrow(d) == 0) return(NULL)
  ann <- before_after_p %>% filter(parameter == param_name)
  ggplot(d, aes(visit, value, fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(width = 0.08, size = 1.1, alpha = 0.5, show.legend = FALSE) +
    stat_pvalue_manual(ann, label = "label", tip.length = 0.01, size = 2.8) +
    facet_wrap(vars(diet), nrow = 1, scales = "free_y", labeller = labeller(diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(title = sprintf("%s: baseline vs intervention within each diet arm (paired Wilcoxon)", PARAM_LABELS[[param_name]]),
         x = NULL, y = NULL) +
    theme_Publication()
}

for (param_name in c("ka", "kel", "F_12C", "F_13C6")) {
  p <- plot_before_after(param_name)
  if (!is.null(p)) {
    out_path <- sprintf("results/diet_before_after_boxplot_%s.pdf", param_name)
    ggsave(out_path, p, width = 7, height = 5, dpi = 150)
    cat("\nSaved", out_path, "\n")
  }
}
cat("Saved results/diet_before_after_wilcoxon.csv\n")

# ---- Two-wave ("second peak") characteristics by diet arm -----------------
# Exploratory only. See "Second-wave characteristics"
# in docs/diet-summary.md.

two_wave_rate <- fits %>% filter(reliable) %>%
  group_by(diet, visit) %>%
  summarise(n = n(), n_two_wave = sum(model == "two_wave"), .groups = "drop") %>%
  mutate(pct_two_wave = round(100 * n_two_wave / n, 1))

write_csv(two_wave_rate, "results/diet_two_wave_rate.csv")
cat("\n=== two_wave selection rate by diet x visit (reliable fits only) ===\n")
print(as.data.frame(two_wave_rate), digits = 3)

two_wave_tbl <- table(fits$diet[fits$reliable], fits$model[fits$reliable] == "two_wave")
if (all(dim(two_wave_tbl) == c(2, 2))) {
  fisher_p <- fisher.test(two_wave_tbl)$p.value
  cat("\nFisher's exact test, diet x two_wave selection: p =", round(fisher_p, 3), "\n")
}

two_wave_params <- fits %>% filter(reliable, model == "two_wave") %>%
  select(subject_id, diet, visit, t_lag, f_delayed, f_delayed_13C6)
write_csv(two_wave_params, "results/diet_two_wave_params.csv")

if (nrow(two_wave_params) >= 4) {
  p_two_wave <- two_wave_params %>%
    pivot_longer(c(t_lag, f_delayed, f_delayed_13C6), names_to = "parameter", values_to = "value") %>%
    filter(!is.na(value)) %>%
    ggplot(aes(diet, value, fill = diet)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(aes(color = diet), width = 0.08, size = 1.4, alpha = 0.7, show.legend = FALSE) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
    facet_wrap(vars(parameter), scales = "free_y",
               labeller = labeller(parameter = c(t_lag = "t_lag (min)", f_delayed = "f_delayed (12C)",
                                                  f_delayed_13C6 = "f_delayed (13C6)"))) +
    scale_x_discrete(labels = DIET_LABELS) +
    scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
    scale_color_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
    labs(title = sprintf("Two-wave characteristics by diet arm (exploratory, n=%d curves)", nrow(two_wave_params)),
         x = NULL, y = NULL) +
    theme_Publication()
  ggsave("results/diet_two_wave_boxplot.pdf", p_two_wave, width = 10, height = 4.5, dpi = 150)
  cat("\nSaved results/diet_two_wave_boxplot.pdf, results/diet_two_wave_rate.csv, results/diet_two_wave_params.csv\n")
} else {
  cat("\nToo few two_wave fits to plot a diet comparison (n =", nrow(two_wave_params), ")\n")
}

# ---- Average concentration-time curves per diet x visit x isotope ---------
# Must mirror simulate_fit() in 02_fit_erie_model.R exactly: the mechanism that
# won for each subject x visit determines the curve, and a mismatch (e.g. feeding
# k_release = NA to simulate_delayed_release()) blanks a whole group's mean. See
# "Median and mean concentration-time curves" in docs/diet-summary.md.
curve_matrix <- function(diet_val, vis, isotope) {
  sub <- fits %>% filter(diet == diet_val, visit == vis, reliable)
  if (nrow(sub) == 0) return(NULL)

  conc_list <- sub %>% pmap(function(...) {
    r <- list(...)
    if (isotope == "12C") {
      if (r$model == "two_wave") {
        simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_12C, dose, r$Vd)
        simulate_lagged_dose(simA, FINE_T_DIET, r$dose_12C_mg, r$f_delayed, r$t_lag)
      } else {
        bateman_conc(FINE_T_DIET, r$ka, r$kel, r$F_12C, r$dose_12C_mg, r$Vd)
      }
      # t_lag1_13C6/t_lag2_13C6 are checked before r$model: those 13C6 mechanisms
      # can win under either 12C model.
    } else if (!is.na(r$t_lag1_13C6)) {
      simB <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_13C6, dose, r$Vd)
      simulate_two_lag_dose(simB, FINE_T_DIET, dose_13C6_mg, f_delayed = 0, t_lag1 = r$t_lag1_13C6, gap = 0)
    } else if (!is.na(r$t_lag2_13C6)) {
      simB <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_13C6, dose, r$Vd)
      simulate_lagged_dose(simB, FINE_T_DIET, dose_13C6_mg, r$f_delayed2_13C6, r$t_lag2_13C6)
    } else if (r$model == "two_wave") {
      simB <- function(t, dose) simulate_delayed_release(t, r$k_release, r$ka, r$kel, r$F_13C6, dose, r$Vd)$conc
      simulate_lagged_dose(simB, FINE_T_DIET, dose_13C6_mg, r$f_delayed_13C6, r$t_lag)
    } else {
      simulate_delayed_release(FINE_T_DIET, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, r$Vd)$conc
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

average_curve_isotope <- function(diet_val, vis, isotope) {
  conc_mat <- curve_matrix(diet_val, vis, isotope)
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
# statistic (median with bootstrap CI band, mean with SEM band) and isotope (12C
# and 13C6 differ ~1000x in concentration). The y-axis is shared across the two
# diet panels.
CURVE_STATS <- list(
  median = list(cols = c("median_conc", "median_lo", "median_hi"), title = "Median %s concentration (95%% bootstrap CI)"),
  mean   = list(cols = c("mean_conc", "mean_lo", "mean_hi"),       title = "Mean %s concentration (SEM band)")
)

plot_diet_curves <- function(isotope_val, stat) {
  curves <- expand_grid(diet = c("low_fructose", "high_fructose"), visit = c("baseline", "intervention")) %>%
    pmap_dfr(function(diet, visit) {
      ac <- average_curve_isotope(diet, visit, isotope_val)
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
    labs(title = sprintf(paste0(CURVE_STATS[[stat]]$title, ": baseline vs intervention"), ISOTOPE_LABELS[[isotope_val]]),
         x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
}

for (stat in names(CURVE_STATS)) {
  for (isotope_val in c("12C", "13C6")) {
    p <- plot_diet_curves(isotope_val, stat)
    if (!is.null(p)) {
      out_path <- sprintf("results/diet_summary_curves_%s_%s.pdf", stat, isotope_val)
      ggsave(out_path, p, width = 9, height = 5, dpi = 150)
      cat("\nSaved", out_path, "\n")
    }
  }
}
# ---- Curve-derived PK metrics (AUC, Cmax, Tmax) per subject x visit -------
# Cmax and Tmax come from each subject's simulated fitted curve on FINE_T_DIET. See "Curve-derived
# metrics" in docs/diet-summary.md.
curve_metrics <- expand_grid(diet = c("low_fructose", "high_fructose"), visit = c("baseline", "intervention"),
                             isotope = c("12C", "13C6")) %>%
  pmap_dfr(function(diet, visit, isotope) {
    m <- curve_matrix(diet, visit, isotope)
    if (is.null(m)) return(NULL)
    f <- fits[match(paste(colnames(m), visit), paste(fits$subject_id, fits$visit)), ]
    # AUC to infinity: every mechanism delivers F * dose to the blood and kel clears it.
    dose_iso <- if (isotope == "12C") f$dose_12C_mg else dose_13C6_mg
    F_iso    <- if (isotope == "12C") f$F_12C else f$F_13C6
    tibble(
      subject_id = colnames(m), diet = diet, visit = visit, isotope = isotope,
      AUC  = F_iso * dose_iso / (f$Vd * f$kel),
      Cmax = apply(m, 2, max),
      Tmax = FINE_T_DIET[apply(m, 2, which.max)]
    )
  })
write_csv(curve_metrics, "results/diet_curve_metrics.csv")

# Extends the parameter LMMs above to the curve-derived metrics, same formula and
# random intercept, one model per metric x isotope. AUC and Cmax stay on the log
# scale like the fitted parameters; Tmax stays on its raw minutes scale, as in the
# paired test below. See "Linear mixed models" in docs/diet-summary.md.
metric_covariates <- fits %>%
  filter(reliable) %>%
  distinct(subject_id, visit, sex, age_years, bw_kg)

fit_metric_lmm <- function(metric_name, isotope_val, log_scale) {
  d <- curve_metrics %>% filter(isotope == isotope_val) %>%
    inner_join(metric_covariates, by = c("subject_id", "visit")) %>%
    transmute(subject_id, diet, visit, sex, age_years, bw_kg, value = .data[[metric_name]])
  run_lmm(d, sprintf("%s_%s", metric_name, isotope_val), log_scale = log_scale)
}

lmm_results_metrics <- bind_rows(
  fit_metric_lmm("AUC", "12C", log_scale = TRUE),
  fit_metric_lmm("AUC", "13C6", log_scale = TRUE),
  fit_metric_lmm("Cmax", "12C", log_scale = TRUE),
  fit_metric_lmm("Cmax", "13C6", log_scale = TRUE),
  fit_metric_lmm("Tmax", "12C", log_scale = FALSE),
  fit_metric_lmm("Tmax", "13C6", log_scale = FALSE)
)

lmm_results <- bind_rows(lmm_results_params, lmm_results_metrics) %>% arrange(parameter, term)
write_csv(lmm_results, "results/diet_lmm_results.csv")
cat("\n=== LMM (diet x visit, subject random intercept): F-tests ===\n")
print(as.data.frame(lmm_results), digits = 3)

curve_metrics_summary <- curve_metrics %>%
  pivot_longer(c(AUC, Cmax, Tmax), names_to = "metric", values_to = "value") %>%
  group_by(isotope, diet, visit, metric) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, 0.25), q75 = quantile(value, 0.75), .groups = "drop") %>%
  arrange(metric, isotope, diet, visit)
write_csv(curve_metrics_summary, "results/diet_curve_metrics_summary.csv")
cat("\n=== Curve-derived metrics: median (IQR) per diet x visit ===\n")
print(as.data.frame(curve_metrics_summary), digits = 3)

# Paired Wilcoxon signed-rank test, baseline vs intervention within each arm and isotope
# (subjects with both visits; needs >= 3 pairs). The signed differences are on the log scale
# for AUC and Cmax and in minutes for Tmax, and no multiplicity correction is applied.
# See "Curve-derived metrics" in docs/diet-summary.md.
metric_deltas <- curve_metrics %>%
  pivot_longer(c(AUC, Cmax, Tmax), names_to = "metric", values_to = "value") %>%
  pivot_wider(names_from = visit, values_from = value) %>%
  filter(!is.na(baseline), !is.na(intervention)) %>%
  mutate(delta = if_else(metric == "Tmax", intervention - baseline, log(intervention / baseline)))

metric_wilcoxon <- metric_deltas %>%
  group_by(metric, isotope, diet) %>%
  summarise(n_paired = n(), median_delta = median(delta),
            p = if (n() >= 3) tryCatch(wilcox.test(delta, mu = 0, exact = if (first(metric) == "Tmax") FALSE else NULL)$p.value,
                                       error = function(e) NA_real_) else NA_real_,
            .groups = "drop") %>%
  mutate(label = case_when(is.na(p) ~ sprintf("n=%d", n_paired),
                           p < 0.001 ~ sprintf("p<0.001 (n=%d)", n_paired),
                           TRUE ~ sprintf("p=%.3f (n=%d)", p, n_paired)))
write_csv(metric_wilcoxon, "results/diet_curve_metrics_wilcoxon.csv")
cat("\n=== Paired Wilcoxon, intervention vs baseline within each arm (curve-derived metrics) ===\n")
print(as.data.frame(metric_wilcoxon %>% select(-label)), digits = 3)

METRIC_LABELS <- c(AUC = "AUC 0-inf (mg*min/L)", Cmax = "Cmax (mg/L)", Tmax = "Tmax (min)")

for (metric in names(METRIC_LABELS)) {
  ann <- metric_wilcoxon %>% filter(metric == !!metric)
  p <- ggplot(curve_metrics, aes(visit, .data[[metric]], fill = visit)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
    geom_jitter(aes(color = visit), width = 0.08, size = 1.2, alpha = 0.6, show.legend = FALSE) +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label), vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
    facet_grid(rows = vars(isotope), cols = vars(diet), scales = "free_y",
               labeller = labeller(isotope = ISOTOPE_LABELS, diet = DIET_LABELS)) +
    scale_x_discrete(labels = VISIT_LABELS) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    scale_color_manual(values = VISIT_COLORS, labels = VISIT_LABELS, name = NULL) +
    labs(title = sprintf("%s by diet arm and visit (from fitted curves, paired Wilcoxon)", METRIC_LABELS[[metric]]), x = NULL, y = NULL) +
    theme_Publication()
  out_path <- sprintf("results/diet_curve_%s_boxplot.pdf", tolower(metric))
  ggsave(out_path, p, width = 8, height = 7, dpi = 150)
  cat("\nSaved", out_path, "\n")
}

cat("\nSaved results/diet_parameter_summary.csv\n")
