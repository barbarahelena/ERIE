# Diet-arm summary plots
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
VISIT_LABELS   <- c(FCT1 = "FCT1 (before diet)", FCT2 = "FCT2 (after diet)")
DIET_LABELS    <- c(low_fructose = "Low fructose diet", high_fructose = "High fructose diet")
DIET_COLORS    <- c(low_fructose = "#1b9e77", high_fructose = "#d95f02")
FINE_T_DIET <- seq(0, 400, by = 2)

# Open data
results      <- read_csv("results/fit_results.csv", show_col_types = FALSE)
covariates   <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants    <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
dose_13C6_mg <- constants$value[constants$constant == "tracer_13C6_dose_mg"]

fits <- results %>%
  left_join(covariates %>% select(subject_id, visit, bw_kg, height_cm, sex, diet), by = c("subject_id", "visit")) %>%
  filter(!is.na(ka))

fits <- fits %>%
  filter(!is.na(diet)) %>%
  mutate(
    Vd = nadler_blood_volume(bw_kg, height_cm, sex),
    dose_12C_mg = 1000 * bw_kg,
    # Reliability gated on the 12C fit only: it isn't stuck at the shared
    # kel bound, and its own R2 isn't flagged. Applied uniformly to every
    # parameter (including F_13C6/capsule dissolution) rather than
    # additionally requiring 13C6's own R2/k_release bound to pass - 13C6's
    # much smaller, noisier signal fails its own bar far more often even
    # when the underlying shared kinetics (from the same joint fit) are
    # trustworthy, which excluded a lot of otherwise-fine subjects. See
    # "Fit quality and what to trust" in docs/pk-model.md.
    #
    # `converged` deliberately excluded: it only reflects whether optim()'s
    # L-BFGS-B hit its strict internal stopping criterion vs. its maxit cap,
    # not whether the winning fit is actually good - checked directly on
    # this cohort's own fit_results.csv, r2_12C for converged=FALSE rows
    # (mean 0.902, n=22) was statistically indistinguishable from, if
    # anything slightly better than, converged=TRUE rows (mean 0.890,
    # n=46), and its minimum (0.762) was far ABOVE the converged group's
    # minimum (0.249, the known-hard subjects). Requiring it anyway was
    # needlessly discarding good fits, especially here where a delta/
    # two-wave comparison requires BOTH visits reliable simultaneously - a
    # 32% per-curve "not converged" rate compounds to exclude the majority
    # of subjects from any paired comparison even though the fits it drops
    # are just as trustworthy by R2 as the ones it keeps.
    reliable = !kel_at_bound & !r2_12C_low
  )

# ---- Volume of distribution (Vd) by diet arm -------------------------------
# Vd isn't a fitted parameter - it's a direct, deterministic function of
# each subject's own weight/height/sex via Nadler's equation (see "Volume
# of distribution (Vd)" in docs/pk-model.md) - so it isn't gated on
# `reliable` (a PK fit's own convergence/R2 has no bearing on it), and isn't
# split by visit either: it's essentially fixed per subject (weight rarely
# changes meaningfully within the study), so one value per subject
# (averaged across whichever visits are present) compared across diet arms
# is more useful than a visit-split view that just adds noise. Useful as a
# covariate-balance check: it should cluster by sex (built into the
# formula) and be reasonably similar between diet arms if randomization
# worked as intended.

vd_data <- fits %>% group_by(subject_id, diet, sex) %>%
  summarise(Vd = mean(Vd), .groups = "drop")

p_vd <- ggplot(vd_data, aes(diet, Vd, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(aes(shape = sex), width = 0.08, size = 1.8, alpha = 0.8) +
  scale_x_discrete(labels = DIET_LABELS) +
  scale_fill_manual(values = DIET_COLORS, guide = "none") +
  labs(title = "Nadler-estimated blood volume (Vd) by diet arm",
       x = NULL, y = "Vd (L)", shape = "Sex") +
  theme_Publication()

ggsave("results/diet_vd_boxplot.pdf", p_vd, width = 6, height = 5, dpi = 150)
cat("\nSaved results/diet_vd_boxplot.pdf\n")

# ---- Fit quality (R2) by diet arm and visit --------------------------------
# Deliberately NOT gated on `reliable` - reliable excludes on r2_12C_low,
# which is DERIVED from r2_12C itself, so filtering by it here would hide
# exactly the poor fits this plot exists to surface. A QC/diagnostic view of
# fit quality across the whole cohort, not a trustworthy-subset comparison
# like the parameter plots below.
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
  summarise_param("F_13C6"),
  summarise_param("capsule_dissolution_halflife_min")
) %>% arrange(parameter, diet, visit)

dir.create("results", showWarnings = FALSE)
write_csv(param_table, "results/diet_parameter_summary.csv")
cat("=== Parameter summary by diet x visit (reliable fits only) ===\n")
print(as.data.frame(param_table), digits = 3)

# ---- Statistical comparison: diet x time linear mixed models --------------
# One LMM per parameter, testing whether diet arm, visit (before/after the
# diet), or their interaction (the actual "did the diet change this
# differently by arm" question) explains variation - subject_id as a random
# intercept, since each subject contributes a paired FCT1/FCT2 observation
# (repeated measures), not two independent ones. Uses the same reliability
# filter as the descriptive summary above; lmer handles the resulting
# unbalanced design (not every subject has both visits reliable) without
# needing complete pairs, unlike a paired t-test.

# log_transform defaults TRUE: ka/kel/k_release-derived halflife are rate
# constants and F_12C/F_13C6 are bioavailable fractions - all positive,
# multiplicative-scale quantities standardly treated as log-normal in PK
# work (not normal on their raw scale). Log-transforming makes the LMM's
# normal-residuals assumption more defensible and turns a diet effect into
# a fold-change rather than an absolute difference, which is the more
# natural scale for a rate constant.
fit_lmm <- function(param, log_transform = TRUE) {
  d <- fits %>% filter(reliable) %>%
    select(subject_id, diet, visit, value = all_of(param))
  if (log_transform) d <- d %>% mutate(value = log(value))
  n_subjects_both <- d %>% count(subject_id) %>% filter(n == 2) %>% nrow()
  if (n_subjects_both < 3) {
    warning(param, ": fewer than 3 subjects with both visits reliable - skipping LMM")
    return(NULL)
  }
  model <- lmer(value ~ diet * visit + (1 | subject_id), data = d)
  a <- anova(model)  # Type III, Satterthwaite df (lmerTest default)
  tibble(parameter = param, term = rownames(a), `F` = a$`F value`, df1 = a$NumDF, df2 = a$DenDF, p = a$`Pr(>F)`)
}

lmm_results <- bind_rows(
  fit_lmm("ka"),
  fit_lmm("kel"),
  fit_lmm("F_12C"),
  fit_lmm("F_13C6"),
  fit_lmm("capsule_dissolution_halflife_min")
)

write_csv(lmm_results, "results/diet_lmm_results.csv")
cat("\n=== LMM (diet x visit, subject random intercept): F-tests ===\n")
print(as.data.frame(lmm_results), digits = 3)

# ---- Boxplot: diet x visit distributions for each parameter ---------------
PARAM_LABELS <- c(ka = "ka (1/min)", kel = "kel (1/min)", F_12C = "F[12C]", F_13C6 = "F[13C6]",
                   capsule_dissolution_halflife_min = "Capsule t1/2 (min)")

box_data <- bind_rows(
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "ka", value = ka),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "kel", value = kel),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "F_12C", value = F_12C),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "F_13C6", value = F_13C6),
  fits %>% filter(reliable) %>% transmute(subject_id, diet, visit, parameter = "capsule_dissolution_halflife_min",
                                           value = capsule_dissolution_halflife_min)
)

p_box <- ggplot(box_data, aes(visit, value, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(aes(color = diet), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.6),
             size = 1.2, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(vars(parameter), scales = "free_y", nrow = 2,
             labeller = labeller(parameter = PARAM_LABELS)) +
  scale_x_discrete(labels = VISIT_LABELS) +
  scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  scale_color_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  labs(title = "Fitted PK parameters by diet arm and visit", x = NULL, y = NULL) +
  theme_Publication()

ggsave("results/diet_parameter_boxplot.pdf", p_box, width = 11, height = 7, dpi = 150)
cat("\nSaved results/diet_parameter_boxplot.pdf and results/diet_lmm_results.csv\n")

# ---- Delta plots: within-subject before/after diet change, by diet arm ----
# The boxplot above compares FCT1 and FCT2 as separate distributions; the
# actual "did the diet change this parameter" question is the per-subject
# FCT2-vs-FCT1 change, only defined for subjects with BOTH visits reliable
# (so every delta is a complete pair, not a mix of paired and unpaired
# values). Reported as a log fold-change (log(FCT2) - log(FCT1)), matching
# the log-transformed LMMs above and for the same reason - these are
# rate/fraction-like quantities better compared multiplicatively. Compared
# between diet arms with a Wilcoxon rank-sum test (ggpubr::stat_compare_means)
# rather than a t-test, since n per arm is small and not assumed normal.

delta_param <- function(param) {
  fits %>% filter(reliable) %>%
    select(subject_id, diet, visit, value = all_of(param)) %>%
    pivot_wider(names_from = visit, values_from = value) %>%
    filter(!is.na(FCT1), !is.na(FCT2)) %>%
    transmute(subject_id, diet, parameter = param, delta = log(FCT2) - log(FCT1))
}

delta_data <- bind_rows(
  delta_param("ka"), delta_param("kel"), delta_param("F_12C"),
  delta_param("F_13C6"), delta_param("capsule_dissolution_halflife_min")
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
  labs(title = "Within-subject diet change (log fold-change, FCT2 vs FCT1) by diet arm",
       x = NULL, y = "log(FCT2 / FCT1)") +
  theme_Publication()

ggsave("results/diet_parameter_delta_boxplot.pdf", p_delta, width = 11, height = 7, dpi = 150)
cat("\nSaved results/diet_parameter_delta_boxplot.pdf and results/diet_parameter_deltas.csv\n")

# ---- FCT1 vs FCT2 ("before/after") WITHIN each diet arm --------------------
# A different question from the delta plot above (which compares the SIZE
# of the FCT2-FCT1 change BETWEEN diet arms): this asks whether FCT1 and
# FCT2 differ at all WITHIN each diet arm on its own. Computed as a paired
# Wilcoxon signed-rank test - equivalent to a one-sample Wilcoxon test of
# each subject's own delta (already computed above) against 0 - rather
# than ggpubr's automatic pairing detection across facets, which silently
# breaks if the same subject doesn't land in the same row order in every
# facet. box_data (below) still shows every reliable single-visit value,
# including subjects who only have one visit reliable - a larger, more
# representative sample than the paired test itself can use (pairing
# necessarily requires both visits), so the boxplot's own n and the test's
# n legitimately differ; that's expected, not a bug.
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
  mutate(group1 = "FCT1", group2 = "FCT2",
         label = if_else(is.na(p), sprintf("n=%d", n_paired), sprintf("p=%.3f (n=%d)", p, n_paired)))

write_csv(before_after_p, "results/diet_before_after_wilcoxon.csv")

p_before_after <- ggplot(box_data, aes(visit, value, fill = diet)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.08, size = 1.1, alpha = 0.5, show.legend = FALSE) +
  stat_pvalue_manual(before_after_p, label = "label", tip.length = 0.01, size = 2.8) +
  facet_grid(rows = vars(parameter), cols = vars(diet), scales = "free_y",
             labeller = labeller(parameter = PARAM_LABELS, diet = DIET_LABELS)) +
  scale_x_discrete(labels = VISIT_LABELS) +
  scale_fill_manual(values = DIET_COLORS, guide = "none") +
  labs(title = "FCT1 vs FCT2 within each diet arm (paired Wilcoxon signed-rank test)",
       x = NULL, y = NULL) +
  theme_Publication()

ggsave("results/diet_before_after_boxplot.pdf", p_before_after, width = 9, height = 11, dpi = 150)
cat("\nSaved results/diet_before_after_boxplot.pdf and results/diet_before_after_wilcoxon.csv\n")

# ---- Two-wave ("second peak") characteristics by diet arm -----------------
# Whether a genuine second wave was detected at all (two_wave selected, via
# has_peak_dip_rise() in 02_fit_erie_model.R) is itself a diet-relevant
# outcome, and - among the subjects who show one - so is its timing
# (t_lag, "time between peaks") and how much of the dose rode it
# (f_delayed for 12C, f_delayed_13C6 for 13C6's own independent choice).
# Exploratory only: just 12/68 curves are two_wave, so this is underpowered
# and reported for completeness, not as a confirmed effect.

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
# Mirrors simulate_fit() in 02_fit_erie_model.R exactly - both curves' shape
# depends on which model/mechanism won for that subject x visit, not just
# on the plain single-compartment equations. Getting this wrong isn't just
# a shape mismatch: two_wave's 13C6 stage always carries k_release, but
# single_wave's newer 13C6 mechanisms (t_lag1_13C6 onset lag,
# f_delayed2_13C6/t_lag2_13C6 second wave) both explicitly set
# k_release = NA once they win (see fit_subject_visit_single_wave()) -
# feeding that NA into simulate_delayed_release() silently returns NA for
# every t>0, and since rowMeans()/sd() below aren't na.rm, ONE such subject
# in a diet x visit group is enough to blank out that group's entire mean
# curve. Confirmed: 14 of 21 single_wave rows in the full cohort have
# k_release = NA (11 onset-lag + 3 second-wave), enough to hit virtually
# every group - this is what silently broke every 13C6 panel.
average_curve_isotope <- function(diet_val, vis, isotope) {
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
      # t_lag1_13C6/t_lag2_13C6 checked FIRST, before dispatching on
      # r$model - 13C6's own onset-lag/independent-second-wave mechanisms
      # can now win under EITHER 12C model (see fit_subject_visit_two_wave()
      # stage 2), so model alone no longer determines which 13C6 mechanism
      # is actually in play.
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

  tibble(
    time_min  = FINE_T_DIET,
    mean_conc = rowMeans(conc_mat),
    sem_conc  = apply(conc_mat, 1, sd) / sqrt(ncol(conc_mat)),
    n         = ncol(conc_mat)
  )
}

diet_curves <- expand_grid(diet = c("low_fructose", "high_fructose"), visit = c("FCT1", "FCT2"), isotope = c("12C", "13C6")) %>%
  pmap_dfr(function(diet, visit, isotope) {
    ac <- average_curve_isotope(diet, visit, isotope)
    if (is.null(ac)) return(NULL)
    ac %>% mutate(diet = diet, visit = visit, isotope = isotope)
  })

# facet_grid (isotope x visit), not facet_wrap, so each ISOTOPE ROW gets its
# own free y-scale shared across the two visit columns - the 12C and 13C6
# curves differ by ~1000x in concentration (same reason 02_fit_erie_model.R's
# per-subject plots use free scales), but FCT1 vs FCT2 for the same isotope
# are on a comparable scale and are usefully left directly comparable.
p <- ggplot(diet_curves, aes(time_min, mean_conc, color = diet, fill = diet)) +
  geom_ribbon(aes(ymin = mean_conc - sem_conc, ymax = mean_conc + sem_conc), alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(rows = vars(isotope), cols = vars(visit), scales = "free_y",
             labeller = labeller(isotope = ISOTOPE_LABELS, visit = VISIT_LABELS)) +
  scale_color_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  scale_fill_manual(values = DIET_COLORS, labels = DIET_LABELS, name = NULL) +
  labs(title = "Mean fructose concentration by dietary arm",
       x = "Time (min)", y = "Concentration (mg/L)") +
  theme_Publication()

ggsave("results/diet_summary_curves.pdf", p, width = 10, height = 7, dpi = 150)
cat("\nSaved results/diet_summary_curves.pdf and results/diet_parameter_summary.csv\n")
