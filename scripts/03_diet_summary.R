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
DIET_LABELS    <- c(low_fructose = "Diet A: low fructose", high_fructose = "Diet B: high fructose")
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
    # Reliability gated on the 12C fit only: it converged, isn't stuck at
    # the shared kel bound, and its own R2 isn't flagged. Applied uniformly
    # to every parameter (including F_13C6/capsule dissolution) rather than
    # additionally requiring 13C6's own R2/k_release bound to pass - 13C6's
    # much smaller, noisier signal fails its own bar far more often even
    # when the underlying shared kinetics (from the same joint fit) are
    # trustworthy, which excluded a lot of otherwise-fine subjects. See
    # "Fit quality and what to trust" in docs/pk-model.md.
    reliable = converged & !kel_at_bound & !r2_12C_low
  )

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

fit_lmm <- function(param) {
  d <- fits %>% filter(reliable) %>%
    select(subject_id, diet, visit, value = all_of(param))
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
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_point(aes(color = diet), position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
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

# ---- Average concentration-time curves per diet x visit x isotope ---------

average_curve_isotope <- function(diet_val, vis, isotope) {
  sub <- fits %>% filter(diet == diet_val, visit == vis, reliable)
  if (nrow(sub) == 0) return(NULL)

  conc_list <- sub %>% pmap(function(...) {
    r <- list(...)
    if (isotope == "12C") {
      bateman_conc(FINE_T_DIET, r$ka, r$kel, r$F_12C, r$dose_12C_mg, r$Vd)
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
