# Compare this model's per-subject R2 against the former (validated) joint
# model's saved results (former_models/MixedModel/Scripts/
# fructose_joint_model_final.R output). Same subjects, same visits, same R2
# definition (1 - SSE/SS_tot on the observed points) - different Vd,
# objective, and fitting procedure (see "Origin of the model" and "Volume of
# distribution" in docs/pk-model.md and docs/problems-and-fixes.md for what
# changed and why).
# This is a check on how those changes affected fit quality, not a claim
# that either model's R2 is "the truth."
# Barbara Verhaar

suppressMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(ggthemes)
})

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
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))

}

ISOTOPE_LABELS <- c("12C" = "Fructose 12C (unlabelled)", "13C6" = "Fructose 13C6 (labelled)")

FORMER_RESULTS_CSV <- "former_models/MixedModel/Results/fit_results_joint.csv"
CURRENT_RESULTS_CSV <- "results/fit_results.csv"

if (!file.exists(CURRENT_RESULTS_CSV)) {
  stop(CURRENT_RESULTS_CSV, " not found - run scripts/02_fit_erie_model.R first.")
}

# Normalizes either visit vocabulary (raw FCT1/FCT2, or this project's coded
# baseline/intervention) to baseline/intervention, so this script doesn't
# care whether 01_clean_data.R's visit recode has been run yet.
normalize_visit <- function(v) {
  dplyr::recode_values(v,
    from = c("FCT1", "FCT2", "baseline", "intervention"),
    to   = c("baseline", "intervention", "baseline", "intervention"))
}

former <- read_csv(FORMER_RESULTS_CSV, show_col_types = FALSE) %>%
  transmute(subject_id, visit = normalize_visit(visit),
            r2_12C_former = r2_12C, r2_13C6_former = r2_13C6)

current <- read_csv(CURRENT_RESULTS_CSV, show_col_types = FALSE) %>%
  transmute(subject_id, visit = normalize_visit(visit),
            r2_12C_current = r2_12C, r2_13C6_current = r2_13C6)

# inner_join: only subject x visit pairs present (and fitted) in both models
# are comparable - e.g. a subject dropped from one model's cohort just isn't
# part of this comparison, rather than showing up as a misleading NA "loss".
cmp <- former %>%
  inner_join(current, by = c("subject_id", "visit")) %>%
  filter(!is.na(r2_12C_former), !is.na(r2_12C_current),
         !is.na(r2_13C6_former), !is.na(r2_13C6_current)) %>%
  mutate(
    delta_12C  = r2_12C_current  - r2_12C_former,
    delta_13C6 = r2_13C6_current - r2_13C6_former
  )

dir.create("results", showWarnings = FALSE)
write_csv(cmp, "results/r2_comparison_vs_former_model.csv")

summary_stats <- tibble(
  isotope           = c("12C", "13C6"),
  n                 = c(sum(!is.na(cmp$delta_12C)), sum(!is.na(cmp$delta_13C6))),
  median_former     = c(median(cmp$r2_12C_former), median(cmp$r2_13C6_former)),
  median_current    = c(median(cmp$r2_12C_current), median(cmp$r2_13C6_current)),
  median_delta      = c(median(cmp$delta_12C), median(cmp$delta_13C6)),
  n_improved        = c(sum(cmp$delta_12C > 0), sum(cmp$delta_13C6 > 0)),
  n_worse           = c(sum(cmp$delta_12C < 0), sum(cmp$delta_13C6 < 0)),
  correlation       = c(cor(cmp$r2_12C_former, cmp$r2_12C_current),
                         cor(cmp$r2_13C6_former, cmp$r2_13C6_current))
)
write_csv(summary_stats, "results/r2_comparison_summary.csv")

# ---- Scatter plot: former (x) vs current (y) R2, per isotope --------------

plot_df <- bind_rows(
  cmp %>% transmute(subject_id, visit, isotope = "12C",
                     former = r2_12C_former, current = r2_12C_current),
  cmp %>% transmute(subject_id, visit, isotope = "13C6",
                     former = r2_13C6_former, current = r2_13C6_current)
)

p <- ggplot(plot_df, aes(former, current)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.7) +
  facet_wrap(vars(isotope), labeller = labeller(isotope = ISOTOPE_LABELS)) +
  coord_equal() +
  labs(title = "Per-subject R²: former model vs. current model",
       subtitle = "Above the dashed line = current model fits better; below = former model fits better",
       x = "Former model R²", y = "Current model R²") +
  theme_Publication()

ggsave("results/r2_comparison_plot.pdf", p, width = 10, height = 5.5)
