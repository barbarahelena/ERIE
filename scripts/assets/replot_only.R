# Quick standalone re-plot from an already-written results/fit_results.csv -
# re-simulates and re-draws results/plots_individual/*.pdf using the
# CURRENT plot_subject_fit()/simulate_fit() logic in 02_fit_erie_model.R,
# without re-running the (slow) fitting itself. Useful for verifying a
# plotting-only fix (e.g. the t_lag1_13C6 rendering bug) against results
# that were already fit under the old plotting code.
#
# Not part of the regular pipeline - a one-off verification tool, not
# sourced or called by 02_fit_erie_model.R or 03_diet_summary.R.

suppressMessages({
  library(readr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(grid)
  library(ggthemes)
})

source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_diagnostics.R")

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
                strip.text.y = element_text(size = rel(1.05)),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

ISOTOPE_LABELS <- c("12C" = "Fructose 12C", "13C6" = "Fructose 13C6")
ISOTOPE_COLORS <- c("12C" = "steelblue", "13C6" = "firebrick")
VISIT_LABELS   <- c(baseline = "Baseline (FCT1)", intervention = "Intervention (FCT2)")
CLEARANCE_FRAC <- 0.01

concentrations <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE)
covariates     <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants      <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
MW_12C  <- constants$value[constants$constant == "MW_12C"]
MW_13C6 <- constants$value[constants$constant == "MW_13C6"]
dose_13C6_mg <- constants$value[constants$constant == "tracer_13C6_dose_mg"]

data <- concentrations %>%
  left_join(covariates, by = c("subject_id", "visit")) %>%
  mutate(conc_mgL = conc_umol_L * if_else(isotope == "12C", MW_12C, MW_13C6) / 1000) %>%
  filter(!is.na(bw_kg), !is.na(conc_mgL))
baseline_12C <- data %>% filter(isotope == "12C", time_min == 0) %>%
  select(subject_id, visit, baseline_mgL = conc_mgL)
data_corrected <- data %>%
  left_join(baseline_12C, by = c("subject_id", "visit")) %>%
  mutate(conc_mgL = if_else(isotope == "12C", conc_mgL - baseline_mgL, conc_mgL)) %>%
  select(-baseline_mgL)

results <- read_csv("results/fit_results.csv", show_col_types = FALSE)

simulate_fit <- function(sid, vis, r) {
  cov <- covariates %>% filter(subject_id == sid, visit == vis)
  Vd <- nadler_blood_volume(cov$bw_kg, cov$height_cm, cov$sex)
  dose_12C <- cov$dose_12C_mg
  fine <- seq(0, 400, length.out = 400)

  if (r$model == "two_wave") {
    simA <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_12C, dose, Vd)
    sim12 <- simulate_lagged_dose(simA, fine, dose_12C, r$f_delayed, r$t_lag)
  } else {
    sim12 <- bateman_conc(fine, r$ka, r$kel, r$F_12C, dose_12C, Vd)
  }

  simB_instant <- function(t, dose) bateman_conc(t, r$ka, r$kel, r$F_13C6, dose, Vd)
  sim13 <- if (!is.na(r$t_lag1_13C6)) {
    simulate_two_lag_dose(simB_instant, fine, dose_13C6_mg, f_delayed = 0, t_lag1 = r$t_lag1_13C6, gap = 0)
  } else if (!is.na(r$t_lag2_13C6)) {
    simulate_lagged_dose(simB_instant, fine, dose_13C6_mg, r$f_delayed2_13C6, r$t_lag2_13C6)
  } else if (r$model == "two_wave") {
    simB <- function(t, dose) simulate_delayed_release(t, r$k_release, r$ka, r$kel, r$F_13C6, dose, Vd)$conc
    simulate_lagged_dose(simB, fine, dose_13C6_mg, r$f_delayed_13C6, r$t_lag)
  } else {
    simulate_delayed_release(fine, r$k_release, r$ka, r$kel, r$F_13C6, dose_13C6_mg, Vd)$conc
  }

  t_end <- max(time_to_clearance(fine, sim12, CLEARANCE_FRAC),
               time_to_clearance(fine, sim13, CLEARANCE_FRAC))
  keep <- fine <= t_end

  bind_rows(
    tibble(time_min = fine[keep], conc_mgL = sim12[keep], isotope = "12C"),
    tibble(time_min = fine[keep], conc_mgL = sim13[keep], isotope = "13C6")
  )
}

fmt_aic <- function(x) if_else(is.na(x), "NA", sprintf("%.0f", x))

plot_subject_fit <- function(sid) {
  obs <- data_corrected %>% filter(subject_id == sid) %>% select(visit, isotope, time_min, conc_mgL)
  fits <- results %>% filter(subject_id == sid, !is.na(ka))
  if (nrow(fits) == 0) return(NULL)

  sim <- fits %>% pmap_dfr(function(...) {
    row <- tibble(...)
    bind_cols(visit = row$visit, simulate_fit(sid, row$visit, row))
  })

  ann <- bind_rows(
    fits %>% transmute(visit, isotope = "12C",
                        label = sprintf("R2=%.3f\nAIC sw/tw=%s/%s", r2_12C,
                                        fmt_aic(aic_single_wave), fmt_aic(aic_two_wave))),
    fits %>% transmute(visit, isotope = "13C6", label = sprintf("R2=%.3f", r2_13C6))
  )

  ggplot(obs, aes(time_min, conc_mgL)) +
    geom_point(color = "black", size = 1.4) +
    geom_line(data = sim, aes(color = isotope), linewidth = 0.8, show.legend = FALSE) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
              hjust = 1.05, vjust = 1.2, size = 2.5, lineheight = 0.9) +
    facet_grid(rows = vars(isotope), cols = vars(visit), scales = "free",
               labeller = labeller(visit = VISIT_LABELS, isotope = ISOTOPE_LABELS)) +
    scale_color_manual(values = ISOTOPE_COLORS) +
    labs(title = sid, x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
}

dir.create("results/plots_individual", showWarnings = FALSE, recursive = TRUE)
n <- 0
for (sid in unique(results$subject_id)) {
  p <- plot_subject_fit(sid)
  if (!is.null(p)) {
    ggsave(file.path("results/plots_individual", sprintf("%s_joint.pdf", sid)),
           p, width = 8, height = 5, dpi = 120)
    n <- n + 1
  }
}
cat("Re-plotted", n, "subjects to results/plots_individual/*.pdf\n")
