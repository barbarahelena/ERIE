# Diet-arm summary: average simulated concentration-time curves and a
# parameter comparison table, split by dietary arm (A vs. B) and visit
# (FCT1 = before diet, FCT2 = after diet). Reuses the fits already produced
# by 02_fit_erie_model.R - does not refit anything.
#
# Diet A: low fructose, calories matched with glucose supplementation.
# Diet B: high fructose.
# Barbara Verhaar

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

source("scripts/assets/pk_curves.R")

# House plotting theme
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
# NOTE: also duplicated in 02_fit_erie_model.R - keep both in sync if changed.
ISOTOPE_LABELS <- c("12C" = "Fructose 12C (unlabelled)", "13C6" = "Fructose 13C6 (labelled)")
VISIT_LABELS   <- c(FCT1 = "FCT1 (before diet)", FCT2 = "FCT2 (after diet)")
# diet values are recoded from the raw "A"/"B" at the source, in
# 01_clean_data.R - see erie_covariates.csv.
DIET_LABELS    <- c(low_fructose = "Diet A: low fructose", high_fructose = "Diet B: high fructose")
DIET_COLORS    <- c(low_fructose = "#1b9e77", high_fructose = "#d95f02")

# Vd = estimated total blood volume (Nadler 1962, scripts/assets/pk_curves.R),
# from each subject's own weight/height/sex - see "Volume of distribution"
# in docs/pk-model.md.
FINE_T_DIET <- seq(0, 400, by = 2)

if (!file.exists("results/fit_results_joint.csv")) {
  stop("results/fit_results_joint.csv not found - run scripts/02_fit_erie_model.R first.")
}

results      <- read_csv("results/fit_results_joint.csv", show_col_types = FALSE)
covariates   <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
constants    <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
dose_13C6_mg <- constants$value[constants$constant == "tracer_13C6_dose_mg"]

fits <- results %>%
  left_join(covariates %>% select(subject_id, visit, bw_kg, height_cm, sex, diet), by = c("subject_id", "visit")) %>%
  filter(!is.na(ka))

if (any(is.na(fits$diet))) {
  warning(sum(is.na(fits$diet)), " fitted subject x visit rows have no diet assignment and are excluded from the diet summary")
}

fits <- fits %>%
  filter(!is.na(diet)) %>%
  mutate(
    Vd = nadler_blood_volume(bw_kg, height_cm, sex),
    dose_12C_mg = 1000 * bw_kg,
    # Per-curve reliability: a joint fit is only as trustworthy as (a) it
    # converged, (b) the shared kel isn't stuck at its bound (which taints
    # both curves, since ka/kel are fit jointly), and (c) that curve's own
    # R2 - and for 13C6, k_release - aren't flagged. See "Fit quality and
    # what to trust" in docs/pk-model.md.
    reliable_12C    = converged & !kel_at_bound & !r2_12C_low,
    reliable_13C6   = converged & !kel_at_bound & !k_release_at_bound & !r2_13C6_low,
    reliable_shared = reliable_12C & reliable_13C6   # for ka/kel, shared across both curves
  )

# ---- Parameter summary table (mean/SD/SEM/n per diet x visit) -------------

summarise_param <- function(param, reliable_col) {
  fits %>%
    filter(.data[[reliable_col]]) %>%
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
  summarise_param("ka",  "reliable_shared"),
  summarise_param("kel", "reliable_shared"),
  summarise_param("F_12C", "reliable_12C"),
  summarise_param("F_13C6", "reliable_13C6"),
  summarise_param("capsule_dissolution_halflife_min", "reliable_13C6")
) %>% arrange(parameter, diet, visit)

dir.create("results", showWarnings = FALSE)
write_csv(param_table, "results/diet_parameter_summary.csv")
cat("=== Parameter summary by diet x visit (reliable fits only) ===\n")
print(as.data.frame(param_table), digits = 3)

# ---- Average concentration-time curves per diet x visit x isotope ---------

average_curve_isotope <- function(diet_val, vis, isotope) {
  reliable_col <- if (isotope == "12C") "reliable_12C" else "reliable_13C6"
  sub <- fits %>% filter(diet == diet_val, visit == vis, .data[[reliable_col]])
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

ggsave("results/diet_summary_curves.png", p, width = 10, height = 7, dpi = 150)
cat("\nSaved results/diet_summary_curves.png and results/diet_parameter_summary.csv\n")
