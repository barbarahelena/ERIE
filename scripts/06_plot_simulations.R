# =============================================================================
# 06_plot_simulations.R - supplementary figure from the replicates of 04 and 05 (no refitting)
#   A  share of replicates with p < 0.05 on ka (95% Wilson CI)
#   B  error in the fitted log ka change (median, 2.5-97.5 percentiles)
# Output: results/02_diagnostics_ka_kel/simulation/simulation_figure.pdf and simulation_figure_ka.pdf
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(ggplot2); library(grid); library(ggthemes); library(ggpubr) })

# House plotting theme (same as 03_diet_summary_shared_kel.R)
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

SIM_DIR <- "results/02_diagnostics_ka_kel/simulation"
FIT_COLORS <- c(independent = "#2a78d6", shared = "#eb6834")
FIT_LABELS <- c(independent = "Each visit fitted alone", shared = "One kel per subject")
FIT_SHAPES <- c(independent = 16, shared = 17)

# Scenario labels: truth at baseline (ka / kel, 1/min) and what changes at the intervention visit
SCENARIOS <- tribble(
  ~scenario,     ~group,                                ~label,
  "null_sep",    "No change",                           "ka 0.040, kel 0.025",
  "null_diag",   "No change",                           "ka 0.025, kel 0.025",
  "ka_up_sep",   "ka x 1.5",                            "ka 0.040, kel 0.025",
  "ka_up_diag",  "ka x 1.5",                            "ka 0.025, kel 0.025",
  "ka_up_flip",  "ka x 1.5",                            "ka 0.020, kel 0.050",
  "kel_up_sep",  "kel x 1.5\n(assumption violated)",    "ka 0.040, kel 0.025",
  "kel_up_diag", "kel x 1.5\n(assumption violated)",    "ka 0.025, kel 0.025"
) %>% mutate(group = factor(group, unique(group)), label = factor(paste(scenario, label)))
# Per panel: title, x-axis title and reference line (the 5% test level in A, no error in B)
PANELS <- tribble(
  ~panel, ~title,                       ~xlab,                                                         ~ref,
  "rej",  "ka change detected",         "Share of replicates with p < 0.05",                           0.05,
  "bias", "Error in fitted ka change",  "Fitted minus true log(ka intervention / ka baseline)",         0
)

indep  <- read_csv(file.path(SIM_DIR, "simulation_replicates.csv"), show_col_types = FALSE) %>% mutate(fit = "independent")
shared <- read_csv(file.path(SIM_DIR, "simulation_shared_kel_replicates.csv"), show_col_types = FALSE) %>% mutate(fit = "shared")

wilson <- function(x, n) {   # 95% Wilson interval for x successes out of n
  z <- 1.96; p <- x / n; mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  tibble(est = p, lo = mid - half, hi = mid + half)
}

rej <- bind_rows(indep, shared) %>% filter(!is.na(p_ka)) %>% group_by(scenario, fit) %>%
  summarise(x = sum(p_ka < 0.05), n = n(), .groups = "drop") %>%
  mutate(wilson(x, n), panel = "rej")
bias <- bind_rows(indep, shared) %>% group_by(scenario, fit) %>%
  summarise(est = median(bias_ka_change, na.rm = TRUE), lo = quantile(bias_ka_change, 0.025, na.rm = TRUE),
            hi = quantile(bias_ka_change, 0.975, na.rm = TRUE), .groups = "drop") %>% mutate(panel = "bias")

d <- bind_rows(rej, bias) %>% inner_join(SCENARIOS, by = "scenario")

dodge <- position_dodge(width = 0.6)
plot_panel <- function(d, which, show_y) {
  pn <- PANELS %>% filter(panel == which)
  p <- ggplot(d %>% filter(panel == which), aes(est, label, color = fit, shape = fit)) +
    geom_vline(xintercept = pn$ref, color = "grey55", linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, position = dodge, na.rm = TRUE) +
    geom_point(size = 2.6, position = dodge) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_y_discrete(labels = function(x) sub("^\\S+ ", "", x), limits = rev) +
    scale_color_manual(values = FIT_COLORS, labels = FIT_LABELS, name = NULL) +
    scale_shape_manual(values = FIT_SHAPES, labels = FIT_LABELS, name = NULL) +
    labs(title = pn$title, x = pn$xlab, y = NULL) +
    theme_Publication(base_size = 12) +
    # Facet rows: scenario group as a horizontal label outside the y axis
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
          legend.key.size = unit(0.5, "cm"), plot.margin = unit(c(5, 5, 2, 2), "mm"))
  # Panel B shares A's rows, so it drops the row labels
  if (!show_y) p <- p + theme(strip.text.y.left = element_blank(), strip.background = element_blank(),
                              axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p
}
plot_sim <- function(d, width_a = 1.45) {
  fig <- ggarrange(plot_panel(d, "rej", TRUE), plot_panel(d, "bias", FALSE), ncol = 2, widths = c(width_a, 1),
                   labels = c("A", "B"), align = "h", common.legend = TRUE, legend = "bottom")
  annotate_figure(fig, bottom = text_grob(paste0(
    "Simulated study design (33 subjects, two visits, single-wave truth; ka and kel are the baseline values in 1/min).\n",
    "A: share of replicates (95% Wilson CI). B: median and 2.5-97.5 percentiles over replicates; ",
    "the true log change in ka is 0.41 for ka x 1.5 and 0 otherwise."), size = 7, face = "italic", hjust = 1, x = 0.99))
}
save_sim <- function(p, name, height) {
  ggsave(file.path(SIM_DIR, paste0(name, ".pdf")), p, width = 9, height = height, bg = "white")
  ggsave(file.path(SIM_DIR, paste0(name, ".png")), p, width = 9, height = height, dpi = 200, bg = "white")
  cat("Saved", file.path(SIM_DIR, paste0(name, ".pdf")), "\n")
}

# All scenarios, including a kel change that violates the shared-kel assumption
save_sim(plot_sim(d, width_a = 1.7), "simulation_figure", 5.6)
# Only the scenarios the shared-kel model is built for: no change and a ka change
save_sim(plot_sim(d %>% filter(!grepl("^kel", scenario)) %>% mutate(group = droplevels(group))), "simulation_figure_ka", 4.4)
write_csv(d %>% select(panel, group, scenario, fit, est, lo, hi), file.path(SIM_DIR, "simulation_figure_data.csv"))
