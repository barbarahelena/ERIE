# =============================================================================
# simulations/05_plot_structure_examples.R - example simulated curves, generated as in 04_structure_selection.R,
# each fitted alone with single_wave and two_wave as in stage 1 of 02. One row per outcome of the structure
# choice (single and two waves recovered, a wave added, a wave missed), to judge by eye what the data can show.
# Env: ERIE_N_CORES. Output: results/02_diagnostics_ka_kel/structure_selection/structure_examples.pdf
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2); library(grid); library(ggthemes) })
source("scripts/assets/model_12C.R")
source("scripts/assets/theme_publication.R")

OUT_DIR <- "results/02_diagnostics_ka_kel/structure_selection"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = max(1, parallel::detectCores() - 2)))
N_CURVES <- 160; N_SHOW <- 4; SEED <- 21
KA_TRUE <- 0.040; KEL_TRUE <- 0.025   # the "null" scenario of 04

cov <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE) %>% filter(!is.na(bw_kg))
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", time_min > 0, !is.na(conc_umol_L))
design <- conc %>% arrange(time_min) %>% group_by(subject_id, visit) %>% summarise(times = list(time_min), .groups = "drop") %>%
  inner_join(cov %>% select(subject_id, visit, vd_L, dose_12C_mg), by = c("subject_id", "visit"))
visit_fits <- read_csv("results/01_shared_kel/visit_fits.csv", show_col_types = FALSE)
shared <- read_csv("results/01_shared_kel/shared_kel_fits.csv", show_col_types = FALSE)
NOISE_SD <- median(sqrt(c(shared$rss_b / shared$n_b, shared$rss_i / shared$n_i)), na.rm = TRUE)
P_TWO_WAVE <- mean(visit_fits$model == "two_wave")
WAVES <- visit_fits %>% filter(model == "two_wave") %>% select(f_delayed, t_lag)

one_curve <- function(k) {
  set.seed(string_seed(paste("example", SEED, k)))
  dd <- design[sample.int(nrow(design), 1), ]
  two <- runif(1) < P_TWO_WAVE; w <- WAVES[sample.int(nrow(WAVES), 1), ]
  th <- c(ka_b = KA_TRUE * exp(rnorm(1, 0, 0.2)), kel = KEL_TRUE * exp(rnorm(1, 0, 0.2)), F_b = 1,
          f_delayed_b = w$f_delayed, t_lag_b = w$t_lag)
  v0 <- list(Vd = dd$vd_L, dose = dd$dose_12C_mg); times <- dd$times[[1]]
  th[["F_b"]] <- 80 * exp(rnorm(1, 0, 0.2)) / max(pred_visit(v0, "b", two, th, seq(0, 400, by = 1)))
  vs <- list(b = make_visit(times, pred_visit(v0, "b", two, th, times) + rnorm(length(times), 0, NOISE_SD), dd$dose_12C_mg, dd$vd_L))
  # both candidates as in fit_alone(), with its seeds
  run <- function(tw, nested = list()) {
    set.seed(string_seed(paste(k, "b", tw, "alone")))
    fit_visits(vs, list(b = tw), make_seeds(vs, list(b = tw), nested))
  }
  sw <- run(FALSE)
  tried <- vs$b$t30_present && (vs$b$evidence || sw$r2[["b"]] < R2_ALWAYS_TRY_TWO_WAVE)
  tw <- if (tried) run(TRUE, list(sw$par)) else NULL
  aic_sw <- aic_of(sw); aic_tw <- if (is.null(tw)) NA_real_ else aic_of(tw)
  chosen_two <- !is.null(tw) && (vs$b$definite || aic_tw < aic_sw)
  t_plot <- seq(0, 360, by = 2)
  lines <- bind_rows(tibble(curve = "truth", time_min = t_plot, conc_mgL = pred_visit(v0, "b", two, th, t_plot)),
                     tibble(curve = "single_wave fit", time_min = t_plot, conc_mgL = pred_visit(vs$b, "b", FALSE, sw$par, t_plot)),
                     if (!is.null(tw)) tibble(curve = "two_wave fit", time_min = t_plot, conc_mgL = pred_visit(vs$b, "b", TRUE, tw$par, t_plot)))
  list(k = k, true_two = two, chosen_two = chosen_two, tried = tried, aic_sw = aic_sw, aic_tw = aic_tw,
       f_delayed = if (two) w$f_delayed else NA_real_, t_lag = if (two) w$t_lag else NA_real_,
       obs = tibble(time_min = c(0, times), conc_mgL = c(0, vs$b$conc)), lines = lines)
}

cat("Simulating and fitting", N_CURVES, "curves on", N_CORES, "cores\n")
res <- parallel::mclapply(seq_len(N_CURVES), one_curve, mc.cores = N_CORES, mc.preschedule = FALSE)
info <- map_dfr(res, ~ tibble(k = .x$k, true_two = .x$true_two, chosen_two = .x$chosen_two, tried = .x$tried,
                              aic_sw = .x$aic_sw, aic_tw = .x$aic_tw, f_delayed = .x$f_delayed, t_lag = .x$t_lag)) %>%
  mutate(outcome = case_when(true_two & chosen_two ~ "two waves, recovered", !true_two & !chosen_two ~ "one wave, recovered",
                             !true_two & chosen_two ~ "wave added", TRUE ~ "wave missed"),
         outcome = factor(outcome, c("one wave, recovered", "two waves, recovered", "wave added", "wave missed")))
cat("\nOutcomes:\n"); print(table(info$outcome))
cat("\nTrue waves: f_delayed and t_lag, recovered vs missed\n")
print(info %>% filter(true_two) %>% group_by(outcome) %>%
        summarise(n = n(), median_f_delayed = median(f_delayed), median_t_lag = median(t_lag), .groups = "drop"))

shown <- info %>% group_by(outcome) %>% slice_head(n = N_SHOW) %>% mutate(col = row_number()) %>% ungroup()
lab <- shown %>% transmute(k, outcome, col, label = sprintf("AIC one/two %.1f/%s%s", aic_sw, if_else(is.na(aic_tw), "-", sprintf("%.1f", aic_tw)),
                                                         if_else(true_two, sprintf("\ntrue f_delayed %.2f, t_lag %.0f", f_delayed, t_lag), "")))
obs <- map_dfr(res[shown$k], ~ mutate(.x$obs, k = .x$k)) %>% left_join(shown, by = "k")
lines <- map_dfr(res[shown$k], ~ mutate(.x$lines, k = .x$k)) %>% left_join(shown, by = "k")
p <- ggplot(obs, aes(time_min, conc_mgL)) +
  geom_line(data = lines, aes(color = curve, linetype = curve), linewidth = 0.6) + geom_point(size = 1.1) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = label), hjust = 1.03, vjust = 1.2, size = 2.1, lineheight = 0.9) +
  facet_grid(outcome ~ col, scales = "free_y", labeller = labeller(col = function(x) paste("example", x))) +
  scale_color_manual(values = c("truth" = "black", "single_wave fit" = "steelblue", "two_wave fit" = "#d95f02"), name = NULL) +
  scale_linetype_manual(values = c("truth" = "dashed", "single_wave fit" = "solid", "two_wave fit" = "solid"), name = NULL) +
  guides(color = guide_legend(keywidth = unit(1.2, "cm"))) +
  labs(title = "Simulated curves: structure chosen per visit (stage 1 of 02)", x = "Time (min)", y = "Concentration (mg/L)") +
  theme_Publication() + theme(strip.text.y = element_text(size = rel(0.6)))
ggsave(file.path(OUT_DIR, "structure_examples.pdf"), p, width = 11, height = 9)
cat("\nSaved", file.path(OUT_DIR, "structure_examples.pdf"), "\n")
