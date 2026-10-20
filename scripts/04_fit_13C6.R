# =============================================================================
# 04_fit_13C6.R - 13C6 fitted with ka and kel fixed at the 12C shared-kel fit
# The model (mechanisms, bounds, selection) is in assets/model_13C6.R; docs/pk-model.md, "The 13C6 tracer with fixed rates".
# Output (results/03_13C6_fixed_rates/): fits_13C6.csv, plots_individual/
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2); library(grid); library(ggthemes) })
source("scripts/assets/model_13C6.R")   # pred_13C6() and fit_13C6()
source("scripts/assets/theme_publication.R")

OUT_DIR <- "results/03_13C6_fixed_rates"; dir.create(file.path(OUT_DIR, "plots_individual"), showWarnings = FALSE, recursive = TRUE)
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = 1))

cov   <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
const <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
MW_12C <- const$value[const$constant == "MW_12C"]; MW_13C6 <- const$value[const$constant == "MW_13C6"]
conc_all <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(!is.na(conc_umol_L)) %>%
  mutate(conc_mgL = conc_umol_L * if_else(isotope == "12C", MW_12C, MW_13C6) / 1000) %>%
  group_by(subject_id, visit, isotope) %>%
  mutate(conc_mgL = if_else(isotope == "12C", conc_mgL - conc_mgL[time_min == 0][1], conc_mgL)) %>% ungroup()
conc13 <- conc_all %>% filter(isotope == "13C6", time_min > 0)

shared <- read_csv("results/01_shared_kel/shared_kel_fits.csv", show_col_types = FALSE)
rates <- bind_rows(
  shared %>% transmute(subject_id, visit = "baseline", model_12C = model_b, ka = ka_b, kel = kel_shared, F_12C = F_b,
                       f_delayed_12C = f_delayed_b, t_lag_12C = t_lag_b, r2_12C = r2_b),
  shared %>% transmute(subject_id, visit = "intervention", model_12C = model_i, ka = ka_i, kel = kel_shared, F_12C = F_i,
                       f_delayed_12C = f_delayed_i, t_lag_12C = t_lag_i, r2_12C = r2_i)) %>%
  filter(!is.na(ka)) %>% arrange(subject_id, visit)

fit_curve_13C6 <- function(i) {
  r <- rates[i, ]; sid <- r$subject_id; vis <- r$visit
  o <- conc13 %>% filter(subject_id == sid, visit == vis) %>% arrange(time_min)
  cv <- cov %>% filter(subject_id == sid, visit == vis)
  f <- fit_13C6(paste(sid, vis, "13C6"), o$time_min, o$conc_mgL, r, cv$dose_13C6_mg, cv$vd_L)
  if (is.null(f)) return(tibble(subject_id = sid, visit = vis))
  best <- f$best
  aic_mech <- function(m) if (m %in% names(f$cands)) f$cands[[m]]$aic else NA_real_

  tibble(subject_id = sid, visit = vis, model_12C = r$model_12C, ka = r$ka, kel = r$kel, t_lag_12C = r$t_lag_12C,
         F_12C = r$F_12C, r2_12C = r$r2_12C, mechanism_13C6 = best$mech,
         F_13C6 = par_or_na(best, "F_13C6"), k_release = par_or_na(best, "k_release"), f_delayed_13C6 = par_or_na(best, "f_delayed_13C6"),
         t_lag1_13C6 = par_or_na(best, "t_lag1_13C6"), f_delayed2_13C6 = par_or_na(best, "f_delayed2_13C6"),
         t_lag2_13C6 = par_or_na(best, "t_lag2_13C6"), r2_13C6 = best$r2, rss_13C6 = best$rss, n_13C6 = nrow(o),
         k_release_at_bound = isTRUE(par_or_na(best, "k_release") > BOUNDS_RELEASE$k_release[2] - 1e-4), converged = best$converged,
         aic_release = aic_mech("release"), aic_onset_lag = aic_mech("onset_lag"), aic_second_wave = aic_mech("second_wave"))
}

res <- bind_rows(parallel::mclapply(seq_len(nrow(rates)), fit_curve_13C6, mc.cores = N_CORES)) %>%
  left_join(cov %>% distinct(subject_id, diet), by = "subject_id")
write_csv(res, file.path(OUT_DIR, "fits_13C6.csv"))

cat("\nCurves:", nrow(res), "| mechanism:\n"); print(table(res$model_12C, res$mechanism_13C6, dnn = c("12C", "13C6")))
cat("13C6 R2 median", round(median(res$r2_13C6, na.rm = TRUE), 3), "| below 0.70:", sum(res$r2_13C6 < 0.70, na.rm = TRUE),
    "| median F_13C6", signif(median(res$F_13C6, na.rm = TRUE), 3), "vs F_12C", signif(median(res$F_12C, na.rm = TRUE), 3), "\n")

# Per-subject plots: observed points (t=0 included) and both fitted curves.
FINE_PLOT <- seq(0, 360, by = 2)
curve_rows <- function(r) {
  cv <- cov %>% filter(subject_id == r$subject_id, visit == r$visit)
  bolus12 <- function(tt, dose) bateman_conc(tt, r$ka, r$kel, r$F_12C, dose, cv$vd_L)
  c12 <- if (r$model_12C == "two_wave") {
    f <- rates %>% filter(subject_id == r$subject_id, visit == r$visit)
    simulate_lagged_dose(bolus12, FINE_PLOT, cv$dose_12C_mg, f$f_delayed_12C, f$t_lag_12C)
  } else bolus12(FINE_PLOT, cv$dose_12C_mg)
  th <- unlist(r[c("F_13C6", "k_release", "f_delayed_13C6", "t_lag1_13C6", "f_delayed2_13C6", "t_lag2_13C6")])
  c13 <- pred_13C6(r$mechanism_13C6, FINE_PLOT, th, r, cv$dose_13C6_mg, cv$vd_L)
  tibble(visit = r$visit, time_min = rep(FINE_PLOT, 2), conc_mgL = c(c12, c13), isotope = rep(c("12C", "13C6"), each = length(FINE_PLOT)))
}
for (sid in unique(res$subject_id)) {
  rs <- res %>% filter(subject_id == sid, !is.na(mechanism_13C6))
  if (nrow(rs) == 0) next
  sim <- bind_rows(lapply(seq_len(nrow(rs)), function(j) curve_rows(rs[j, ])))
  ann <- bind_rows(rs %>% transmute(visit, isotope = "12C", label = sprintf("%s, R2=%.3f", model_12C, r2_12C)),
                   rs %>% transmute(visit, isotope = "13C6", label = sprintf("%s, R2=%.3f", mechanism_13C6, r2_13C6)))
  p <- ggplot(conc_all %>% filter(subject_id == sid, visit %in% rs$visit), aes(time_min, conc_mgL)) +
    geom_point(size = 1.4) + geom_line(data = sim, aes(color = isotope), linewidth = 0.8, show.legend = FALSE) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3, size = 2.5) +
    facet_grid(isotope ~ visit, scales = "free_y") +
    scale_color_manual(values = c("12C" = "steelblue", "13C6" = "firebrick")) +
    labs(title = sid, subtitle = "ka and kel from the 12C shared-kel fit, held fixed for 13C6", x = "Time (min)", y = "Concentration (mg/L)") +
    theme_Publication()
  ggsave(file.path(OUT_DIR, "plots_individual", sprintf("%s_13C6.pdf", sid)), p, width = 8, height = 5)
}
cat("Done:", file.path(OUT_DIR, "fits_13C6.csv"), "and plots_individual/\n")
