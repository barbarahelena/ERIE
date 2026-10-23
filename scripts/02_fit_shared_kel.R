# =============================================================================
# 02_fit_shared_kel.R - 12C fitted with one kel per subject, and ka, F and the wave structure per visit
# Model and assumptions: docs/pk-model.md, "Why one kel per subject" and "The post-peak dip and the lagged-dose model".
# Stage 1 fits each visit alone and chooses its structure; stage 2 fits one kel per subject with that structure.
# The model (settings, stage 1, stage 2) is in assets/model_12C.R; this script reads the data, runs it and writes the output.
# Output (results/01_shared_kel/): visit_fits.csv, shared_kel_fits.csv, shared_kel_wilcoxon.csv, plots_individual/
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2); library(grid); library(ggthemes) })
source("scripts/assets/model_12C.R")   # make_visit(), pred_visit(), fit_alone() and fit_shared()
source("scripts/assets/theme_publication.R")

OUT_DIR <- "results/01_shared_kel"; dir.create(file.path(OUT_DIR, "plots_individual"), showWarnings = FALSE, recursive = TRUE)
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = 1))

cov  <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
MW   <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE) %>% filter(constant == "MW_12C") %>% pull(value)
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", !is.na(conc_umol_L)) %>% mutate(conc_mgL = conc_umol_L * MW / 1000) %>%
  group_by(subject_id, visit) %>% mutate(conc_mgL = conc_mgL - conc_mgL[time_min == 0][1]) %>% ungroup() %>%
  filter(time_min > 0)
completers <- cov %>% filter(!is.na(bw_kg)) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)

visit_data <- function(sid, vis) {
  o <- conc %>% filter(subject_id == sid, visit == vis) %>% arrange(time_min)
  cv <- cov %>% filter(subject_id == sid, visit == vis)
  make_visit(o$time_min, o$conc_mgL, cv$dose_12C_mg, cv$vd_L)
}

# One figure per subject: observed points (t=0 included), the shared fit (solid) and each visit fitted alone (dashed)
plot_pair <- function(vs, best, alone, row) {
  t_plot <- seq(0, 360, by = 2)
  line <- function(f, s, fit) if (is.null(f)) NULL else
    tibble(visit = VISITS[[s]], fit = fit, time_min = t_plot, conc_mgL = pred_visit(vs[[s]], s, f$waves[[s]], f$par, t_plot))
  lines <- bind_rows(lapply(names(VISITS), function(s) bind_rows(line(best, s, "one kel per subject"), line(alone[[s]], s, "visit fitted alone"))))
  obs <- bind_rows(lapply(names(VISITS), function(s) tibble(visit = VISITS[[s]], time_min = c(0, vs[[s]]$times), conc_mgL = c(0, vs[[s]]$conc))))
  ann <- tibble(visit = VISITS, label = c(sprintf("%s, R2=%.3f\nka=%.4f", row$model_b, row$r2_b, row$ka_b),
                                          sprintf("%s, R2=%.3f\nka=%.4f", row$model_i, row$r2_i, row$ka_i)))
  aic_txt <- paste(sprintf("%s %.1f/%s", VISITS, map_dbl(alone, ~ .x$aic[["single_wave"]]),
                           map_chr(alone, ~ sprintf("%.1f", .x$aic[["two_wave"]]))), collapse = ", ")
  p <- ggplot(obs, aes(time_min, conc_mgL)) +
    geom_line(data = lines, aes(linetype = fit, color = fit), linewidth = 0.7) + geom_point(size = 1.4) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3, size = 2.5) +
    facet_wrap(vars(visit)) + scale_linetype_manual(values = c("one kel per subject" = "solid", "visit fitted alone" = "dashed"), name = NULL) +
    scale_color_manual(values = c("one kel per subject" = "steelblue", "visit fitted alone" = "grey50"), name = NULL) +
    guides(linetype = guide_legend(keywidth = unit(1.2, "cm")), color = guide_legend(keywidth = unit(1.2, "cm"))) +
    labs(title = row$subject_id, x = "Time (min)", y = "Concentration (mg/L, baseline-corrected)",
         subtitle = sprintf("kel = %.4f /min (t1/2 %.0f min) | AIC alone single/two wave: %s", row$kel_shared, log(2) / row$kel_shared, aic_txt)) +
    theme_Publication() + theme(plot.subtitle = element_text(size = rel(0.55)))
  ggsave(file.path(OUT_DIR, "plots_individual", sprintf("%s_shared_kel.pdf", row$subject_id)), p, width = 8, height = 4.5)
}

fit_pair <- function(sid) {
  vs <- list(b = visit_data(sid, "baseline"), i = visit_data(sid, "intervention"))
  alone <- imap(vs, function(v, s) fit_alone(sid, vs, s))
  visit_rows <- bind_rows(imap(compact(alone), function(f, s) tibble(subject_id = sid, visit = VISITS[[s]],
    model = if (f$waves[[s]]) "two_wave" else "single_wave", ka = f$par[[paste0("ka_", s)]], kel = f$par[["kel"]],
    F = f$par[[paste0("F_", s)]], f_delayed = par_or_na(f, paste0("f_delayed_", s)), t_lag = par_or_na(f, paste0("t_lag_", s)),
    r2 = f$r2[[s]], kel_at_bound = f$par[["kel"]] < KEL[1] + 1e-4 | f$par[["kel"]] > KEL[2] - 1e-4,
    converged = f$converged, retried = f$retried, retry_improved = f$retry_improved, two_wave_tried = f$tried, definite_two_wave = vs[[s]]$definite,
    aic_single_wave = f$aic[["single_wave"]], aic_two_wave = f$aic[["two_wave"]], rss = f$rss[[s]], n = f$n[[s]])))
  if (length(compact(alone)) < 2) return(list(pair = tibble(subject_id = sid), visits = visit_rows))

  best <- fit_shared(sid, vs, alone)
  if (is.null(best)) return(list(pair = tibble(subject_id = sid), visits = visit_rows))
  p <- best$par
  row <- tibble(subject_id = sid, model_b = visit_rows$model[1], model_i = visit_rows$model[2],
         kel_shared = p[["kel"]], ka_b = p[["ka_b"]], ka_i = p[["ka_i"]], F_b = p[["F_b"]], F_i = p[["F_i"]],
         f_delayed_b = par_or_na(best, "f_delayed_b"), t_lag_b = par_or_na(best, "t_lag_b"),
         f_delayed_i = par_or_na(best, "f_delayed_i"), t_lag_i = par_or_na(best, "t_lag_i"),
         r2_b = best$r2[["b"]], r2_i = best$r2[["i"]],
         kel_at_bound = p[["kel"]] < KEL[1] + 1e-4 | p[["kel"]] > KEL[2] - 1e-4, converged = best$converged,
         retried = best$retried, retry_improved = best$retry_improved,
         ka_b_indep = visit_rows$ka[1], ka_i_indep = visit_rows$ka[2], kel_b_indep = visit_rows$kel[1], kel_i_indep = visit_rows$kel[2],
         n_b = best$n[["b"]], n_i = best$n[["i"]], rss_b = best$rss[["b"]], rss_i = best$rss[["i"]],
         rss_shared = sum(best$rss), rss_indep = sum(visit_rows$rss), n_obs = sum(best$n),
         delta_aic = aic_of(best) - aic_of(alone$b) - aic_of(alone$i),   # the separate fits have one parameter more
         ka_change_log = log(p[["ka_i"]] / p[["ka_b"]]))
  plot_pair(vs, best, alone, row)
  list(pair = row, visits = visit_rows)
}

out <- parallel::mclapply(completers, fit_pair, mc.cores = N_CORES)
diet <- cov %>% distinct(subject_id, diet)
res <- bind_rows(map(out, "pair")) %>% left_join(diet, by = "subject_id") %>% mutate(structure_rule = STRUCTURE_RULE)
visit_fits <- bind_rows(map(out, "visits")) %>% left_join(diet, by = "subject_id") %>% mutate(structure_rule = STRUCTURE_RULE)
write_csv(res, file.path(OUT_DIR, "shared_kel_fits.csv"))
write_csv(visit_fits, file.path(OUT_DIR, "visit_fits.csv"))

cat("\nStructure rule:", STRUCTURE_RULE, "| subjects:", nrow(res), "| two_wave curves:", sum(visit_fits$model == "two_wave"), "of", nrow(visit_fits),
    "(tried in", sum(visit_fits$two_wave_tried), ")\n")
cat("delta_AIC (shared - independent) median", round(median(res$delta_aic, na.rm = TRUE), 2),
    "| shared kel preferred (delta_AIC < 0):", sum(res$delta_aic < 0, na.rm = TRUE), "| clearly worse (> 2):", sum(res$delta_aic > 2, na.rm = TRUE), "\n")
cat("Retried: stage 1", sum(visit_fits$retried), "curves (improved", sum(visit_fits$retry_improved, na.rm = TRUE),
    ") | stage 2", sum(res$retried, na.rm = TRUE), "subjects (improved", sum(res$retry_improved, na.rm = TRUE), ")\n")
cat("Subjects with |log ka change| < 0.05 (kel/ka label unresolved):", sum(abs(res$ka_change_log) < 0.05, na.rm = TRUE), "\n")

wil <- bind_rows(
  res %>% mutate(delta = log(ka_i / ka_b), what = "ka (shared kel)"),
  res %>% mutate(delta = log(ka_i_indep / ka_b_indep), what = "ka (independent fits)"),
  res %>% mutate(delta = log(kel_i_indep / kel_b_indep), what = "kel (independent fits)")) %>%
  filter(is.finite(delta)) %>% group_by(what, diet) %>%
  summarise(n_paired = n(), median_delta = median(delta),
            p = tryCatch(wilcox.test(delta, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_), .groups = "drop")
write_csv(wil, file.path(OUT_DIR, "shared_kel_wilcoxon.csv"))
cat("\n=== Paired Wilcoxon on the log ratio, intervention vs baseline ===\n"); print(as.data.frame(wil), digits = 3)
