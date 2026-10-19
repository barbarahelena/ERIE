# =============================================================================
# 02_fit_shared_kel.R - 12C fitted with one kel per subject, and ka, F and the wave structure per visit
# Model and assumptions: docs/pk-model.md, "Why one kel per subject" and "The post-peak dip and the lagged-dose model".
# Output (results/01_shared_kel/): shared_kel_fits.csv, shared_kel_wilcoxon.csv, plots_individual/
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2); library(grid); library(ggthemes) })
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_fit.R")
source("scripts/assets/theme_publication.R")
source("scripts/assets/pk_diagnostics.R")

OUT_DIR <- "results/01_shared_kel"; dir.create(file.path(OUT_DIR, "plots_individual"), showWarnings = FALSE, recursive = TRUE)
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = 1))
MIN_TMAX <- 10; CMAX_TOL <- 0.10; CMAX_LAMBDA <- 20; TMAX_LAMBDA <- 50
KEL <- log(2) / c(70, 7); KA <- c(1e-4, 1); FBIO <- c(1e-4, 1); F_DELAYED <- c(0.001, 0.999); T_LAG <- c(5, 150)
EARLY_LAG_OK_MGL <- 5             # t=30 at or below this: t_lag may go down to 5 min, otherwise 30
R2_ALWAYS_TRY_TWO_WAVE <- 0.50    # single_wave R2 below this: two_wave is tried without raw-data evidence
EXCESS_DEFINITE_TWO_WAVE <- 0.90  # peak_dip_rise_info()$excess at or above this: two_wave without the AIC comparison
N_RANDOM <- 30; MAXIT <- 80
FINE_T <- seq(0, 400, length.out = 50)
RATE_SEEDS <- c(list(c(ka = 0.03, kel = 0.02, F = 0.3), c(ka = 0.08, kel = 0.05, F = 0.15),
                     c(ka = 0.01, kel = 0.01, F = 0.5), c(ka = 0.05, kel = 0.08, F = 0.2)),
                lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) c(ka = v, kel = v, F = 0.03)))
WAVE_SEEDS <- list(c(0.01, 60), c(0.01, 90), c(0.3, 60), c(0.4, 45), c(0.5, 60), c(0.35, 75), c(0.45, 90))   # (f_delayed, t_lag)

cov  <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
MW   <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE) %>% filter(constant == "MW_12C") %>% pull(value)
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", !is.na(conc_umol_L)) %>% mutate(conc_mgL = conc_umol_L * MW / 1000) %>%
  group_by(subject_id, visit) %>% mutate(conc_mgL = conc_mgL - conc_mgL[time_min == 0][1]) %>% ungroup() %>%
  filter(time_min > 0)
completers <- cov %>% filter(!is.na(bw_kg)) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)

# One visit: data, dose, Vd and the raw-data evidence for a second wave
visit_data <- function(sid, vis) {
  o <- conc %>% filter(subject_id == sid, visit == vis) %>% arrange(time_min)
  cv <- cov %>% filter(subject_id == sid, visit == vis)
  dip <- peak_dip_rise_info(o$time_min, o$conc_mgL)
  t30 <- o$conc_mgL[o$time_min == 30]
  list(times = o$time_min, conc = o$conc_mgL, dose = cv$dose_12C_mg, Vd = cv$vd_L, t30_present = length(t30) == 1,
       evidence = dip$detected || has_near_peak_neighbor(o$time_min, o$conc_mgL),
       definite = isTRUE(dip$excess >= EXCESS_DEFINITE_TWO_WAVE), trigger_time = dip$trigger_time,
       t_lag_lower = if (length(t30) == 1 && t30 > EARLY_LAG_OK_MGL) 30 else T_LAG[1])
}

# Parameters are named by visit suffix (ka_b, F_i, t_lag_b, ...); kel has no suffix and is shared by every visit fitted together
pred_visit <- function(v, s, two_wave, th, t) {
  simA <- function(tt, dose) bateman_conc(tt, th[[paste0("ka_", s)]], th[["kel"]], th[[paste0("F_", s)]], dose, v$Vd)
  if (two_wave) simulate_lagged_dose(simA, t, v$dose, th[[paste0("f_delayed_", s)]], th[[paste0("t_lag_", s)]]) else simA(t, v$dose)
}

param_bounds <- function(vs, waves) {
  b <- list(kel = KEL)
  for (s in names(vs)) {
    b[[paste0("ka_", s)]] <- KA; b[[paste0("F_", s)]] <- FBIO
    if (waves[[s]]) { b[[paste0("f_delayed_", s)]] <- F_DELAYED; b[[paste0("t_lag_", s)]] <- c(vs[[s]]$t_lag_lower, T_LAG[2]) }
  }
  b
}

# Starting points: the fixed and ka = kel rate seeds, each with one wave seed (cycled), every `nested` fit with
# each wave seed, and random seeds. The wave seeds include the trigger time of the dip and 20 and 40 min before it.
make_seeds <- function(vs, waves, nested = list()) {
  wave_opts <- lapply(names(vs), function(s) {
    tt <- vs[[s]]$trigger_time
    c(WAVE_SEEDS, if (!is.na(tt)) list(c(0.4, tt), c(0.3, tt - 20), c(0.35, tt - 40)))
  }) %>% setNames(names(vs))
  add_waves <- function(p, j) {
    for (s in names(vs)) if (waves[[s]] && !(paste0("t_lag_", s) %in% names(p))) {
      w <- wave_opts[[s]]; p[paste0(c("f_delayed_", "t_lag_"), s)] <- w[[(j - 1) %% length(w) + 1]]
    }
    p
  }
  rates <- lapply(RATE_SEEDS, function(r) {
    p <- c(kel = r[["kel"]])
    for (s in names(vs)) p[paste0(c("ka_", "F_"), s)] <- r[c("ka", "F")]
    p
  })
  n_wave <- max(lengths(wave_opts))
  bounds <- param_bounds(vs, waves)
  seeds <- c(imap(rates, add_waves),
             unlist(lapply(nested, function(p) lapply(seq_len(n_wave), function(j) add_waves(p, j))), recursive = FALSE),
             random_seeds(N_RANDOM, bounds, log_scale = c("kel", paste0("ka_", names(vs)))))
  unique(lapply(seeds, function(p) p[names(bounds)]))
}

# Fit the visits in `vs` (a list named by suffix) together; `waves` says per visit whether it is two_wave
fit_visits <- function(vs, waves, seeds) {
  curves <- imap(vs, function(v, s) list(times = v$times, conc = v$conc,
    simulate = function(th) pred_visit(v, s, waves[[s]], th, v$times),
    fine_simulate = function(th) list(time = FINE_T, conc = pred_visit(v, s, waves[[s]], th, FINE_T))))
  joint <- build_joint_objective(curves, min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                 tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
  # The Tmax floor of the combined curve only reaches the taller wave. Both waves have the shape of the
  # plain Bateman curve, so the floor on that shape applies it to each wave.
  two <- names(waves)[unlist(waves)]
  objective <- function(th) {
    total <- joint(th)
    for (s in two) {
      shape <- bateman_conc(FINE_T, th[[paste0("ka_", s)]], th[["kel"]], 1, 1, 1)
      total <- total + TMAX_LAMBDA * max(0, MIN_TMAX - FINE_T[which.max(shape)])^2
    }
    total
  }
  bounds <- param_bounds(vs, waves)
  fit <- fit_multistart(objective, map_dbl(bounds, 1), map_dbl(bounds, 2), seeds = seeds, control = list(maxit = MAXIT))
  if (is.null(fit)) return(NULL)
  list(par = fit$par, waves = waves, converged = fit$convergence == 0,
       rss = map_dbl(curves, ~ sum((.x$conc - .x$simulate(fit$par))^2)), n = map_int(curves, ~ length(.x$conc)))
}

# AIC with a residual variance per visit: each visit is judged on its own residuals, so the noisier visit cannot
# hide a second wave at the other one
aic_of <- function(f) sum(f$n * log(f$rss / f$n)) + 2 * length(f$par)

# One figure per subject: observed points (t=0 included), the shared fit (solid) and each visit fitted alone (dashed)
plot_pair <- function(vs, best, indep, aic, row) {
  vis <- c(b = "baseline", i = "intervention"); t_plot <- seq(0, 360, by = 2)
  line <- function(f, s, fit) if (is.null(f)) NULL else
    tibble(visit = vis[[s]], fit = fit, time_min = t_plot, conc_mgL = pred_visit(vs[[s]], s, f$waves[[s]], f$par, t_plot))
  lines <- bind_rows(lapply(names(vis), function(s) bind_rows(line(best, s, "one kel per subject"), line(indep[[s]], s, "visit fitted alone"))))
  obs <- bind_rows(lapply(names(vis), function(s) tibble(visit = vis[[s]], time_min = c(0, vs[[s]]$times), conc_mgL = c(0, vs[[s]]$conc))))
  ann <- tibble(visit = vis, label = c(sprintf("%s, R2=%.3f\nka=%.4f", row$model_b, row$r2_b, row$ka_b),
                                       sprintf("%s, R2=%.3f\nka=%.4f", row$model_i, row$r2_i, row$ka_i)))
  p <- ggplot(obs, aes(time_min, conc_mgL)) +
    geom_line(data = lines, aes(linetype = fit, color = fit), linewidth = 0.7) + geom_point(size = 1.4) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = label), hjust = 1.05, vjust = 1.3, size = 2.5) +
    facet_wrap(vars(visit)) + scale_linetype_manual(values = c("one kel per subject" = "solid", "visit fitted alone" = "dashed"), name = NULL) +
    scale_color_manual(values = c("one kel per subject" = "steelblue", "visit fitted alone" = "grey50"), name = NULL) +
    guides(linetype = guide_legend(keywidth = unit(1.2, "cm")), color = guide_legend(keywidth = unit(1.2, "cm"))) +
    labs(title = row$subject_id, x = "Time (min)", y = "Concentration (mg/L, baseline-corrected)",
         subtitle = sprintf("kel = %.4f /min (t1/2 %.0f min) | AIC %s", row$kel_shared, log(2) / row$kel_shared,
                            paste(names(aic), sprintf("%.1f", aic), collapse = ", "))) +
    theme_Publication() + theme(plot.subtitle = element_text(size = rel(0.55)))
  ggsave(file.path(OUT_DIR, "plots_individual", sprintf("%s_shared_kel.pdf", row$subject_id)), p, width = 8, height = 4.5)
}

fit_pair <- function(sid) {
  vs <- list(b = visit_data(sid, "baseline"), i = visit_data(sid, "intervention"))
  run <- function(b, i, nested = list()) {
    set.seed(string_seed(paste(sid, b, i)))
    fit_visits(vs, list(b = b, i = i), make_seeds(vs, list(b = b, i = i), nested))
  }
  sw_sw <- run(FALSE, FALSE)
  if (is.null(sw_sw)) return(tibble(subject_id = sid))
  # two_wave is tried where t=30 exists and the raw data show a second wave or single_wave fits very poorly
  r2_sw <- map2_dbl(vs, sw_sw$rss, ~ 1 - .y / sum((.x$conc - mean(.x$conc))^2))
  try2 <- map2_lgl(vs, r2_sw, ~ .x$t30_present && (.x$evidence || .y < R2_ALWAYS_TRY_TWO_WAVE))
  fits <- list(sw_sw = sw_sw)
  if (try2[["b"]]) fits$tw_sw <- run(TRUE, FALSE, list(sw_sw$par))
  if (try2[["i"]]) fits$sw_tw <- run(FALSE, TRUE, list(sw_sw$par))
  if (all(try2)) fits$tw_tw <- run(TRUE, TRUE, compact(list(sw_sw$par, fits$tw_sw$par, fits$sw_tw$par)))
  fits <- compact(fits)
  aic <- map_dbl(fits, aic_of)
  # A visit with definite raw-data evidence keeps two_wave, unless no such fit succeeded
  definite <- map_lgl(vs, "definite")
  ok <- map_lgl(fits, ~ all(unlist(.x$waves)[definite]))
  if (!any(ok)) ok[] <- TRUE
  best <- fits[ok][[which.min(aic[ok])]]
  p <- best$par; w <- best$waves

  # Each visit fitted alone with the selected structure, for the AIC comparison. Seeds include the shared
  # solution and the same with ka and kel swapped (same curve, F scaled by ka / kel).
  indep <- imap(vs, function(v, s) {
    swapped <- p; ka <- p[[paste0("ka_", s)]]
    swapped[c("kel", paste0("ka_", s), paste0("F_", s))] <- c(ka, p[["kel"]], p[[paste0("F_", s)]] * ka / p[["kel"]])
    set.seed(string_seed(paste(sid, s, "alone")))
    fit_visits(vs[s], w[s], make_seeds(vs[s], w[s], list(p, swapped)))
  })
  get <- function(f, nm) if (is.null(f) || !(nm %in% names(f$par))) NA_real_ else f$par[[nm]]
  r2 <- function(s) 1 - best$rss[[s]] / sum((vs[[s]]$conc - mean(vs[[s]]$conc))^2)
  N <- sum(best$n)
  rss_indep <- map_dbl(indep, ~ if (is.null(.x)) NA_real_ else sum(.x$rss))
  row <- tibble(subject_id = sid, model_b = if (w$b) "two_wave" else "single_wave", model_i = if (w$i) "two_wave" else "single_wave",
         kel_shared = p[["kel"]], ka_b = p[["ka_b"]], ka_i = p[["ka_i"]], F_b = p[["F_b"]], F_i = p[["F_i"]],
         f_delayed_b = get(best, "f_delayed_b"), t_lag_b = get(best, "t_lag_b"),
         f_delayed_i = get(best, "f_delayed_i"), t_lag_i = get(best, "t_lag_i"),
         r2_b = r2("b"), r2_i = r2("i"), kel_at_bound = p[["kel"]] < KEL[1] + 1e-4 | p[["kel"]] > KEL[2] - 1e-4,
         converged = best$converged, two_wave_tried_b = try2[["b"]], two_wave_tried_i = try2[["i"]],
         ka_b_indep = get(indep$b, "ka_b"), ka_i_indep = get(indep$i, "ka_i"),
         kel_b_indep = get(indep$b, "kel"), kel_i_indep = get(indep$i, "kel"),
         n_b = best$n[["b"]], n_i = best$n[["i"]], rss_b = best$rss[["b"]], rss_i = best$rss[["i"]],
         rss_shared = sum(best$rss), rss_indep = sum(rss_indep), n_obs = N,
         delta_aic = aic_of(best) - (sum(best$n * log(rss_indep / best$n)) + 2 * (length(p) + 1)),
         ka_change_log = log(p[["ka_i"]] / p[["ka_b"]]),
         aic_sw_sw = aic[["sw_sw"]], aic_tw_sw = unname(aic["tw_sw"]), aic_sw_tw = unname(aic["sw_tw"]), aic_tw_tw = unname(aic["tw_tw"]))
  plot_pair(vs, best, indep, aic, row)
  row
}

res <- bind_rows(parallel::mclapply(completers, fit_pair, mc.cores = N_CORES)) %>%
  left_join(cov %>% distinct(subject_id, diet), by = "subject_id")
write_csv(res, file.path(OUT_DIR, "shared_kel_fits.csv"))

cat("\nSubjects:", nrow(res), "| two_wave curves:", sum(c(res$model_b, res$model_i) == "two_wave", na.rm = TRUE), "of", 2 * nrow(res),
    "(tried in", sum(c(res$two_wave_tried_b, res$two_wave_tried_i), na.rm = TRUE), ")\n")
cat("delta_AIC (shared - independent) median", round(median(res$delta_aic, na.rm = TRUE), 2),
    "| shared kel preferred (delta_AIC < 0):", sum(res$delta_aic < 0, na.rm = TRUE), "| clearly worse (> 2):", sum(res$delta_aic > 2, na.rm = TRUE), "\n")
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
