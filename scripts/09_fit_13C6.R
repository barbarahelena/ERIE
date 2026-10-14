# =============================================================================
# 09_fit_13C6.R - the 13C6 tracer curve with the rates of the 12C fit held fixed
#
# For each completer x visit, ka (that visit) and kel (shared per subject) come
# from 02_fit_shared_kel.R and are not refitted; 12C's mechanism and t_lag come
# from the independent fit in results/01_independent_fit/. Only the parameters
# specific to the tracer capsule are fitted. Three 13C6 mechanisms, as in the
# archived joint model (scripts/archive/joint_12C_13C6/02_fit_erie_model.R):
#   release      capsule dissolution at rate k_release (K=2). When 12C is two_wave,
#                13C6 gets its own f_delayed_13C6 at 12C's t_lag (K=3).
#   onset_lag    instant bolus starting at t_lag1_13C6 in [0, 90] (K=2); tried
#                when the first sample is below 25% of the peak.
#   second_wave  instant bolus with its own late wave, f_delayed2_13C6 in
#                [0.05, 0.95] and t_lag2_13C6 in [30, 150] (K=3); tried when 13C6's
#                own data show a dip or plateau, rejected at the f lower bound.
# The mechanism with the lowest AIC on 13C6's raw residuals is kept. Same
# objective and penalties as 02. 13C6 has no endogenous level, so no baseline
# correction. Every result is conditional on the stable-kel assumption of 02.
#
# Outputs (results/04_13C6_fixed_rates/):
#   fits_13C6.csv            per subject x visit: fixed rates, mechanism, F_13C6, R2, AIC per candidate
#   plots_individual/*.pdf   12C (shared-kel fit) and 13C6 per visit
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2) })
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_fit.R")
source("scripts/assets/pk_diagnostics.R")
source("scripts/assets/pred_13C6.R")   # pred_13C6(): the prediction per mechanism

OUT_DIR <- "results/04_13C6_fixed_rates"; dir.create(file.path(OUT_DIR, "plots_individual"), showWarnings = FALSE, recursive = TRUE)
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = 1))
MIN_TMAX <- 10; CMAX_TOL <- 0.10; CMAX_LAMBDA <- 20; TMAX_LAMBDA <- 50
BOUNDS_RELEASE <- list(F_13C6 = c(1e-4, 1), k_release = c(0.001, 1), f_delayed_13C6 = c(0.001, 0.999))
BOUNDS_ONSET   <- list(F_13C6 = c(1e-4, 1), t_lag1_13C6 = c(0, 90))
BOUNDS_WAVE    <- list(F_13C6 = c(1e-4, 1), f_delayed2_13C6 = c(0.05, 0.95), t_lag2_13C6 = c(30, 150))
N_RANDOM <- 20
FINE_T <- seq(0, 400, length.out = 50)

cov   <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
const <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
MW_12C <- const$value[const$constant == "MW_12C"]; MW_13C6 <- const$value[const$constant == "MW_13C6"]
conc_all <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(!is.na(conc_umol_L)) %>%
  mutate(conc_mgL = conc_umol_L * if_else(isotope == "12C", MW_12C, MW_13C6) / 1000) %>%
  group_by(subject_id, visit, isotope) %>%
  mutate(conc_mgL = if_else(isotope == "12C", conc_mgL - conc_mgL[time_min == 0][1], conc_mgL)) %>% ungroup()
conc13 <- conc_all %>% filter(isotope == "13C6", time_min > 0)

shared <- read_csv("results/02_shared_kel/shared_kel_fits.csv", show_col_types = FALSE)
indep  <- read_csv("results/01_independent_fit/fit_results.csv", show_col_types = FALSE) %>%
  select(subject_id, visit, model_12C = model, f_delayed_12C = f_delayed, t_lag_12C = t_lag)
rates <- bind_rows(
  shared %>% transmute(subject_id, visit = "baseline", ka = ka_b, kel = kel_shared, F_12C = F_b, r2_12C = r2_b),
  shared %>% transmute(subject_id, visit = "intervention", ka = ka_i, kel = kel_shared, F_12C = F_i, r2_12C = r2_i)) %>%
  filter(!is.na(ka)) %>% left_join(indep, by = c("subject_id", "visit")) %>% arrange(subject_id, visit)

fit_curve_13C6 <- function(i) {
  r <- rates[i, ]; sid <- r$subject_id; vis <- r$visit
  set.seed(string_seed(paste(sid, vis, "13C6")))
  o <- conc13 %>% filter(subject_id == sid, visit == vis) %>% arrange(time_min)
  cv <- cov %>% filter(subject_id == sid, visit == vis)
  dose <- cv$dose_13C6_mg; Vd <- cv$vd_L; n <- nrow(o)
  if (n < 3) return(tibble(subject_id = sid, visit = vis))

  run <- function(mech, bounds, seeds) {
    curve <- list(times = o$time_min, conc = o$conc_mgL,
                  simulate = function(th) pred_13C6(mech, o$time_min, th, r, dose, Vd),
                  fine_simulate = function(th) list(time = FINE_T, conc = pred_13C6(mech, FINE_T, th, r, dose, Vd)))
    objective <- build_joint_objective(list(curve), min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                       tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
    fit <- fit_multistart(objective, vapply(bounds, `[`, numeric(1), 1), vapply(bounds, `[`, numeric(1), 2),
                          seeds = seeds, control = list(maxit = 60))
    if (is.null(fit)) return(NULL)
    rss <- sum((o$conc_mgL - curve$simulate(fit$par))^2)
    list(mech = mech, par = fit$par, rss = rss, aic = n * log(rss / n) + 2 * length(bounds), converged = fit$convergence == 0)
  }

  F0 <- r$F_12C
  two_wave <- r$model_12C == "two_wave"
  b_rel <- if (two_wave) BOUNDS_RELEASE else BOUNDS_RELEASE[c("F_13C6", "k_release")]
  seeds_rel <- c(
    unlist(lapply(c(0.01, 0.03, 0.1, 0.3), function(k) lapply(c(F0, 0.02), function(f)
      c(F_13C6 = f, k_release = k, f_delayed_13C6 = 0.3)[names(b_rel)])), recursive = FALSE),
    if (two_wave) list(c(F_13C6 = F0, k_release = 0.05, f_delayed_13C6 = 0.01), c(F_13C6 = F0, k_release = 0.05, f_delayed_13C6 = 0.99)),
    random_seeds(N_RANDOM, b_rel, log_scale = "k_release"))
  cands <- list(release = run("release", b_rel, seeds_rel))
  if (has_onset_lag_evidence(o$time_min, o$conc_mgL))
    cands$onset_lag <- run("onset_lag", BOUNDS_ONSET, lapply(c(10, 20, 30, 40, 50), function(tl) c(F_13C6 = F0, t_lag1_13C6 = tl)))
  if (has_peak_dip_rise(o$time_min, o$conc_mgL) || has_near_peak_neighbor(o$time_min, o$conc_mgL)) {
    sw <- run("second_wave", BOUNDS_WAVE, lapply(c(40, 60, 90, 120), function(tl) c(F_13C6 = F0, f_delayed2_13C6 = 0.5, t_lag2_13C6 = tl)))
    # pinned at the lower bound of f_delayed2_13C6 means no real second wave
    if (!is.null(sw) && sw$par[["f_delayed2_13C6"]] > BOUNDS_WAVE$f_delayed2_13C6[1] + 1e-3) cands$second_wave <- sw
  }
  cands <- compact(cands)
  if (length(cands) == 0) return(tibble(subject_id = sid, visit = vis))
  best <- cands[[which.min(map_dbl(cands, "aic"))]]
  p <- function(nm) if (nm %in% names(best$par)) best$par[[nm]] else NA_real_
  aic_of <- function(m) if (m %in% names(cands)) cands[[m]]$aic else NA_real_

  tibble(subject_id = sid, visit = vis, model_12C = r$model_12C, ka = r$ka, kel = r$kel, t_lag_12C = r$t_lag_12C,
         F_12C = r$F_12C, r2_12C = r$r2_12C, mechanism_13C6 = best$mech,
         F_13C6 = p("F_13C6"), k_release = p("k_release"), f_delayed_13C6 = p("f_delayed_13C6"), t_lag1_13C6 = p("t_lag1_13C6"),
         f_delayed2_13C6 = p("f_delayed2_13C6"), t_lag2_13C6 = p("t_lag2_13C6"),
         r2_13C6 = 1 - best$rss / sum((o$conc_mgL - mean(o$conc_mgL))^2), rss_13C6 = best$rss, n_13C6 = n,
         k_release_at_bound = isTRUE(p("k_release") > BOUNDS_RELEASE$k_release[2] - 1e-4), converged = best$converged,
         aic_release = aic_of("release"), aic_onset_lag = aic_of("onset_lag"), aic_second_wave = aic_of("second_wave"))
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
    f <- indep %>% filter(subject_id == r$subject_id, visit == r$visit)
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
    theme_bw()
  ggsave(file.path(OUT_DIR, "plots_individual", sprintf("%s_13C6.pdf", sid)), p, width = 8, height = 5)
}
cat("Done:", file.path(OUT_DIR, "fits_13C6.csv"), "and plots_individual/\n")
