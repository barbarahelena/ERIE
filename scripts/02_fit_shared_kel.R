# =============================================================================
# 02_fit_shared_kel.R - the real cohort fitted with one kel per subject
#
# CONDITIONAL analysis: assumes kel does not change between a subject's two
# visits. Each completer's visits are fitted together with a shared kel and a
# separate ka and F per visit. Each visit keeps the mechanism that
# scripts/archive/fructose_independent/02_fit_erie_model_fructose_independent.R selected (single_wave, or two_wave
# with its f_delayed and t_lag held at the independent estimates), so only
# ka, kel and F are re-estimated. Same objective and penalties as the independent fit; kel bounded to t½ 7-70 min.
#
# Outputs (results/02_shared_kel/):
#   shared_kel_fits.csv      per subject: shared kel, ka and F by visit, R2, RSS, AIC vs the independent fits
#   (03_diet_summary_shared_kel.R summarises these by diet arm)
#   shared_kel_wilcoxon.csv  paired Wilcoxon on ka (log ratio) and the independent-fit ka, per diet
# AIC = N ln(RSS/N) + 2k pooled over both visits (raw mg/L residuals). The shared
# model has one parameter fewer than the independent one; negative delta_aic
# favours the shared kel. Compatibility with a stable kel is not proof of it.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
source("scripts/assets/pk_curves.R")
source("scripts/assets/pk_fit.R")

RES <- "results/01_independent_fit"   # the independent fits this script starts from
OUT_DIR <- "results/02_shared_kel"; dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
N_CORES <- as.integer(Sys.getenv("ERIE_N_CORES", unset = 1))
MIN_TMAX <- 10; CMAX_TOL <- 0.10; CMAX_LAMBDA <- 20; TMAX_LAMBDA <- 50
BOUNDS <- list(kel = log(2) / c(70, 7), ka_b = c(1e-4, 1), ka_i = c(1e-4, 1), F_b = c(1e-4, 1), F_i = c(1e-4, 1))
N_RANDOM <- 25
FINE_T <- seq(0, 400, length.out = 50)

cov  <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
MW   <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE) %>% filter(constant == "MW_12C") %>% pull(value)
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>%
  filter(isotope == "12C", !is.na(conc_umol_L)) %>% mutate(conc_mgL = conc_umol_L * MW / 1000) %>%
  group_by(subject_id, visit) %>% mutate(conc_mgL = conc_mgL - conc_mgL[time_min == 0][1]) %>% ungroup() %>%
  filter(time_min > 0)
fits <- read_csv(file.path(RES, "fit_results.csv"), show_col_types = FALSE)
rss_tab <- bind_rows(single_wave = read_csv(file.path(RES, "fit_results_single_wave.csv"), show_col_types = FALSE) %>% select(subject_id, visit, rss, n),
                     two_wave    = read_csv(file.path(RES, "fit_results_two_wave.csv"), show_col_types = FALSE) %>% select(subject_id, visit, rss, n),
                     .id = "model")
fits <- fits %>% left_join(rss_tab, by = c("subject_id", "visit", "model"))
completers <- cov %>% filter(!is.na(bw_kg)) %>% count(subject_id) %>% filter(n == 2) %>% pull(subject_id)

visit_curve <- function(sid, vis, ka_name, F_name) {
  o <- conc %>% filter(subject_id == sid, visit == vis) %>% arrange(time_min)
  cv <- cov %>% filter(subject_id == sid, visit == vis)
  f <- fits %>% filter(subject_id == sid, visit == vis)
  dose <- cv$dose_12C_mg; Vd <- cv$vd_L
  pred <- function(t, th) {
    simA <- function(tt, dose) bateman_conc(tt, th[[ka_name]], th[["kel"]], th[[F_name]], dose, Vd)
    if (f$model == "two_wave") simulate_lagged_dose(simA, t, dose, f$f_delayed, f$t_lag) else simA(t, dose)
  }
  list(times = o$time_min, conc = o$conc_mgL,
       simulate = function(th) pred(o$time_min, th),
       fine_simulate = function(th) list(time = FINE_T, conc = pred(FINE_T, th)))
}

fit_pair <- function(sid) {
  set.seed(sum(utf8ToInt(sid)) + 7)
  cb <- visit_curve(sid, "baseline", "ka_b", "F_b"); ci <- visit_curve(sid, "intervention", "ka_i", "F_i")
  objective <- build_joint_objective(list(cb, ci), min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                     tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
  fb <- fits %>% filter(subject_id == sid, visit == "baseline"); fi <- fits %>% filter(subject_id == sid, visit == "intervention")
  sd <- function(kel, ka_b, ka_i, Fb = 0.1, Fi = Fb) c(kel = kel, ka_b = ka_b, ka_i = ka_i, F_b = Fb, F_i = Fi)
  seeds <- c(
    list(sd(sqrt(fb$kel * fi$kel), fb$ka, fi$ka, fb$F, fi$F),                 # independent fits, kel averaged
         sd(min(fb$kel, fb$ka), max(fb$kel, fb$ka), max(fi$kel, fi$ka), fb$F, fi$F),
         sd(min(fi$kel, fi$ka), max(fb$kel, fb$ka), max(fi$kel, fi$ka), fb$F, fi$F)),
    list(sd(0.02, 0.03, 0.03, 0.3), sd(0.05, 0.08, 0.08, 0.15), sd(0.01, 0.01, 0.01, 0.5), sd(0.08, 0.05, 0.05, 0.2)),
    lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) sd(v, v, v, 0.03)),
    random_seeds(N_RANDOM, BOUNDS, log_scale = c("kel", "ka_b", "ka_i")))
  fit <- fit_multistart(objective, vapply(BOUNDS, `[`, numeric(1), 1), vapply(BOUNDS, `[`, numeric(1), 2), seeds = seeds,
                        control = list(maxit = 60))
  if (is.null(fit)) return(tibble(subject_id = sid))
  p <- fit$par
  rss_b <- sum((cb$conc - cb$simulate(p))^2); rss_i <- sum((ci$conc - ci$simulate(p))^2)
  N <- length(cb$conc) + length(ci$conc)
  r2 <- function(cv, rss) 1 - rss / sum((cv$conc - mean(cv$conc))^2)
  k_ind <- sum(if_else(c(fb$model, fi$model) == "two_wave", 5, 3))
  aic_shared <- N * log((rss_b + rss_i) / N) + 2 * (k_ind - 1)
  aic_indep  <- N * log((fb$rss + fi$rss) / N) + 2 * k_ind
  tibble(subject_id = sid, model_b = fb$model, model_i = fi$model, kel_shared = p[["kel"]],
         ka_b = p[["ka_b"]], ka_i = p[["ka_i"]], F_b = p[["F_b"]], F_i = p[["F_i"]],
         r2_b = r2(cb, rss_b), r2_i = r2(ci, rss_i), kel_at_bound = p[["kel"]] < BOUNDS$kel[1] + 1e-4 | p[["kel"]] > BOUNDS$kel[2] - 1e-4,
         ka_b_indep = fb$ka, ka_i_indep = fi$ka, kel_b_indep = fb$kel, kel_i_indep = fi$kel,
         rss_shared = rss_b + rss_i, rss_indep = fb$rss + fi$rss, n_obs = N,
         delta_aic = aic_shared - aic_indep, ka_change_log = log(p[["ka_i"]] / p[["ka_b"]]))
}

res <- bind_rows(parallel::mclapply(completers, fit_pair, mc.cores = N_CORES)) %>%
  left_join(cov %>% distinct(subject_id, diet), by = "subject_id")
write_csv(res, file.path(OUT_DIR, "shared_kel_fits.csv"))

cat("\nSubjects:", nrow(res), "| delta_AIC (shared - independent) median", round(median(res$delta_aic, na.rm = TRUE), 2),
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
