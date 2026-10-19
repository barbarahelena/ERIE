# =============================================================================
# 05_simulate_shared_kel.R - simulated study design as in 04, fitted with one kel per subject
# Design, scenarios and results: docs/pk-model.md, "Why one kel per subject".
# Env: ERIE_SIM_REPS (default 50), ERIE_N_CORES. Output: results/02_diagnostics_ka_kel/simulation/
# =============================================================================
options(erie.sim.lib = TRUE)
source("scripts/04_simulate_ka_kel_recovery.R")

BOUNDS_SHARED <- list(kel = log(2) / c(70, 7), ka_b = c(1e-4, 1), ka_i = c(1e-4, 1), F_b = c(1e-4, 1), F_i = c(1e-4, 1))
N_RANDOM_SHARED <- 25   # more than the 15 per curve in 04: 5 parameters instead of 3

fit_shared <- function(times_b, y_b, times_i, y_i, dose_b, dose_i, Vd_b, Vd_i, seed) {
  set.seed(seed)
  mk <- function(times, y, ka_name, F_name, dose, Vd) list(times = times, conc = y,
    simulate      = function(th) bateman_conc(times, th[[ka_name]], th[["kel"]], th[[F_name]], dose, Vd),
    fine_simulate = function(th) list(time = FINE_T, conc = bateman_conc(FINE_T, th[[ka_name]], th[["kel"]], th[[F_name]], dose, Vd)))
  objective <- build_joint_objective(list(mk(times_b, y_b, "ka_b", "F_b", dose_b, Vd_b), mk(times_i, y_i, "ka_i", "F_i", dose_i, Vd_i)),
                                     min_tmax = MIN_TMAX, cmax_tol = CMAX_TOL, cmax_lambda = CMAX_LAMBDA,
                                     tmax_lambda = TMAX_LAMBDA, proportional_weighting = FALSE)
  sd <- function(kel, ka_b, ka_i, F = 0.1) c(kel = kel, ka_b = ka_b, ka_i = ka_i, F_b = F, F_i = F)
  seeds <- c(list(sd(0.02, 0.03, 0.03, 0.3), sd(0.05, 0.08, 0.08, 0.15), sd(0.01, 0.01, 0.01, 0.5),
                  sd(0.08, 0.05, 0.05, 0.2), sd(0.04, 0.02, 0.04, 0.2), sd(0.025, 0.04, 0.06, 0.1)),
             lapply(c(0.008, 0.012, 0.016, 0.02, 0.025, 0.03, 0.04, 0.05, 0.07), function(v) sd(v, v, v, 0.03)),
             random_seeds(N_RANDOM_SHARED, BOUNDS_SHARED, log_scale = c("kel", "ka_b", "ka_i")))
  fit <- fit_multistart(objective, vapply(BOUNDS_SHARED, `[`, numeric(1), 1), vapply(BOUNDS_SHARED, `[`, numeric(1), 2),
                        seeds = seeds, control = list(maxit = 60))
  if (is.null(fit)) return(c(kel = NA_real_, ka_b = NA_real_, ka_i = NA_real_))
  fit$par[c("kel", "ka_b", "ka_i")]
}

run_rep_shared <- function(sc, rep) {
  set.seed(SIM_SEED * 100000 + rep * 100 + match(sc$scenario, scenarios$scenario) + 50)
  subj <- tibble(subject_id = completers,
                 ka_s  = sc$ka  * exp(rnorm(length(completers), 0, 0.2)),
                 kel_s = sc$kel * exp(rnorm(length(completers), 0, 0.2)),
                 cmax_s = 80 * exp(rnorm(length(completers), 0, 0.2)))
  fits <- map_dfr(completers, function(sid) {
    s <- subj %>% filter(subject_id == sid)
    one <- function(vis) {
      dd <- design %>% filter(subject_id == sid, visit == vis)
      ka_t  <- s$ka_s  * if (vis == "intervention") sc$ka_x else 1
      kel_t <- s$kel_s * if (vis == "intervention") sc$kel_x else 1
      times <- dd$times[[1]]
      F_t <- s$cmax_s / max(bateman_conc(FINE_T, ka_t, kel_t, 1, dd$dose_12C_mg, dd$vd_L))
      list(times = times, dose = dd$dose_12C_mg, Vd = dd$vd_L, ka = ka_t, kel = kel_t,
           y = bateman_conc(times, ka_t, kel_t, F_t, dd$dose_12C_mg, dd$vd_L) + rnorm(length(times), 0, NOISE_SD))
    }
    b <- one("baseline"); i <- one("intervention")
    f <- fit_shared(b$times, b$y, i$times, i$y, b$dose, i$dose, b$Vd, i$Vd, seed = sample.int(1e6, 1))
    tibble(subject_id = sid, kel = f[["kel"]], ka_b = f[["ka_b"]], ka_i = f[["ka_i"]],
           ka_b_true = b$ka, ka_i_true = i$ka, kel_b_true = b$kel)
  })
  w <- fits %>% filter(is.finite(ka_b), is.finite(ka_i))
  p <- tryCatch(wilcox.test(log(w$ka_i / w$ka_b), mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
  sep <- abs(fits$ka_b_true - fits$kel_b_true) > 1e-9
  tibble(scenario = sc$scenario, rep = rep, p_ka = p,
         bias_log_ka = median(c(log(w$ka_b / w$ka_b_true), log(w$ka_i / w$ka_i_true))),
         bias_log_kel = median(log(w$kel / w$kel_b_true)),
         bias_ka_change = median(log(w$ka_i / w$ka_b)) - median(log(w$ka_i_true / w$ka_b_true)),
         swapped = if (any(sep)) mean(((fits$ka_b < fits$kel) != (fits$ka_b_true < fits$kel_b_true))[sep], na.rm = TRUE) else NA_real_,
         on_diagonal = mean(abs(w$ka_b - w$kel) / pmax(w$ka_b, w$kel) < 1e-3 | abs(w$ka_i - w$kel) / pmax(w$ka_i, w$kel) < 1e-3))
}

jobs <- expand_grid(s = seq_len(nrow(scenarios)), rep = seq_len(N_REPS))
cat("Running", nrow(jobs), "shared-kel replicates on", N_CORES, "cores (", length(completers), "subject fits each )\n")
t0 <- Sys.time()
out <- bind_rows(parallel::mclapply(seq_len(nrow(jobs)), function(i) run_rep_shared(scenarios[jobs$s[i], ], jobs$rep[i]),
                                    mc.cores = N_CORES))
cat("Elapsed:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
write_csv(out, file.path(OUT_DIR, "simulation_shared_kel_replicates.csv"))

summ <- out %>% group_by(scenario) %>%
  summarise(reps = n(), rej_ka = mean(p_ka < 0.05, na.rm = TRUE),
            bias_ka = median(bias_log_ka), bias_kel = median(bias_log_kel), bias_ka_change = median(bias_ka_change),
            swapped = mean(swapped, na.rm = TRUE), on_diagonal = mean(on_diagonal), .groups = "drop") %>%
  arrange(match(scenario, scenarios$scenario))
write_csv(summ, file.path(OUT_DIR, "simulation_shared_kel_summary.csv"))
cat("\n=== Shared-kel fit: ka rejection rate (p < 0.05), bias (median log fitted/true), label swaps ===\n")
print(as.data.frame(summ), digits = 2)
