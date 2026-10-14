# =============================================================================
# 04_profile_ka_kel.R - diagnostic for flip-flop kinetics in the single_wave fits
#
# For every subject x visit that the independent fit selected single_wave for, profile the 12C
# residual sum of squares over a log-spaced (ka, kel) grid with F solved by
# least squares at every grid point (the model is linear in F). Compares the best
# point on the ka = kel diagonal with the best point off it, to show whether the
# diagonal solutions the optimizer returns are real minima or sit on a ridge.
#
# Because the Bateman curve is symmetric in (ka, kel) up to the scale of F, the
# profile is symmetric about the diagonal: an off-diagonal minimum comes as a
# mirrored pair. Output: results/03_diagnostics_ka_kel/profile_ka_kel.csv
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr) })
source("scripts/assets/pk_curves.R")

FIT_DIR <- "results/01_independent_fit"
OUT <- "results/03_diagnostics_ka_kel"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fits <- read_csv(file.path(FIT_DIR, "fit_results.csv"), show_col_types = FALSE) %>%
  filter(model == "single_wave")
conc <- read_csv("data/processed/erie_concentrations.csv", show_col_types = FALSE) %>% filter(isotope == "12C")
cov  <- read_csv("data/processed/erie_covariates.csv", show_col_types = FALSE)
MW   <- read_csv("data/processed/erie_constants.csv", show_col_types = FALSE)
MW   <- MW$value[MW$constant == "MW_12C"]

grid_k <- 10^seq(log10(0.003), log10(0.3), length.out = 120)   # same span for ka and kel

profile_one <- function(sid, vis) {
  o <- conc %>% filter(subject_id == sid, visit == vis) %>% arrange(time_min)
  y <- o$conc_umol_L * MW / 1000
  y <- y - y[o$time_min == 0]
  keep <- o$time_min > 0 & is.finite(y)
  t <- o$time_min[keep]; y <- y[keep]
  cv <- cov %>% filter(subject_id == sid, visit == vis)
  dose <- cv$dose_12C_mg; Vd <- cv$vd_L
  rss <- outer(grid_k, grid_k, Vectorize(function(ka, kel) {
    s <- bateman_conc(t, ka, kel, 1, dose, Vd)
    sum((y - s * sum(s * y) / sum(s * s))^2)     # F by least squares
  }))
  if (!any(is.finite(rss))) stop("non-finite profile for ", sid, " ", vis, ": n=", length(t), " Vd=", Vd, " dose=", dose)
  idx <- which(rss == min(rss, na.rm = TRUE), arr.ind = TRUE)[1, ]
  diag_rss <- min(diag(rss), na.rm = TRUE)
  tibble(subject_id = sid, visit = vis,
         rss_grid_min = min(rss), ka_grid = grid_k[idx[1]], kel_grid = grid_k[idx[2]],
         rss_diag_min = diag_rss,
         diag_excess = diag_rss / min(rss) - 1,          # 0 = the diagonal holds the global minimum
         sst = sum((y - mean(y))^2))
}

res <- bind_rows(lapply(seq_len(nrow(fits)), function(i) profile_one(fits$subject_id[i], fits$visit[i])))
res <- res %>% left_join(fits %>% select(subject_id, visit, ka, kel, r2), by = c("subject_id", "visit")) %>%
  mutate(r2_grid = 1 - rss_grid_min / sst, r2_diag = 1 - rss_diag_min / sst)
write_csv(res, file.path(OUT, "profile_ka_kel.csv"))
print(res %>% transmute(subject_id, visit, ka, kel, r2 = round(r2, 3),
                        ka_grid = round(ka_grid, 4), kel_grid = round(kel_grid, 4),
                        r2_grid = round(r2_grid, 3), r2_diag = round(r2_diag, 3)), n = Inf)
