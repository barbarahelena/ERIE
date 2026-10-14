# =============================================================================
# pk_stats.R - paired-visit statistics for per-subject PK summaries
#
# Generic: knows subjects, two visits and two groups, nothing about fructose.
# =============================================================================

#' Random-intercept mixed model for a visit effect within and between groups
#'
#' Fits `y ~ group * visit + (1 | subject)` by REML with [nlme::lme()], on the
#' log scale or the raw scale. Reports, per group, the change from baseline to
#' intervention with a 95% CI and p-value (from the nested parameterisation
#' `0 + group + group:visit`), and the group x visit interaction (the difference
#' between the groups' changes) from the full model. Unlike the paired Wilcoxon,
#' the model keeps subjects with only one visit.
#'
#' @param d Data frame with `subject_id`, `diet` (the group), `visit`
#'   (`baseline`/`intervention`) and `value`.
#' @param log_scale TRUE: the model is on log(value) and estimates are log ratios.
#' @param group_levels Order of the groups; the interaction is the change of the
#'   second group minus the change of the first.
#' @return One row per group: `n_subjects`, `lmm_estimate`, `lmm_lo`,
#'   `lmm_hi`, `p_lmm`, plus `interaction_estimate` and `p_interaction`
#'   (repeated on both rows). All NA when the model cannot be fitted (for
#'   example a value that is identical at both visits for every subject).
mixed_visit_tests <- function(d, log_scale = TRUE, group_levels = c("low_fructose", "high_fructose")) {
  empty <- data.frame(diet = group_levels, n_subjects = NA_integer_, lmm_estimate = NA_real_, lmm_lo = NA_real_,
                      lmm_hi = NA_real_, p_lmm = NA_real_, interaction_estimate = NA_real_, p_interaction = NA_real_)
  d <- d[is.finite(d$value) & (!log_scale | d$value > 0) & d$diet %in% group_levels, , drop = FALSE]
  d$y <- if (log_scale) log(d$value) else d$value
  d$diet <- factor(d$diet, levels = group_levels)
  d$visit <- factor(d$visit, levels = c("baseline", "intervention"))
  d$subject_id <- factor(d$subject_id)
  if (nlevels(droplevels(d$diet)) < 2 || length(unique(d$subject_id)) < 6) return(empty)
  n_subj <- tapply(as.character(d$subject_id), d$diet, function(x) length(unique(x)))
  empty$n_subjects <- as.integer(n_subj[group_levels])

  # Nothing varies within a subject (for example one kel shared by both visits): the visit
  # effect is zero by construction and the model would return a meaningless p-value.
  wide <- tapply(d$y, list(d$subject_id, d$visit), mean)
  within <- wide[, "intervention"] - wide[, "baseline"]
  if (all(is.na(within)) || max(abs(within), na.rm = TRUE) < 1e-8) return(empty)

  fit <- function(formula) tryCatch(
    suppressWarnings(nlme::lme(formula, random = ~ 1 | subject_id, data = d, method = "REML",
                               control = nlme::lmeControl(opt = "optim"))),
    error = function(e) NULL)
  nested <- fit(y ~ 0 + diet + diet:visit)
  full   <- fit(y ~ diet * visit)
  if (is.null(nested) || is.null(full)) return(empty)

  tt <- summary(nested)$tTable
  out <- empty
  for (g in group_levels) {
    row <- paste0("diet", g, ":visitintervention")
    if (!row %in% rownames(tt)) next
    est <- tt[row, "Value"]; se <- tt[row, "Std.Error"]; df <- tt[row, "DF"]
    out$lmm_estimate[out$diet == g] <- est
    out$lmm_lo[out$diet == g] <- est - qt(0.975, df) * se
    out$lmm_hi[out$diet == g] <- est + qt(0.975, df) * se
    out$p_lmm[out$diet == g] <- tt[row, "p-value"]
  }
  tf <- summary(full)$tTable
  irow <- paste0("diet", group_levels[2], ":visitintervention")
  if (irow %in% rownames(tf)) {
    out$interaction_estimate <- tf[irow, "Value"]
    out$p_interaction <- tf[irow, "p-value"]
  }
  out
}

#' Two-line annotation text: paired Wilcoxon and linear mixed model (LMM) p-values side by side
#' @param n_paired,p_wilcoxon,p_lmm Vectors, one element per annotation.
test_label <- function(n_paired, p_wilcoxon, p_lmm) {
  fmt <- function(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
  sprintf("Wilcoxon p=%s (n=%d)\nLMM p=%s", fmt(p_wilcoxon), n_paired, fmt(p_lmm))
}

#' Subtitle with the group x visit interaction p-value of the LMM
#' @param ann Rows of the test table for one parameter or metric (needs `p_interaction`).
interaction_subtitle <- function(ann) {
  p <- unique(ann$p_interaction[!is.na(ann$p_interaction)])
  if (length(p) == 0) return("Diet x visit (LMM): n/a")
  sprintf("Diet x visit (LMM): p=%s", if (p[1] < 0.001) "<0.001" else sprintf("%.3f", p[1]))
}
