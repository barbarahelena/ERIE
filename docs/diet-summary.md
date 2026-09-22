# Diet-arm summary

This documents `scripts/03_diet_summary.R`: how the fitted PK parameters are compared between the low- and high-fructose diet arms, which fits are included, and why the analysis choices were made. It reads `results/fit_results.csv` (see "Output files" in `docs/pk-model.md`) and `data/processed/erie_covariates.csv`, and writes everything to `results/`.

Read `docs/pk-model.md`, "Fit quality and what to trust", first: it explains which of the parameters below can be interpreted as absolute numbers (`ka`, `kel`) and which only within-subject (`F_12C`, `F_13C6`).

## Which fits are included

A subject x visit counts as `reliable` when 12C's fit is not stuck at the shared `kel` bound (`kel_at_bound`) and its own R² is at least `R2_INCLUDE_MIN` (0.90 in the script). All parameter summaries, the LMMs, the boxplots, delta plots, two-wave summaries and average curves use only `reliable` fits.

- **`R2_INCLUDE_MIN` is separate from, and stricter than, `R2_RELIABLE_MIN` (0.70)** in `02_fit_erie_model.R`. That one only flags a fit for the fitting pipeline's own retry logic; this one decides whether a fit belongs in a diet-arm comparison. Change it in the script and re-run; the number of fits dropped at each cutoff is visible in the `n` column of `results/diet_parameter_summary.csv`.
- **Gated on the 12C fit only**, and applied to every parameter, including `F_13C6`. 13C6's much smaller, noisier signal fails its own R² bar far more often even when the shared kinetics from the same joint fit (`ka`, `kel`) are trustworthy; requiring 13C6 to pass too excluded many otherwise-fine subjects.
- **`converged` is left out of the gate.** It records whether `optim()`'s L-BFGS-B reached its strict stopping criterion before the `maxit` cap, which is a different question from whether the winning fit is good. Checked on this cohort's own results: `r2_12C` for `converged = FALSE` rows (mean 0.953, minimum 0.759, n=26) was essentially the same as for `converged = TRUE` rows (mean 0.959, minimum 0.812, n=42), so the flag does not separate good fits from bad ones. Requiring it anyway would discard good fits - and because the paired comparisons need *both* visits reliable at once, a ~38% per-curve "not converged" rate (26/68) compounds to exclude a large share of subjects from any paired test. (Counted on the run that produced the current `results/fit_results.csv`; recompute if the fits change.)

## Outputs

Two plots use all fitted subjects, without the `reliable` filter:
- **Vd by diet arm** (`diet_vd_boxplot.pdf`). Vd is 0.4 x each visit's measured total body water (see "Volume of distribution" in `docs/pk-model.md`) and has no fit of its own, so a PK fit's R² has no bearing on it. It is shown once per subject (mean over the visits present). It works as a covariate-balance check: it should cluster by sex and be similar between arms if randomization worked.
- **Fit quality (R²) by diet arm and visit** (`diet_r2_boxplot.pdf`). This is a QC view of the whole cohort; `reliable` excludes on R², so filtering by it would hide exactly the poor fits the plot exists to show. The dashed line is `R2_RELIABLE_MIN`.

Everything else covers `ka`, `kel`, `F_12C` and `F_13C6`:

| Output | What it is |
|---|---|
| `diet_parameter_summary.csv` | mean, SD, SEM and n per diet x visit |
| `diet_lmm_results.csv` | F-tests from the LMMs (below) |
| `diet_parameter_boxplot_<param>.pdf` | one PDF per parameter, faceted by diet arm, coloured by visit |
| `diet_parameter_deltas.csv`, `diet_parameter_delta_boxplot.pdf` | within-subject change (below) |
| `diet_before_after_wilcoxon.csv`, `diet_before_after_boxplot_<param>.pdf` | baseline vs intervention within each arm (below) |
| `diet_two_wave_rate.csv`, `diet_two_wave_params.csv`, `diet_two_wave_boxplot.pdf` | second-wave characteristics (below) |
| `diet_summary_curves_median_<isotope>.pdf`, `diet_summary_curves_mean_<isotope>.pdf` | median and mean fitted curves per isotope (below) |
| `diet_curve_metrics.csv`, `diet_curve_metrics_summary.csv`, `diet_curve_metrics_wilcoxon.csv`, `diet_curve_auc_boxplot.pdf`, `diet_curve_cmax_boxplot.pdf`, `diet_curve_tmax_boxplot.pdf` | AUC, Cmax and Tmax from the fitted curves (below) |

The per-parameter PDFs are faceted by diet arm and coloured by visit, so the baseline-vs-intervention comparison - the actual before/after-diet question - is the primary visual read within each panel.

### Linear mixed models

One LMM per parameter: `value ~ diet * visit + sex + age_years + bw_kg + (1 | subject_id)`, with Type III F-tests and Satterthwaite degrees of freedom (the `lmerTest` default). The random intercept is there because each subject contributes a paired baseline/intervention observation. `lmer` handles the unbalanced design (some subjects have only one reliable visit), which a paired t-test cannot. A parameter is skipped, with a warning, if fewer than 3 subjects have both visits reliable. All results are in one `diet_lmm_results.csv`.

Ten models in total: `ka`, `kel`, `F_12C` and `F_13C6` from the fitted parameters, and `AUC`, `Cmax` and `Tmax` per isotope (`AUC_12C`, `AUC_13C6`, ...) from the curve-derived metrics ("Curve-derived metrics" below) - the metric models use the same `reliable`-filtered fits, since `curve_metrics` is itself built only from reliable fits.

- **Log scale, except Tmax.** `ka`, `kel`, the `F` values, `AUC` and `Cmax` are positive, multiplicative-scale quantities, conventionally treated as log-normal in PK work. Log-transforming makes the normal-residuals assumption more defensible and turns a diet effect into a fold-change, the more natural scale for a rate constant or an exposure. The raw scale was run alongside it and log was consistently the better-behaved fit, so the raw-scale model was dropped for those five. `Tmax` stays on its raw minutes scale, matching the paired Wilcoxon test on curve metrics below, since it is not a multiplicative quantity and can be 0.
- **Sex, age and body weight as covariates.** Sex, age and weight plausibly affect absorption, clearance and `F` independent of diet. Vd follows the measured TBW, which correlates with sex and weight, but `ka`, `kel` and `F` are fit independent of Vd, so including these covariates lets the diet/visit effect be read net of that variation. The 12C dose is 1000 mg per kg body weight, so `bw_kg` also adjusts the 12C-derived metrics (`F_12C`, `AUC_12C`, `Cmax_12C`) for that dose-weight link directly, on top of any effect through Vd. Age and weight are continuous and untransformed, both in years and kg as recorded. BMI was tried as a covariate as well and dropped.
- **`ka`'s random intercept is at the boundary (variance 0, `isSingular()` is `TRUE`).** Between-subject variance in log(`ka`) is close to zero once the fixed effects (including the two new covariates) are accounted for, so `lmer` estimates no extra subject-level spread on top of them. `lmerTest`'s F-tests are still valid in this case; they reduce to what a fixed-effects-only model on the same data would give. The other nine models are not singular.

### Within-subject change (delta)

The boxplots compare baseline and intervention as separate distributions. The delta plot answers "did the diet change this parameter, and differently by arm": for each subject with *both* visits reliable (so every delta is a complete pair), `log(intervention) - log(baseline)`, matching the log-scale LMMs and for the same reason. The two arms are compared with a Wilcoxon rank-sum test (`ggpubr::stat_compare_means`), since n per arm is small and the test needs no normality assumption.

### baseline vs intervention within each arm

A different question from the delta plot (which compares the *size* of the change *between* arms): does baseline differ from intervention at all *within* an arm. This is a paired Wilcoxon signed-rank test, computed as a one-sample Wilcoxon test of each subject's own delta against 0 (NA if fewer than 3 pairs). It is computed explicitly, because ggpubr's automatic pairing detection silently breaks if a subject lands in a different row order in some facet.

The boxplot itself shows every reliable single-visit value, including subjects with only one reliable visit - a larger, more representative sample than the paired test can use, since pairing needs both visits. So the boxplot's n and the test's n legitimately differ; the annotation gives the test's `n`.

### Second-wave characteristics

Whether `two_wave` was selected at all is itself a diet-relevant outcome, and among subjects who show one, so are its timing (`t_lag`, the time between peaks) and how much of the dose rode it (`f_delayed` for 12C, `f_delayed_13C6` for 13C6's own independent choice; see `docs/pk-model.md`). This is **exploratory only**: `two_wave` is selected for most fits (49/68 in the current results), so whether it was selected barely separates subjects, and `t_lag`/`f_delayed` depend on that model choice. It is reported for completeness as an exploratory result. The script prints Fisher's exact test of diet x `two_wave` selection without saving it, and draws the parameter boxplot when at least 4 `two_wave` fits are available.

### Median and mean concentration-time curves

Each subject's fitted curve is simulated on a 2-min grid, then summarised at every time point per diet x visit. The median with a 95% bootstrap confidence band and the mean with a SEM band are separate plots, one PDF per statistic and isotope, with diet arms as panels and baseline/intervention overlaid by colour. Both diet panels share one y-axis that starts at 0, so the arms can be compared directly. The median is a Harrell-Davis median: a weighted average of the sorted subject values that changes smoothly where subjects swap rank, which the plain sample median does not (it kinks at every swap). It smooths across subjects and leaves the time axis alone, so peak height and time stay close to the plain median's. The band is the 2.5th to 97.5th percentile of 1000 bootstrap medians, resampling subjects with the same estimator (seed `BOOT_SEED`). It is a pointwise band, and with 14 to 17 subjects per group it is approximate. The mean plot shows the group average and the precision of that average. Both are pointwise summaries, so neither is the curve of a single subject or of median parameters. Two things matter here:

- **The simulation must mirror `simulate_fit()` in `02_fit_erie_model.R` exactly.** Each curve's shape depends on which model and 13C6 mechanism won for that subject x visit, and 13C6's onset-lag and second-wave mechanisms can win under *either* 12C model, so `t_lag1_13C6` and `t_lag2_13C6` are checked before `model`. Those mechanisms set `k_release = NA`, and feeding that into `simulate_delayed_release()` silently returns `NA` for every t>0. The group summaries (`hd_median_rows()`, `rowMeans()`, `sd`) have no `na.rm`, so one such subject would blank out the whole diet x visit group's curve.
- **One PDF per isotope.** 12C and 13C6 differ ~1000x in concentration (the same reason the per-subject plots in `02_fit_erie_model.R` use free y-scales), so a combined figure would either hide 13C6's shape on a shared axis or need free scales that make the isotopes hard to compare or caption as one figure.

### Curve-derived metrics (AUC, Cmax, Tmax)

Each subject x visit x isotope curve (the same simulated curves as above) gives three non-compartmental metrics: AUC from time 0 until the curve has cleared (AUC to infinity), Cmax (the highest concentration) and Tmax (the time of Cmax, on the 2-min grid). The per-subject values are in `diet_curve_metrics.csv`, and `diet_curve_metrics_summary.csv` has the median, quartiles and n per diet x visit x isotope. One boxplot PDF per metric has isotopes as rows and diet arms as columns, coloured by visit, annotated with a paired Wilcoxon signed-rank test of intervention vs baseline within each arm and isotope (`diet_curve_metrics_wilcoxon.csv`). The test uses subjects with both visits reliable (at least 3 pairs), on the log ratio for AUC and Cmax and on the difference in minutes for Tmax, which has many ties on its 2-min grid, so its p-value uses the normal approximation. No correction is applied for the 12 tests.

- **They come from the fitted curves.** AUC, Cmax and Tmax are derived from the model, so they follow from `ka`, `kel`, `F`, Vd and the second-wave parameters. With Vd fixed by the measured TBW, AUC scales with `F / kel`, and Cmax and Tmax depend on `ka`, `kel` and any second wave. They are a view of the same fits as the parameter boxplots and give a link between those boxplots and the curve plots.
- **Median of AUCs differs from AUC of the median curve.** The per-subject median is the value to compare with the boxplots. The pointwise median curve has an area that need not equal it.
- **Tmax is coarse and unstable.** It has 2-min resolution, and on flat or two-peaked curves it can jump between the two peaks with a small change in the fit. Cmax and AUC are more stable.
- **AUC is computed analytically.** Every mechanism delivers `F * dose` to the blood and `kel` clears it, so the area to infinity is `F * dose / (Vd * kel)` for each isotope, with no cut-off time. It agrees with trapezoid integration of the simulated curves over 0 to 400 min (ratio 0.95 to 1.00, lowest for subjects with a small `ka` or `kel`, whose curves have not cleared by 400 min).
