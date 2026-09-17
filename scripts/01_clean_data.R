# Data cleaning of ERIE data
# Barbara Verhaar
#
# Rationale for the non-obvious decisions below (column-to-subject mapping,
# per-file decimal marks, subject ID convention, tracer dose default) is in
# docs/data-cleaning-notes.md.

# Libraries
suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(tibble)
})

# Paths
raw_dir <- "data"
out_dir <- "data/processed"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
eu_locale <- locale(decimal_mark = ",", grouping_mark = ".")

# =============================================================================
# 1. Fructose conc (12C liquid dose + 13C6 tracer)
# =============================================================================

read_fructose_file <- function(path, isotope_name, locale) {
  df <- read_delim(path, delim = ";", locale = locale, na = c("nan", "NaN", ""),
             col_types = cols(.default = "d"), show_col_types = FALSE) %>%
    rename(time_min = 1) %>%
    pivot_longer(-time_min, names_to = "col_label", values_to = "conc_umol_L") %>%
    mutate(
      index          = as.integer(str_extract(col_label, "\\d+$")),
      subject_id     = sprintf("ER%02d", if_else(index <= 35, index, index - 35L)),
      visit_label    = str_extract(col_label, "FCT[12]"),
      visit_position = if_else(index <= 35, "FCT1", "FCT2"),
      visit          = coalesce(visit_label, visit_position),
      isotope        = isotope_name
    )
  df %>% select(subject_id, visit, isotope, time_min, conc_umol_L)
}

# 12C file uses comma decimals, 13C6 file uses period decimals
conc_12C <- read_fructose_file(file.path(raw_dir, "ERIE_fructose_12C.csv"), "12C", locale = eu_locale)
conc_13C6 <- read_fructose_file(file.path(raw_dir, "ERIE_fructose_13C.csv"), "13C6", locale = default_locale())

# Both isotopes are drawn from the same blood samples, double check if matching
time_check <- full_join(
  conc_12C %>% distinct(subject_id, visit, time_min) %>% mutate(has_12C = TRUE),
  conc_13C6 %>% distinct(subject_id, visit, time_min) %>% mutate(has_13C6 = TRUE),
  by = c("subject_id", "visit", "time_min")
) %>% filter(is.na(has_12C) | is.na(has_13C6))
if (nrow(time_check) > 0) {
  warning(nrow(time_check), " (subject, visit, time) combinations present in one isotope file but not the other")
}

concentrations <- bind_rows(conc_12C, conc_13C6) %>%
  arrange(subject_id, visit, isotope, time_min)
write_csv(concentrations, file.path(out_dir, "erie_concentrations.csv"))

# =============================================================================
# 2. Weight and height
# =============================================================================

wh_raw <- read_csv(file.path(raw_dir, "weight_height_data_ERIE.csv"), show_col_types = FALSE)

bodyweights <- wh_raw %>%
  transmute(
    subject_id = sprintf("ER%02d", as.integer(str_extract(Subject_ID, "\\d+$"))),
    FCT1 = fct1_gewicht,
    FCT2 = fct2_gewicht
  ) %>%
  pivot_longer(c(FCT1, FCT2), names_to = "visit", values_to = "bw_kg")
# ER33/ER34 have no FCT2 body weight, dropped out of study

heights <- wh_raw %>%
  transmute(
    subject_id = sprintf("ER%02d", as.integer(str_extract(Subject_ID, "\\d+$"))),
    height_cm  = dem_height
  )

# =============================================================================
# 3. Sex
# =============================================================================
sex <- read_csv(file.path(raw_dir, "ERIE_metadata_sex.csv"), show_col_types = FALSE) %>%
  transmute(
    subject_id = sprintf("ER%02d", as.integer(str_extract(Subject_ID, "\\d+$"))),
    sex = Sex
  )

# =============================================================================
# 4. Diets: Diet A -> low_fructose (calorie suppl w/ gluc), B -> high_fructose
# =============================================================================
diet_raw <- read_xlsx(file.path(raw_dir, "ERIE_Diets.xlsx"))
diet <- tibble(
  subject_id = diet_raw[[1]],
  diet = recode_values(diet_raw$Diet, from = c("A", "B"), to = c("low_fructose", "high_fructose"))
)

# =============================================================================
# 4. Constants (dose/molecular-weight provenance)
# =============================================================================

const_xlsx <- read_xlsx(file.path(raw_dir, "ERIE_constants.xlsx"), col_names = c("label", "value", "unit"))

mw_13C6_g_per_mol <- const_xlsx$value[const_xlsx$label == "molecular weight 13C fructose"]

# The molar dose recalculated from the 120 mg mass
administered_dose_mg   <- 120
administered_dose_umol <- administered_dose_mg * 1000 / mw_13C6_g_per_mol   # mg * 1000 / (g/mol) -> umol

MW_12C <- 180.16   # g/mol, unlabeled fructose - standard value, not subject/visit-specific

constants <- tribble(
  ~constant,                ~value,                  ~unit,
  "tracer_13C6_dose_mg",    administered_dose_mg,    "mg",
  "tracer_13C6_dose_umol",  administered_dose_umol,  "umol",
  "MW_13C6",                mw_13C6_g_per_mol,       "g/mol",
  "MW_12C",                 MW_12C,                  "g/mol",
  "dose_12C_per_kg",        1000,                    "mg/kg"
)
write_csv(constants, file.path(out_dir, "erie_constants.csv"))

# =============================================================================
# 5. Covariates (one row per subject x visit)
# =============================================================================

covariates <- bodyweights %>%
  left_join(sex, by = "subject_id") %>%
  left_join(diet, by = "subject_id") %>%
  left_join(heights, by = "subject_id") %>%
  mutate(
    dose_12C_mg    = 1000 * bw_kg,
    dose_13C6_umol = administered_dose_umol,
    dose_13C6_mg   = administered_dose_mg
  ) %>%
  arrange(subject_id, visit)
write_csv(covariates, file.path(out_dir, "erie_covariates.csv"))
