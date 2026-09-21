# Data cleaning notes

This documents the decisions made in `scripts/01_clean_data.R`, why they were necessary, and what still needs a human decision. Read this before trusting `data/processed/*.csv`, and update it if you change the cleaning script.

## What the raw files actually look like

All raw files live in `data/` and are untouched by the cleaning script - it only reads them. They are exports from whatever system originally recorded the ERIE fructose challenge test (FCT) data, and have several rough edges:

| File | Format problem |
|---|---|
| `ERIE_fructose_12C.csv` | Wide, semicolon-delimited, comma decimals |
| `ERIE_fructose_13C.csv` | Wide, semicolon-delimited, **period** decimals, one typo'd column header |
| `ERIE_weight.csv` | Wide (single row!), semicolon-delimited, comma decimals; replaced by `weight_height_data_ERIE.csv`, see below |
| `weight_height_data_ERIE.csv` | Tidy-ish, 3-digit subject IDs, one column quoted for no obvious reason |
| `ERIE_metadata_sex.csv` | Tidy already, but a different subject-ID format |
| `ERIE_constants.csv` / `.xlsx` | Two files, disagree on one value |
| `ERIE_Diets.xlsx` | Tidy already, 2-digit subject IDs - no reformatting needed |

## Decoding the wide column labels

Columns are labelled `FCT1 - 1` ... `FCT1 - 35`, `FCT2 - 36` ... `FCT2 - 70`. The number is **not** a subject ID - it's a serial column position from 1 to
70. Columns 1-35 are FCT1 (baseline visit) for subjects 1-35, in order, and columns 36-70 are FCT2 (post-diet visit) for the *same* 35 subjects in the *same* order (column 36 = subject 1's FCT2, ..., column 70 = subject 35's FCT2).

This mapping is not stated anywhere in the raw files. It was reverse-engineered by cross-referencing values against the older tidy spreadsheets in `former_models/My12C6FructoseModel/Data/ERIE_12C6Fructose_tidy.xlsx`, which label columns with real subject IDs. For example, raw column "FCT1 - 1", t=0 = 23.9687... matches "ER01" FCT1 t=0 in the tidy spreadsheet exactly, and column position 36's single value in `ERIE_weight.csv` (101 kg) matches ER01's separately-recorded FCT2 weight. Several more values were checked the same way before trusting the mapping for the full cohort.

**Subject number** is always derived from the numeric column position, and label text is ignored for it. **Visit** is read from the "FCT1"/"FCT2" text in the label when it is present and parses cleanly, and falls back to the position rule (`<= 35` -> FCT1, `> 35` -> FCT2) otherwise. Label and position are used as given, and a mismatch between them would silently change which visit a column's data is assigned to. The only automated consistency check is that every (subject, visit, time) point in one isotope file also exists in the other (both come from the same blood samples), which raises a `warning()` otherwise.

Once resolved to `FCT1`/`FCT2`, the value is recoded to `baseline`/ `intervention` - that's the coded `visit` value used everywhere downstream (`erie_concentrations.csv`, `erie_covariates.csv`, and every script that reads them). `FCT1`/`FCT2` stays raw-file vocabulary only.

## Inconsistent decimal marks across files

`ERIE_fructose_12C.csv` and `ERIE_constants.csv` use comma decimals (e.g. `23,96874488`). `ERIE_fructose_13C.csv` and `weight_height_data_ERIE.csv` both use period decimals instead (e.g. `0.876411032`) - verified by grepping the fructose files for both patterns: zero comma-decimals in the 13C6 file and zero period-decimals in the 12C file, so each file follows one consistent convention.

**This matters:** parsing all four files with the same (comma-decimal) locale does not error loudly - it silently fails to parse most of the 13C6 tracer curve, and `readr` reports this only as a "parsing issues" warning that's easy to scroll past. The cleaning script parses each file with the locale it actually uses. This was found by inspecting `problems()` on the parsed 13C6 table during development, since no automatic check exists, so check the decimal mark of any new raw file explicitly.

## Subject ID convention

`ERIE_metadata_sex.csv` and `weight_height_data_ERIE.csv` use 3-digit IDs (`ER001`...`ER035`). Everything else (raw column labels, once decoded, and all of `former_models/`) uses 2-digit IDs (`ER01`...`ER35`). The cleaning script standardizes everything to the 2-digit form, for consistency with all prior modeling work.

## Weight and height: `weight_height_data_ERIE.csv`

This file replaces `ERIE_weight.csv` (which stays in `data/` and is read by no script). It has one row per subject with `dem_height` (cm, constant across visits), `fct1_gewicht` and `fct2_gewicht` (kg, one column per visit), and the height is needed for the weight+height-based Vd formula in `scripts/assets/pk_curves.R` (see "Volume of distribution" in `docs/pk-model.md`).

`fct2_gewicht` is double-quoted in the raw CSV (`fct1_gewicht` is unquoted) for no apparent reason, including `""` for ER33/ER34's missing FCT2 weight. `readr::read_csv()` infers the column as numeric (quoting in CSV only escapes delimiters) and turns `""` into `NA`, so no special handling is needed.

## The 13C6 tracer dose discrepancy

`ERIE_constants.csv` and `ERIE_constants.xlsx` agree exactly on two values:

- Molar amount: **537.3166407 µmol**
- Molecular weight of 13C6-fructose: **186.11 g/mol**

They disagree on a third: the csv lists the weighed mass as `0,12` (0.12 g / 120 mg), while the xlsx lists it as `0.1` (100 mg).

The arithmetic settles which one is consistent with the molar amount both files agree on:

```
537.3166407 µmol x 186.11 g/mol = 99,999.99997 µg = 99.99999997 mg ≈ 100.000 mg
```

That matches the xlsx's 100 mg to 7 significant figures, while 120 mg would require 644.78 µmol instead of 537.32 µmol. The csv's "0,12" looks like a data-entry typo (an extra digit) unrelated to the molar figure.

Separately, the manuscript's Methods section and every prior model script in `former_models/` state the capsule as **120 mg** - almost certainly the *nominal* protocol target dose (what the capsule was designed to contain), which can legitimately differ from what was actually weighed into a given batch.

**Decision (2026-09-17):** whoever weighed the capsules confirmed that 120 mg was the dose administered, so the molar amount stated in the constants files (537.3166407 µmol) is wrong, as is the xlsx's mass of 100 mg. The cleaning script takes 120 mg as ground truth and **recalculates** the molar dose from it (120 mg x 1000 / 186.11 g/mol = **644.7799688 µmol**), ignoring the xlsx's stated molar amount. `data/processed/erie_constants.csv` carries the confirmed and recalculated pair (`tracer_13C6_dose_mg` / `tracer_13C6_dose_umol`) and omits the xlsx-derived 100 mg / 537.3166407 µmol values, which are known to be wrong; the arithmetic above records how that was established.

## Known missingness

- **ER33 and ER34 have no intervention data at all** (body weight, both fructose curves) - both are documented study dropouts, consistent with the manuscript ("35 participants completing the study compared with the intended 40" and per-subject attrition described in Fig. 1). Their baseline weight is left out of the missing intervention `bw_kg`: both subjects also have no intervention concentration data, so a filled-in weight would never be used, and leaving it `NA` reflects what is known about their intervention visit.
- Two 13C6 concentration values are slightly negative (ER06 intervention t=360: -0.0011 µmol/L; ER27 intervention t=240: -0.0172 µmol/L). Both are small in magnitude, at late timepoints where the true tracer concentration is near zero, and consistent with ordinary assay noise near the limit of detection. They stay as they are in the cleaned data; whether to clip them to zero before fitting is a modeling decision (see `docs/pk-model.md`).

## Diet group assignment

`data/ERIE_Diets.xlsx` (copied from `former_models/MixedModel/Data/`) gives one row per subject: `Diet` = `A` or `B`. Diet A is low fructose, with calories matched by glucose supplementation; Diet B is high fructose. This is a subject-level assignment that stays fixed across baseline and intervention (the diet intervention happens between the two visits), so it is joined into `erie_covariates.csv` by `subject_id` alone. Used by `scripts/03_diet_summary.R` for the diet-arm comparisons (see `docs/diet-summary.md`); the PK fitting itself (`02_fit_erie_model.R`) doesn't need it.

## Output files

All in `data/processed/`, one row per (subject, visit, isotope, time) unless noted:

- `erie_concentrations.csv` - `subject_id, visit, isotope, time_min, conc_umol_L`
- `erie_covariates.csv` - one row per subject x visit: `bw_kg, sex, diet, height_cm, dose_12C_mg, dose_13C6_umol, dose_13C6_mg`
- `erie_constants.csv` - `constant, value, unit`: physical/dosing constants (13C6 dose in mg and µmol, both molecular weights, and the 12C dose per kg). `MW_12C` (180.16 g/mol, unlabeled fructose) is a standard value; `MW_13C6` is read from `ERIE_constants.xlsx`.
