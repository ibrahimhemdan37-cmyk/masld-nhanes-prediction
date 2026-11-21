# NHANES MASLD cross-sectional cohort

Built by `build_cohort_masld.py` from the NHANES public-use XPT files in
`../raw/` (downloaded from the CDC on 2026-11-21). Protocol v1.2.

## Design

* Development sample: NHANES 2017-March 2020 pre-pandemic (P_ files, 3.2 years,
  MEC weight WTMECPRP).
* Temporal validation sample: NHANES August 2021-August 2023 (L files, 2.0
  years, MEC weight WTMEC2YR).
* Pooled weights are rescaled by release years / 5.2 (3.2/5.2 and 2/5.2).
* No NHANES data collected after August 2023 had been publicly released.

## Why the development sample is not split by cycle

The combined 2017-March 2020 files cannot be split into 2017-2018 and
2019-March 2020 participants: CDC reassigned SEQN values in the combined
release (P_DEMO SEQN 109263-124822 versus DEMO_J 93703-102956, zero overlap),
regenerated the masked variance strata/PSU, and set SDDSRVYR to 66 for every
row. Verified empirically on 2026-11-21: all 9,254 DEMO_J participants match a
P_DEMO row on age, sex, race, education and poverty-income ratio, so the
combined files do contain the 2017-2018 participants, but under new SEQN
values with no cycle indicator. Protocol v1.1's SEQN-subsetting instruction
was therefore wrong and is superseded.

## Outcome

MASLD (Rinella et al. 2023): CAP >= 248 dB/m plus at least one of five
cardiometabolic criteria, alcohol within limits (<=2 drinks/day women,
<=3 men), and no viral hepatitis B (HBsAg, LBDHBG = 1) or C (HCV RNA,
LBXHCR = 1). A reliable exam is LUAXSTAT = 1 (fasted >= 3 h, >= 10 complete
stiffness measures, liver stiffness IQR/median < 30%).

Metabolic criteria are each recorded as met (1), not met (0) or undetermined
(missing). A participant is excluded for missing metabolic data only if no
criterion is met AND at least one criterion is undetermined: anyone with
>= 1 criterion met is classifiable regardless of missing components. This
matters because triglycerides (and fasting glucose) are measured only in the
morning fasting subsample, so requiring all five criteria to be known would
discard about half the cohort.

## Predictors (primary set, 10)

Age, sex, race/Hispanic origin, smoking status, education, leisure-time
moderate-to-vigorous physical activity (minutes/week), total cholesterol,
ALT, AST, GGT. Physical activity is leisure-time only, for comparability
across cycles: 2017-March 2020 uses the GPAQ recreational items (vigorous
days x minutes/day + moderate days x minutes/day); 2021-2023 uses the
leisure-time items (frequency per day/week/month/year x minutes per session,
separately for moderate and vigorous activity). No outcome-defining variable
(BMI, waist, blood pressure, glucose, HbA1c, triglycerides, HDL, diabetes
status) is among the predictors.

## Data fixes applied

* The value 5.397605347e-79 (= 2^-260) in these files is a SAS zero misread
  by pandas' XPORT reader, not a missing-value code: codebook counts prove
  it (ALQ121 code 0, "Never in the last year", n = 1,049, exactly equals the
  number of 2^-260 values in ALQ121). Such values are restored to 0.0; see
  `sas_zero_log.csv`. Genuine missing values already read as NaN.
* Blood pressure: the mean is taken over readings > 0 for both systolic and
  diastolic (symmetric rule, disclosed here).
* Alcohol: ALQ130 is drinks per drinking day, so the average is ALQ130 x
  drinking days per year / 365 (frequency from ALQ121 category midpoints);
  never-drinkers are coded 0. Using ALQ130 directly as drinks/day would
  wrongly exclude about 40% of the sample; the corrected rule excludes
  about 3-5%.

## Strata

`0` masked variance strata contain participants from both
releases (stratum numbers were reused across releases). Analyses must nest
the stratum within the release (use survey_release + sdmvstra as the
stratum), or use the release-specific design.

## Files

* `nhanes_masld_cross_sectional.csv` — analysis cohort (one row per participant)
* `cohort_flow.csv` — flow of exclusions by release
* `source_manifest.csv` — CDC source file for every component and release
* `outcome_summary.csv` — counts and weighted/unweighted prevalence
* `variable_completeness.csv` — missingness per variable and release
* `data_dictionary.csv` — variable definitions and sources
* `sas_zero_log.csv` — SAS zeros restored from the 2^-260 reading (if any)

## Cohort flow

                                           stage  all_releases  2017-March2020  2021-2023
                    All interviewed participants         27493           15560      11933
Age 20 or older with positive examination weight         14608            8544       6064
                           Valid CAP measurement         13536            7922       5614
       Reliable elastography exam (LUAXSTAT = 1)         12671            7395       5276
     After excluding survey-identified pregnancy         12671            7395       5276
             After excluding viral hepatitis B/C         11764            6870       4894
         After excluding significant alcohol use         10447            6344       4103
   After excluding undetermined metabolic status         10423            6320       4103

## Outcome summary

survey_release     n  masld  masld_unweighted_pct  masld_weighted_pct  steatosis  steatosis_weighted_pct  masld_with_significant_fibrosis  masld_with_significant_fibrosis_weighted_pct  n_complete_core  masld_complete_core
  All releases 10423   5974                 57.32               56.07       6039                   56.81                              912                                          8.29            10096                 5796
2017-March2020  6320   3709                 58.69               56.77       3743                   57.40                              522                                          7.86             6198                 3645
     2021-2023  4103   2265                 55.20               54.83       2296                   55.77                              390                                          9.06             3898                 2151
