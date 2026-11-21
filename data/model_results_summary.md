# MASLD predictive modelling — results summary

Protocol v1.2. Development: NHANES 2017–March 2020 (n = 6,320; 3,709 MASLD
cases). Temporal validation: NHANES August 2021–August 2023 (n = 4,103;
2,265 cases). Outcome: MASLD = CAP ≥248 dB/m + ≥1 cardiometabolic
criterion + alcohol within limits + no viral hepatitis. Same 10 predictors
in both models. Survey-weighted analysis throughout; 95% CIs from 1,500
rescaled Rao–Wu bootstrap replicates (seed 1729). Odds-ratio CIs from 1,500
nonparametric bootstrap refits. Software: Python 3.12.3, scikit-learn 1.9.1,
SciPy 1.18.1, pandas 3.1.1. Analysis script: `model_masld.py`.

## Tuning (5-fold cross-validation within development, weighted Brier)

- Logistic regression: every C in {0.01, 0.1, 1, 10} gave CV Brier 0.2090 —
  the penalty is negligible relative to the loss at n ≈ 6,300 (fitted
  coefficient norms 1.262–1.265), so C = 0.01 was taken by tie-break.
- Gradient boosting: best configuration 7 leaves, learning rate 0.05,
  100 iterations, min_samples_leaf 30, L2 = 1.0 — CV Brier 0.2002.
  Larger/deeper configurations overfit (CV Brier up to 0.225).

## Primary outcome — temporal validation

| Model | AUROC (95% CI) | Average precision | Brier (95% CI) | Calibration intercept | Calibration slope (95% CI) | Mean predicted vs observed |
|---|---|---|---|---|---|---|
| Logistic regression | 0.698 (0.690–0.706) | 0.729 (0.722–0.737) | 0.218 (0.216–0.221) | 0.027 (0.003–0.051) | 0.838 (0.792–0.884) | 54.3% vs 54.8% |
| Gradient boosting | 0.721 (0.713–0.729) | 0.743 (0.735–0.751) | 0.211 (0.209–0.214) | 0.044 (0.018–0.070) | 0.891 (0.854–0.928) | 53.9% vs 54.8% |
| Benchmark (constant development prevalence) | — | — | 0.248 | — | — | — |

Paired differences (GBM − LR): AUROC +0.023 (95% CI 0.020–0.026);
Brier −0.0072 (−0.0083 to −0.0062). Both CIs exclude zero: gradient
boosting is genuinely, if modestly, better — unlike the CKD report, where
the two models were tied (paired AUROC difference −0.004).

## Threshold metrics (validation, survey-weighted)

| Model | Rule | Threshold | Sensitivity | Specificity | PPV | NPV |
|---|---|---|---|---|---|---|
| Logistic | Youden's J | 0.502 | 72.3% (71.0–73.6) | 56.8% (55.3–58.3) | 67.0% | 62.8% |
| Logistic | 90% sensitivity | 0.357 | 90.2% (89.4–90.9) | 31.3% (29.7–32.8) | 61.4% | 72.4% |
| GBM | Youden's J | 0.511 | 71.9% (70.9–73.0) | 60.6% (59.3–61.9) | 68.9% | 64.0% |
| GBM | 90% sensitivity | 0.342 | 90.0% (89.3–90.8) | 36.7% (35.3–38.0) | 63.3% | 75.2% |

## Odds ratios — logistic model (development fit, bootstrap 95% CI)

| Predictor | OR (95% CI) |
|---|---|
| Age, per 10 years | 1.32 (1.26–1.40) |
| MVPA, per 60 min/week | 0.96 (0.94–0.98) |
| Total cholesterol, per 20 mg/dL | 1.03 (0.99–1.08) |
| ALT, per 10 U/L | 1.90 (1.58–2.34) |
| AST, per 10 U/L | 0.60 (0.47–0.71) — see note below |
| GGT, per 10 U/L | 1.07 (1.01–1.14) |
| Male (vs female) | 1.18 (0.97–1.44) |
| Mexican American (vs non-Hispanic White) | 1.90 (1.48–2.45) |
| Non-Hispanic Asian | 1.11 (0.89–1.40) |
| Non-Hispanic Black | 0.89 (0.75–1.07) |
| Other Hispanic | 1.00 (0.77–1.26) |
| Other or multiracial | 1.54 (1.08–2.18) |
| Current smoking (vs never) | 0.87 (0.69–1.10) |
| Former smoking (vs never) | 1.18 (0.96–1.46) |
| Some college (vs college graduate or above) | 1.28 (1.03–1.58) |
| High school graduate (vs college graduate) | 1.22 (0.98–1.55) |
| 9th–11th grade (vs college graduate) | 1.16 (0.87–1.57) |
| Less than 9th grade (vs college graduate) | 1.10 (0.81–1.50) |

**AST note.** The AST odds ratio below 1 is a collinearity artefact, not a
protective effect: ALT and AST correlate at r = 0.76 in development, and
the two coefficients must be read jointly (ALT strongly positive, AST
negative once ALT is held fixed). Permutation importance confirms AST still
contributes to the gradient boosting model. Do not report the AST odds
ratio on its own as a protective association.

## Permutation importance — gradient boosting (validation, Brier increase)

ALT 0.026, GGT 0.024, age 0.011, AST 0.009, MVPA 0.003, race/ethnicity
0.001, smoking 0.0005, total cholesterol 0.0004, sex 0.0001,
education 0.00004. The liver enzymes dominate; demographics and lifestyle
contribute little once the enzymes are in the model.

## Secondary outcomes and sensitivity analyses (validation)

| Analysis | LR AUROC | GBM AUROC | LR Brier | GBM Brier | LR slope | GBM slope |
|---|---|---|---|---|---|---|
| Steatosis alone (CAP ≥248) | 0.696 | 0.719 | 0.218 | 0.211 | 0.86 | 0.90 |
| MASLD + significant fibrosis (LSM ≥8 kPa) | 0.675 | 0.697 | 0.080 | 0.079 | 0.65 | 0.68 |
| MASLD at CAP ≥285 dB/m | 0.702 | 0.737 | 0.202 | 0.191 | 0.73 | 0.96 |
| MASLD by FLI ≥60 (fasting subsample, n = 2,150) | 0.792 | 0.808 | 0.187 | 0.179 | 0.94 | 1.07 |
| Criterion-containing predictor set | 0.853 | 0.853 | 0.154 | 0.154 | 0.96 | 1.01 |
| Routine-clinical predictor set (core + BMI + SBP) | 0.831 | 0.833 | 0.167 | 0.165 | 0.94 | 0.99 |
| Complete-case analysis (n = 3,898) | 0.698 | 0.723 | 0.218 | 0.210 | 0.83 | 0.90 |
| Unweighted analysis | 0.693 | 0.725 | 0.221 | 0.209 | 1.02 | 1.08 |

## Interpretation for the manuscript

1. **Gradient boosting wins, modestly.** +0.023 AUROC and −0.007 Brier,
   both with CIs excluding zero. Report the paired differences, not just
   the two AUROCs side by side.
2. **Discrimination is modest (0.70–0.72).** Expected: the predictors
   deliberately exclude everything that defines the outcome, and MASLD at
   CAP ≥248 is very common (55% weighted), which caps discrimination.
   Comparable to the CKD report (0.67).
3. **Calibration is reasonable and better than the CKD report.**
   Slopes 0.84 (LR) and 0.89 (GBM) — mild over-dispersion, versus 0.63 in
   the CKD report. Intercepts 0.03–0.04 — good calibration-in-the-large.
4. **The criterion-containing set (AUROC 0.853) is the leakage
   demonstration**, not a better model: adding BMI, waist, HbA1c,
   triglycerides, HDL and systolic BP — variables that define the outcome
   — inflates AUROC by 0.13. Report it as a sensitivity analysis with the
   caution, exactly as the protocol planned.
5. **The routine-clinical set (AUROC 0.831/0.833, slopes 0.94/0.99) is
   the deployable model**: core predictors + BMI + systolic BP only — both
   measured free at every clinic visit, no fasting sample. It answers the
   clinical objection to excluding anthropometrics, and it belongs in the
   Discussion as the answer to "why not use what's already in the chart".
5. **FLI ≥60 (AUROC 0.79–0.81) is partly definitional**: FLI is built
   from triglycerides, BMI, waist and GGT, which correlate with the
   predictors; and it uses the fasting subsample only. Do not present it
   as better prediction.
6. **MASLD + fibrosis is the weakest analysis**: AUROC 0.68–0.70 and
   slopes 0.65–0.68, on only 390 validation events — report it, do not
   over-interpret it.
7. **Robustness is good**: complete-case and unweighted analyses give
   essentially identical AUROCs; imputation and weighting are not driving
   the results.
8. **Thresholds are honest but not clinically decisive**: at Youden's J,
   ~72% sensitivity with ~57–61% specificity; at 90% sensitivity,
   specificity falls to ~31–37% — a screening-style threshold would flag
   most of the adult population. The model is a risk-stratification aid,
   not a diagnostic rule.

## Files

- `model_masld.py` — analysis script (reproducible from the cohort CSV)
- `model_results.json` — all numbers above, machine-readable
- `metrics_table.csv` — validation metrics, primary + variants
- `odds_ratios.csv` — logistic odds ratios with bootstrap CIs
- `permutation_importance.csv` — GBM permutation importance
- `threshold_metrics.csv` — threshold metrics with bootstrap CIs
- `tuning_results.csv` — 5-fold CV Brier for every tuned configuration
- `validation_predictions.csv` — per-participant validation predictions
  (SEQN, outcome, both model probabilities, weight)
