#!/usr/bin/env python3
"""Table 1 (survey-weighted baseline characteristics by MASLD status) and
Figures 1-5 for the MASLD manuscript. Reads the built cohort CSV and the
modelling outputs. Writes table1.csv and figures/*.png (300 dpi).
Run with the masld-venv python (matplotlib, pandas, numpy, scipy)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False})

df = pd.read_csv(HERE / "nhanes_masld_cross_sectional.csv")
pred = pd.read_csv(HERE / "validation_predictions.csv")
ors = pd.read_csv(HERE / "odds_ratios.csv")
imp = pd.read_csv(HERE / "permutation_importance.csv")
flow = pd.read_csv(HERE / "cohort_flow.csv")


def wmean(x, w):
    x = np.asarray(x, dtype=float); w = np.asarray(w, dtype=float)
    ok = ~np.isnan(x)
    return float(np.average(x[ok], weights=w[ok]))


def wsd(x, w):
    x = np.asarray(x, dtype=float); w = np.asarray(w, dtype=float)
    ok = ~np.isnan(x)
    x, w = x[ok], w[ok]
    m = np.average(x, weights=w)
    return float(np.sqrt(np.average((x - m) ** 2, weights=w)))


def wmedian(x, w):
    x = np.asarray(x, dtype=float); w = np.asarray(w, dtype=float)
    ok = ~np.isnan(x)
    x, w = x[ok], w[ok]
    o = np.argsort(x); x, w = x[o], w[o]
    cw = np.cumsum(w)
    return float(x[np.searchsorted(cw, 0.5 * cw[-1])])


def wq(x, w, q):
    x = np.asarray(x, dtype=float); w = np.asarray(w, dtype=float)
    ok = ~np.isnan(x)
    x, w = x[ok], w[ok]
    o = np.argsort(x); x, w = x[o], w[o]
    cw = np.cumsum(w)
    return float(x[np.searchsorted(cw, q * cw[-1])])


def wpct(mask, w):
    mask = np.asarray(mask, dtype=bool); w = np.asarray(w, dtype=float)
    return 100.0 * float(np.sum(w[mask]) / np.sum(w))


def fmt_mean_sd(m, s):
    return f"{m:.1f} ({s:.1f})"


def fmt_med_iqr(med, q1, q3):
    return f"{med:.1f} [{q1:.1f}–{q3:.1f}]"


def fmt_pct(p):
    return f"{p:.1f}"


# ---------------------------------------------------------------- Table 1
def col_stats(sub, col, kind):
    w = sub["mec_weight_pooled"].values.astype(float)
    x = sub[col]
    if kind == "mean_sd":
        return fmt_mean_sd(wmean(x.values, w), wsd(x.values, w))
    if kind == "median_iqr":
        return fmt_med_iqr(wmedian(x.values, w), wq(x.values, w, 0.25),
                           wq(x.values, w, 0.75))
    raise ValueError(kind)


def cat_rows(sub, col, levels):
    w = sub["mec_weight_pooled"].values.astype(float)
    out = {}
    for lv in levels:
        out[lv] = fmt_pct(wpct(sub[col].values == lv, w))
    return out


DEV = "2017-March2020"
VAL = "2021-2023"
groups = {
    "Dev MASLD": df[(df.survey_release == DEV) & (df.masld == 1)],
    "Dev no MASLD": df[(df.survey_release == DEV) & (df.masld == 0)],
    "Val MASLD": df[(df.survey_release == VAL) & (df.masld == 1)],
    "Val no MASLD": df[(df.survey_release == VAL) & (df.masld == 0)],
}
for k, g in groups.items():
    print(k, "n =", len(g), flush=True)

rows = []
CONT_ROWS = [
    ("Age, years", "age_years", "mean_sd"),
    ("Body mass index, kg/m²", "bmi_kg_m2", "mean_sd"),
    ("Waist circumference, cm", "waist_cm", "mean_sd"),
    ("Systolic blood pressure, mmHg", "mean_systolic_bp_mm_hg", "mean_sd"),
    ("Diastolic blood pressure, mmHg", "mean_diastolic_bp_mm_hg", "mean_sd"),
    ("Total cholesterol, mg/dL", "total_cholesterol_mg_dl", "mean_sd"),
    ("HDL cholesterol, mg/dL", "hdl_cholesterol_mg_dl", "mean_sd"),
    ("HbA1c, %", "hba1c_pct", "mean_sd"),
    ("Physical activity, min/week", "mvpa_min_per_week", "median_iqr"),
    ("Triglycerides, mg/dL *", "triglycerides_mg_dl", "median_iqr"),
    ("Fasting glucose, mg/dL *", "fasting_glucose_mg_dl", "median_iqr"),
    ("ALT, U/L", "alt_u_l", "median_iqr"),
    ("AST, U/L", "ast_u_l", "median_iqr"),
    ("GGT, U/L", "ggt_u_l", "median_iqr"),
    ("CAP, dB/m", "cap_db_m", "median_iqr"),
    ("Liver stiffness, kPa", "liver_stiffness_kpa", "median_iqr"),
    ("Number of metabolic criteria met", "n_metabolic_criteria_met",
     "mean_sd"),
]
for label, col, kind in CONT_ROWS:
    rows.append({"Characteristic": label,
                 **{g: col_stats(sub, col, kind)
                    for g, sub in groups.items()}})

CAT_ROWS = [
    ("Female sex, %", "sex", ["Female"]),
    ("Race and Hispanic origin, %", "race_ethnicity",
     ["Mexican American", "Other Hispanic", "Non-Hispanic White",
      "Non-Hispanic Black", "Non-Hispanic Asian", "Other or multiracial"]),
    ("Education, %", "education",
     ["Less than 9th grade", "9th to 11th grade",
      "High school graduate or GED", "Some college or associate degree",
      "College graduate or above"]),
    ("Smoking status, %", "smoking_status", ["Never", "Former", "Current"]),
    ("Overweight or obesity criterion, %", "criterion_adiposity", [1]),
    ("Hyperglycaemia criterion, %", "criterion_hyperglycaemia", [1]),
    ("Hypertension criterion, %", "criterion_hypertension", [1]),
    ("Hypertriglyceridaemia criterion, %",
     "criterion_hypertriglyceridaemia", [1]),
    ("Low HDL criterion, %", "criterion_low_hdl", [1]),
    ("On antihypertensive treatment, %", "on_antihypertensive", [1]),
    ("On lipid-lowering treatment, %", "on_lipid_lowering", [1]),
]
for label, col, levels in CAT_ROWS:
    if len(levels) == 1 and levels[0] == 1:
        # binary indicator: single row, % with value 1
        rows.append({"Characteristic": label,
                     **{g: fmt_pct(wpct(sub[col].values == 1,
                                        sub["mec_weight_pooled"].values.astype(float)))
                        for g, sub in groups.items()}})
    else:
        rows.append({"Characteristic": label,
                     **{g: "" for g in groups}})
        for lv in levels:
            rows.append({"Characteristic": f"  {lv}",
                         **{g: fmt_pct(wpct(sub[col].values == lv,
                                            sub["mec_weight_pooled"].values.astype(float)))
                            for g, sub in groups.items()}})

t1 = pd.DataFrame(rows, columns=["Characteristic"] + list(groups))
t1.to_csv(HERE / "table1.csv", index=False)
print("table1 rows:", len(t1), flush=True)
print(t1.to_string(index=False), flush=True)

# ---------------------------------------------------------------- Figure 1
stages = [
    ("All interviewed participants", flow["all_releases"].iloc[0]),
    ("Age ≥20 years with positive\nMEC weight",
     flow["all_releases"].iloc[1]),
    ("Valid CAP measurement", flow["all_releases"].iloc[2]),
    ("Reliable elastography exam\n(LUAXSTAT = 1)",
     flow["all_releases"].iloc[3]),
    ("No viral hepatitis B or C", flow["all_releases"].iloc[5]),
    ("Alcohol within MASLD limits", flow["all_releases"].iloc[6]),
    ("Metabolic status determined\n(analysed cohort)",
     flow["all_releases"].iloc[7]),
]
fig, ax = plt.subplots(figsize=(7.2, 8.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
y = 9.4
for i, (label, n) in enumerate(stages):
    ax.add_patch(FancyBboxPatch((1.2, y - 0.55), 7.6, 0.85,
                                boxstyle="round,pad=0.08",
                                fc="#eef3f8" if i < len(stages) - 1
                                else "#dbe9d6",
                                ec="#1f4e79"))
    ax.text(5.0, y - 0.12, f"{label}\n n = {int(n):,}",
            ha="center", va="center", fontsize=10)
    if i < len(stages) - 1:
        ax.annotate("", xy=(5.0, y - 0.80), xytext=(5.0, y - 0.60),
                    arrowprops=dict(arrowstyle="->", color="#1f4e79"))
    y -= 1.15
ax.text(5.0, 0.35, "Development (2017–March 2020): n = 6,320 "
        "(3,709 MASLD)\nValidation (Aug 2021–Aug 2023): n = 4,103 "
        "(2,265 MASLD)", ha="center", va="center", fontsize=10,
        style="italic")
fig.tight_layout()
fig.savefig(FIG / "fig1_flow.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
from sklearn.metrics import roc_curve, roc_auc_score
y = pred["masld"].values
w = pred["mec_weight_pooled"].values.astype(float)
fig, ax = plt.subplots(figsize=(6.2, 5.4))
for col, label, color in [("p_logistic", "Logistic regression", "#1f4e79"),
                          ("p_gbm", "Gradient boosting", "#e67e22")]:
    p = pred[col].values
    fpr, tpr, _ = roc_curve(y, p, sample_weight=w)
    auc = roc_auc_score(y, p, sample_weight=w)
    ax.plot(fpr, tpr, color=color, lw=2, label=f"{label} (AUROC {auc:.3f})")
ax.plot([0, 1], [0, 1], color="grey", lw=1, ls="--")
ax.set_xlabel("1 − specificity")
ax.set_ylabel("Sensitivity")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.legend(loc="lower right", frameon=False)
ax.set_title("Receiver operating characteristic curves\n"
             "temporal validation (Aug 2021–Aug 2023, n = 4,103)")
fig.tight_layout()
fig.savefig(FIG / "fig2_roc.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3
fig, ax = plt.subplots(figsize=(6.2, 5.4))
for col, label, color, marker in [
        ("p_logistic", "Logistic regression", "#1f4e79", "o"),
        ("p_gbm", "Gradient boosting", "#e67e22", "s")]:
    p = pred[col].values
    qs = np.quantile(p, np.linspace(0, 1, 11))
    xs, ys = [], []
    for a, b in zip(qs[:-1], qs[1:]):
        m = (p >= a) & (p <= b if b == qs[-1] else p < b)
        xs.append(float(np.average(p[m], weights=w[m])))
        ys.append(float(np.average(y[m], weights=w[m])))
    ax.plot(xs, ys, color=color, lw=1.5, marker=marker, ms=5,
            label=label)
ax.plot([0, 1], [0, 1], color="grey", lw=1, ls="--",
        label="Perfect calibration")
ax.set_xlabel("Mean predicted probability (decile)")
ax.set_ylabel("Observed proportion with MASLD")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(loc="upper left", frameon=False)
ax.set_title("Calibration plots\ntemporal validation (Aug 2021–Aug 2023)")
fig.tight_layout()
fig.savefig(FIG / "fig3_calibration.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 4
o = ors.copy()
o = o.iloc[::-1]
fig, ax = plt.subplots(figsize=(7.6, 6.4))
ypos = np.arange(len(o))
log_or = np.log(o["odds_ratio"].values)
lo = np.log(o["or_ci95_lo"].values)
hi = np.log(o["or_ci95_hi"].values)
ax.errorbar(log_or, ypos, xerr=[log_or - lo, hi - log_or],
            fmt="o", color="#1f4e79", ecolor="#1f4e79", capsize=3, ms=5)
ax.axvline(0, color="grey", lw=1, ls="--")
ax.set_yticks(ypos)
ax.set_yticklabels(o["predictor"], fontsize=9)
ax.set_xticks([np.log(0.5), 0.0, np.log(2.0)])
ax.set_xticklabels(["0.5", "1", "2"])
ax.set_xlabel("Odds ratio (log scale)")
ax.set_xlim(np.log(0.4), np.log(3.0))
ax.set_title("Adjusted odds ratios from the logistic regression model\n"
             "(development sample; bars are 95% CIs from 1,500 bootstraps)")
fig.tight_layout()
fig.savefig(FIG / "fig4_forest.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 5
im = imp.sort_values("mean_brier_increase")
DISPLAY = {"alt_u_l": "ALT", "ggt_u_l": "GGT", "age_years": "Age",
           "ast_u_l": "AST", "mvpa_min_per_week": "Physical activity (MVPA)",
           "race_ethnicity": "Race/ethnicity", "smoking_status": "Smoking",
           "total_cholesterol_mg_dl": "Total cholesterol", "sex": "Sex",
           "education": "Education"}
im = im.assign(label=im["predictor"].map(DISPLAY))
fig, ax = plt.subplots(figsize=(7.0, 4.6))
ax.barh(im["label"], im["mean_brier_increase"], color="#1f4e79")
ax.set_xlabel("Mean increase in weighted Brier score when permuted\n"
              "(validation sample; 10 repeats)")
ax.set_title("Permutation variable importance — gradient boosting model")
fig.tight_layout()
fig.savefig(FIG / "fig5_importance.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 6
# Decision curve analysis (validation sample, survey-weighted net benefit)
def net_benefit(yv, pv, wv, pt):
    pred_pos = pv >= pt
    tp = float(np.sum(wv[pred_pos & (yv == 1)]))
    fp = float(np.sum(wv[pred_pos & (yv == 0)]))
    n = float(np.sum(wv))
    return (tp - fp * pt / (1 - pt)) / n


prev = float(np.average(y, weights=w))
ths = np.arange(0.05, 0.96, 0.01)
nb_all = prev - (1 - prev) * ths / (1 - ths)
fig, ax = plt.subplots(figsize=(6.2, 5.4))
ax.plot(ths, nb_all, color="grey", lw=1.5, ls="--", label="Treat all")
ax.axhline(0, color="black", lw=1, label="Treat none")
for col, label, color in [("p_logistic", "Logistic regression", "#1f4e79"),
                          ("p_gbm", "Gradient boosting", "#e67e22")]:
    nb = [net_benefit(y, pred[col].values, w, t) for t in ths]
    ax.plot(ths, nb, color=color, lw=2, label=label)
ax.set_xlabel("Threshold probability")
ax.set_ylabel("Net benefit")
ax.set_xlim(0.05, 0.95)
ax.set_ylim(-0.05, 0.6)
ax.legend(loc="lower left", frameon=False)
ax.set_title("Decision curve analysis\ntemporal validation "
             "(Aug 2021–Aug 2023, n = 4,103)")
fig.tight_layout()
fig.savefig(FIG / "fig6_dca.png", dpi=300, bbox_inches="tight")
plt.close(fig)

print("figures done:", sorted(p.name for p in FIG.glob("*.png")), flush=True)
