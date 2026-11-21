#!/usr/bin/env python3
"""
MASLD predictive modelling — NHANES 2017-March 2020 (development) vs
August 2021-August 2023 (temporal validation).

Protocol v1.2. Models: L2-penalised logistic regression vs histogram
gradient boosting, same 10 predictors. Tuning by 5-fold CV within
development on survey-weighted Brier. CIs from 1,500 rescaled Rao-Wu
bootstrap replicates (seed 1729). Odds-ratio CIs from 1,500 nonparametric
bootstrap refits of the logistic model (seed 1729).

Outputs (same directory as this script):
  model_results.json, metrics_table.csv, odds_ratios.csv,
  permutation_importance.csv, tuning_results.csv,
  validation_predictions.csv, model_results_summary.md
"""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import expit, logit
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, brier_score_loss,
                             roc_auc_score, roc_curve)
from sklearn.model_selection import StratifiedKFold

warnings.filterwarnings("ignore")

HERE = Path(__file__).resolve().parent
DATA = HERE / "nhanes_masld_cross_sectional.csv"
SEED = 1729
N_BOOT = 1500
DEV = "2017-March2020"
VAL = "2021-2023"

# ---- predictor specification -------------------------------------------
# (column, scale divisor, label) — continuous terms centred at the
# training median and scaled to clinically interpretable units.
CONT = [
    ("age_years", 10.0, "Age, per 10 years"),
    ("mvpa_min_per_week", 60.0, "MVPA, per 60 min/week"),
    ("total_cholesterol_mg_dl", 20.0, "Total cholesterol, per 20 mg/dL"),
    ("alt_u_l", 10.0, "ALT, per 10 U/L"),
    ("ast_u_l", 10.0, "AST, per 10 U/L"),
    ("ggt_u_l", 10.0, "GGT, per 10 U/L"),
]
# (column, reference category, label prefix)
CAT = [
    ("sex", "Female", "Sex"),
    ("race_ethnicity", "Non-Hispanic White", "Race/ethnicity"),
    ("smoking_status", "Never", "Smoking"),
    ("education", "College graduate or above", "Education"),
]
# extra predictors for the routine-clinical sensitivity analysis:
# BMI and systolic blood pressure only — measured free at every clinic
# visit, no fasting blood sample required.
CONT_ROUTINE = [
    ("bmi_kg_m2", 5.0, "BMI, per 5 kg/m2"),
    ("mean_systolic_bp_mm_hg", 10.0, "Systolic BP, per 10 mmHg"),
]
# extra predictors for the criterion-containing sensitivity analysis
CONT_EXTRA = [
    ("bmi_kg_m2", 5.0, "BMI, per 5 kg/m2"),
    ("waist_cm", 10.0, "Waist, per 10 cm"),
    ("hba1c_pct", 1.0, "HbA1c, per 1%"),
    ("triglycerides_mg_dl", 50.0, "Triglycerides, per 50 mg/dL"),
    ("hdl_cholesterol_mg_dl", 10.0, "HDL cholesterol, per 10 mg/dL"),
    ("mean_systolic_bp_mm_hg", 10.0, "Systolic BP, per 10 mmHg"),
]


def weighted_median(values, weights):
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = ~np.isnan(v)
    v, w = v[ok], w[ok]
    order = np.argsort(v)
    v, w = v[order], w[order]
    cw = np.cumsum(w)
    return float(v[np.searchsorted(cw, 0.5 * cw[-1])])


def weighted_mode(values, weights):
    s = pd.Series(values)
    w = pd.Series(weights, index=s.index)
    ok = s.notna()
    return s[ok].groupby(s[ok]).apply(lambda idx: w.loc[idx.index].sum()).idxmax()


class Preprocessor:
    """Weighted-median/mode imputation, centring at the training median,
    one-hot encoding with a fixed reference category."""

    def __init__(self, cont, cat):
        self.cont = cont
        self.cat = cat
        self.names_ = None

    def fit(self, df, w):
        self.med_ = {}
        self.centre_ = {}
        for col, scale, _ in self.cont:
            self.med_[col] = weighted_median(df[col].values, w)
            self.centre_[col] = weighted_median(
                (df[col].fillna(self.med_[col]) / scale).values, w)
        self.mode_ = {}
        self.levels_ = {}
        for col, ref, _ in self.cat:
            self.mode_[col] = weighted_mode(df[col].values, w)
            levels = sorted(df[col].dropna().unique().tolist())
            self.levels_[col] = [lv for lv in levels if lv != ref]
        names = [lbl for _, _, lbl in self.cont]
        for col, ref, prefix in self.cat:
            for lv in self.levels_[col]:
                names.append(f"{prefix}: {lv} (vs {ref})")
        self.names_ = names
        return self

    def transform(self, df):
        cols = []
        for col, scale, _ in self.cont:
            cols.append(df[col].fillna(self.med_[col]).values / scale
                        - self.centre_[col])
        for col, ref, _ in self.cat:
            filled = df[col].fillna(self.mode_[col])
            for lv in self.levels_[col]:
                cols.append((filled == lv).astype(float).values)
        return np.column_stack(cols)


# ---- metrics ------------------------------------------------------------
def weighted_brier(y, p, w):
    return float(np.average((p - y) ** 2, weights=w))


def calibration_intercept(y, p, w):
    """Intercept of y ~ offset(logit p), weighted MLE (brentq root)."""
    lp = logit(np.clip(p, 1e-6, 1 - 1e-6))
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    def score(a):
        return float(np.sum(w * (y - expit(a + lp))))

    lo, hi = -10.0, 10.0
    if score(lo) * score(hi) > 0:
        return float(np.nan)
    return float(brentq(score, lo, hi))


def calibration_slope(y, p, w):
    """Slope of y ~ a + b*logit(p), weighted MLE."""
    lp = logit(np.clip(p, 1e-6, 1 - 1e-6))
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)

    def nll(ab):
        a, b = ab
        eta = a + b * lp
        ll = y * np.log(expit(eta) + 1e-300) + (1 - y) * np.log(
            expit(-eta) + 1e-300)
        return -float(np.sum(w * ll))

    res = minimize(nll, x0=np.array([0.0, 1.0]), method="Nelder-Mead",
                   options={"maxiter": 500, "xatol": 1e-4, "fatol": 1e-4})
    return float(res.x[1])


def metrics(y, p, w):
    return {
        "auroc": float(roc_auc_score(y, p, sample_weight=w)),
        "average_precision": float(average_precision_score(
            y, p, sample_weight=w)),
        "brier": weighted_brier(y, p, w),
        "calibration_intercept": calibration_intercept(y, p, w),
        "calibration_slope": calibration_slope(y, p, w),
        "mean_predicted": float(np.average(p, weights=w)),
        "observed": float(np.average(y, weights=w)),
    }


def threshold_metrics(y, p, w, thr):
    pred = (p >= thr).astype(int)
    tp = float(np.sum(w * (pred == 1) * (y == 1)))
    fp = float(np.sum(w * (pred == 1) * (y == 0)))
    tn = float(np.sum(w * (pred == 0) * (y == 0)))
    fn = float(np.sum(w * (pred == 0) * (y == 1)))
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    ppv = tp / (tp + fp) if tp + fp else float("nan")
    npv = tn / (tn + fn) if tn + fn else float("nan")
    return {"threshold": float(thr), "sensitivity": sens, "specificity": spec,
            "ppv": ppv, "npv": npv}


def youden_threshold(y, p, w):
    fpr, tpr, thr = roc_curve(y, p, sample_weight=w)
    j = tpr - fpr
    i = int(np.argmax(j))
    return float(thr[i]), float(tpr[i]), float(1 - fpr[i])


def sens90_threshold(y, p, w):
    fpr, tpr, thr = roc_curve(y, p, sample_weight=w)
    ok = np.where(tpr >= 0.90)[0]
    if len(ok) == 0:
        return float(thr[-1]), float(tpr[-1]), float(1 - fpr[-1])
    # most specific threshold still achieving >=90% sensitivity
    i = int(ok[np.argmax(thr[ok])])
    return float(thr[i]), float(tpr[i]), float(1 - fpr[i])


# ---- Rao-Wu bootstrap ----------------------------------------------------
def rao_wu_multipliers(df_val, rng):
    """Per-row weight multipliers for one rescaled Rao-Wu replicate:
    within each stratum, draw (n_h - 1) PSUs with replacement."""
    mult = np.ones(len(df_val))
    strata = df_val["sdmvstra"].values
    psu = df_val["sdmvpsu"].values
    for h in np.unique(strata):
        idx = np.where(strata == h)[0]
        psus = np.unique(psu[idx])
        n_h = len(psus)
        if n_h <= 1:
            continue
        drawn, counts = np.unique(
            rng.choice(psus, size=n_h - 1, replace=True),
            return_counts=True)
        factor = n_h / (n_h - 1)
        for d, c in zip(drawn, counts):
            mult[idx[psu[idx] == d]] = factor * c
    return mult


def bootstrap_metrics(y, p_lr, p_gbm, w, df_val, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    out = {"lr": [], "gbm": []}
    for _ in range(n_boot):
        wb = w * rao_wu_multipliers(df_val, rng)
        out["lr"].append(metrics(y, p_lr, wb))
        out["gbm"].append(metrics(y, p_gbm, wb))
    return out


def ci(values):
    lo, hi = np.percentile(values, [2.5, 97.5])
    return float(lo), float(hi)


# ---- main ----------------------------------------------------------------
def main():
    df = pd.read_csv(DATA)
    dev = df[df["survey_release"] == DEV].copy()
    val = df[df["survey_release"] == VAL].copy()
    w_dev = dev["mec_weight_pooled"].values.astype(float)
    w_val = val["mec_weight_pooled"].values.astype(float)
    y_dev = dev["masld"].values.astype(int)
    y_val = val["masld"].values.astype(int)
    print(f"dev n={len(dev)} events={y_dev.sum()} | "
          f"val n={len(val)} events={y_val.sum()}", flush=True)

    # ---- tuning: 5-fold CV within development, weighted Brier -----------
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    tuning_rows = []

    def cv_brier(make_model, pre):
        scores = []
        for tr, te in skf.split(dev, y_dev):
            pp = Preprocessor(*pre).fit(dev.iloc[tr], w_dev[tr])
            Xtr = pp.transform(dev.iloc[tr])
            Xte = pp.transform(dev.iloc[te])
            m = make_model()
            m.fit(Xtr, y_dev[tr], sample_weight=w_dev[tr])
            p = m.predict_proba(Xte)[:, 1]
            scores.append(weighted_brier(y_dev[te], p, w_dev[te]))
        return float(np.mean(scores))

    pre10 = (CONT, CAT)
    lr_grid = [0.01, 0.1, 1.0, 10.0]
    best_lr = (None, np.inf)
    for c in lr_grid:
        s = cv_brier(lambda c=c: LogisticRegression(C=c, max_iter=2000), pre10)
        tuning_rows.append({"model": "logistic", "params": f"C={c}",
                            "cv_brier": s})
        print(f"LR C={c}: CV Brier {s:.4f}", flush=True)
        if s < best_lr[1]:
            best_lr = (c, s)
    gbm_grid = [(l, r, it, l2) for l in (7, 15, 31) for r in (0.05, 0.1)
                for it in (100, 200) for l2 in (0.0, 1.0)]
    best_gbm = (None, np.inf)
    for l, r, it, l2 in gbm_grid:
        s = cv_brier(lambda l=l, r=r, it=it, l2=l2:
                     HistGradientBoostingClassifier(
                         max_leaf_nodes=l, learning_rate=r, max_iter=it,
                         min_samples_leaf=30, l2_regularization=l2,
                         random_state=SEED), pre10)
        tuning_rows.append({"model": "gbm",
                            "params": f"leaves={l}, lr={r}, iter={it}, l2={l2}",
                            "cv_brier": s})
        print(f"GBM leaves={l} lr={r} iter={it} l2={l2}: CV Brier {s:.4f}",
              flush=True)
        if s < best_gbm[1]:
            best_gbm = ((l, r, it, l2), s)
    print(f"BEST LR C={best_lr[0]} (CV Brier {best_lr[1]:.4f}) | "
          f"BEST GBM {best_gbm[0]} (CV Brier {best_gbm[1]:.4f})", flush=True)

    # ---- refit on all development, assess on validation -----------------
    pre = Preprocessor(CONT, CAT).fit(dev, w_dev)
    X_dev, X_val = pre.transform(dev), pre.transform(val)
    lr = LogisticRegression(C=best_lr[0], max_iter=2000)
    lr.fit(X_dev, y_dev, sample_weight=w_dev)
    l, r, it, l2 = best_gbm[0]
    gbm = HistGradientBoostingClassifier(
        max_leaf_nodes=l, learning_rate=r, max_iter=it, min_samples_leaf=30,
        l2_regularization=l2, random_state=SEED)
    gbm.fit(X_dev, y_dev, sample_weight=w_dev)
    p_lr = lr.predict_proba(X_val)[:, 1]
    p_gbm = gbm.predict_proba(X_val)[:, 1]

    # benchmark: constant development prevalence
    p_const = np.full(len(val), float(np.average(y_dev, weights=w_dev)))
    bench_brier = weighted_brier(y_val, p_const, w_val)

    res = {"n_dev": int(len(dev)), "n_val": int(len(val)),
           "events_dev": int(y_dev.sum()), "events_val": int(y_val.sum()),
           "best_lr_C": best_lr[0], "best_lr_cv_brier": best_lr[1],
           "best_gbm": {"max_leaf_nodes": l, "learning_rate": r,
                        "max_iter": it, "l2_regularization": l2,
                        "cv_brier": best_gbm[1]},
           "benchmark_brier": bench_brier,
           "seed": SEED, "n_boot": N_BOOT}

    m_lr, m_gbm = metrics(y_val, p_lr, w_val), metrics(y_val, p_gbm, w_val)
    res["validation"] = {"logistic": m_lr, "gbm": m_gbm}
    print("LR val:", {k: round(v, 4) for k, v in m_lr.items()}, flush=True)
    print("GBM val:", {k: round(v, 4) for k, v in m_gbm.items()}, flush=True)

    # threshold metrics (point estimates, thresholds from full validation)
    thr_rows = []
    for name, p in [("logistic", p_lr), ("gbm", p_gbm)]:
        t_y, s_y, sp_y = youden_threshold(y_val, p, w_val)
        t_90, s_90, sp_90 = sens90_threshold(y_val, p, w_val)
        tm_y = threshold_metrics(y_val, p, w_val, t_y)
        tm_y["rule"] = "Youden J"
        tm_90 = threshold_metrics(y_val, p, w_val, t_90)
        tm_90["rule"] = "90% sensitivity"
        thr_rows += [{"model": name, **tm_y}, {"model": name, **tm_90}]
        res.setdefault("thresholds", {})[name] = {
            "youden": tm_y, "sens90": tm_90}
    print("thresholds done", flush=True)

    # ---- Rao-Wu bootstrap CIs for validation metrics --------------------
    boots = bootstrap_metrics(y_val, p_lr, p_gbm, w_val, val)
    for name in ("lr", "gbm"):
        key = "logistic" if name == "lr" else "gbm"
        cis = {}
        for k in ("auroc", "average_precision", "brier",
                  "calibration_intercept", "calibration_slope"):
            cis[k] = ci([b[k] for b in boots[name]])
        res["validation"][key]["ci95"] = cis
        # threshold-metric CIs at the fixed point-estimate thresholds
        t_y = res["thresholds"][key]["youden"]["threshold"]
        t_90 = res["thresholds"][key]["sens90"]["threshold"]
        p = p_lr if name == "lr" else p_gbm
        rng = np.random.default_rng(SEED + 1)
        for tlabel, t in [("youden", t_y), ("sens90", t_90)]:
            sens_b, spec_b, ppv_b, npv_b = [], [], [], []
            for _ in range(N_BOOT):
                wb = w_val * rao_wu_multipliers(val, rng)
                tm = threshold_metrics(y_val, p, wb, t)
                sens_b.append(tm["sensitivity"]); spec_b.append(tm["specificity"])
                ppv_b.append(tm["ppv"]); npv_b.append(tm["npv"])
            res["thresholds"][key][tlabel]["ci95"] = {
                "sensitivity": ci(sens_b), "specificity": ci(spec_b),
                "ppv": ci(ppv_b), "npv": ci(npv_b)}
    print("bootstrap CIs done", flush=True)

    # ---- odds ratios with bootstrap CIs ---------------------------------
    rng = np.random.default_rng(SEED + 2)
    coefs = []
    for b in range(N_BOOT):
        idx = rng.integers(0, len(dev), len(dev))
        d_b, y_b, w_b = dev.iloc[idx], y_dev[idx], w_dev[idx]
        pp = Preprocessor(CONT, CAT).fit(d_b, w_b)
        m = LogisticRegression(C=best_lr[0], max_iter=2000)
        m.fit(pp.transform(d_b), y_b, sample_weight=w_b)
        coefs.append(m.coef_[0])
        if (b + 1) % 300 == 0:
            print(f"OR bootstrap {b + 1}/{N_BOOT}", flush=True)
    coefs = np.array(coefs)
    or_rows = []
    for j, name in enumerate(pre.names_):
        lo, hi = np.percentile(coefs[:, j], [2.5, 97.5])
        or_rows.append({"predictor": name,
                        "log_odds": float(lr.coef_[0][j]),
                        "odds_ratio": float(np.exp(lr.coef_[0][j])),
                        "or_ci95_lo": float(np.exp(lo)),
                        "or_ci95_hi": float(np.exp(hi))})
    res["logistic_intercept"] = float(lr.intercept_[0])
    print("OR bootstrap done", flush=True)

    # ---- GBM permutation importance (validation, weighted Brier) ---------
    base_brier = weighted_brier(y_val, p_gbm, w_val)
    rng = np.random.default_rng(SEED + 3)
    imp_rows = []
    raw_cols = [c for c, _, _ in CONT] + [c for c, _, _ in CAT]
    for col in raw_cols:
        diffs = []
        for _ in range(10):
            v2 = val.copy()
            v2[col] = rng.permutation(v2[col].values)
            p2 = gbm.predict_proba(pre.transform(v2))[:, 1]
            diffs.append(weighted_brier(y_val, p2, w_val) - base_brier)
        imp_rows.append({"predictor": col,
                         "mean_brier_increase": float(np.mean(diffs)),
                         "sd_brier_increase": float(np.std(diffs))})
    imp_rows.sort(key=lambda d: -d["mean_brier_increase"])
    print("permutation importance done", flush=True)

    # ---- secondary outcomes & sensitivity analyses -----------------------
    def run_variant(y_d, y_v, d_dev, d_val, w_d, w_v, pre_spec,
                    label, note=""):
        pp = Preprocessor(*pre_spec).fit(d_dev, w_d)
        Xd, Xv = pp.transform(d_dev), pp.transform(d_val)
        m1 = LogisticRegression(C=best_lr[0], max_iter=2000)
        m1.fit(Xd, y_d, sample_weight=w_d)
        l_, r_, it_, l2_ = best_gbm[0]
        m2 = HistGradientBoostingClassifier(
            max_leaf_nodes=l_, learning_rate=r_, max_iter=it_,
            min_samples_leaf=30, l2_regularization=l2_, random_state=SEED)
        m2.fit(Xd, y_d, sample_weight=w_d)
        out = {"label": label, "note": note,
               "n_dev": int(len(d_dev)), "n_val": int(len(d_val)),
               "events_dev": int(y_d.sum()), "events_val": int(y_v.sum()),
               "logistic": metrics(y_v, m1.predict_proba(Xv)[:, 1], w_v),
               "gbm": metrics(y_v, m2.predict_proba(Xv)[:, 1], w_v)}
        print(label, "LR AUROC",
              round(out["logistic"]["auroc"], 3), "GBM AUROC",
              round(out["gbm"]["auroc"], 3), flush=True)
        return out

    variants = []
    variants.append(run_variant(
        dev["steatosis"].values.astype(int),
        val["steatosis"].values.astype(int), dev, val, w_dev, w_val, pre10,
        "Steatosis alone (CAP >=248)"))
    variants.append(run_variant(
        dev["masld_with_significant_fibrosis"].values.astype(int),
        val["masld_with_significant_fibrosis"].values.astype(int),
        dev, val, w_dev, w_val, pre10,
        "MASLD with significant fibrosis (LSM >=8 kPa)"))
    masld285_dev = ((dev["cap_db_m"] >= 285)
                    & (dev["n_metabolic_criteria_met"] >= 1)).astype(int).values
    masld285_val = ((val["cap_db_m"] >= 285)
                    & (val["n_metabolic_criteria_met"] >= 1)).astype(int).values
    variants.append(run_variant(masld285_dev, masld285_val, dev, val,
                                w_dev, w_val, pre10,
                                "MASLD at CAP >=285 dB/m"))
    fast_dev = dev[dev["fatty_liver_index"].notna()]
    fast_val = val[val["fatty_liver_index"].notna()]
    fli_dev = ((fast_dev["fli_ge_60"] == 1)
               & (fast_dev["n_metabolic_criteria_met"] >= 1)).astype(int).values
    fli_val = ((fast_val["fli_ge_60"] == 1)
               & (fast_val["n_metabolic_criteria_met"] >= 1)).astype(int).values
    variants.append(run_variant(
        fli_dev, fli_val, fast_dev, fast_val,
        fast_dev["mec_weight_pooled"].values.astype(float),
        fast_val["mec_weight_pooled"].values.astype(float), pre10,
        "MASLD by Fatty Liver Index >=60",
        note="fasting subsample only"))
    pre_crit = (CONT + CONT_EXTRA, CAT)
    variants.append(run_variant(
        y_dev, y_val, dev, val, w_dev, w_val, pre_crit,
        "Criterion-containing predictor set",
        note="adds BMI, waist, HbA1c, triglycerides, HDL, systolic BP"))
    pre_routine = (CONT + CONT_ROUTINE, CAT)
    variants.append(run_variant(
        y_dev, y_val, dev, val, w_dev, w_val, pre_routine,
        "Routine-clinical predictor set",
        note="adds BMI and systolic BP only (no fasting labs)"))
    cc = [c for c, _, _ in CONT] + [c for c, _, _ in CAT]
    cc_dev = dev.dropna(subset=cc)
    cc_val = val.dropna(subset=cc)
    variants.append(run_variant(
        cc_dev["masld"].values.astype(int),
        cc_val["masld"].values.astype(int), cc_dev, cc_val,
        cc_dev["mec_weight_pooled"].values.astype(float),
        cc_val["mec_weight_pooled"].values.astype(float), pre10,
        "Complete-case analysis"))
    variants.append(run_variant(
        y_dev, y_val, dev, val,
        np.ones(len(dev)), np.ones(len(val)), pre10,
        "Unweighted analysis"))
    res["variants"] = variants

    # ---- save outputs ----------------------------------------------------
    pd.DataFrame(tuning_rows).to_csv(HERE / "tuning_results.csv", index=False)
    pd.DataFrame(or_rows).to_csv(HERE / "odds_ratios.csv", index=False)
    pd.DataFrame(imp_rows).to_csv(HERE / "permutation_importance.csv",
                                  index=False)
    pd.DataFrame(thr_rows).to_csv(HERE / "threshold_metrics.csv", index=False)
    pred = val[["seqn", "survey_release", "sdmvstra", "sdmvpsu",
                "mec_weight_pooled"]].copy()
    pred["masld"] = y_val
    pred["p_logistic"] = p_lr
    pred["p_gbm"] = p_gbm
    pred.to_csv(HERE / "validation_predictions.csv", index=False)
    with open(HERE / "model_results.json", "w") as f:
        json.dump(res, f, indent=2, default=float)
    print("saved outputs", flush=True)
    return res


if __name__ == "__main__":
    main()
