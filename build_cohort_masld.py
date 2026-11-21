"""Rebuild the NHANES cross-sectional MASLD dataset and feasibility summaries.

Design (protocol v1.2): cross-sectional analysis of NHANES, two releases.
    Development sample: NHANES 2017-March 2020 pre-pandemic (P_ files, 3.2 years,
    MEC weight WTMECPRP). Temporal validation sample: NHANES August 2021-August
    2023 (L files, 2.0 years, MEC weight WTMEC2YR).

Why there is no 2017-2018 / 2019-March 2020 split: CDC reassigned SEQN values
in the combined 2017-March 2020 files (P_DEMO SEQN range 109263-124822, while
DEMO_J runs 93703-102956, zero overlap), regenerated the masked variance
strata/PSU, and set SDDSRVYR to 66 for every row. The 2019-March 2020
participants therefore cannot be identified inside the P_ files, so the
development sample is the combined pre-pandemic sample exactly as released.
Protocol v1.1's SEQN-subsetting instruction was wrong and is superseded here.

Outcome: MASLD = CAP >= 248 dB/m plus at least one of five cardiometabolic
criteria (Rinella et al. 2023), alcohol within limits (<=2 drinks/day women,
<=3 men), and no viral hepatitis B/C. A reliable elastography exam is
LUAXSTAT == 1 (fasted >= 3 h, >= 10 complete stiffness measures, liver
stiffness IQR/median < 30%).

Data fixes carried over from the CKD review, corrected 2026-11-21:
  * The value 5.397605347e-79 (= 2^-260) that appears in these files is NOT a
    missing-value code: it is how pandas' XPORT reader returns a SAS zero.
    Proved by exact count matches against the CDC codebooks (ALQ121 code 0,
    "Never in the last year", n = 1,049, equals the number of 2^-260 values in
    ALQ121; PAD790Q/PAD810Q zeros are the "never" answers). Such values are
    therefore restored to 0.0, not to missing. (The CKD review had called the
    same value a missing-value artifact; in the CKD data the DBP zeros are
    recorded 0 readings and the PIR zeros are genuine zero income, so the
    same correction applies there.)
  * Blood pressure: the mean is taken over readings > 0 for BOTH systolic and
    diastolic (the CKD build kept zero diastolic readings; here the rule is
    symmetric and disclosed).
  * Alcohol: ALQ130 is drinks per DRINKING day, not per calendar day, so the
    average is ALQ130 x drinking days per year / 365, with drinking frequency
    from ALQ121 category midpoints (every day 365, nearly every day 312,
    3-4 times/week 182, 2 times/week 104, once/week 52, 2-3 times/month 30,
    once/month 12, 7-11 times/year 9, 3-6 times/year 4.5, 1-2 times/year 1.5).
    Never-drinkers (ALQ111 = 2, or ALQ121 = 0) are coded 0. Using ALQ130
    directly as drinks/day would wrongly exclude about 40% of the sample;
    the corrected rule excludes about 3-5%.

Run from any directory with Python 3, pandas and numpy installed:
    python build_cohort_masld.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RAW = ROOT.parent / "raw"
CYCLES = {
    "2017-March2020": ("2017", "P", 3.2),
    "2021-2023": ("2021", "L", 2.0),
}
TOTAL_YEARS = sum(c[2] for c in CYCLES.values())  # 5.2

CORE = [
    "age_years", "sex", "race_ethnicity", "smoking_status", "education",
    "mvpa_min_per_week", "total_cholesterol_mg_dl",
    "alt_u_l", "ast_u_l", "ggt_u_l",
]
OPTIONAL = [
    "waist_cm", "bmi_kg_m2", "hba1c_pct", "triglycerides_mg_dl",
    "hdl_cholesterol_mg_dl", "mean_systolic_bp_mm_hg",
    "mean_diastolic_bp_mm_hg", "poverty_income_ratio", "uric_acid_mg_dl",
    "fatty_liver_index",
]
OUTCOMES = [
    "cap_db_m", "liver_stiffness_kpa", "steatosis",
    "masld", "masld_with_significant_fibrosis", "fli_ge_60",
]
CRITERIA = [
    "criterion_adiposity", "criterion_hyperglycaemia",
    "criterion_hypertension", "criterion_hypertriglyceridaemia",
    "criterion_low_hdl",
]
YES_NO = {1: 1, 2: 0}
ARTIFACT_LIMIT = 1e-30

# component -> variables to read (SEQN is always added). BP and PAQ differ by cycle.
FILES = {
    "DEMO": "RIDAGEYR RIAGENDR RIDRETH3 RIDEXPRG DMDEDUC2 INDFMPIR SDMVSTRA SDMVPSU",
    "LUX": "LUAXSTAT LUXCAPM LUXSMED LUXSIQRM LUANMVGP",
    "BMX": "BMXBMI BMXWAIST",
    "BIOPRO": "LBXSATSI LBXSASSI LBXSGTSI LBXSUA",
    "TRIGLY": "LBXTR",
    "HDL": "LBDHDD",
    "TCHOL": "LBXTC",
    "GLU": "LBXGLU",
    "GHB": "LBXGH",
    "DIQ": "DIQ010 DIQ050 DIQ070",
    "SMQ": "SMQ020 SMQ040",
    "ALQ": "ALQ111 ALQ121 ALQ130",
    "HEPBD": "LBDHBG",
    "HEPC": "LBXHCR",
}
BPQ_VARS = {
    # BPQ050A = BP medication; BPQ100D = cholesterol medication; BPQ080 identifies
    # the skip pattern (not told to take the medication -> not taking it).
    "P": "BPQ020 BPQ050A BPQ100D BPQ080",
    # BPQ150 / BPQ101D are the 2021-2023 renames (BPQ_L analytic notes)
    "L": "BPQ020 BPQ150 BPQ101D BPQ080",
}
PAQ_VARS = {
    # GPAQ leisure-time recreational activity (2017-March 2020)
    "P": "PAQ650 PAQ655 PAD660 PAQ665 PAQ670 PAD675",
    # Leisure-time physical activity, frequency + minutes per session (2021-2023)
    "L": "PAD790Q PAD790U PAD800 PAD810Q PAD810U PAD820",
}


def source_name(component, tag):
    return f"P_{component}" if tag == "P" else f"{component}_{tag}"


def recode_sas_zero(data, cycle, component, log):
    """Restore SAS zeros misread by pandas' XPORT reader as 2^-260.

    The value 5.397605347e-79 (= 2^-260) is the only sub-1e-30 value in these
    files, and codebook counts prove it is a zero, not a missing code (e.g.
    ALQ121 code 0, "Never in the last year", n = 1,049, exactly equals the
    number of 2^-260 values in ALQ121). It is therefore restored to 0.0.
    Genuine SAS missing values already read as NaN and are untouched.
    """
    num = data.select_dtypes(include=[np.number])
    for col in num.columns:
        bad = num[col].notna() & (num[col] != 0) & (num[col].abs() < ARTIFACT_LIMIT)
        if bad.any():
            log.append([cycle, component, col, int(bad.sum())])
            data.loc[bad, col] = 0.0
    return data


def weekly_minutes(days, minutes):
    """Minutes per week from days/week (1-7) and minutes/day (10-840)."""
    days = days.where(days.between(1, 7))
    minutes = minutes.where(minutes.between(10, 1440))
    return days * minutes


def read_cycle(cycle, year, tag, years, artifact_log):
    folder = RAW / ("2019-March2020" if tag == "P" else "2021-2023")
    weight = "WTMECPRP" if tag == "P" else "WTMEC2YR"
    fields = dict(FILES)
    fields["DEMO"] = FILES["DEMO"] + " " + weight
    fields["BPQ"] = BPQ_VARS[tag]
    fields["PAQ"] = PAQ_VARS[tag]
    bp_component = "BPXO"
    fields[bp_component] = "BPXOSY1 BPXOSY2 BPXOSY3 BPXODI1 BPXODI2 BPXODI3"
    parts = []
    sources = []
    for component, names in fields.items():
        stem = source_name(component, tag)
        path = folder / f"{stem}.xpt"
        # Triglycerides were renamed LBXTR -> LBXTLG in the 2021-2023 cycle.
        raw_names = names.replace("LBXTR", "LBXTLG") if (component == "TRIGLY" and tag == "L") else names
        data = pd.read_sas(path, format="xport")[["SEQN", *raw_names.split()]]
        if component == "TRIGLY" and tag == "L":
            data = data.rename(columns={"LBXTLG": "LBXTR"})
        assert data["SEQN"].notna().all() and data["SEQN"].is_unique, path
        data = recode_sas_zero(data, cycle, component, artifact_log)
        parts.append(data)
        url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{year}/DataFiles/{stem}"
        sources.append([cycle, component, f"{url}.xpt", f"{url}.htm", len(data)])
    merged = parts[0]
    for part in parts[1:]:
        assert part["SEQN"].isin(merged["SEQN"]).all(), (cycle, part.columns.tolist()[:2])
        merged = merged.merge(part, how="left", on="SEQN", validate="one_to_one")
    merged["survey_release"] = cycle
    merged["cycle_tag"] = tag
    merged["mec_weight_original"] = merged[weight]
    merged["mec_weight_pooled"] = merged[weight] * years / TOTAL_YEARS
    sbp = ["BPXOSY1", "BPXOSY2", "BPXOSY3"]
    dbp = ["BPXODI1", "BPXODI2", "BPXODI3"]
    # Symmetric rule: mean over readings > 0 for both systolic and diastolic.
    merged["mean_systolic_bp_mm_hg"] = merged[sbp].where(merged[sbp].gt(0)).mean(axis=1)
    merged["mean_diastolic_bp_mm_hg"] = merged[dbp].where(merged[dbp].gt(0)).mean(axis=1)
    merged["bp_protocol"] = "oscillometric"
    merged["bp_med_var"] = "BPQ050A" if tag == "P" else "BPQ150"
    merged["chol_med_var"] = "BPQ100D" if tag == "P" else "BPQ101D"
    return merged, sources


def mvpa_minutes(data):
    """Leisure-time moderate-to-vigorous physical activity, minutes per week.

    2017-March 2020 (GPAQ): vigorous recreational days x minutes/day plus
    moderate recreational days x minutes/day. A screener answer of No gives 0
    for that domain; a Yes with missing days/minutes leaves the domain missing.
    2021-2023 (LTPA): frequency (times per day/week/month/year) x minutes per
    session, separately for moderate (PAD790Q/U, PAD800) and vigorous
    (PAD810Q/U, PAD820) activity. A frequency of 0 gives 0 for that domain.
    """
    out = pd.Series(np.nan, index=data.index, dtype=float)
    is_p = data["cycle_tag"].eq("P")
    if is_p.any():
        sub = data.loc[is_p]
        out.loc[is_p] = _mvpa_gpaq(sub)
    if (~is_p).any():
        sub = data.loc[~is_p]
        out.loc[~is_p] = _mvpa_ltpa(sub)
    return out


def _domain_gpaq(screener, days, minutes):
    wk = weekly_minutes(days, minutes)
    result = pd.Series(np.nan, index=screener.index, dtype=float)
    result = result.mask(screener.eq(2), 0.0)   # No -> 0 minutes
    result = result.mask(screener.eq(1), wk)    # Yes -> days x minutes (NaN if either missing)
    return result


def _mvpa_gpaq(data):
    vig = _domain_gpaq(data["PAQ650"], data["PAQ655"], data["PAD660"])
    mod = _domain_gpaq(data["PAQ665"], data["PAQ670"], data["PAD675"])
    return vig + mod


def _domain_ltpa(freq, unit, minutes):
    unit = unit.map(lambda x: x.decode("utf-8", "ignore") if isinstance(x, bytes) else x)
    unit = unit.astype(str).str.strip().str.upper()
    factor = unit.map({"D": 7.0, "W": 1.0, "M": 12.0 / 52.0, "Y": 1.0 / 52.0})
    freq = freq.where(freq.between(0, 365))
    minutes = minutes.where(minutes.between(1, 1440))
    out = freq * factor * minutes
    out = out.where(freq.notna() & factor.notna() & minutes.notna())
    return out.mask(freq.eq(0), 0.0)


def _mvpa_ltpa(data):
    mod = _domain_ltpa(data["PAD790Q"], data["PAD790U"], data["PAD800"])
    vig = _domain_ltpa(data["PAD810Q"], data["PAD810U"], data["PAD820"])
    return mod + vig


def select_cohort(data):
    counts = []

    def record(stage, rows):
        counts.append([stage, len(rows), *(rows["survey_release"].eq(k).sum() for k in CYCLES)])

    record("All interviewed participants", data)
    data = data.loc[data["RIDAGEYR"].ge(20) & data["mec_weight_original"].gt(0)].copy()
    record("Age 20 or older with positive examination weight", data)
    data = data.loc[data["LUXCAPM"].notna()].copy()
    record("Valid CAP measurement", data)
    data = data.loc[data["LUAXSTAT"].eq(1)].copy()
    record("Reliable elastography exam (LUAXSTAT = 1)", data)
    data = data.loc[~data["RIDEXPRG"].eq(1)].copy()
    record("After excluding survey-identified pregnancy", data)
    hep_b = data["LBDHBG"].eq(1)
    hep_c = data["LBXHCR"].eq(1)
    hep_unknown = data["LBDHBG"].isna() & data["LBXHCR"].isna()
    data = data.loc[~(hep_b | hep_c | hep_unknown)].copy()
    record("After excluding viral hepatitis B/C", data)
    # ALQ130 is drinks per DRINKING day; average drinks/day = ALQ130 x
    # drinking days per year / 365, with frequency from ALQ121 midpoints.
    days_per_year = data["ALQ121"].map({
        0: 0.0, 1: 365.0, 2: 312.0, 3: 182.0, 4: 104.0, 5: 52.0,
        6: 30.0, 7: 12.0, 8: 9.0, 9: 4.5, 10: 1.5,
    })
    data["drinks_per_day"] = (
        data["ALQ130"].where(data["ALQ130"].lt(777)) * days_per_year / 365.0
    )
    data.loc[data["ALQ121"].eq(0) | data["ALQ111"].eq(2), "drinks_per_day"] = 0.0
    limit = np.where(data["RIAGENDR"].eq(2), 2.0, 3.0)
    alcohol_ok = data["drinks_per_day"].notna() & data["drinks_per_day"].le(limit)
    data = data.loc[alcohol_ok].copy()
    record("After excluding significant alcohol use", data)
    flow = pd.DataFrame(counts, columns=["stage", "all_releases", *CYCLES])
    return data, flow


def derive_fields(data):
    data = data.rename(columns={
        "SEQN": "seqn", "RIDAGEYR": "age_years",
        "LUXCAPM": "cap_db_m", "LUXSMED": "liver_stiffness_kpa",
        "LBXSATSI": "alt_u_l", "LBXSASSI": "ast_u_l", "LBXSGTSI": "ggt_u_l",
        "LBXSUA": "uric_acid_mg_dl", "LBXTR": "triglycerides_mg_dl",
        "LBDHDD": "hdl_cholesterol_mg_dl", "LBXTC": "total_cholesterol_mg_dl",
        "LBXGLU": "fasting_glucose_mg_dl", "LBXGH": "hba1c_pct",
        "BMXBMI": "bmi_kg_m2", "BMXWAIST": "waist_cm",
        "INDFMPIR": "poverty_income_ratio",
        "SDMVSTRA": "sdmvstra", "SDMVPSU": "sdmvpsu",
    })
    data["sex"] = data["RIAGENDR"].map({1: "Male", 2: "Female"})
    data["race_ethnicity"] = data["RIDRETH3"].map({
        1: "Mexican American", 2: "Other Hispanic", 3: "Non-Hispanic White",
        4: "Non-Hispanic Black", 6: "Non-Hispanic Asian", 7: "Other or multiracial",
    })
    data["education"] = data["DMDEDUC2"].map({
        1: "Less than 9th grade", 2: "9th to 11th grade",
        3: "High school graduate or GED", 4: "Some college or associate degree",
        5: "College graduate or above",
    })
    smoking = data["SMQ040"].map({1: "Current", 2: "Current", 3: "Former"})
    data["smoking_status"] = smoking.where(data["SMQ020"].eq(1)).mask(data["SMQ020"].eq(2), "Never")
    data["age_top_coded_80_plus"] = data["age_years"].eq(80).astype(int)
    data["mvpa_min_per_week"] = mvpa_minutes(data)

    female = data["sex"].eq("Female")
    asian = data["race_ethnicity"].eq("Non-Hispanic Asian")
    bmi_limit = np.where(asian, 23.0, 25.0)
    waist_limit = np.where(female, 80.0, 94.0)
    crit1_met = data["bmi_kg_m2"].ge(bmi_limit) | data["waist_cm"].gt(waist_limit)
    crit1_known = data["bmi_kg_m2"].notna() | data["waist_cm"].notna()

    diabetes = (
        data["DIQ010"].eq(1) | data["hba1c_pct"].ge(6.5)
        | data["DIQ050"].eq(1) | data["DIQ070"].eq(1)
    )
    crit2_met = data["fasting_glucose_mg_dl"].ge(100) | data["hba1c_pct"].ge(5.7) | diabetes
    crit2_known = (
        data["fasting_glucose_mg_dl"].notna() | data["hba1c_pct"].notna()
        | data["DIQ010"].notna() | data["DIQ050"].notna() | data["DIQ070"].notna()
    )

    bp_med_raw = pd.Series(
        np.where(data["cycle_tag"].eq("P"), data["BPQ050A"], data["BPQ150"]),
        index=data.index,
    )
    chol_med_raw = pd.Series(
        np.where(data["cycle_tag"].eq("P"), data["BPQ100D"], data["BPQ101D"]),
        index=data.index,
    )
    # The medication questions are skipped for people never told they have the
    # condition (BPQ020 = 2 for blood pressure, BPQ080 = 2 for cholesterol):
    # for them "not taking the medication" is known, not missing. An actual
    # answer is never overridden; only a skipped (missing) answer is filled.
    on_htn = bp_med_raw.map(YES_NO)
    on_htn = on_htn.mask(on_htn.isna() & data["BPQ020"].eq(2), 0.0)
    on_chol = chol_med_raw.map(YES_NO)
    on_chol = on_chol.mask(on_chol.isna() & data["BPQ080"].eq(2), 0.0)
    data["on_antihypertensive"] = on_htn.astype("Int64")
    data["on_lipid_lowering"] = on_chol.astype("Int64")
    crit3_met = (
        data["mean_systolic_bp_mm_hg"].ge(130) | data["mean_diastolic_bp_mm_hg"].ge(85)
        | data["on_antihypertensive"].eq(1)
    )
    crit3_known = (
        data["mean_systolic_bp_mm_hg"].notna() | data["mean_diastolic_bp_mm_hg"].notna()
        | data["on_antihypertensive"].notna()
    )
    crit4_met = data["triglycerides_mg_dl"].ge(150) | data["on_lipid_lowering"].eq(1)
    crit4_known = data["triglycerides_mg_dl"].notna() | data["on_lipid_lowering"].notna()
    hdl_limit = np.where(female, 50.0, 40.0)
    crit5_met = data["hdl_cholesterol_mg_dl"].le(hdl_limit) | data["on_lipid_lowering"].eq(1)
    crit5_known = data["hdl_cholesterol_mg_dl"].notna() | data["on_lipid_lowering"].notna()

    mets = [crit1_met, crit2_met, crit3_met, crit4_met, crit5_met]
    knowns = [crit1_known, crit2_known, crit3_known, crit4_known, crit5_known]
    for name, met, known in zip(CRITERIA, mets, knowns):
        data[name] = met.astype("Int64").mask(~known)
    data["n_metabolic_criteria_met"] = pd.DataFrame(
        {n: m.astype("Int64") for n, m in zip(CRITERIA, mets)}
    ).sum(axis=1, min_count=1)
    data["n_metabolic_criteria_known"] = pd.DataFrame(
        {n: k.astype(int) for n, k in zip(CRITERIA, knowns)}
    ).sum(axis=1)

    # Only participants whose metabolic status is undetermined AND who have no
    # criterion met are excluded: anyone with >= 1 criterion met is classifiable
    # as MASLD-eligible regardless of missing components (e.g. triglycerides,
    # measured only in the fasting subsample).
    met0 = data["n_metabolic_criteria_met"].fillna(0).eq(0)
    undetermined = met0 & data["n_metabolic_criteria_known"].lt(5)
    n_excluded_metabolic = int(undetermined.sum())
    data = data.loc[~undetermined].copy()

    data["steatosis"] = data["cap_db_m"].ge(248).astype(int)
    data["masld"] = (
        data["steatosis"].eq(1) & data["n_metabolic_criteria_met"].ge(1)
    ).astype(int)
    data["masld_with_significant_fibrosis"] = (
        data["masld"].eq(1) & data["liver_stiffness_kpa"].ge(8.0)
    ).astype(int)
    y = (
        0.953 * np.log(data["triglycerides_mg_dl"]) + 0.139 * data["bmi_kg_m2"]
        + 0.718 * np.log(data["ggt_u_l"]) + 0.053 * data["waist_cm"] - 15.745
    )
    data["fatty_liver_index"] = 100 * np.exp(y) / (1 + np.exp(y))
    fli = data["fatty_liver_index"]
    data["fli_ge_60"] = fli.ge(60).astype("Int64").mask(fli.isna())
    return data, n_excluded_metabolic


def weighted_prevalence(subset, outcome):
    w = subset["mec_weight_pooled"]
    return float((w * subset[outcome]).sum() / w.sum())


def write_outputs(data, flow, sources, artifact_log):
    design = ["sdmvstra", "sdmvpsu", "mec_weight_original", "mec_weight_pooled"]
    fields = [
        "seqn", "survey_release", *design, *OUTCOMES, *CRITERIA,
        "n_metabolic_criteria_met", "n_metabolic_criteria_known",
        *CORE, *OPTIONAL,
        "age_top_coded_80_plus", "bp_protocol", "on_antihypertensive",
        "on_lipid_lowering", "drinks_per_day", "fasting_glucose_mg_dl",
    ]
    cohort = data[fields].sort_values("seqn").copy()
    for name in ["seqn", "sdmvstra", "sdmvpsu", "age_years"]:
        cohort[name] = cohort[name].astype("Int64")
    assert cohort["seqn"].is_unique
    assert cohort[[o for o in OUTCOMES if o != "fli_ge_60"] + design].notna().all().all()
    assert cohort["masld"].eq(1).sum() == (
        (cohort["steatosis"].eq(1)) & (cohort["n_metabolic_criteria_met"].ge(1))
    ).sum()
    strata_collision = (
        cohort.groupby("sdmvstra")["survey_release"].nunique().gt(1).sum()
    )

    cohort.to_csv(ROOT / "nhanes_masld_cross_sectional.csv", index=False, float_format="%.10g")
    flow.to_csv(ROOT / "cohort_flow.csv", index=False)
    pd.DataFrame(
        sources, columns=["survey_release", "component", "cdc_xpt_url", "cdc_codebook_url", "source_rows"]
    ).to_csv(ROOT / "source_manifest.csv", index=False)
    if artifact_log:
        pd.DataFrame(
            artifact_log, columns=["survey_release", "component", "variable", "n_zero_restored"]
        ).to_csv(ROOT / "sas_zero_log.csv", index=False)

    summaries = []
    missingness = []
    for label, subset in [("All releases", cohort), *cohort.groupby("survey_release", sort=False)]:
        has_core = subset[CORE].notna().all(axis=1)
        summaries.append({
            "survey_release": label, "n": len(subset),
            "masld": int(subset["masld"].sum()),
            "masld_unweighted_pct": round(100 * subset["masld"].mean(), 2),
            "masld_weighted_pct": round(100 * weighted_prevalence(subset, "masld"), 2),
            "steatosis": int(subset["steatosis"].sum()),
            "steatosis_weighted_pct": round(100 * weighted_prevalence(subset, "steatosis"), 2),
            "masld_with_significant_fibrosis": int(subset["masld_with_significant_fibrosis"].sum()),
            "masld_with_significant_fibrosis_weighted_pct": round(
                100 * weighted_prevalence(subset, "masld_with_significant_fibrosis"), 2),
            "n_complete_core": int(has_core.sum()),
            "masld_complete_core": int(subset.loc[has_core, "masld"].sum()),
        })
        for name in CORE + OPTIONAL + OUTCOMES + CRITERIA:
            available = subset[name].notna().sum()
            missingness.append([label, name, len(subset), int(available),
                                int(len(subset) - available),
                                round(100 * (1 - available / len(subset)), 2)])
    pd.DataFrame(summaries).to_csv(ROOT / "outcome_summary.csv", index=False)
    pd.DataFrame(
        missingness, columns=["survey_release", "variable", "n", "n_available", "n_missing", "pct_missing"]
    ).to_csv(ROOT / "variable_completeness.csv", index=False)

    dictionary = [
        ["seqn", "Respondent sequence number", "All files"],
        ["survey_release", "NHANES release (development 2017-March2020; validation 2021-2023)", "DEMO"],
        ["sdmvstra", "Masked variance pseudo-stratum", "DEMO"],
        ["sdmvpsu", "Masked variance pseudo-PSU", "DEMO"],
        ["mec_weight_original", "MEC examination weight as released (WTMECPRP or WTMEC2YR)", "DEMO"],
        ["mec_weight_pooled", "MEC weight rescaled by release years / 5.2", "Derived"],
        ["cap_db_m", "Controlled attenuation parameter, median (dB/m)", "LUX"],
        ["liver_stiffness_kpa", "Liver stiffness, median (kPa)", "LUX"],
        ["steatosis", "CAP >= 248 dB/m", "Derived"],
        ["masld", "MASLD: steatosis + >=1 cardiometabolic criterion, alcohol within limits, no viral hepatitis", "Derived"],
        ["masld_with_significant_fibrosis", "MASLD and liver stiffness >= 8.0 kPa", "Derived"],
        ["fli_ge_60", "Fatty Liver Index >= 60 (sensitivity outcome; fasting subsample only)", "Derived"],
        ["criterion_adiposity", "BMI >=25 (>=23 Asian) or waist >94 cm men / >80 cm women", "Derived"],
        ["criterion_hyperglycaemia", "Fasting glucose >=100 mg/dL or HbA1c >=5.7% or diabetes", "Derived"],
        ["criterion_hypertension", "BP >=130/85 mmHg or antihypertensive treatment", "Derived"],
        ["criterion_hypertriglyceridaemia", "Triglycerides >=150 mg/dL or lipid-lowering treatment", "Derived"],
        ["criterion_low_hdl", "HDL <=40 mg/dL men / <=50 mg/dL women or lipid-lowering treatment", "Derived"],
        ["n_metabolic_criteria_met", "Number of the five criteria met (0-5)", "Derived"],
        ["n_metabolic_criteria_known", "Number of the five criteria determinable (0-5)", "Derived"],
        ["age_years", "Age at examination (top-coded at 80)", "DEMO"],
        ["sex", "Male / Female", "DEMO"],
        ["race_ethnicity", "Race/Hispanic origin (RIDRETH3)", "DEMO"],
        ["smoking_status", "Never / Former / Current (SMQ020, SMQ040)", "SMQ"],
        ["education", "Highest education (DMDEDUC2)", "DEMO"],
        ["mvpa_min_per_week", "Leisure-time moderate-to-vigorous activity, min/week (GPAQ recreational in 2017-March2020; LTPA in 2021-2023)", "PAQ"],
        ["total_cholesterol_mg_dl", "Total cholesterol (LBXTC)", "TCHOL"],
        ["alt_u_l", "Alanine aminotransferase (LBXSATSI)", "BIOPRO"],
        ["ast_u_l", "Aspartate aminotransferase (LBXSASSI)", "BIOPRO"],
        ["ggt_u_l", "Gamma-glutamyl transferase (LBXSGTSI)", "BIOPRO"],
        ["waist_cm", "Waist circumference (BMXWAIST)", "BMX"],
        ["bmi_kg_m2", "Body mass index (BMXBMI)", "BMX"],
        ["hba1c_pct", "Glycated haemoglobin (LBXGH)", "GHB"],
        ["triglycerides_mg_dl", "Triglycerides (LBXTR; fasting subsample only)", "TRIGLY"],
        ["hdl_cholesterol_mg_dl", "HDL cholesterol (LBDHDD)", "HDL"],
        ["mean_systolic_bp_mm_hg", "Mean systolic BP over readings > 0 mmHg (oscillometric)", "BPXO"],
        ["mean_diastolic_bp_mm_hg", "Mean diastolic BP over readings > 0 mmHg (oscillometric)", "BPXO"],
        ["poverty_income_ratio", "Family poverty income ratio (INDFMPIR)", "DEMO"],
        ["uric_acid_mg_dl", "Uric acid (LBXSUA)", "BIOPRO"],
        ["fatty_liver_index", "Fatty Liver Index, 0-100 (Bedogni et al. 2006)", "Derived"],
        ["age_top_coded_80_plus", "Age top-coded at 80 (1) or not (0)", "Derived"],
        ["bp_protocol", "Blood pressure measurement protocol", "Derived"],
        ["on_antihypertensive", "Taking prescribed antihypertensive medication (BPQ050A / BPQ150; 0 if never told high BP, BPQ020 = 2)", "BPQ"],
        ["on_lipid_lowering", "Taking prescribed lipid-lowering medication (BPQ100D / BPQ101D; 0 if not told to take it, BPQ080 = 2)", "BPQ"],
        ["drinks_per_day", "Average alcoholic drinks per day, past 12 months (ALQ130 x ALQ121 frequency / 365; 0 if never drank)", "ALQ"],
        ["fasting_glucose_mg_dl", "Fasting plasma glucose (LBXGLU; fasting subsample only)", "GLU"],
    ]
    pd.DataFrame(dictionary, columns=["variable", "description", "source"]).to_csv(
        ROOT / "data_dictionary.csv", index=False
    )

    readme = f"""# NHANES MASLD cross-sectional cohort

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

`{strata_collision}` masked variance strata contain participants from both
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

{flow.to_string(index=False)}

## Outcome summary

{pd.DataFrame(summaries).to_string(index=False)}
"""
    (ROOT / "README.md").write_text(readme, encoding="utf-8")
    print(flow.to_string(index=False))
    print(pd.DataFrame(summaries).to_string(index=False))
    print(f"Strata shared across releases: {strata_collision}")
    if artifact_log:
        print("SAS zeros restored (2^-260 -> 0.0):")
        print(pd.DataFrame(artifact_log).to_string(index=False))


if __name__ == "__main__":
    frames = []
    manifest = []
    artifact_log = []
    for cycle, args in CYCLES.items():
        frame, sources = read_cycle(cycle, *args, artifact_log)
        frames.append(frame)
        manifest.extend(sources)
    data = pd.concat(frames, ignore_index=True)
    assert data["SEQN"].is_unique
    cohort, flow = select_cohort(data)
    cohort, n_excluded_metabolic = derive_fields(cohort)
    flow = pd.concat([flow, pd.DataFrame([[
        "After excluding undetermined metabolic status",
        len(cohort),
        *(cohort["survey_release"].eq(k).sum() for k in CYCLES),
    ]], columns=flow.columns)], ignore_index=True)
    print(f"Excluded for undetermined metabolic status: {n_excluded_metabolic}")
    write_outputs(cohort, flow, manifest, artifact_log)
