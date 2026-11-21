#!/usr/bin/env python3
"""Build the MASLD manuscript docx (schema-safe python-docx).

Reads the verified analysis outputs in this directory:
  table1.csv, model_results.json, metrics_table.csv, threshold_metrics.csv,
  odds_ratios.csv, permutation_importance.csv, validation_predictions.csv
and the bioinformatics results in bioinformatics/ (degs.csv, enrichment.csv).
Writes MASLD_NHANES_Manuscript.docx in this directory.
"""
import csv
import json
from pathlib import Path

import docx
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt, Inches, RGBColor

HERE = Path(__file__).resolve().parent
OUT = HERE / "MASLD_NHANES_Manuscript.docx"

TITLE = ("Identifying Metabolic Dysfunction-Associated Steatotic Liver Disease "
         "in US Adults: A Cross-Sectional NHANES Study with Temporal "
         "Validation and Transcriptomic Context")
RUNNING = "Identifying MASLD in US adults: NHANES prediction models"


def load_json(name):
    with open(HERE / name) as f:
        return json.load(f)


RESULTS = load_json("model_results.json")


def shade_cell(cell, fill):
    """Schema-safe cell shading: shd appended to an otherwise empty tcPr
    (or after tcW/gridSpan, both of which precede shd in the tcPr sequence)."""
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), fill)
    tcPr.append(shd)


def set_cell_text(cell, text, bold=False, size=9, color=None, align=None):
    cell.text = ""
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(1)
    p.paragraph_format.space_after = Pt(1)
    if align is not None:
        p.alignment = align
    r = p.add_run(str(text))
    r.font.name = "Times New Roman"
    r.font.size = Pt(size)
    r.bold = bold
    if color is not None:
        r.font.color.rgb = RGBColor.from_string(color)


def add_table(doc, title, headers, rows, widths=None, size=9):
    doc.add_paragraph(title, style="CaptionTitle")
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    if widths:
        for i, w in enumerate(widths):
            table.columns[i].width = Inches(w)
    hdr = table.rows[0].cells
    for i, h in enumerate(headers):
        shade_cell(hdr[i], "1F4E79")
        set_cell_text(hdr[i], h, bold=True, size=size, color="FFFFFF")
    for row in rows:
        cells = table.add_row().cells
        for i, val in enumerate(row):
            bold = str(val).startswith("**")
            set_cell_text(cells[i], str(val).replace("**", ""), bold=bold,
                          size=size)
    doc.add_paragraph("", style="Normal")
    return table


def para(doc, text, style="Normal", bold_prefix=None, italic=False):
    p = doc.add_paragraph(style=style)
    if bold_prefix:
        r = p.add_run(bold_prefix)
        r.bold = True
    r = p.add_run(text)
    if italic:
        r.italic = True
    return p


def build():
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    st.paragraph_format.space_after = Pt(6)
    st.paragraph_format.line_spacing = 2.0
    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1)
        section.right_margin = Inches(1)
    for name, size, before in [("Heading 1", 14, 14), ("Heading 2", 12, 10)]:
        s = doc.styles[name]
        s.font.name = "Times New Roman"
        s.font.size = Pt(size)
        s.font.bold = True
        s.paragraph_format.space_before = Pt(before)
        s.paragraph_format.space_after = Pt(4)
        s.paragraph_format.line_spacing = 1.0
    cap = doc.styles.add_style("CaptionTitle",
                                 docx.enum.style.WD_STYLE_TYPE.PARAGRAPH)
    cap.font.name = "Times New Roman"
    cap.font.size = Pt(10)
    cap.font.bold = True
    cap.paragraph_format.space_before = Pt(10)
    cap.paragraph_format.space_after = Pt(2)
    cap.paragraph_format.line_spacing = 1.0
    fn = doc.styles.add_style("Footnote",
                              docx.enum.style.WD_STYLE_TYPE.PARAGRAPH)
    fn.font.name = "Times New Roman"
    fn.font.size = Pt(9)
    fn.font.italic = True
    fn.paragraph_format.space_after = Pt(6)
    fn.paragraph_format.line_spacing = 1.0

    # ---------------------------------------------------------- title page
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(TITLE)
    r.bold = True
    r.font.size = Pt(16)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Ibrahim H. Hemdan")
    r.font.size = Pt(12)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("[Affiliation and institution — to be completed by the author]")
    r.italic = True
    r.font.size = Pt(11)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Corresponding author: Ibrahim H. Hemdan — "
                  "[email address to be added]")
    r.font.size = Pt(11)
    for line in [
        "Running title: " + RUNNING,
        "Article type: Original research",
        "Keywords: metabolic dysfunction-associated steatotic liver disease; "
        "MASLD; NHANES; prediction model; gradient boosting; calibration; "
        "transient elastography; transcriptomics",
        "Word count: abstract 243 words; main text 4,492 words; "
        "6 tables; 8 figures.",
        "Target journal: to be confirmed by the author (formatted generically "
        "for a general medical journal; reference style Vancouver).",
    ]:
        para(doc, line, style="Normal")
    doc.add_page_break()

    # ---------------------------------------------------------- abstract
    doc.add_heading("Abstract", level=1)
    abstract = [
        ("Background: ", "Metabolic dysfunction-associated steatotic liver "
         "disease (MASLD) affects about one in three US adults, yet most "
         "cases are undetected because diagnosis requires imaging or biopsy. "
         "Routine data could focus further testing on those most likely to "
         "have the disease."),
        ("Methods: ", "We analysed the National Health and Nutrition "
         "Examination Survey (NHANES) 2017–March 2020 (development, n = "
         "6,320) and August 2021–August 2023 (temporal validation, n = "
         "4,103). MASLD was defined as a controlled attenuation parameter "
         "≥248 dB/m with at least one cardiometabolic risk factor, alcohol "
         "intake within MASLD limits, and no viral hepatitis. Logistic "
         "regression and histogram gradient boosting were compared using the "
         "same ten predictors, excluding every outcome-defining variable. "
         "Analyses were survey-weighted, with bootstrap confidence "
         "intervals. Public liver RNA-seq data (GSE135251; 216 biopsies) "
         "provided transcriptomic context."),
        ("Results: ", "Weighted MASLD prevalence was 56.8% and 54.8% in "
         "development and validation. On validation, gradient boosting achieved "
         "an AUROC of 0.721 (95% CI 0.713–0.729) versus 0.698 (0.690–0.706) "
         "for logistic regression (paired difference +0.023, 0.020–0.026); "
         "Brier scores were 0.211 versus 0.218 and calibration slopes 0.89 "
         "versus 0.84. At Youden's threshold, sensitivity was ~72% and "
         "specificity 57–61%. Alanine aminotransferase was the strongest "
         "predictor (odds ratio 1.90 per 10 U/L). In liver tissue, 606 genes "
         "were differentially expressed in steatohepatitis, including "
         "up-regulation of FASN, GGT1, PNPLA3 and COL1A1."),
        ("Conclusions: ", "Routine-variable models identify MASLD with "
         "modest discrimination; gradient boosting is modestly better than "
         "treating all patients from thresholds of about 0.25 upward. The "
         "models are a risk-stratification aid, not a diagnostic rule."),
    ]
    for bold_prefix, text in abstract:
        para(doc, text, bold_prefix=bold_prefix)

    # ---------------------------------------------------------- introduction
    doc.add_heading("Introduction", level=1)
    for text in [
        "Metabolic dysfunction-associated steatotic liver disease (MASLD) is "
        "the most common chronic liver disease worldwide. A 2023 systematic "
        "review estimated that about 38% of adults globally had non-alcoholic "
        "fatty liver disease in 2016–2019, and the burden continues to rise "
        "with obesity and type 2 diabetes [1, 2]. In the United States, "
        "recent National Health and Nutrition Examination Survey (NHANES) "
        "analyses using vibration-controlled transient elastography (VCTE) "
        "found an age-adjusted MASLD prevalence of 31.9% at a controlled "
        "attenuation parameter (CAP) threshold of 285 dB/m [3], and 25.6% at "
        "a threshold of 300 dB/m in the 2017–2023 cycles [4]. MASLD is the "
        "leading cause of liver-related morbidity in many countries, and it "
        "is strongly associated with cardiovascular disease, type 2 diabetes "
        "and chronic kidney disease [2, 5].",
        "MASLD is defined by hepatic steatosis together with at least one of "
        "five cardiometabolic risk factors, in a person drinking alcohol "
        "below the thresholds that define alcohol-related disease and "
        "without another cause of steatosis [1]. The definition was adopted "
        "by a multisociety Delphi consensus in 2023 and has been incorporated "
        "into the 2024 EASL–EASD–EASO clinical practice guidelines, which "
        "recommend case-finding in people with metabolic risk factors rather "
        "than population screening [1, 5]. In practice, however, most cases "
        "remain undetected: diagnosis requires imaging or biopsy, and liver "
        "biopsy is invasive while elastography is not universally available. "
        "A model that identifies people likely to have MASLD from data "
        "already collected in routine care could focus elastography and "
        "further assessment where they are most useful — a need that has "
        "grown more urgent since resmetirom became the first approved "
        "pharmacological treatment for metabolic dysfunction-associated "
        "steatohepatitis (MASH) with fibrosis [6].",
        "Several non-invasive scores exist, notably the Fatty Liver Index "
        "(FLI), which combines triglycerides, body mass index (BMI), waist "
        "circumference and gamma-glutamyl transferase (GGT) [7]. Machine "
        "learning models for MASLD have also been reported: an international "
        "validation study of gradient boosting reported areas under the "
        "receiver operating characteristic curve (AUROCs) above 0.80 across "
        "five cohorts [8], and a hospital-based study externally validated "
        "in NHANES reported AUROCs of 0.78–0.83 [9]. These models, however, "
        "include predictors that themselves define the outcome — BMI, waist "
        "circumference, triglycerides and glucose — so their discrimination "
        "is inflated by construction and they cannot answer the more "
        "informative question: how well can MASLD be identified from "
        "variables that do not define it? Few studies have evaluated "
        "prediction models with temporal validation in a nationally "
        "representative sample, reported calibration and decision-curve "
        "analysis alongside discrimination, and linked the predictors to "
        "liver biology.",
        "We therefore developed and temporally validated two models for "
        "identifying MASLD in US adults — L2-penalised logistic regression "
        "and histogram gradient boosting — using NHANES 2017–March 2020 for "
        "development and August 2021–August 2023 for validation, with the "
        "same ten predictors in both models and every outcome-defining "
        "variable excluded. We report discrimination, calibration, threshold "
        "performance and clinical utility, test the sensitivity of the "
        "results to the CAP threshold, the outcome definition and the "
        "predictor set, and place the clinical predictors in biological "
        "context by re-analysing public human liver transcriptomic data. "
        "Reporting follows the TRIPOD+AI statement [10].",
    ]:
        para(doc, text)

    # ---------------------------------------------------------- methods
    doc.add_heading("Methods", level=1)
    doc.add_heading("Study design and data source", level=2)
    for text in [
        "This is a cross-sectional analysis of NHANES, a programme of the US "
        "National Center for Health Statistics (NCHS) that examines a "
        "nationally representative sample of the civilian "
        "non-institutionalised population using a complex, multistage, "
        "stratified probability design [11]. Data come from two releases: "
        "NHANES 2017–March 2020 pre-pandemic (development sample) and "
        "NHANES August 2021–August 2023 (temporal validation sample). The "
        "2017–March 2020 files were released as a single combined "
        "pre-pandemic sample; the 2017–2018 and 2019–March 2020 participants "
        "cannot be separated within them, because the Centers for Disease "
        "Control and Prevention (CDC) reassigned participant identifiers and "
        "masked design variables in the combined release. Both releases "
        "include VCTE (FibroScan 502 V2 Touch, M/XL probe) with CAP and "
        "liver stiffness measurement (LSM) [12, 13]. NHANES is publicly "
        "available de-identified data collected with informed consent under "
        "NCHS Research Ethics Review Board approval; no additional ethics "
        "approval was required for this secondary analysis.",
        "Participants were adults aged 20 years or older with a positive "
        "mobile examination centre (MEC) weight and a valid CAP measurement "
        "from a reliable elastography examination (LUAXSTAT = 1, defined by "
        "the CDC as fasting for at least 3 hours, at least 10 complete "
        "stiffness measurements, and an interquartile range to median LSM "
        "ratio below 30%). Pregnant participants are ineligible for "
        "elastography and none remained after that step. We then excluded "
        "participants with viral hepatitis B (hepatitis B surface antigen "
        "positive) or C (hepatitis C virus RNA positive) or unknown viral "
        "hepatitis status, participants whose average alcohol intake exceeded "
        "the MASLD limits (see below) or whose alcohol data were missing, "
        "and participants whose metabolic status could not be determined "
        "(no cardiometabolic criterion met and at least one criterion "
        "undetermined).",
    ]:
        para(doc, text)
    doc.add_heading("Outcome definition", level=2)
    for text in [
        "The outcome was MASLD as defined by the 2023 multisociety consensus "
        "[1]: hepatic steatosis (CAP ≥248 dB/m, the threshold for steatosis "
        "grade S1 or higher in an individual-patient meta-analysis [14]) "
        "plus at least one of five cardiometabolic risk factors — (i) BMI "
        "≥25 kg/m² (≥23 kg/m² for Asian ancestry) or waist circumference "
        ">94 cm in men or >80 cm in women; (ii) fasting glucose ≥100 mg/dL "
        "or 2-hour glucose ≥140 mg/dL or glycated haemoglobin (HbA1c) ≥5.7% "
        "or type 2 diabetes or glucose-lowering treatment; (iii) blood "
        "pressure ≥130/85 mmHg or antihypertensive treatment; (iv) "
        "triglycerides ≥150 mg/dL or lipid-lowering treatment; (v) "
        "high-density lipoprotein (HDL) cholesterol ≤40 mg/dL in men or "
        "≤50 mg/dL in women or lipid-lowering treatment — with average "
        "alcohol intake of no more than 2 drinks per day in women and 3 in "
        "men, and no viral hepatitis B or C. Average alcohol intake was "
        "calculated from the number of drinks per drinking day multiplied by "
        "drinking frequency and divided by 365, because the questionnaire "
        "item records drinks per drinking day, not per calendar day; using "
        "it directly would have wrongly excluded about 40% of the sample. "
        "Secondary outcomes were steatosis alone (CAP ≥248 dB/m) and MASLD "
        "with significant fibrosis (LSM ≥8.0 kPa).",
        "Two sensitivity definitions were pre-specified: MASLD at a higher "
        "CAP threshold (≥285 dB/m, the threshold that optimises sensitivity "
        "and specificity and the one used in recent US prevalence estimates "
        "[3, 14]) and MASLD defined by a Fatty Liver Index ≥60 [7] in the "
        "morning fasting subsample. Each metabolic criterion was recorded as "
        "met, not met, or undetermined; a participant was excluded for "
        "undetermined metabolic status only if no criterion was met and at "
        "least one was undetermined, so that missing triglycerides in the "
        "non-fasting half of the sample did not discard participants who "
        "were already classifiable.",
    ]:
        para(doc, text)
    doc.add_heading("Predictors", level=2)
    for text in [
        "Both models used the same ten predictors: age, sex, race and "
        "Hispanic origin, smoking status, education, leisure-time "
        "moderate-to-vigorous physical activity (MVPA), total cholesterol, "
        "alanine aminotransferase (ALT), aspartate aminotransferase (AST) and "
        "GGT. The predictor set was chosen by a single rule: no variable "
        "that defines the outcome, or that is a direct component of it, may "
        "be a predictor. "
        "BMI, waist circumference, blood pressure, triglycerides, HDL "
        "cholesterol, fasting glucose and HbA1c are therefore excluded from "
        "the primary models; a sensitivity analysis adds them back to show "
        "the size of the resulting inflation. Continuous predictors were "
        "centred at the training weighted median and scaled (age per 10 "
        "years, MVPA per 60 minutes per week, total cholesterol per 20 "
        "mg/dL, ALT, AST and GGT per 10 U/L); categorical predictors were "
        "one-hot encoded with reference categories female sex, non-Hispanic "
        "White race, never smoking, and college graduate or above. Missing "
        "predictor values (0–2.5% for the ten core predictors) were imputed "
        "with the training weighted median (continuous) or weighted mode "
        "(categorical), estimated on training data only. Physical activity "
        "was harmonised as leisure-time MVPA only, from the Global Physical "
        "Activity Questionnaire recreational items in 2017–March 2020 and "
        "the leisure-time items in 2021–2023; the instrument changed between "
        "cycles, so MVPA values are not directly comparable across cycles, "
        "and this is treated as a limitation rather than a modelling "
        "correction.",
    ]:
        para(doc, text)
    doc.add_heading("Statistical analysis", level=2)
    for text in [
        "All analyses used the MEC examination weights, rescaled by years of "
        "data collection (3.2/5.2 for 2017–March 2020 and 2/5.2 for August "
        "2021–August 2023) so that the pooled weights represent the US "
        "population over the 5.2 years covered. Hyperparameters were tuned "
        "by five-fold stratified cross-validation within the development "
        "sample, optimising the survey-weighted Brier score: logistic "
        "regression over C in {0.01, 0.1, 1, 10}, and histogram gradient "
        "boosting over maximum leaves {7, 15, 31}, learning rate {0.05, 0.1}, "
        "iterations {100, 200} and L2 regularisation {0, 1}, with a minimum "
        "of 30 samples per leaf. Each model was then refitted on all "
        "development data and evaluated, without recalibration, on the "
        "validation sample. We report the AUROC, average precision, Brier "
        "score, calibration intercept and slope, and mean predicted versus "
        "observed probability, with 95% confidence intervals from 1,500 "
        "rescaled Rao–Wu bootstrap replicates of the validation sample "
        "(strata resampled at the primary sampling unit level; single-PSU "
        "strata retained with multiplier 1). Paired differences in AUROC and "
        "Brier score between the two models were bootstrapped the same way. "
        "The constant-prevalence benchmark Brier score (everyone assigned "
        "the development prevalence) is reported for reference.",
        "Threshold performance is reported at Youden's J and at 90% "
        "sensitivity, with survey-weighted sensitivity, specificity, "
        "positive and negative predictive values and bootstrap confidence "
        "intervals at the fixed point-estimate thresholds. Clinical utility "
        "was assessed by decision curve analysis [15, 16], computing the "
        "survey-weighted net benefit of each model across threshold "
        "probabilities from 0.05 to 0.95 against the default strategies of "
        "treating all patients and treating none. Logistic-regression odds "
        "ratios with 95% confidence intervals come from 1,500 nonparametric "
        "bootstrap refits on the development sample. Gradient boosting "
        "variable importance is permutation importance on the validation "
        "sample (mean increase in weighted Brier score over 10 repeats per "
        "predictor). Sensitivity analyses refit both models with the same "
        "tuned hyperparameters for each secondary outcome and definition, "
        "for a criterion-containing predictor set (adding BMI per 5 kg/m², "
        "waist per 10 cm, HbA1c per 1%, triglycerides per 50 mg/dL, HDL per "
        "10 mg/dL and systolic blood pressure per 10 mmHg), for a "
        "routine-clinical predictor set (adding BMI per 5 kg/m² and systolic "
        "blood pressure per 10 mmHg only — both measured without a blood "
        "sample), for a complete-case analysis, and for an unweighted "
        "analysis.",
        "Analyses were performed in Python 3.12.3 (pandas 3.1.1, NumPy "
        "2.5.5, SciPy 1.18.1, scikit-learn 1.9.1). The cohort build, the "
        "modelling script and all outputs are deposited in a public code "
        "repository (see Declarations); NHANES source files are publicly "
        "available from the CDC.",
    ]:
        para(doc, text)
    doc.add_heading("Transcriptomic re-analysis", level=2)
    for text in [
        "To place the clinical predictors in biological context, we "
        "re-analysed publicly available human liver RNA sequencing data "
        "(Gene Expression Omnibus accession GSE135251 [17]), comprising 216 "
        "snap-frozen liver biopsies: 155 with non-alcoholic steatohepatitis "
        "(NASH), 51 with simple steatosis (NAFL) and 10 histologically "
        "normal controls. The dataset predates the 2023 MASLD rename; its "
        "NAFL and NASH groups are the same disease spectrum as MASLD and "
        "MASH under the current nomenclature. Gene-level HTSeq count files "
        "were downloaded from GEO. Genes were kept if counts per million "
        "exceeded 1 in at least 10 samples (15,619 Ensembl identifiers), "
        "mapped to HUGO Gene Nomenclature Committee symbols (14,971 genes "
        "tested), and compared between NASH and controls by Welch's t-test "
        "on log2 counts per million, with Benjamini–Hochberg correction. "
        "Differentially expressed genes were defined at a false discovery "
        "rate below 0.05 and an absolute log2 fold-change above 1. Pathway "
        "enrichment used the hypergeometric test against the Kyoto "
        "Encyclopedia of Genes and Genomes (KEGG) human pathways (373 "
        "pathways retrieved from the KEGG REST API), with Benjamini–Hochberg "
        "correction, run separately for all, up-regulated and down-regulated "
        "differentially expressed genes. Ten MASLD-relevant genes were "
        "examined individually: PNPLA3, TM6SF2, GGT1, GPT, GOT1, FASN, "
        "SREBF1, PPARA, COL1A1 and ACTA2. This analysis provides biological "
        "context only; the biopsy-selected cohort cannot validate the "
        "NHANES prediction model.",
    ]:
        para(doc, text)

    # ---------------------------------------------------------- results
    doc.add_heading("Results", level=1)
    doc.add_heading("Study population", level=2)
    for text in [
        "Of 27,493 interviewed participants, 14,608 were aged 20 years or "
        "older with a positive MEC weight, 13,536 had a valid CAP "
        "measurement, and 12,671 had a reliable elastography examination. "
        "After excluding viral hepatitis (907 participants), significant "
        "alcohol use (1,317) and undetermined metabolic status (24), 10,423 "
        "participants remained for analysis: 6,320 in development (3,709 "
        "with MASLD) and 4,103 in validation (2,265 with MASLD) (Figure 1). "
        "The weighted prevalence of MASLD (CAP ≥248 dB/m) was 56.8% in "
        "development and 54.8% in validation; weighted steatosis prevalence "
        "was 57.4% and 55.8%, and MASLD with significant fibrosis (LSM ≥8.0 "
        "kPa) was 7.9% and 9.1%. At the higher CAP threshold of 285 dB/m, "
        "weighted MASLD prevalence was 36.6% and 34.2%, close to published "
        "US estimates at that threshold (31.9% in Kim et al. [3] and "
        "age-standardised 26.8% falling to 23.6% at a 300 dB/m threshold in "
        "Unalp-Arida et al. [4]).",
        "Participants with MASLD were older and had substantially higher BMI "
        "and waist circumference, higher blood pressure, total cholesterol "
        "and HbA1c, lower HDL cholesterol, and higher ALT, AST and GGT than "
        "participants without MASLD, in both samples (Table 1). Median ALT "
        "was 21 U/L in MASLD versus 16 U/L without MASLD; median GGT was 23 "
        "versus 16 U/L; median CAP was 301 versus 213 dB/m. Participants "
        "with MASLD met a mean of 3.0 of the five metabolic criteria versus "
        "1.7 in those without MASLD. The adiposity criterion was met by "
        "96.6% of participants with MASLD and 65.3% of those without it; "
        "the hyperglycaemia criterion by 64.3% versus 34.7%. Median "
        "leisure-time MVPA was lower in MASLD in both samples, but the "
        "absolute values differed between cycles (30 versus 148 minutes per "
        "week among participants with MASLD), reflecting the change of "
        "physical activity instrument between NHANES cycles.",
    ]:
        para(doc, text)
    doc.add_heading("Model tuning and primary validation", level=2)
    for text in [
        "Five-fold cross-validation within development gave a weighted Brier "
        "score of 0.2090 for logistic regression at every value of C tested "
        "(0.01 to 10): at this sample size the penalty is negligible "
        "relative to the loss, so C = 0.01 was taken by tie-break. The best "
        "gradient boosting configuration was 7 leaves, learning rate 0.05, "
        "100 iterations and L2 = 1.0 (cross-validated Brier 0.2002); larger "
        "and deeper configurations overfitted (Brier up to 0.225).",
        "On temporal validation (Table 2, Figures 2 and 3), gradient boosting "
        "achieved an AUROC of 0.721 (95% CI 0.713–0.729) versus 0.698 "
        "(0.690–0.706) for logistic regression, and an average precision of "
        "0.743 (0.735–0.751) versus 0.729 (0.722–0.737). The paired "
        "difference in AUROC was +0.023 (95% CI 0.020–0.026) and in Brier "
        "score −0.0072 (95% CI −0.0083 to −0.0062); both intervals exclude "
        "zero, so gradient boosting is genuinely, if modestly, better. Both "
        "models beat the constant-prevalence benchmark Brier score of 0.248 "
        "(0.218 and 0.211). Calibration was reasonable, with slopes of 0.838 "
        "(0.792–0.884) for logistic regression and 0.891 (0.854–0.928) for "
        "gradient boosting — mild over-dispersion — and intercepts of 0.027 "
        "(0.003–0.051) and 0.044 (0.018–0.070), indicating good "
        "calibration-in-the-large. Mean predicted probabilities (54.3% and "
        "53.9%) were close to the observed prevalence (54.8%).",
    ]:
        para(doc, text)
    doc.add_heading("Threshold performance and clinical utility", level=2)
    for text in [
        "At Youden's threshold (Table 3), logistic regression (threshold "
        "0.502) had a sensitivity of 72.3% (95% CI 71.0–73.6) and a "
        "specificity of 56.8% (55.3–58.3), with a positive predictive value "
        "of 67.0% and a negative predictive value of 62.8%. Gradient boosting "
        "(threshold 0.511) had a sensitivity of 71.9% (70.9–73.0) and a "
        "specificity of 60.6% (59.3–61.9), with predictive values of 68.9% "
        "and 64.0%. At a threshold chosen for 90% sensitivity, specificity "
        "fell to 31.3% for logistic regression and 36.7% for gradient "
        "boosting — a screening-style threshold would flag most of the adult "
        "population. Decision curve analysis (Figure 6) showed that both "
        "models offer net benefit over treating all patients from threshold "
        "probabilities of about 0.25 upward: at a threshold of 0.50, net "
        "benefit was 0.213 for gradient boosting and 0.201 for logistic "
        "regression versus 0.097 for treating all; at 0.70 it was 0.072 and "
        "0.058 versus −0.506. Both models remained above treating none up "
        "to a threshold of about 0.90.",
    ]:
        para(doc, text)
    doc.add_heading("Predictor effects and variable importance", level=2)
    for text in [
        "In the logistic model (Table 4, Figure 4), the strongest predictor "
        "was ALT (odds ratio 1.90 per 10 U/L, 95% CI 1.58–2.34), followed by "
        "Mexican American ethnicity versus non-Hispanic White (1.90, "
        "1.48–2.45) and age (1.32 per 10 years, 1.26–1.40). GGT (1.07 per 10 "
        "U/L, 1.01–1.14) and MVPA (0.96 per 60 minutes per week, 0.94–0.98) "
        "were also significant. Total cholesterol (1.03, 0.99–1.08), male sex "
        "(1.18, 0.97–1.44), smoking and education were mostly not "
        "significant. AST had an odds ratio below 1 (0.60 per 10 U/L, "
        "0.47–0.71). This is a collinearity artefact, not a protective "
        "effect: ALT and AST correlate at r = 0.76 in development, and the "
        "two coefficients must be read jointly — ALT strongly positive, AST "
        "negative once ALT is held fixed. Permutation importance confirms "
        "that AST still contributes to the gradient boosting model (Figure "
        "5). The AST odds ratio must not be reported on its own as a "
        "protective association.",
        "Permutation importance for gradient boosting was dominated by the "
        "liver enzymes: ALT (mean Brier increase 0.026), GGT (0.024), age "
        "(0.011) and AST (0.009), followed by MVPA (0.003) and race and "
        "Hispanic origin (0.001); smoking, total cholesterol, sex and "
        "education each contributed less than 0.0005.",
    ]:
        para(doc, text)
    doc.add_heading("Secondary outcomes and sensitivity analyses", level=2)
    for text in [
        "Results were similar for steatosis alone (AUROC 0.696 for logistic "
        "regression and 0.719 for gradient boosting) and slightly better at "
        "the higher CAP threshold of 285 dB/m (0.702 and 0.737, with a "
        "gradient boosting calibration slope of 0.96) (Table 5). MASLD with "
        "significant fibrosis was the weakest analysis (AUROC 0.675 and "
        "0.697; calibration slopes 0.65 and 0.68), on only 390 validation "
        "events. MASLD defined by a Fatty Liver Index ≥60 gave higher "
        "AUROCs (0.792 and 0.808) in the fasting subsample, but this is "
        "partly definitional: the index is built from triglycerides, BMI, "
        "waist circumference and GGT, which correlate with the predictors. "
        "Complete-case analysis (0.698 and 0.723) and unweighted analysis "
        "(0.693 and 0.725) gave essentially identical AUROCs, so imputation "
        "and weighting are not driving the results.",
        "Adding the outcome-defining variables back into the predictor set — "
        "BMI, waist circumference, HbA1c, triglycerides, HDL cholesterol "
        "and systolic blood pressure — raised the AUROC to 0.853 for both "
        "models and lowered the Brier score to 0.154. This is the leakage "
        "demonstration, not a better model: it quantifies how much "
        "discrimination is bought by predicting an outcome from its own "
        "definition, and it is the reason those variables were excluded "
        "from the primary models.",
    ]:
        para(doc, text)
    doc.add_heading("Transcriptomic context", level=2)
    for text in [
        "In the liver RNA-seq re-analysis, 606 of 14,971 tested genes were "
        "differentially expressed between NASH (n = 155) and controls (n = "
        "10) at a false discovery rate below 0.05 and an absolute log2 "
        "fold-change above 1 (415 up-regulated, 191 down-regulated) "
        "(Figure 7). The largest fold-changes were IGFBP1 (−5.14), FOS "
        "(−4.77) and EGR1 (−3.85) on the down-regulated side, and AKR1B10 "
        "(+3.52, adjusted p = 4.5 × 10⁻⁹) and FADS2 (+3.08) on the "
        "up-regulated side; AKR1B10 was also the top biomarker gene reported "
        "by Govaere et al. in the source publication, so this independent "
        "re-analysis reproduces that finding [17].",
        "The lipogenesis axis was transcriptionally activated in NASH liver "
        "(Table 6): FASN (+2.56, adjusted p = 2.0 × 10⁻⁴), SREBF1 (+1.00, "
        "0.019) and PNPLA3 (+1.25, 6.0 × 10⁻⁴), as was the fibrosis marker "
        "COL1A1 (+1.10, 1.1 × 10⁻³). GGT1 (+1.45, 1.3 × 10⁻³) and GPT "
        "(+0.89, 0.011) — the genes encoding serum GGT and ALT, two "
        "predictors in the NHANES model — were up-regulated in NASH liver, "
        "whereas GOT1 (AST), TM6SF2, PPARA and ACTA2 were not differentially "
        "expressed. KEGG enrichment of the differentially expressed genes "
        "(Figure 8) highlighted carbon metabolism, bile secretion and "
        "extracellular-matrix and focal-adhesion programmes among "
        "up-regulated genes (adjusted p = 0.0021 for carbon metabolism, bile "
        "secretion, integrin signalling and extracellular-matrix–receptor "
        "interaction), and mineral absorption (0.00036), interleukin-17 "
        "signalling (0.0021), tumour necrosis factor signalling (0.026) and "
        "adipocytokine signalling (0.035) among down-regulated genes. The "
        "down-regulation of interleukin-17 and tumour necrosis factor "
        "signalling at bulk-tissue mRNA level is unexpected for an "
        "inflammatory disease and most likely reflects cell-composition "
        "differences and the small control group; it is reported, not "
        "interpreted further.",
    ]:
        para(doc, text)

    # ---------------------------------------------------------- discussion
    doc.add_heading("Discussion", level=1)
    for text in [
        "In a nationally representative sample of 10,423 US adults, two "
        "models built only from routine demographic and laboratory variables "
        "identified MASLD with modest discrimination on temporal validation "
        "(AUROC 0.698 for logistic regression and 0.721 for gradient "
        "boosting), with reasonable calibration and with net benefit over "
        "treating all patients from threshold probabilities of about 0.25 "
        "upward. Gradient boosting was genuinely, if modestly, better than "
        "logistic regression on both discrimination and overall accuracy. "
        "The liver enzymes dominated both models, and liver transcriptomic "
        "data support the biological plausibility of that signal: the genes "
        "encoding ALT and GGT are up-regulated in steatohepatitis, as are "
        "the lipogenesis and fibrosis programmes.",
        "Our AUROCs are lower than those reported by recent machine learning "
        "studies of MASLD — above 0.80 in the international validation study "
        "of Ye et al. [8] and 0.78–0.83 in the hospital-based study of Zhu "
        "et al. [9]. The difference is by design, not by failure: those "
        "models include BMI, waist circumference, triglycerides, glucose and "
        "other variables that define the outcome, and our criterion-containing "
        "sensitivity analysis shows exactly what that buys — an AUROC of "
        "0.853. A model that predicts MASLD from its own defining variables "
        "is, in effect, re-deriving the definition. Our stricter design "
        "answers the harder and more useful question of how much signal "
        "remains once the definition is removed, and the answer — an AUROC "
        "of about 0.72 with good calibration. Discrimination is also capped "
        "by the outcome itself: "
        "at a CAP threshold of 248 dB/m, MASLD affects more than half the "
        "adult population, and distinguishing cases from non-cases among "
        "people who are mostly metabolically unhealthy is intrinsically "
        "harder than detecting a rare disease.",
        "The modest advantage of gradient boosting over logistic regression "
        "is consistent with non-linear relationships and interactions among "
        "the predictors — for example, the joint behaviour of ALT and AST — "
        "that a linear model cannot capture. That the advantage is small "
        "(+0.023 AUROC) is itself informative: most of the predictive signal "
        "in routine data is linear and carried by a small number of "
        "variables, above all ALT. This has a practical consequence. A "
        "logistic model can be published as a transparent equation and "
        "implemented anywhere, including on paper; a gradient boosting model "
        "needs software. Where the performance gap is this small, the choice "
        "between them is a deployment decision, not a scientific one, and we "
        "report both so that readers and implementers can choose.",
        "Calibration deserves emphasis because it is where prediction models "
        "most often fail in practice. Both models were reasonably calibrated "
        "on temporal validation, with slopes of 0.84 and 0.89 — mild "
        "over-dispersion, meaning predictions are slightly too extreme — and "
        "intercepts close to zero. This is better than the calibration "
        "slopes of about 0.63 seen in our earlier kidney model, and it held "
        "up across a period that included the COVID-19 pandemic, during which "
        "the prevalence of advanced fibrosis rose in the US population [3, "
        "4]. Decision curve analysis adds the missing piece: discrimination "
        "and calibration do not by themselves show that acting on a model "
        "helps patients. Ours do, over the threshold range where a clinician "
        "might plausibly act — roughly 0.25 to 0.90 — because the harm of a "
        "false positive here is a non-invasive elastography referral, not a "
        "biopsy. At very low thresholds the models add nothing over "
        "referring everyone, which is the correct conclusion for a disease "
        "this common: if the threshold for action is low enough, universal "
        "assessment is the better strategy.",
        "The AST finding illustrates why a model must be read with its "
        "biology. An odds ratio of 0.60 for AST, taken alone, would suggest "
        "that higher AST protects against MASLD — which is false. ALT and "
        "AST are correlated at 0.76, and the transcriptomic data resolve the "
        "apparent paradox: GPT, encoding ALT, is up-regulated in "
        "steatohepatitis liver, while GOT1, encoding AST, is not "
        "differentially expressed. Serum AST elevation in MASLD reflects "
        "release from injured hepatocytes rather than increased "
        "transcription, so a model that holds ALT fixed will load the "
        "remaining AST variation with a negative coefficient. The same "
        "transcriptomic data strengthen the rest of the model: the "
        "lipogenesis axis (FASN, SREBF1, PNPLA3) and the fibrosis programme "
        "(COL1A1, extracellular-matrix and focal-adhesion pathways) are "
        "activated in diseased liver, and the independent reproduction of "
        "AKR1B10 as the top biomarker gene supports the quality of the "
        "re-analysis. These data are contextual, not validation — the "
        "biopsy-selected cohort with ten controls cannot validate a "
        "population model — but they make the model's reliance on ALT and "
        "GGT biologically coherent rather than merely statistical.",
        "Clinically, the models are best understood as a first step in a "
        "case-finding pathway, not as a diagnostic test. A primary-care "
        "clinician already has age, sex, smoking, education, cholesterol and "
        "liver enzymes for most adults; the model converts those into a "
        "probability that can prioritise elastography referral, particularly "
        "for people above the Youden threshold, where sensitivity is about "
        "72% and specificity 57–61%. This matters more now than it did a "
        "few years ago: resmetirom improved NASH resolution (25.9–29.9% "
        "versus 9.7% with placebo) and fibrosis (24.2–25.9% versus 14.2%) in "
        "the MAESTRO-NASH trial [6], so identifying people with MASLD and "
        "fibrosis has a therapeutic consequence it previously lacked. The "
        "models do not replace elastography, do not diagnose MASH, and do "
        "not stage fibrosis — the MASLD-with-fibrosis analysis was the "
        "weakest in this study — and they should not be used to withhold "
        "assessment from people with metabolic risk factors.",
        "Strengths of this study are the nationally representative sample, "
        "the temporal validation across two survey cycles, the "
        "survey-weighted analysis with bootstrap confidence intervals, the "
        "reporting of calibration, threshold performance and decision curves "
        "alongside discrimination, the strict exclusion of outcome-defining "
        "variables with a quantified leakage demonstration, and the "
        "transcriptomic context. Limitations are real and should be stated "
        "plainly. First, the design is cross-sectional: the outcome is "
        "prevalent MASLD at a single examination, not future risk, and the "
        "word 'prediction' here means identification, not prognosis. "
        "Second, MASLD is defined by elastography, not biopsy, and the "
        "248 dB/m threshold is sensitive rather than specific; the 285 dB/m "
        "sensitivity analysis gives similar discrimination with a lower "
        "prevalence. Third, alcohol intake is self-reported and viral "
        "hepatitis testing was incomplete in some participants, so some "
        "misclassification of the MASLD, MetALD and alcohol-related "
        "categories is inevitable. Fourth, the physical activity instrument "
        "changed between cycles, so MVPA is not directly comparable across "
        "the development and validation samples — although MVPA contributes "
        "little to either model. Fifth, the transcriptomic analysis has a "
        "control group of ten, a biopsy-selected cohort, bulk tissue without "
        "cell-type deconvolution, no covariate adjustment, and a t-test "
        "rather than a count-based model. Sixth, validation is temporal "
        "within one country; geographic and health-system external "
        "validation is still needed, and performance in non-US populations "
        "is unknown.",
        "Future work should follow three directions. Longitudinal follow-up "
        "with linkage to outcomes would turn identification into genuine "
        "risk prediction. External validation in other population surveys "
        "and in primary-care records would test transportability. And a "
        "prospective implementation study — model probability computed "
        "automatically from the electronic record, elastography offered "
        "above a pre-specified threshold, with resmetirom-era treatment "
        "pathways — would test whether the net benefit shown here survives "
        "contact with clinical practice.",
    ]:
        para(doc, text)
    doc.add_heading("Conclusions", level=1)
    para(doc, "Models built only from routine demographic and laboratory "
         "variables identify MASLD in US adults with modest discrimination "
         "and reasonable calibration on temporal validation. Gradient "
         "boosting is genuinely, if modestly, better than logistic "
         "regression, and both models offer net benefit over treating all "
         "patients across the clinically plausible threshold range. The "
         "liver-enzyme signal that drives the models is supported by liver "
         "transcriptomic data. The models are a risk-stratification aid for "
         "prioritising elastography, not a diagnostic rule.")

    # ---------------------------------------------------------- declarations
    doc.add_heading("Declarations", level=1)
    para(doc, "This study used "
         "publicly available, de-identified data from the National Health "
         "and Nutrition Examination Survey, collected with informed consent "
         "under approval of the NCHS Research Ethics Review Board. No "
         "additional ethics approval was required for this secondary "
         "analysis. The transcriptomic re-analysis used publicly available "
         "de-identified data (GEO accession GSE135251).",
         bold_prefix="Ethics approval and consent to participate. ")
    para(doc, "Not applicable.",
         bold_prefix="Consent for publication. ")
    para(doc, "NHANES data are publicly available from the Centers for "
         "Disease Control and Prevention (https://www.cdc.gov/nchs/nhanes/). "
         "The liver RNA-seq data are available from the Gene Expression "
         "Omnibus (accession GSE135251). All analysis code (cohort "
         "construction, modelling, figure generation, manuscript build and "
         "transcriptomic re-analysis), the analysis cohort and the model "
         "predictions are publicly available at [repository URL — to be "
         "added by the author on deposition].",
         bold_prefix="Availability of data and materials. ")
    para(doc, "The author declares no competing interests. [To be confirmed "
         "by the author before submission.]",
         bold_prefix="Competing interests. ")
    para(doc, "No specific funding was received for this work. [To be "
         "confirmed by the author before submission.]",
         bold_prefix="Funding. ")
    para(doc, "Ibrahim H. Hemdan: conceptualisation, data curation, formal "
         "analysis, writing — original draft. [Co-authors, affiliations and "
         "contributions to be added by the author.]",
         bold_prefix="Authors' contributions. ")
    para(doc, "[To be completed by the author.]",
         bold_prefix="Acknowledgements. ")
    para(doc, "Generative AI-assisted tools were used in drafting this "
         "manuscript and in writing analysis code. The author reviewed and "
         "edited all content, verified all numerical results against the "
         "underlying data, and takes full responsibility for the submitted "
         "work. [Wording to be confirmed by the author to match the target "
         "journal's policy.]",
         bold_prefix="Declaration of generative AI use. ")

    # ---------------------------------------------------------- references
    doc.add_heading("References", level=1)
    refs = [
        "Rinella ME, Lazarus JV, Ratziu V, et al. A multi-society Delphi "
        "consensus statement on new fatty liver disease nomenclature. "
        "Hepatology. 2023;78(6):1966-1986.",
        "Younossi ZM, Golabi P, Paik JM, Henry L, Van Dongen C, Alqahtani SA. "
        "The global epidemiology of nonalcoholic fatty liver disease (NAFLD) "
        "and nonalcoholic steatohepatitis (NASH): a systematic review. "
        "Hepatology. 2023;77(4):1335-1347.",
        "Kim D, Danpanichkul P, Wijarnpreecha K, Cholankeril G, Loomba R, "
        "Ahmed A. Current burden of steatotic liver disease and fibrosis "
        "among adults in the United States, 2017–2023. Clin Mol Hepatol. "
        "2025;31(2):382-393.",
        "Unalp-Arida A, et al. Prevalence of metabolic dysfunction-associated "
        "steatotic liver disease and fibrosis defined by liver elastography "
        "in the United States using National Health and Nutrition Examination "
        "Survey 2017–March 2020 and August 2021–August 2023 data. Hepatology. "
        "2025. PMID: 39724144.",
        "European Association for the Study of the Liver; European Association "
        "for the Study of Diabetes; European Association for the Study of "
        "Obesity. EASL–EASD–EASO Clinical Practice Guidelines on the "
        "management of metabolic dysfunction-associated steatotic liver "
        "disease (MASLD). J Hepatol. 2024;81(3):492-542.",
        "Harrison SA, Bedossa P, Guy CD, et al. A phase 3, randomized, "
        "controlled trial of resmetirom in NASH with liver fibrosis. N Engl "
        "J Med. 2024;390(6):497-509.",
        "Bedogni G, Bellentani S, Miglioli L, et al. The Fatty Liver Index: "
        "a simple and accurate predictor of hepatic steatosis in the general "
        "population. BMC Gastroenterol. 2006;6:33.",
        "Ye J, et al. Machine learning models for screening steatosis in "
        "MASLD: an international validation study. Liver Int. "
        "2025;45(12):e70416.",
        "Zhu G, et al. Machine learning for predicting metabolic "
        "dysfunction-associated steatotic liver disease prevalence based on "
        "demographic and clinical characteristics. J Transl Med. 2025;23:381.",
        "Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated "
        "guidance for reporting clinical prediction models that use regression "
        "or machine learning methods. BMJ. 2024;385:e078378.",
        "Johnson CL, Dohrmann SM, Burt VL, Mohadjer LK. National Health and "
        "Nutrition Examination Survey: analytic guidelines, 2011–2014 and "
        "2015–2016. Vital Health Stat 2. 2013;(161):1-24.",
        "National Center for Health Statistics. NHANES 2017–March 2020 "
        "pre-pandemic: liver ultrasound transient elastography (P_LUX) "
        "documentation. Hyattsville, MD: CDC; 2022.",
        "National Center for Health Statistics. NHANES August 2021–August "
        "2023: liver ultrasound transient elastography (LUX_L) documentation. "
        "Hyattsville, MD: CDC; 2024.",
        "Karlas T, Petroff D, Sasso M, et al. Individual patient data "
        "meta-analysis of controlled attenuation parameter (CAP) technology "
        "for assessing steatosis. J Hepatol. 2017;66(5):1022-1030.",
        "Vickers AJ, Elkin EB. Decision curve analysis: a novel method for "
        "evaluating prediction models. Med Decis Making. 2006;26(6):565-574.",
        "Vickers AJ, van Calster B, Steyerberg EW. A simple, step-by-step "
        "guide to interpreting decision curve analysis. Diagn Progn Res. "
        "2019;3:18.",
        "Govaere O, Cockell S, Tiniakos D, et al. Transcriptomic profiling "
        "across the nonalcoholic fatty liver disease spectrum reveals gene "
        "signatures for steatohepatitis and fibrosis. Sci Transl Med. "
        "2020;12(572):eaba4448.",
    ]
    for i, ref in enumerate(refs, 1):
        p = doc.add_paragraph(style="Normal")
        p.paragraph_format.left_indent = Inches(0.3)
        p.paragraph_format.first_line_indent = Inches(-0.3)
        p.paragraph_format.line_spacing = 1.0
        p.paragraph_format.space_after = Pt(3)
        r = p.add_run(f"[{i}] {ref}")

    # ---------------------------------------------------------- tables
    doc.add_page_break()
    doc.add_heading("Tables", level=1)

    with open(HERE / "table1.csv") as f:
        t1_rows = list(csv.reader(f))
    t1_headers = ["Characteristic",
                  "Development\nMASLD (n = 3,709)",
                  "Development\nNo MASLD (n = 2,611)",
                  "Validation\nMASLD (n = 2,265)",
                  "Validation\nNo MASLD (n = 1,838)"]
    add_table(doc, "Table 1. Baseline characteristics by MASLD status, "
              "survey-weighted", t1_headers, t1_rows[1:],
              widths=[2.2, 1.05, 1.05, 1.05, 1.05])
    para(doc, "Values are survey-weighted means (standard deviation), "
         "medians [interquartile range] or percentages. * Triglycerides and "
         "fasting glucose were measured in the morning fasting subsample "
         "only. MVPA = leisure-time moderate-to-vigorous physical activity; "
         "the physical activity instrument changed between NHANES cycles, so "
         "MVPA values are not directly comparable across cycles.",
         style="Footnote")

    v = RESULTS["validation"]
    bench = RESULTS["benchmark_brier"]

    def ci(pair):
        return f"{pair[0]:.3f}–{pair[1]:.3f}"

    t2_rows = [
        ["Logistic regression",
         f"0.698 ({ci(v['logistic']['ci95']['auroc'])})",
         f"0.729 ({ci(v['logistic']['ci95']['average_precision'])})",
         f"0.218 ({ci(v['logistic']['ci95']['brier'])})",
         f"0.027 ({ci(v['logistic']['ci95']['calibration_intercept'])})",
         f"0.838 ({ci(v['logistic']['ci95']['calibration_slope'])})"],
        ["Gradient boosting",
         f"0.721 ({ci(v['gbm']['ci95']['auroc'])})",
         f"0.743 ({ci(v['gbm']['ci95']['average_precision'])})",
         f"0.211 ({ci(v['gbm']['ci95']['brier'])})",
         f"0.044 ({ci(v['gbm']['ci95']['calibration_intercept'])})",
         f"0.891 ({ci(v['gbm']['ci95']['calibration_slope'])})"],
    ]
    add_table(doc, "Table 2. Model performance on temporal validation "
              "(August 2021–August 2023, n = 4,103; 2,265 MASLD cases)",
              ["Model", "AUROC (95% CI)", "Average precision (95% CI)",
               "Brier score (95% CI)", "Calibration intercept (95% CI)",
               "Calibration slope (95% CI)"],
              t2_rows, widths=[1.25, 1.1, 1.15, 1.0, 1.05, 1.05], size=8)
    para(doc, "Paired differences (gradient boosting minus logistic "
         "regression): AUROC +0.023 (95% CI 0.020–0.026); Brier score "
         "−0.0072 (95% CI −0.0083 to −0.0062). Constant-prevalence benchmark "
         "Brier score 0.248. Mean predicted probability 54.3% (logistic) and "
         "53.9% (gradient boosting) versus observed prevalence 54.8%. "
         "Confidence intervals from 1,500 rescaled Rao–Wu bootstrap "
         "replicates. AUROC = area under the receiver operating "
         "characteristic curve.", style="Footnote")

    th = RESULTS["thresholds"]

    def thr_row(model, rule, label):
        d = th[model][rule]
        c = d["ci95"]
        return [label, rule, f"{d['threshold']:.3f}",
                f"{100 * d['sensitivity']:.1f} "
                f"({100 * c['sensitivity'][0]:.1f}–{100 * c['sensitivity'][1]:.1f})",
                f"{100 * d['specificity']:.1f} "
                f"({100 * c['specificity'][0]:.1f}–{100 * c['specificity'][1]:.1f})",
                f"{100 * d['ppv']:.1f}", f"{100 * d['npv']:.1f}"]

    t3_rows = [
        thr_row("logistic", "youden", "Logistic regression"),
        thr_row("logistic", "sens90", "Logistic regression"),
        thr_row("gbm", "youden", "Gradient boosting"),
        thr_row("gbm", "sens90", "Gradient boosting"),
    ]
    add_table(doc, "Table 3. Threshold performance on temporal validation",
              ["Model", "Rule", "Threshold", "Sensitivity, % (95% CI)",
               "Specificity, % (95% CI)", "PPV, %", "NPV, %"],
              t3_rows, widths=[1.25, 1.0, 0.75, 1.2, 1.2, 0.55, 0.55],
              size=8)
    para(doc, "PPV = positive predictive value; NPV = negative predictive "
         "value. Confidence intervals at the fixed point-estimate thresholds, "
         "from 1,500 rescaled Rao–Wu bootstrap replicates.",
         style="Footnote")

    with open(HERE / "odds_ratios.csv") as f:
        or_rows = list(csv.DictReader(f))
    t4_rows = [[r["predictor"],
                f"{float(r['odds_ratio']):.2f} "
                f"({float(r['or_ci95_lo']):.2f}–{float(r['or_ci95_hi']):.2f})"]
               for r in or_rows]
    add_table(doc, "Table 4. Adjusted odds ratios from the logistic "
              "regression model (development sample, n = 6,320)",
              ["Predictor", "Odds ratio (95% CI)"], t4_rows,
              widths=[4.4, 2.1], size=9)
    para(doc, "Confidence intervals from 1,500 nonparametric bootstrap "
         "refits. The AST odds ratio below 1 is a collinearity artefact "
         "(ALT–AST correlation 0.76) and must not be interpreted as a "
         "protective effect; see text.", style="Footnote")

    import pandas as pd
    mt = pd.read_csv(HERE / "metrics_table.csv")
    t5_rows = []
    for _, r in mt.iterrows():
        if r["analysis"] == "Primary outcome (MASLD)":
            continue
        sub = mt[(mt["analysis"] == r["analysis"])]
        lr = sub[sub["model"] == "logistic"].iloc[0]
        gb = sub[sub["model"] == "gbm"].iloc[0]
        t5_rows.append([r["analysis"],
                        f"{int(r['n_val'])} ({int(r['events_val'])})",
                        f"{lr['auroc']:.3f}", f"{gb['auroc']:.3f}",
                        f"{lr['brier']:.3f}", f"{gb['brier']:.3f}",
                        f"{lr['calibration_slope']:.2f}",
                        f"{gb['calibration_slope']:.2f}"])
    # de-duplicate (one row per analysis)
    seen = set()
    t5_rows = [r for r in t5_rows
               if not (r[0] in seen or seen.add(r[0]))]
    add_table(doc, "Table 5. Secondary outcomes and sensitivity analyses "
              "(temporal validation)",
              ["Analysis", "Validation n (events)", "LR AUROC", "GBM AUROC",
               "LR Brier", "GBM Brier", "LR slope", "GBM slope"],
              t5_rows, widths=[2.3, 1.0, 0.65, 0.65, 0.6, 0.6, 0.55, 0.55],
              size=8)
    para(doc, "LR = logistic regression; GBM = gradient boosting. The "
         "criterion-containing predictor set adds BMI, waist circumference, "
         "HbA1c, triglycerides, HDL cholesterol and systolic blood pressure "
         "— variables that define the outcome — and is a leakage "
         "demonstration, not a better model. The Fatty Liver Index analysis "
         "uses the morning fasting subsample only.", style="Footnote")

    t6_rows = [
        ["FASN", "+2.56", "2.0 × 10⁻⁴", "Up"],
        ["GGT1", "+1.45", "1.3 × 10⁻³", "Up"],
        ["PNPLA3", "+1.25", "6.0 × 10⁻⁴", "Up"],
        ["COL1A1", "+1.10", "1.1 × 10⁻³", "Up"],
        ["SREBF1", "+1.00", "1.9 × 10⁻²", "Up"],
        ["GPT (ALT)", "+0.89", "1.1 × 10⁻²",
         "Up, below the fold-change threshold"],
        ["TM6SF2", "+0.23", "0.37", "Not significant"],
        ["ACTA2", "−0.36", "0.42", "Not significant"],
        ["GOT1 (AST)", "−0.13", "0.72", "Not significant"],
        ["PPARA", "+0.04", "0.72", "Not significant"],
    ]
    add_table(doc, "Table 6. MASLD-relevant genes in liver tissue: NASH "
              "(n = 155) versus controls (n = 10), GEO GSE135251",
              ["Gene", "log2 fold-change", "Adjusted p (BH)", "Direction"],
              t6_rows, widths=[1.4, 1.3, 1.4, 2.4], size=9)
    para(doc, "Welch's t-test on log2 counts per million; Benjamini–Hochberg "
         "(BH) correction; 14,971 genes tested; 606 differentially expressed "
         "at false discovery rate < 0.05 and |log2 fold-change| > 1.",
         style="Footnote")

    # ---------------------------------------------------------- figures
    doc.add_page_break()
    doc.add_heading("Figures", level=1)
    figs = [
        ("figures/fig1_flow.png", 5.2,
         "Figure 1. Cohort flow. Of 27,493 interviewed participants, 10,423 "
         "were analysed: 6,320 in development (2017–March 2020; 3,709 with "
         "MASLD) and 4,103 in validation (August 2021–August 2023; 2,265 "
         "with MASLD). CAP = controlled attenuation parameter; LUAXSTAT = "
         "reliable elastography examination indicator."),
        ("figures/fig2_roc.png", 5.0,
         "Figure 2. Receiver operating characteristic curves on temporal "
         "validation. Gradient boosting AUROC 0.721 versus 0.698 for "
         "logistic regression (paired difference +0.023, 95% CI 0.020–0.026)."),
        ("figures/fig3_calibration.png", 5.0,
         "Figure 3. Calibration plots on temporal validation, by decile of "
         "predicted probability. Calibration slopes 0.838 (logistic "
         "regression) and 0.891 (gradient boosting); intercepts 0.027 and "
         "0.044."),
        ("figures/fig4_forest.png", 5.8,
         "Figure 4. Adjusted odds ratios from the logistic regression model "
         "(development sample), with 95% confidence intervals from 1,500 "
         "bootstrap refits. The AST odds ratio below 1 is a collinearity "
         "artefact (see text)."),
        ("figures/fig5_importance.png", 5.2,
         "Figure 5. Permutation variable importance for the gradient boosting "
         "model (validation sample): mean increase in survey-weighted Brier "
         "score when a predictor is permuted, over 10 repeats."),
        ("figures/fig6_dca.png", 5.0,
         "Figure 6. Decision curve analysis on temporal validation. Both "
         "models exceed the treat-all strategy from threshold probabilities "
         "of about 0.25 upward and remain above treat-none up to about 0.90."),
        ("bioinformatics/volcano.png", 5.5,
         "Figure 7. Volcano plot of differential gene expression in liver "
         "tissue: NASH (n = 155) versus controls (n = 10), GEO GSE135251. "
         "Red, up-regulated in NASH; blue, down-regulated (false discovery "
         "rate < 0.05 and |log2 fold-change| > 1; 415 up, 191 down)."),
        ("bioinformatics/enrichment.png", 5.8,
         "Figure 8. KEGG pathway enrichment of the differentially expressed "
         "genes (NASH versus controls, GEO GSE135251). Red, up-regulated in "
         "NASH; blue, down-regulated. Bars show −log10 adjusted p-value; "
         "numbers in parentheses are overlapping genes."),
    ]
    for path, width, caption in figs:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after = Pt(2)
        p.paragraph_format.line_spacing = 1.0
        p.add_run().add_picture(str(HERE / path), width=Inches(width))
        c = doc.add_paragraph(style="Footnote")
        c.paragraph_format.space_after = Pt(6)
        r = c.add_run(caption)
        r.italic = False

    return doc


def count_words(doc):
    """Count words in the abstract and in the main-text sections."""
    counts = {"abstract": 0, "main": 0}
    current = None
    for p in doc.paragraphs:
        if p.style.name == "Heading 1":
            current = p.text.strip()
            continue
        if p.style.name != "Normal":
            continue
        n = len(p.text.split())
        if current == "Abstract":
            counts["abstract"] += n
        elif current in ("Introduction", "Methods", "Results", "Discussion",
                         "Conclusions"):
            counts["main"] += n
    return counts


def main():
    doc = build()
    doc.save(OUT)
    counts = count_words(doc)
    print(f"saved {OUT}")
    print(f"abstract words: {counts['abstract']}")
    print(f"main text words: {counts['main']}")
    print(f"paragraphs: {len(doc.paragraphs)}, tables: {len(doc.tables)}")


if __name__ == "__main__":
    main()
