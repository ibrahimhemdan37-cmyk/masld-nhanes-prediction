#!/usr/bin/env python3
"""Build the MASLD supplementary material docx (schema-safe python-docx).

Reads tuning_results.csv, model_results.json and data_dictionary.csv from
this directory. Writes MASLD_NHANES_Supplementary.docx.
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
OUT = HERE / "MASLD_NHANES_Supplementary.docx"


def shade_cell(cell, fill):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), fill)
    tcPr.append(shd)


def set_cell_text(cell, text, bold=False, size=9, color=None):
    cell.text = ""
    p = cell.paragraphs[0]
    p.paragraph_format.space_before = Pt(1)
    p.paragraph_format.space_after = Pt(1)
    r = p.add_run(str(text))
    r.font.name = "Times New Roman"
    r.font.size = Pt(size)
    r.bold = bold
    if color is not None:
        r.font.color.rgb = RGBColor.from_string(color)


def add_table(doc, title, headers, rows, widths=None, size=9):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(10)
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.line_spacing = 1.0
    r = p.add_run(title)
    r.bold = True
    r.font.size = Pt(10)
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
            set_cell_text(cells[i], val, size=size)
    doc.add_paragraph("", style="Normal")
    return table


def para(doc, text, bold_prefix=None, size=12, spacing=2.0):
    p = doc.add_paragraph(style="Normal")
    p.paragraph_format.line_spacing = spacing
    if bold_prefix:
        r = p.add_run(bold_prefix)
        r.bold = True
    p.add_run(text)
    return p


def main():
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
    s = doc.styles["Heading 1"]
    s.font.name = "Times New Roman"
    s.font.size = Pt(14)
    s.font.bold = True

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Supplementary Material")
    r.bold = True
    r.font.size = Pt(16)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Identifying Metabolic Dysfunction-Associated Steatotic "
                  "Liver Disease in US Adults: A Cross-Sectional NHANES "
                  "Study with Temporal Validation and Transcriptomic Context")
    r.italic = True
    r.font.size = Pt(11)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Ibrahim H. Hemdan")
    r.font.size = Pt(12)

    doc.add_heading("Table S1. Hyperparameter tuning (five-fold stratified "
                    "cross-validation within the development sample, "
                    "survey-weighted Brier score)", level=1)
    with open(HERE / "tuning_results.csv") as f:
        tun = list(csv.DictReader(f))
    rows = [[r["model"], r["params"], f"{float(r['cv_brier']):.4f}"]
            for r in tun]
    add_table(doc, "", ["Model", "Hyperparameters", "Cross-validated Brier "
              "score"], rows, widths=[1.2, 3.4, 1.9], size=9)
    para(doc, "The selected configurations were C = 0.01 for logistic "
         "regression (all four values of C gave a cross-validated Brier "
         "score of 0.2090 to four decimal places) and 7 leaves, learning "
         "rate 0.05, 100 iterations and L2 = 1.0 for gradient boosting "
         "(cross-validated Brier 0.2002; the L2 = 0.0 configuration gave the "
         "same value to six decimal places).", size=10, spacing=1.0)

    doc.add_heading("Table S2. Threshold performance on temporal validation, "
                    "with full bootstrap confidence intervals", level=1)
    with open(HERE / "model_results.json") as f:
        res = json.load(f)
    th = res["thresholds"]

    def fmt(x):
        return f"{100 * x:.1f}"

    rows = []
    for model, label in [("logistic", "Logistic regression"),
                         ("gbm", "Gradient boosting")]:
        for rule, rlabel in [("youden", "Youden's J"),
                             ("sens90", "90% sensitivity")]:
            d = th[model][rule]
            c = d["ci95"]
            rows.append([label, rlabel, f"{d['threshold']:.3f}",
                         f"{fmt(d['sensitivity'])} ({fmt(c['sensitivity'][0])}"
                         f"–{fmt(c['sensitivity'][1])})",
                         f"{fmt(d['specificity'])} ({fmt(c['specificity'][0])}"
                         f"–{fmt(c['specificity'][1])})",
                         f"{fmt(d['ppv'])} ({fmt(c['ppv'][0])}"
                         f"–{fmt(c['ppv'][1])})",
                         f"{fmt(d['npv'])} ({fmt(c['npv'][0])}"
                         f"–{fmt(c['npv'][1])})"])
    add_table(doc, "", ["Model", "Rule", "Threshold",
              "Sensitivity, % (95% CI)", "Specificity, % (95% CI)",
              "PPV, % (95% CI)", "NPV, % (95% CI)"], rows,
              widths=[1.15, 0.95, 0.65, 1.1, 1.1, 0.95, 0.95], size=8)
    para(doc, "PPV = positive predictive value; NPV = negative predictive "
         "value. Confidence intervals at the fixed point-estimate thresholds, "
         "from 1,500 rescaled Rao–Wu bootstrap replicates.", size=10,
         spacing=1.0)

    doc.add_heading("Table S3. Data dictionary for the analysis cohort "
                    "(nhanes_masld_cross_sectional.csv)", level=1)
    with open(HERE / "data_dictionary.csv") as f:
        dd = list(csv.DictReader(f))
    rows = [[r["variable"], r["description"], r["source"]] for r in dd]
    add_table(doc, "", ["Variable", "Description", "Source"], rows,
              widths=[1.9, 3.3, 1.3], size=8)

    doc.add_heading("TRIPOD+AI reporting", level=1)
    para(doc, "Reporting follows the TRIPOD+AI statement (Collins et al., "
         "BMJ 2024;385:e078378). The location of each reporting domain in "
         "the main manuscript is:", size=11, spacing=1.0)
    items = [
        ("Title and abstract", "Title page; Abstract"),
        ("Background and objectives", "Introduction"),
        ("Data sources", "Methods — Study design and data source"),
        ("Participants", "Methods — Study design and data source; Results — "
         "Study population; Figure 1"),
        ("Outcome definition", "Methods — Outcome definition"),
        ("Predictors", "Methods — Predictors"),
        ("Sample size", "Results — Study population"),
        ("Missing data", "Methods — Predictors; Results — Study population"),
        ("Statistical analysis methods", "Methods — Statistical analysis"),
        ("Model development and tuning", "Methods — Statistical analysis; "
         "Results — Model tuning and primary validation; Table S1"),
        ("Model evaluation (discrimination, calibration, clinical utility)",
         "Methods — Statistical analysis; Results — Model tuning and primary "
         "validation; Threshold performance and clinical utility; Tables 2, "
         "3; Figures 2, 3, 6"),
        ("Results — participants and model performance",
         "Results; Tables 1–5; Figures 1–6"),
        ("Interpretation", "Discussion"),
        ("Limitations", "Discussion"),
        ("Clinical use and implications", "Discussion"),
        ("Ethics, funding, competing interests, data availability",
         "Declarations"),
        ("Transcriptomic analysis", "Methods — Transcriptomic re-analysis; "
         "Results — Transcriptomic context; Table 6; Figures 7, 8"),
    ]
    add_table(doc, "", ["TRIPOD+AI reporting domain",
              "Location in the manuscript"],
              [[a, b] for a, b in items], widths=[2.6, 3.9], size=9)

    doc.save(OUT)
    print(f"saved {OUT}")


if __name__ == "__main__":
    main()
