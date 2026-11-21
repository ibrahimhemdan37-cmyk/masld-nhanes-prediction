#!/usr/bin/env python3
"""Stage 1: build count matrix from GSE135251 per-sample HTSeq count files,
parse phenotype from the series matrix, QC diagnostics."""
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
MATRIX = HERE / "GSE135251_series_matrix.txt.gz"
COUNTS_DIR = HERE / "raw_counts"

SUMMARY_ROWS = {"__no_feature", "__ambiguous", "__too_low_aQual",
                "__not_aligned", "__alignment_not_unique"}


def parse_phenotype():
    lines = {}
    with gzip.open(MATRIX, "rt") as f:
        for line in f:
            if line.startswith("!Sample_"):
                key, _, rest = line[1:].partition("\t")
                lines.setdefault(key, []).append(
                    rest.rstrip("\n").split("\t"))

    def clean(v):
        return v.strip().strip('"')

    acc = [clean(x) for x in lines["Sample_geo_accession"][0]]
    chars = {}
    for vals in lines["Sample_characteristics_ch1"]:
        vals = [clean(x) for x in vals]
        field = vals[0].split(":")[0].strip()
        chars[field] = [v.split(":", 1)[1].strip() for v in vals]
    pheno = pd.DataFrame({
        "gsm": acc,
        "source_name": [clean(x) for x in lines["Sample_source_name_ch1"][0]],
        "nas_score": pd.to_numeric(chars["nas score"]),
        "fibrosis_stage": pd.to_numeric(chars["fibrosis stage"]),
        "group": chars["group in paper"],
        "disease": chars["disease"],
        "stage": chars["Stage"],
    })
    supp = [clean(x).split("/")[-1] for x in lines["Sample_supplementary_file_1"][0]]
    pheno["count_file"] = supp
    return pheno


def build_counts(pheno):
    frames = []
    for _, row in pheno.iterrows():
        path = COUNTS_DIR / row["count_file"]
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["gene", row["gsm"]])
        frames.append(df.set_index("gene"))
    mat = pd.concat(frames, axis=1)
    mat = mat.drop(index=[r for r in SUMMARY_ROWS if r in mat.index])
    return mat


if __name__ == "__main__":
    pheno = parse_phenotype()
    print("Phenotype:")
    print(pheno["group"].value_counts().to_string())
    print("\nNAS by group:")
    print(pheno.groupby("group")["nas_score"].describe().round(2).to_string())
    print("\nFibrosis stage by group:")
    print(pd.crosstab(pheno["group"], pheno["fibrosis_stage"]).to_string())
    pheno.to_csv(HERE / "phenotype.csv", index=False)

    mat = build_counts(pheno)
    print(f"\nCount matrix: {mat.shape[0]} genes x {mat.shape[1]} samples")
    lib = mat.sum(axis=0)
    det = (mat > 0).sum(axis=0)
    qc = pd.DataFrame({"gsm": mat.columns, "library_size": lib.values,
                       "detected_genes": det.values})
    qc = qc.merge(pheno[["gsm", "group"]], on="gsm")
    print("\nLibrary size (total counts per sample):")
    print(qc.groupby("group")["library_size"].describe().round(0).to_string())
    print("\nDetected genes per sample:")
    print(qc.groupby("group")["detected_genes"].describe().round(0).to_string())
    print("\nSmallest libraries:")
    print(qc.nsmallest(10, "library_size").to_string(index=False))
    print("\nFewest detected genes:")
    print(qc.nsmallest(10, "detected_genes").to_string(index=False))
    mat.to_csv(HERE / "counts_ensembl.csv")
    qc.to_csv(HERE / "qc_library.csv", index=False)
    print("\nSaved phenotype.csv, counts_ensembl.csv, qc_library.csv")
