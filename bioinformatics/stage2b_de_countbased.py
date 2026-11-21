#!/usr/bin/env python3
"""Stage 2b: GSE135251 differential expression with a COUNT-BASED model.

Replaces the Welch t-test of stage2_de_enrichment.py (kept as
stage2_de_enrichment_ttest.py) with a negative-binomial / count-based
method, as required for Q1 review.

Input: counts_ensembl.csv + phenotype.csv (stage 1 outputs — unchanged).
Filtering and HGNC mapping are IDENTICAL to stage 2:
  - CPM > 1 in >= 10 samples at Ensembl level
  - Ensembl -> HGNC symbol (local hgnc_complete_set.txt, mygene.info
    batches of 100 for the remainder), counts summed per symbol
  - CPM > 1 in >= 10 samples again at symbol level -> genes tested
Method (chosen at run time, in preference order):
  (a) limma-voom in R (edgeR::DGEList + calcNormFactors + voom + lmFit +
      eBayes), one fit per contrast on the relevant sample subset;
  (b) PyDESeq2 (DESeq2 implemented in Python), same three contrasts,
      Wald test, BH correction, independent filtering and Cook's-distance
      filtering OFF so every tested gene receives an adjusted p-value.
Contrasts: NASH (155) vs control (10); NAFL (51) vs control;
all-NAFLD (206) vs control. Thresholds: FDR < 0.05 and |log2FC| > 1.
KEGG enrichment: identical to stage 2 (hypergeometric test, BH, 373 human
pathways from KEGG REST API, background = tested genes, all/up/down).

Outputs (same names/formats as stage 2; old t-test files kept as *_ttest.*):
  degs.csv, degs_secondary_contrasts.csv, enrichment.csv,
  counts_symbol.csv (symbol-level integer count matrix used for DE).
The mean_log2CPM_* columns in degs.csv are descriptive (group means of
log2(CPM+1), computed exactly as in stage 2); log2FC / p_value / padj
come from the count-based model.
"""
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

HERE = Path(__file__).resolve().parent
TARGET_GENES = ["PNPLA3", "TM6SF2", "GGT1", "GPT", "GOT1", "FASN",
                "SREBF1", "PPARA", "COL1A1", "ACTA2", "AKR1B10"]
FLAGGED = ["Small cell lung cancer", "Cornified envelope formation",
           "Synaptic vesicle cycle"]


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1.0)
    return out


def map_symbols(ens_ids):
    hgnc = pd.read_csv(HERE / "hgnc_complete_set.txt", sep="\t",
                       usecols=["ensembl_gene_id", "symbol"], dtype=str)
    hgnc = hgnc.dropna(subset=["ensembl_gene_id", "symbol"])
    mapping = dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"]))
    missing = [e for e in ens_ids if e not in mapping]
    print(f"HGNC mapped {len(ens_ids) - len(missing)} / {len(ens_ids)}; "
          f"{len(missing)} to resolve via mygene.info", flush=True)
    url = "https://mygene.info/v3/query"
    for i in range(0, len(missing), 100):
        chunk = missing[i:i + 100]
        for attempt in range(3):
            try:
                r = requests.post(url, json={"q": chunk, "scopes": "ensembl.gene",
                                             "fields": "symbol", "size": 100},
                                  timeout=30)
                r.raise_for_status()
                for hit in r.json():
                    if "symbol" in hit and "query" in hit:
                        mapping[hit["query"]] = hit["symbol"]
                break
            except Exception as e:
                if attempt == 2:
                    print(f"  mygene batch {i} failed after 3 attempts: {e}",
                          flush=True)
                time.sleep(5)
    return mapping


def build_symbol_counts():
    """Identical filtering/mapping to stage 2; returns (symbol_counts,
    logcpm, phenotype-aligned groups)."""
    pheno = pd.read_csv(HERE / "phenotype.csv")
    mat = pd.read_csv(HERE / "counts_ensembl.csv", index_col=0)
    mat = mat[pheno["gsm"]]
    lib = mat.sum(axis=0)
    cpm = mat.div(lib, axis=1) * 1e6
    keep = (cpm > 1).sum(axis=1) >= 10
    print(f"Ensembl IDs passing CPM>1 in >=10 samples: {keep.sum()} / {len(mat)}",
          flush=True)
    mat_f = mat.loc[keep]
    mapping = map_symbols(mat_f.index)
    n_mapped = sum(1 for e in mat_f.index if e in mapping)
    print(f"Total mapped {n_mapped} / {len(mat_f)}", flush=True)
    mat_m = mat_f.loc[mat_f.index.isin(mapping)].copy()
    mat_m.index = [mapping[i] for i in mat_m.index]
    n_before = len(mat_m)
    sym_counts = mat_m.groupby(level=0).sum()
    print(f"Symbol-level genes after collapsing duplicates: {len(sym_counts)} "
          f"(from {n_before} Ensembl IDs)", flush=True)
    lib2 = sym_counts.sum(axis=0)
    cpm2 = sym_counts.div(lib2, axis=1) * 1e6
    keep2 = (cpm2 > 1).sum(axis=1) >= 10
    sym_counts = sym_counts.loc[keep2]
    cpm2 = cpm2.loc[keep2]
    logcpm = np.log2(cpm2 + 1)
    print(f"Genes tested: {len(sym_counts)}", flush=True)
    return sym_counts, logcpm, pheno


def de_table_from_model(logcpm, groups, g1, g0, label, res):
    """Assemble the degs.csv-style table from a model result frame with
    columns gene, log2FC, p_value, padj (already BH-adjusted)."""
    idx1 = groups == g1
    idx0 = groups == g0
    m1 = logcpm.iloc[:, idx1.values].mean(axis=1)
    m0 = logcpm.iloc[:, idx0.values].mean(axis=1)
    out = res.set_index("gene").reindex(logcpm.index)
    out.index.name = "gene"  # reindex() drops the index name (logcpm index is unnamed)
    out = out.reset_index()
    out[f"mean_log2CPM_{label}_cases"] = m1.values
    out[f"mean_log2CPM_{label}_controls"] = m0.values
    out = out[["gene", "log2FC", "p_value", "padj",
               f"mean_log2CPM_{label}_cases",
               f"mean_log2CPM_{label}_controls"]]
    return out.sort_values("padj")


def run_voom(sym_counts, pheno, outdir):
    """Run limma-voom in R, one fit per contrast subset. Returns dict of
    result DataFrames keyed by contrast label."""
    r_script = outdir / "run_voom.R"
    r_script.write_text('''
suppressMessages({library(edgeR); library(limma)})
counts <- read.csv("counts_symbol.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("phenotype.csv", stringsAsFactors = FALSE)
meta <- meta[match(colnames(counts), meta$gsm), ]
stopifnot(all(meta$gsm == colnames(counts)))
grp3 <- ifelse(grepl("^NASH", meta$group), "NASH",
        ifelse(meta$group == "NAFL", "NAFL", "control"))
meta$grp3 <- factor(grp3, levels = c("control", "NAFL", "NASH"))
run_one <- function(cases, label) {
  sel <- meta$grp3 %in% c("control", cases)
  sub <- counts[, sel]
  g <- droplevels(meta$grp3[sel])
  dge <- DGEList(counts = sub)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~0 + g)
  colnames(design) <- levels(g)
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  cm <- makeContrasts(contrasts = paste0(cases, " - control"),
                      levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, number = Inf, sort.by = "P")
  out <- data.frame(gene = rownames(tt),
                    log2FC = tt$logFC,
                    p_value = tt$P.Value,
                    padj = tt$adj.P.Val,
                    stringsAsFactors = FALSE)
  write.csv(out, paste0("voom_", label, ".csv"), row.names = FALSE)
  cat(label, "done: n genes", nrow(out), "\\n")
}
run_one("NASH", "nash")
run_one("NAFL", "nafl")
# all-NAFLD vs control: pooled two-level design on all samples
sel <- rep(TRUE, nrow(meta))
sub <- counts
g <- factor(ifelse(meta$grp3 == "control", "control", "NAFLD"),
            levels = c("control", "NAFLD"))
dge <- DGEList(counts = sub)
dge <- calcNormFactors(dge)
design <- model.matrix(~0 + g)
colnames(design) <- levels(g)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
cm <- makeContrasts(contrasts = "NAFLD - control", levels = design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, number = Inf, sort.by = "P")
out <- data.frame(gene = rownames(tt), log2FC = tt$logFC,
                  p_value = tt$P.Value, padj = tt$adj.P.Val,
                  stringsAsFactors = FALSE)
write.csv(out, "voom_nafld.csv", row.names = FALSE)
cat("nafld done: n genes", nrow(out), "\\n")
''')
    r = subprocess.run(["Rscript", str(r_script)], cwd=str(outdir),
                       capture_output=True, text=True, timeout=3600)
    print(r.stdout[-3000:])
    if r.returncode != 0:
        print("R STDERR:", r.stderr[-3000:])
        raise RuntimeError(f"Rscript failed (exit {r.returncode})")
    return {
        "NASH": pd.read_csv(outdir / "voom_nash.csv"),
        "NAFL": pd.read_csv(outdir / "voom_nafl.csv"),
        "NAFLD": pd.read_csv(outdir / "voom_nafld.csv"),
    }


def run_pydeseq2(sym_counts, pheno, outdir):
    """PyDESeq2, one fit per contrast subset. Wald test, BH, independent
    filtering and Cook's-distance filtering OFF so every tested gene gets
    an adjusted p-value."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    grp3 = pd.Series(np.where(pheno["group"].str.startswith("NASH"), "NASH",
                              np.where(pheno["group"] == "NAFL", "NAFL",
                                       "control")),
                     index=pheno["gsm"])
    results = {}
    for cases, label in [("NASH", "nash"), ("NAFL", "nafl")]:
        sel = grp3.isin(["control", cases])
        counts = sym_counts.loc[:, sel].T  # samples x genes
        meta = pd.DataFrame({"group": grp3[sel].astype(str)})
        meta.index = counts.index
        dds = DeseqDataSet(counts=counts, metadata=meta,
                           design_factors="group", quiet=True)
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["group", cases, "control"],
                          independent_filter=False, cooks_filter=False,
                          quiet=True)
        stat.summary()
        df = stat.results_df.reset_index().rename(
            columns={"index": "gene", "log2FoldChange": "log2FC",
                     "pvalue": "p_value"})
        df = df[["gene", "log2FC", "p_value", "padj"]]
        print(f"PyDESeq2 {label}: {len(df)} genes", flush=True)
        results[label.upper() if label != "nafl" else "NAFL"] = df
    sel = pd.Series(True, index=grp3.index)
    counts = sym_counts.T
    meta = pd.DataFrame({"group": np.where(grp3 == "control", "control",
                                           "NAFLD")})
    meta.index = counts.index
    dds = DeseqDataSet(counts=counts, metadata=meta,
                       design_factors="group", quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["group", "NAFLD", "control"],
                      independent_filter=False, cooks_filter=False,
                      quiet=True)
    stat.summary()
    df = stat.results_df.reset_index().rename(
        columns={"index": "gene", "log2FoldChange": "log2FC",
                 "pvalue": "p_value"})
    results["NAFLD"] = df[["gene", "log2FC", "p_value", "padj"]]
    print(f"PyDESeq2 nafld: {len(df)} genes", flush=True)
    return results


def kegg_pathways():
    base = "https://rest.kegg.jp"
    r = requests.get(f"{base}/list/pathway/hsa", timeout=120)
    r.raise_for_status()
    names = {}
    for line in r.text.strip().split("\n"):
        pid, name = line.split("\t")
        names[pid.replace("path:", "")] = name
    r = requests.get(f"{base}/link/hsa/pathway", timeout=120)
    r.raise_for_status()
    pw_genes = {}
    for line in r.text.strip().split("\n"):
        pid, gid = line.split("\t")
        pw_genes.setdefault(pid.replace("path:", ""), set()).add(gid)
    r = requests.get(f"{base}/list/hsa", timeout=120)
    r.raise_for_status()
    gid2sym = {}
    for line in r.text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        gid, definition = parts[0], parts[3]
        sym = definition.split(";")[0].split(",")[0].strip()
        if sym:
            gid2sym[gid] = sym
    out = {}
    for pid, gids in pw_genes.items():
        syms = {gid2sym[g] for g in gids if g in gid2sym}
        if pid in names and syms:
            out[pid] = (names[pid], syms)
    return out


def enrich(degs, background, pathways, direction):
    degs = set(degs)
    M = len(background)
    N = len(degs & background)
    rows = []
    for pid, (name, genes) in pathways.items():
        pw_bg = genes & background
        n = len(pw_bg)
        if n < 10 or n > 500:
            continue
        k = len(pw_bg & degs)
        if k == 0:
            continue
        p = stats.hypergeom.sf(k - 1, M, n, N)
        rows.append({"pathway_id": pid, "pathway": name, "direction": direction,
                     "overlap": k, "pathway_size_in_background": n,
                     "deg_count": N, "background_size": M, "p_value": p,
                     "overlap_genes": ";".join(sorted(pw_bg & degs))})
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["padj"] = bh_adjust(df["p_value"].values)
    return df.sort_values("padj")


def main():
    method = sys.argv[1] if len(sys.argv) > 1 else "auto"
    sym_counts, logcpm, pheno = build_symbol_counts()
    sym_counts.to_csv(HERE / "counts_symbol.csv")
    print("Saved counts_symbol.csv", flush=True)

    grp = pd.Series(np.where(pheno["group"].str.startswith("NASH"), "NASH",
                             np.where(pheno["group"] == "NAFL", "NAFL",
                                      "control")),
                    index=pheno["gsm"])

    use_r = False
    if method in ("auto", "voom"):
        chk = subprocess.run(["Rscript", "--version"],
                             capture_output=True, text=True)
        if chk.returncode == 0:
            lib = subprocess.run(
                ["Rscript", "-e",
                 'cat(requireNamespace("limma", quietly=TRUE), requireNamespace("edgeR", quietly=TRUE))'],
                capture_output=True, text=True)
            if lib.stdout.strip().endswith("TRUE TRUE"):
                use_r = True
    if method == "voom" and not use_r:
        raise RuntimeError("R + limma/edgeR not available; cannot run voom")

    if use_r:
        print("METHOD: limma-voom (R)", flush=True)
        res = run_voom(sym_counts, pheno, HERE)
        method_name = "limma-voom"
    else:
        print("METHOD: PyDESeq2", flush=True)
        res = run_pydeseq2(sym_counts, pheno, HERE)
        method_name = "PyDESeq2"

    de_nash = de_table_from_model(logcpm, grp, "NASH", "control", "NASH",
                                  res["NASH"])
    de_nafl = de_table_from_model(logcpm, grp, "NAFL", "control", "NAFL",
                                  res["NAFL"])
    de_all = de_table_from_model(logcpm, grp, "NAFLD", "control", "NAFLD",
                                 res["NAFLD"])
    # fix the NAFLD group labels in de_all (grp has no "NAFLD" level)
    is_nafld = pheno["group"].str.startswith("NASH") | (pheno["group"] == "NAFL")
    is_ctrl = pheno["group"] == "control"
    m1 = logcpm.iloc[:, is_nafld.values].mean(axis=1)
    m0 = logcpm.iloc[:, is_ctrl.values].mean(axis=1)
    de_all["mean_log2CPM_NAFLD_cases"] = m1.values
    de_all["mean_log2CPM_NAFLD_controls"] = m0.values

    for name, de in [("NASH vs control", de_nash),
                     ("NAFL vs control", de_nafl),
                     ("NAFLD vs control", de_all)]:
        sig = de[(de["padj"] < 0.05) & (de["log2FC"].abs() > 1)]
        up = (sig["log2FC"] > 0).sum()
        print(f"{name}: {len(sig)} DEGs at FDR<0.05 & |log2FC|>1 "
              f"({up} up, {len(sig) - up} down); "
              f"{(de['padj'] < 0.05).sum()} at FDR<0.05 alone", flush=True)

    de_nash.to_csv(HERE / "degs.csv", index=False)
    sec = de_nafl.rename(columns={"log2FC": "log2FC_NAFL_vs_control",
                                  "p_value": "p_NAFL_vs_control",
                                  "padj": "padj_NAFL_vs_control"}).drop(
        columns=[c for c in de_nafl.columns if c.startswith("mean_")])
    sec2 = de_all.rename(columns={"log2FC": "log2FC_NAFLD_vs_control",
                                  "p_value": "p_NAFLD_vs_control",
                                  "padj": "padj_NAFLD_vs_control"}).drop(
        columns=[c for c in de_all.columns if c.startswith("mean_")])
    sec = sec.merge(sec2, on="gene", how="outer")
    sec.to_csv(HERE / "degs_secondary_contrasts.csv", index=False)

    print("\nTarget genes (NASH vs control):", flush=True)
    tgt = de_nash[de_nash["gene"].isin(TARGET_GENES)][
        ["gene", "log2FC", "p_value", "padj"]]
    print(tgt.to_string(index=False), flush=True)
    missing = [g for g in TARGET_GENES if g not in de_nash["gene"].values]
    print("Target genes not in tested set:",
          missing if missing else "none", flush=True)

    print("\nFetching KEGG pathways ...", flush=True)
    pathways = kegg_pathways()
    print(f"KEGG human pathways retrieved: {len(pathways)}", flush=True)
    background = set(de_nash["gene"])
    sig = de_nash[(de_nash["padj"] < 0.05) & (de_nash["log2FC"].abs() > 1)]
    up = sig[sig["log2FC"] > 0]["gene"]
    down = sig[sig["log2FC"] < 0]["gene"]
    en_all = enrich(sig["gene"], background, pathways, "all_DEGs")
    en_up = enrich(up, background, pathways, "upregulated")
    en_down = enrich(down, background, pathways, "downregulated")
    frames = [d for d in (en_all, en_up, en_down) if len(d)]
    en = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["pathway_id", "pathway", "direction", "overlap",
                 "pathway_size_in_background", "deg_count", "background_size",
                 "p_value", "padj", "overlap_genes"])
    en = en.sort_values(["direction", "padj"])
    en.to_csv(HERE / "enrichment.csv", index=False)
    print(f"Enrichment rows: {len(en)}", flush=True)
    for d in ["all_DEGs", "upregulated", "downregulated"]:
        sub = en[en["direction"] == d]
        print(f"  {d}: {(sub['padj'] < 0.05).sum()} pathways at padj<0.05",
              flush=True)
    print("\nFlagged pathways (new analysis):", flush=True)
    for p in FLAGGED:
        r = en[en["pathway"].str.contains(p, na=False)]
        if len(r):
            print(" ", p, "->",
                  r[["direction", "padj", "overlap"]].to_dict("records"),
                  flush=True)
        else:
            print(" ", p, "-> absent from enrichment.csv", flush=True)
    de_nash.to_pickle(HERE / "de_nash.pkl")
    en.to_pickle(HERE / "enrichment.pkl")
    with open(HERE / "stage2b_method.txt", "w") as f:
        f.write(method_name + "\n")
    print(f"\nMETHOD USED: {method_name}", flush=True)
    print("Saved degs.csv, degs_secondary_contrasts.csv, enrichment.csv",
          flush=True)


if __name__ == "__main__":
    main()
