#!/usr/bin/env python3
"""Stage 2: GSE135251 differential expression + KEGG enrichment.

Primary contrast: NASH (n=155) vs control (n=10).
Secondary contrasts: NAFL (n=51) vs control; all NAFLD (n=206) vs control.
Test: Welch's t-test on log2(CPM+1). BH correction.
Thresholds: FDR < 0.05 and |log2FC| > 1.
Enrichment: KEGG (REST API), hypergeometric test, BH correction.
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(1729)

TARGET_GENES = ["PNPLA3", "TM6SF2", "GGT1", "GPT", "GOT1", "FASN",
                "SREBF1", "PPARA", "COL1A1", "ACTA2"]


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
    """Ensembl -> approved symbol. Primary source: local HGNC complete set
    (authoritative, no network). Remainder: mygene.info in small batches
    of 100 (large POST batches time out from this environment)."""
    hgnc = pd.read_csv(HERE / "hgnc_complete_set.txt", sep="\t",
                       usecols=["ensembl_gene_id", "symbol"], dtype=str)
    hgnc = hgnc.dropna(subset=["ensembl_gene_id", "symbol"])
    mapping = dict(zip(hgnc["ensembl_gene_id"], hgnc["symbol"]))
    missing = [e for e in ens_ids if e not in mapping]
    print(f"HGNC mapped {len(ens_ids) - len(missing)} / {len(ens_ids)}; "
          f"{len(missing)} to resolve via mygene.info")
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
                    print(f"  mygene batch {i} failed after 3 attempts: {e}")
                time.sleep(5)
    return mapping


def de_table(logcpm, groups, g1, g0, label):
    idx1 = groups == g1
    idx0 = groups == g0
    m1 = logcpm.iloc[:, idx1.values].mean(axis=1)
    m0 = logcpm.iloc[:, idx0.values].mean(axis=1)
    log2fc = (m1 - m0).values
    t, p = stats.ttest_ind(logcpm.iloc[:, idx1.values].values,
                           logcpm.iloc[:, idx0.values].values,
                           axis=1, equal_var=False)
    t = np.nan_to_num(t, nan=0.0)
    p = np.nan_to_num(p, nan=1.0)
    out = pd.DataFrame({
        "gene": logcpm.index,
        "log2FC": log2fc,
        "p_value": p,
        "padj": bh_adjust(p),
        f"mean_log2CPM_{label}_cases": m1.values,
        f"mean_log2CPM_{label}_controls": m0.values,
    }).sort_values("padj")
    return out


def kegg_pathways():
    """Return {pathway_id: (name, set(symbols))} for human KEGG pathways."""
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
        # fields: gene id, type (CDS/ncRNA/...), position, definition.
        # The symbol is the first token of the definition field,
        # e.g. "VN1R5, V1RL5; vomeronasal type-1 receptor 5" -> VN1R5.
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
    pheno = pd.read_csv(HERE / "phenotype.csv")
    mat = pd.read_csv(HERE / "counts_ensembl.csv", index_col=0)
    mat = mat[pheno["gsm"]]
    groups = pheno.set_index("gsm")["group"]
    is_nash = groups.str.startswith("NASH")
    is_nafl = groups == "NAFL"
    is_ctrl = groups == "control"

    # --- filter Ensembl IDs, map to symbols, collapse counts ---
    lib = mat.sum(axis=0)
    cpm = mat.div(lib, axis=1) * 1e6
    keep = (cpm > 1).sum(axis=1) >= 10
    print(f"Ensembl IDs passing CPM>1 in >=10 samples: {keep.sum()} / {len(mat)}")
    mat_f = mat.loc[keep]
    mapping = map_symbols(mat_f.index)
    n_mapped = sum(1 for e in mat_f.index if e in mapping)
    print(f"Total mapped {n_mapped} / {len(mat_f)}")
    mat_m = mat_f.loc[mat_f.index.isin(mapping)].copy()
    mat_m.index = [mapping[i] for i in mat_m.index]
    n_before = len(mat_m)
    sym_counts = mat_m.groupby(level=0).sum()
    print(f"Symbol-level genes after collapsing duplicates: {len(sym_counts)} "
          f"(from {n_before} Ensembl IDs)")

    lib2 = sym_counts.sum(axis=0)
    cpm2 = sym_counts.div(lib2, axis=1) * 1e6
    keep2 = (cpm2 > 1).sum(axis=1) >= 10
    sym_counts = sym_counts.loc[keep2]
    cpm2 = cpm2.loc[keep2]
    logcpm = np.log2(cpm2 + 1)
    print(f"Genes tested: {len(sym_counts)}")

    # --- differential expression ---
    grp = pd.Series(np.where(is_nash, "NASH",
                             np.where(is_nafl, "NAFL", "control")),
                    index=mat.columns)
    grp_ctrl = grp[is_ctrl | is_nash]
    de_nash = de_table(logcpm.loc[:, is_nash | is_ctrl], grp_ctrl,
                       "NASH", "control", "NASH")
    sub_nafl = is_nafl | is_ctrl
    grp2 = pd.Series(np.where(is_nafl[sub_nafl], "NAFL", "control"),
                     index=mat.columns[sub_nafl])
    de_nafl = de_table(logcpm.loc[:, is_nafl | is_ctrl], grp2,
                       "NAFL", "control", "NAFL")
    grp3 = pd.Series(np.where(is_nash | is_nafl, "NAFLD", "control"),
                     index=mat.columns)
    de_all = de_table(logcpm, grp3, "NAFLD", "control", "NAFLD")

    for name, de in [("NASH vs control", de_nash),
                     ("NAFL vs control", de_nafl),
                     ("NAFLD vs control", de_all)]:
        sig = de[(de["padj"] < 0.05) & (de["log2FC"].abs() > 1)]
        up = (sig["log2FC"] > 0).sum()
        print(f"{name}: {len(sig)} DEGs at FDR<0.05 & |log2FC|>1 "
              f"({up} up, {len(sig) - up} down); "
              f"{(de['padj'] < 0.05).sum()} at FDR<0.05 alone")

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

    # --- target genes ---
    print("\nTarget genes (NASH vs control):")
    tgt = de_nash[de_nash["gene"].isin(TARGET_GENES)][
        ["gene", "log2FC", "p_value", "padj"]]
    print(tgt.to_string(index=False))
    missing = [g for g in TARGET_GENES if g not in de_nash["gene"].values]
    print("Target genes not in tested set:", missing if missing else "none")

    # --- KEGG enrichment on primary DEGs ---
    print("\nFetching KEGG pathways ...")
    pathways = kegg_pathways()
    print(f"KEGG human pathways retrieved: {len(pathways)}")
    background = set(de_nash["gene"])
    sig = de_nash[(de_nash["padj"] < 0.05) & (de_nash["log2FC"].abs() > 1)]
    up = sig[sig["log2FC"] > 0]["gene"]
    down = sig[sig["log2FC"] < 0]["gene"]
    en_all = enrich(sig["gene"], background, pathways, "all_DEGs")
    en_up = enrich(up, background, pathways, "upregulated")
    en_down = enrich(down, background, pathways, "downregulated")
    frames = [d for d in (en_all, en_up, en_down) if len(d)]
    if frames:
        en = pd.concat(frames, ignore_index=True)
        en = en.sort_values(["direction", "padj"])
    else:
        en = pd.DataFrame(columns=["pathway_id", "pathway", "direction",
                                   "overlap", "pathway_size_in_background",
                                   "deg_count", "background_size", "p_value",
                                   "padj", "overlap_genes"])
        print("WARNING: no KEGG pathway had any overlap with the DEGs")
    en.to_csv(HERE / "enrichment.csv", index=False)
    print(f"Enrichment rows: {len(en)}")
    print(en[en["direction"] == "all_DEGs"].head(12)[
        ["pathway", "overlap", "pathway_size_in_background",
         "p_value", "padj"]].to_string(index=False))
    # save objects for plotting
    de_nash.to_pickle(HERE / "de_nash.pkl")
    en.to_pickle(HERE / "enrichment.pkl")
    print("\nSaved degs.csv, degs_secondary_contrasts.csv, enrichment.csv")


if __name__ == "__main__":
    main()
