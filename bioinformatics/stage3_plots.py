#!/usr/bin/env python3
"""Stage 3: volcano plot + enrichment bar plot for GSE135251 analysis.
Also strips the ' - Homo sapiens (human)' suffix from KEGG pathway names
in enrichment.csv for readability."""
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent

# --- clean pathway names in enrichment.csv ---
en = pd.read_csv(HERE / "enrichment.csv")
en["pathway"] = en["pathway"].str.replace(" - Homo sapiens (human)", "", regex=False)
en.to_csv(HERE / "enrichment.csv", index=False)

# --- volcano plot: NASH vs control ---
de = pd.read_csv(HERE / "degs.csv")
de["neglog10padj"] = -np.log10(de["padj"].clip(lower=1e-300))
sig = (de["padj"] < 0.05) & (de["log2FC"].abs() > 1)
up = sig & (de["log2FC"] > 0)
down = sig & (de["log2FC"] < 0)

fig, ax = plt.subplots(figsize=(9, 6.5))
ax.scatter(de.loc[~sig, "log2FC"], de.loc[~sig, "neglog10padj"],
           s=6, c="#9aa0a6", alpha=0.55, linewidths=0, label="Not significant")
ax.scatter(de.loc[down, "log2FC"], de.loc[down, "neglog10padj"],
           s=10, c="#2c7fb8", alpha=0.8, linewidths=0,
           label=f"Down in NASH (n={down.sum()})")
ax.scatter(de.loc[up, "log2FC"], de.loc[up, "neglog10padj"],
           s=10, c="#d73027", alpha=0.8, linewidths=0,
           label=f"Up in NASH (n={up.sum()})")
ax.axvline(1, color="grey", lw=0.8, ls="--")
ax.axvline(-1, color="grey", lw=0.8, ls="--")
ax.axhline(-np.log10(0.05), color="grey", lw=0.8, ls="--")
ax.set_xlabel("log$_2$ fold-change (NASH vs control)")
ax.set_ylabel("$-log_{10}$ adjusted p-value (BH)")
ax.set_title("Liver transcriptome: NASH (n=155) vs control (n=10)\nGSE135251",
             fontsize=12)
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.set_xlim(-7.5, 7.5)
ax.spines[["top", "right"]].set_visible(False)

labels = {"AKR1B10": (0.3, 0.3), "FADS2": (0.35, -0.55), "IGFBP1": (-0.3, 0.3),
          "FOS": (-0.3, -0.4), "EGR1": (-0.3, 0.3), "EEF1A2": (0.3, 0.3),
          "SCUBE1": (0.3, 0.3), "LOXL4": (0.3, -0.4), "CYP51A1": (0.3, 0.3),
          "FASN": (0.35, 0.5), "PNPLA3": (-1.5, 0.15), "GGT1": (0.35, 0.55),
          "COL1A1": (0.35, -0.85)}
for gene, (dx, dy) in labels.items():
    row = de[de["gene"] == gene]
    if row.empty:
        continue
    x, y = float(row["log2FC"].iloc[0]), float(row["neglog10padj"].iloc[0])
    ax.annotate(gene, (x, y), xytext=(x + dx, y + dy), fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.6))
fig.tight_layout()
fig.savefig(HERE / "volcano.png", dpi=200)
print("saved volcano.png")

# --- enrichment bar plot: top 10 up + top 5 down ---
up_en = en[(en["direction"] == "upregulated") & (en["padj"] < 0.05)].head(10)
down_en = en[(en["direction"] == "downregulated") & (en["padj"] < 0.05)].head(5)
plot_df = pd.concat([up_en, down_en])
plot_df = plot_df.assign(neglog10padj=-np.log10(plot_df["padj"]))
plot_df = plot_df.sort_values("neglog10padj")
colors = ["#d73027" if d == "upregulated" else "#2c7fb8"
          for d in plot_df["direction"]]

fig, ax = plt.subplots(figsize=(10, 6.5))
y = np.arange(len(plot_df))
ax.barh(y, plot_df["neglog10padj"], color=colors, height=0.65)
ax.set_yticks(y)
ax.set_yticklabels([f"{p} (n={int(o)})" for p, o in
                    zip(plot_df["pathway"], plot_df["overlap"])], fontsize=9)
ax.set_xlabel("$-log_{10}$ adjusted p-value (BH)")
ax.set_title("KEGG enrichment of differentially expressed genes\n"
             "NASH vs control, GSE135251\n(FDR < 0.05, |log$_2$FC| > 1)",
             fontsize=11)
for i, (v, o) in enumerate(zip(plot_df["neglog10padj"], plot_df["overlap"])):
    ax.text(v + 0.05, i, f"{v:.1f}", va="center", fontsize=8)
ax.set_xlim(0, plot_df["neglog10padj"].max() * 1.15)
ax.spines[["top", "right"]].set_visible(False)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#d73027", label="Upregulated in NASH"),
                   Patch(color="#2c7fb8", label="Downregulated in NASH")],
          frameon=False, fontsize=9, loc="lower right")
fig.tight_layout()
fig.savefig(HERE / "enrichment.png", dpi=200)
print("saved enrichment.png")
print(f"up pathways plotted: {len(up_en)}, down: {len(down_en)}")
