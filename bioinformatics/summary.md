# Transcriptomic context — MASLD liver RNA-seq (GSE135251)

## Dataset

- **Accession:** GSE135251 (GEO), "TRANSCRIPTOMIC PROFILING ACROSS THE
  SPECTRUM OF NON-ALCOHOLIC FATTY LIVER DISEASE".
- **Source publication:** Govaere O, et al. Transcriptomic profiling across
  the non-alcoholic fatty liver disease spectrum reveals gene signatures for
  steatohepatitis and fibrosis. *Science Translational Medicine* 2020
  (PMID 33268509).
- **Platform:** Illumina NextSeq 500 (GPL18573), RNA-seq of snap-frozen
  liver biopsies; gene-level HTSeq counts supplied per sample.
- **Samples:** 216 biopsies — 206 NAFLD cases and 10 controls. Histological
  groups as recorded in the GEO metadata ("group in paper"): control n = 10,
  NAFL n = 51, NASH F0–F1 n = 34, NASH F2 n = 53, NASH F3 n = 54,
  NASH F4 n = 14. Mean NAFLD activity score (NAS): controls 0.0, NAFL 2.47,
  NASH F0–F1 4.82, NASH F2 5.21, NASH F3 5.48, NASH F4 4.86.
- **Terminology note:** the dataset predates the 2023 MASLD rename; its
  NAFL/NASH groups are the same disease spectrum as MASLD/MASH under the
  current nomenclature.

## Quality control

- Count files: 64,253 Ensembl genes per sample; the five HTSeq summary rows
  (`__no_feature`, `__ambiguous`, `__too_low_aQual`, `__not_aligned`,
  `__alignment_not_unique`) were removed before analysis.
- Library sizes 14.7–37.8 million counts per sample (median ≈ 24 million);
  detected genes 19,096–25,640 per sample. No sample had aberrantly low
  depth.
- PCA on log2(CPM+1) (15,619 genes): PC1 27.3%, PC2 11.1% of variance; the
  most extreme sample was 3.07 SD from the mean on PC1 and 3.89 SD on PC2.
- **No samples were excluded: all 216 were retained for analysis.**

## Methods

- Genes kept if CPM > 1 in at least 10 samples (the size of the smallest
  group): 15,619 Ensembl IDs.
- Ensembl IDs mapped to HGNC-approved symbols (local HGNC complete set;
  remainder via mygene.info): 14,972 mapped, 647 unmapped and dropped;
  counts summed per symbol, leaving **14,971 genes tested**.
- Differential expression: **PyDESeq2 0.6.0** — the DESeq2 negative-binomial
  model implemented in Python (Muzellec B, Teleńczuk M, Cabeli V, Andreux
  M. *Bioinformatics* 2023;39(9):btad547; implementing Love MI, Huber W,
  Anders S. *Genome Biol* 2014;15:550). One DESeq2 fit per contrast
  subset, median-of-ratios size factors, Wald test, Benjamini–Hochberg
  correction. Independent filtering and Cook's-distance filtering were
  switched OFF so every tested gene receives an adjusted p-value.
  log2 fold-change is the DESeq2 shrunken-free (MLE) estimate on the
  count scale.
- **Method change (21 November 2026):** the original analysis used
  Welch's t-test on log2(CPM+1). That was replaced with the count-based
  model above after reviewer feedback, because RNA-seq counts follow a
  negative-binomial distribution. limma-voom in R was the first choice,
  but R could not be installed in this environment — `apt-get update`
  stalled with no network progress for >10 minutes, so the R route was
  abandoned within the bounded install effort and PyDESeq2 used instead.
  The old t-test results are kept as backups with a `_ttest` suffix
  (`degs_ttest.csv`, `degs_secondary_contrasts_ttest.csv`,
  `enrichment_ttest.csv`, `volcano_ttest.png`, `enrichment_ttest.png`,
  `stage2_de_enrichment_ttest.py`); see the comparison table below.
- **Thresholds: FDR < 0.05 and |log2 fold-change| > 1.**
- Primary contrast: **NASH (n = 155, all fibrosis stages) vs control
  (n = 10)**. Secondary contrasts: NAFL (n = 51) vs control; all NAFLD
  (n = 206) vs control.
- Pathway enrichment: KEGG human pathways (373 retrieved live from the KEGG
  REST API on 21 November 2026), hypergeometric test against the 14,971
  tested genes as background, pathways restricted to 10–500 genes in the
  background, Benjamini–Hochberg correction; run separately for all, up-
  and down-regulated differentially expressed genes.
- Analysis scripts: `stage1_build_matrix.py`, `stage2b_de_countbased.py`
  (replaces `stage2_de_enrichment.py`, kept as
  `stage2_de_enrichment_ttest.py`), `stage3_plots.py` (Python 3.12;
  pandas, SciPy, matplotlib; PyDESeq2 0.6.0 in venv
  `~/workspace/masld-venv-bioinfo`).

## Results — differential expression

| Contrast | Genes at FDR < 0.05 | DEGs at FDR < 0.05 and \|log2FC\| > 1 (up / down) |
|---|---|---|
| NASH (n = 155) vs control (n = 10) | 7,322 | **1,427 (844 up / 583 down)** |
| NAFL (n = 51) vs control (n = 10) | 7,145 | 1,438 (756 up / 682 down) |
| All NAFLD (n = 206) vs control (n = 10) | 7,304 | 1,398 (812 up / 586 down) |

Largest fold-changes in NASH vs control: FOSB log2FC −7.34
(padj 6.6 × 10⁻⁴⁶), FOS −6.07 (padj 8.8 × 10⁻³⁰), IGFBP1 −5.04
(padj 3.3 × 10⁻¹⁶), EEF1A2 +6.94 (padj 1.6 × 10⁻¹¹),
AKR1B10 +5.43 (padj 3.8 × 10⁻¹⁸), SCDP1 +4.99 (padj 4.3 × 10⁻⁶).
Most significant: RGS1 −4.40 (padj 1.1 × 10⁻⁴⁸), FOSB −7.34
(padj 6.6 × 10⁻⁴⁶), TCIM −2.44 (padj 9.9 × 10⁻⁴⁵); most significant
up-regulated: KDM4A +1.07 (padj 2.9 × 10⁻³¹), TBC1D16 +1.43
(padj 1.9 × 10⁻³⁰), UNC119B +1.36 (padj 3.9 × 10⁻³⁰).

**AKR1B10 — still a top biomarker, with a stronger signal than before:**
AKR1B10 is the 2nd most up-regulated DEG by fold-change (+5.43, only
EEF1A2 is higher at +6.94) and the 27th most significant up-regulated
gene by adjusted p (3.8 × 10⁻¹⁸). It was the source paper's own top
biomarker (Govaere et al.), and this count-based re-analysis reproduces
that finding with a larger fold-change and a far smaller adjusted p than
the t-test gave (+3.52, padj 4.5 × 10⁻⁹). It is not the single #1 gene
by either criterion — EEF1A2 exceeds it on fold-change and KDM4A on
adjusted p — so the manuscript should describe it as "among the most
strongly up-regulated genes", not as the top gene.

## Results — MASLD-relevant genes (NASH vs control)

| Gene | log2FC | p | Adjusted p (BH) | Direction |
|---|---|---|---|---|
| FASN | +2.58 | 3.9 × 10⁻¹⁵ | 1.9 × 10⁻¹³ | Up |
| GGT1 | +1.32 | 7.4 × 10⁻¹² | 1.6 × 10⁻¹⁰ | Up |
| PNPLA3 | +1.25 | 1.7 × 10⁻⁸ | 1.7 × 10⁻⁷ | Up |
| COL1A1 | +1.50 | 1.6 × 10⁻⁴ | 6.2 × 10⁻⁴ | Up |
| SREBF1 | +0.97 | 3.7 × 10⁻⁴ | 1.3 × 10⁻³ | Up, significant but just below the \|log2FC\| > 1 threshold |
| GPT (ALT) | +0.77 | 4.7 × 10⁻⁴ | 1.6 × 10⁻³ | Up, significant but below the fold-change threshold |
| GOT1 (AST) | −0.50 | 2.8 × 10⁻⁴ | 1.0 × 10⁻³ | Down, significant but below the fold-change threshold |
| TM6SF2 | +0.13 | 0.50 | 0.60 | Not significant |
| ACTA2 | −0.53 | 0.088 | 0.15 | Not significant |
| PPARA | −0.07 | 0.52 | 0.62 | Not significant |

Interpretation: the lipogenesis axis is transcriptionally activated in NASH
liver (FASN, SREBF1, PNPLA3), the fibrosis marker COL1A1 is up-regulated,
and GGT1 and GPT — the genes encoding serum GGT and ALT, two of the
NHANES model's predictors — are up-regulated in NASH liver. GOT1 (AST)
is significantly **down**-regulated at transcript level (−0.50,
padj 1.0 × 10⁻³, below the fold-change threshold), so serum AST
elevation in NASH does not reflect increased GOT1 transcription; it is
consistent with release from injured hepatocytes. TM6SF2 and PPARA show
no expression change, and ACTA2 (the hepatic stellate cell activation
marker) is not significant at bulk tissue level — see Limitations.

## Results — KEGG enrichment (NASH vs control DEGs)

7 pathways significant at adjusted p < 0.05 for all DEGs; 10 for
up-regulated genes; 2 for down-regulated genes. Top results:

**Upregulated in NASH** — ECM–receptor interaction (hsa04512, 18/71 genes,
padj 1.6 × 10⁻⁵), Cornified envelope formation (hsa04382, 18/93,
padj 4.4 × 10⁻⁴), Cytoskeleton in muscle cells (hsa04820, 27/183,
padj 4.4 × 10⁻⁴), Bile secretion (hsa04976, 15/71, padj 6.1 × 10⁻⁴),
Steroid biosynthesis (hsa00100, 7/17, padj 1.3 × 10⁻³),
ABC transporters (hsa02010, 10/37, padj 1.4 × 10⁻³),
Motor proteins (hsa04814, 21/141, padj 1.9 × 10⁻³),
Small cell lung cancer (hsa05222, 15/84, padj 2.5 × 10⁻³),
Integrin signaling (hsa04518, 20/138, padj 3.2 × 10⁻³),
Focal adhesion (hsa04510, 21/177, padj 0.032).

**Downregulated in NASH** — IL-17 signaling (hsa04657, 11/66 genes,
padj 0.012), NF-kappa B signaling (hsa04064, 12/93, padj 0.036).

The up-regulated signal is coherent: bile secretion and steroid
biosynthesis (lipid/sterol metabolism), and extracellular-matrix /
focal-adhesion programmes (fibrogenesis: ECM–receptor interaction,
integrin signaling, focal adhesion). Carbon metabolism, the top
up-regulated pathway under the t-test (padj 0.0021), is no longer
significant (padj 0.97); mineral absorption, TNF and adipocytokine
signaling are likewise no longer significant among down-regulated genes.

**Flagged pathways — verdict after the count-based re-analysis.**
Two of the three biologically implausible pathways flagged in the
t-test analysis remain significant and still appear prominently in
Figure 8 (enrichment.png); they are reported here, not hidden.

- **"Small cell lung cancer" — STILL significant.** All DEGs
  padj 1.4 × 10⁻³ (overlap 22, the joint-top pathway for all DEGs);
  up-regulated padj 2.5 × 10⁻³ (overlap 15), 8th of the 10 bars in
  Figure 8. Why: the KEGG "Small cell lung cancer" gene set is largely
  a cell-cycle/proliferation plus ECM/focal-adhesion set. Its overlap
  genes here are cell-cycle regulators (CCND1, CCNE2, E2F1, CDKN1A,
  CDKN1B, GADD45B/G) and basement-membrane/ECM genes (COL4A5, LAMA2,
  LAMA3, LAMA5, LAMB2, LAMB3) — i.e. hepatocyte proliferation/regeneration
  and fibrotic ECM remodelling in NASH liver, not lung-cancer biology.
  Most of its overlap is shared with ECM–receptor interaction.
- **"Cornified envelope formation" — STILL significant.** All DEGs
  padj 0.026 (overlap 20); up-regulated padj 4.4 × 10⁻⁴ (overlap 18),
  3rd bar in Figure 8. Why: its overlap genes are keratins (KRT7, KRT8,
  KRT10, KRT23), claudins (CLDN3, CLDN4, CLDN9), collagens (COL4A5,
  COL7A1) and laminins — epithelial/ductular-reaction and ECM genes
  up-regulated in fibrotic liver, not skin cornification. KRT7/KRT8 are
  established ductular-reaction markers in NASH.
- **"Synaptic vesicle cycle" — NO LONGER significant.** All DEGs
  padj 0.95 (overlap 6); up-regulated padj 0.76 (overlap 5). It has
  disappeared from the results entirely under the count-based model,
  consistent with it having been statistical noise under the t-test.

The down-regulation of IL-17 and NF-kappa B signaling at bulk-tissue
mRNA level is unexpected for an inflammatory disease and most likely
reflects cell-composition differences and the small control group, not
a genuine anti-inflammatory state — it is reported, not interpreted
further.

## Limitations

1. **Small control group (n = 10)** — the minimum for a two-group
   comparison. Control means are imprecisely estimated, power is limited
   (especially for modest fold-changes), and one control with fibrosis
   stage 2 is included as classified by the source study.
2. **Biopsy-selected cohort.** These are patients who underwent liver
   biopsy for suspected disease, not a general-population sample. The
   cohort differs fundamentally from the NHANES sample in the main paper,
   so this analysis provides biological context only — it cannot validate
   the NHANES prediction model.
3. **Bulk RNA-seq.** Expression mixes hepatocytes, stellate cells, immune
   and endothelial cells. Cell-type composition confounds the results
   (e.g. ACTA2, a stellate-cell marker, is not significant at bulk level).
   A muscle-related signal is present among the top up-regulated genes
   (EEF1A2 +6.94, MYH7 +4.96, MYL1 +4.85, MYH2 +4.72) and pathways
   ("Cytoskeleton in muscle cells", "Motor proteins"), which may reflect
   biopsy composition rather than hepatocyte biology; it is reported, not
   interpreted as a NASH mechanism.
4. **No covariate adjustment.** The GEO metadata contain only NAS,
   fibrosis stage and group — no age, sex, BMI, diabetes status or batch
   information — so no adjustment was possible.
5. **Count-based model (PyDESeq2), not limma-voom.** R could not be
   installed in this environment, so the preferred limma-voom route was
   not run; PyDESeq2 implements the same DESeq2 negative-binomial model
   in Python. Results would differ slightly under limma-voom or R DESeq2,
   but the headline findings (lipogenesis, fibrosis, AKR1B10) are
   high-count, high-fold-change genes and are robust to the choice of
   count-based method.
6. **Single dataset, no replication cohort** in this analysis.
7. **KEGG membership is release-dependent** (gene sets retrieved live
   from the KEGG REST API on 21 November 2026).

## Files

- `degs.csv` — full differential-expression table, NASH vs control
  (14,971 genes: gene, log2FC, p, adjusted p, group mean log2CPM).
  PyDESeq2 results; the t-test version is kept as `degs_ttest.csv`.
- `degs_secondary_contrasts.csv` — NAFL vs control and all-NAFLD vs
  control (log2FC, p, adjusted p per gene); t-test version kept as
  `degs_secondary_contrasts_ttest.csv`.
- `enrichment.csv` — KEGG enrichment (all / up / down DEGs: pathway,
  overlap, pathway size, p, adjusted p, overlapping genes); t-test
  version kept as `enrichment_ttest.csv`.
- `volcano.png` — volcano plot, NASH vs control (t-test version kept as
  `volcano_ttest.png`).
- `enrichment.png` — top enriched KEGG pathways (10 up, 2 down — only
  two down-regulated pathways are significant now; t-test version kept
  as `enrichment_ttest.png`).
- `counts_symbol.csv` — symbol-level count matrix used by the DESeq2
  fits (14,971 genes × 216 samples).
- `stage2b_de_countbased.py` — count-based DE + enrichment code
  (PyDESeq2; includes the limma-voom R route, unused because R could
  not be installed). `stage2b_log.txt`, `stage2b_method.txt`.
- `phenotype.csv`, `qc_library.csv`, `counts_ensembl.csv` — sample
  metadata, QC metrics, Ensembl-level count matrix.
- `stage1_build_matrix.py`, `stage2_de_enrichment_ttest.py`,
  `stage3_plots.py` — analysis code. Raw GEO files:
  `GSE135251_series_matrix.txt.gz`, `GSE135251_RAW.tar`, `raw_counts/`.

## Draft text for the paper ("Transcriptomic context" section)

> To place the clinical predictors in biological context, we re-analysed
> publicly available human liver RNA-seq data (GEO accession GSE135251;
> Govaere et al., Science Translational Medicine 2020), comprising
> 216 snap-frozen liver biopsies (155 with non-alcoholic steatohepatitis
> [NASH], 51 with simple steatosis [NAFL] and 10 histologically normal
> controls). Comparing NASH with controls using DESeq2 (PyDESeq2;
> negative-binomial model, Wald test, Benjamini–Hochberg correction)
> identified 1,427 differentially expressed genes at FDR < 0.05 and
> |log2 fold-change| > 1 (844 up-regulated, 583 down-regulated) out of
> 14,971 genes tested. The lipogenesis axis was activated in NASH liver
> (FASN log2FC +2.58, adjusted p = 1.9 × 10⁻¹³; SREBF1 +0.97,
> adjusted p = 1.3 × 10⁻³; PNPLA3 +1.25, adjusted p = 1.7 × 10⁻⁷), as
> was the fibrosis marker COL1A1 (+1.50, adjusted p = 6.2 × 10⁻⁴).
> GGT1 (+1.32, adjusted p = 1.6 × 10⁻¹⁰) and GPT (+0.77,
> adjusted p = 1.6 × 10⁻³) — encoding serum GGT and ALT, two predictors
> in our model — were up-regulated in NASH liver, whereas GOT1 (AST) was
> down-regulated (−0.50, adjusted p = 1.0 × 10⁻³) and TM6SF2, PPARA and
> ACTA2 were not differentially expressed. AKR1B10, the top biomarker
> reported by Govaere et al., was the second most strongly up-regulated
> gene (+5.43, adjusted p = 3.8 × 10⁻¹⁸). KEGG enrichment of the
> differentially expressed genes highlighted bile secretion, steroid
> biosynthesis and extracellular-matrix / focal-adhesion programmes among
> up-regulated genes, and IL-17 and NF-kappa B signaling among
> down-regulated genes. These data are from a biopsy-selected cohort with
> only 10 controls and provide biological context only; they do not
> validate the prediction model.

## Comparison — Welch's t-test (original) vs PyDESeq2 (current)

Same count matrix, filtering (CPM > 1 in ≥10 samples), HGNC mapping
(14,971 genes tested), contrasts, thresholds (FDR < 0.05,
|log2FC| > 1) and KEGG method in both; only the differential-expression
model changed. Old values are from the `_ttest` backup files.

| Contrast | t-test DEGs (up / down) | PyDESeq2 DEGs (up / down) |
|---|---|---|
| NASH vs control | 606 (415 / 191) | **1,427 (844 / 583)** |
| NAFL vs control | 667 (415 / 252) | **1,438 (756 / 682)** |
| All NAFLD vs control | 603 (403 / 200) | **1,398 (812 / 586)** |
| Genes at FDR < 0.05, NASH vs control | 6,884 | 7,322 |
| KEGG pathways significant (all / up / down) | 15 / 22 / 5 | **7 / 10 / 2** |

| Gene (NASH vs control) | t-test log2FC (padj) | PyDESeq2 log2FC (padj) |
|---|---|---|
| AKR1B10 | +3.52 (4.5 × 10⁻⁹) | **+5.43 (3.8 × 10⁻¹⁸)** |
| FASN | +2.56 (2.0 × 10⁻⁴) | **+2.58 (1.9 × 10⁻¹³)** |
| GGT1 | +1.45 (1.3 × 10⁻³) | **+1.32 (1.6 × 10⁻¹⁰)** |
| PNPLA3 | +1.25 (6.0 × 10⁻⁴) | **+1.25 (1.7 × 10⁻⁷)** |
| COL1A1 | +1.10 (1.1 × 10⁻³) | **+1.50 (6.2 × 10⁻⁴)** |
| SREBF1 | +1.00 (1.9 × 10⁻²) | **+0.97 (1.3 × 10⁻³)** |
| GPT (ALT) | +0.89 (1.1 × 10⁻²) | **+0.77 (1.6 × 10⁻³)** |
| GOT1 (AST) | −0.13 (0.72, NS) | **−0.50 (1.0 × 10⁻³, significant)** |
| TM6SF2 | +0.23 (0.37, NS) | **+0.13 (0.60, NS)** |
| ACTA2 | −0.36 (0.42, NS) | **−0.53 (0.15, NS)** |
| PPARA | +0.04 (0.72, NS) | **−0.07 (0.62, NS)** |

What changed and what did not: the count-based model finds ~2.4× more
DEGs at the same thresholds, and adjusted p-values are far smaller
(the model borrows strength across genes and handles the 155-vs-10
imbalance better than a t-test on log-CPM). The biological story is
unchanged — same lipogenesis/fibrosis genes up, same non-significant
TM6SF2/PPARA/ACTA2 — with one reversal: GOT1 is now significantly
down-regulated rather than non-significant. Enrichment is more focused
(7/10/2 significant pathways vs 15/22/5); carbon metabolism, lysosome,
mineral absorption, TNF and adipocytokine signaling lost significance,
while ECM–receptor interaction, bile secretion and steroid biosynthesis
remain at the top. Of the three flagged implausible pathways, "Synaptic
vesicle cycle" disappeared (padj 0.95), but "Small cell lung cancer"
(all DEGs padj 1.4 × 10⁻³) and "Cornified envelope formation"
(up-regulated padj 4.4 × 10⁻⁴) remain significant — see the verdicts in
the KEGG section above.
