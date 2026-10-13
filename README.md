# GSE206648 – DESeq2 pipeline (ILC2 in Dengue)

## Run
1. Install R (≥ 4.2) and RStudio.
2. Open `GSE206648_DESeq2_pipeline.R`, check `DATA_DIR` at the top, and click **Source**.
   The first run installs DESeq2, apeglm, fgsea, msigdbr, ggplot2, ggrepel and pheatmap (~5–10 min).
3. The results are written to `Dengue_virus_analysis/results/`.

You don't need to unzip anything, because R reads the `.gz` files directly.

## Design decisions
| Issue | What the script does |
|---|---|
| Each donor was sequenced 3× (`DF1`, `DF1b`, `DF1c`) | Technical replicates are summed per donor with `collapseReplicates()`, which gives **n = 4 per group**. Treating the 34 libraries as independent would inflate significance. Set `COLLAPSE_TECH_REPS <- FALSE` only for exploration. |
| GSM6260025 (DHF2) is missing from NCBI's count matrix | It's excluded automatically. DHF2 stays in the analysis through its other 2 libraries. |
| Low counts | A gene is kept if it has ≥ 10 reads in ≥ 4 donors. |
| Fold changes | apeglm-shrunken LFC is used for ranking, plots and the significance cut-off. The unshrunk LFC and the Wald stat are kept in the tables. |
| Significance | padj < 0.05 and \|log2FC\| ≥ 1 (editable at the top of the script). |

Contrasts: **DF vs Healthy**, **DHF vs Healthy**, **DHF vs DF**.

## Outputs
**tables/**
- metadata
- library QC
- normalized counts
- VST matrix
- for each contrast: all genes, plus a list of significant genes only
- DE summary
- type I IFN gene score per donor and ISG results
- Hallmark GSEA per contrast

**plots/**
- library sizes and library correlation
- PCA (all libraries and per donor)
- dispersion plot
- sample-distance heatmap
- ILC2 marker check
- volcano, MA and top-40 heatmap for each contrast
- ISG heatmap and ISG score
- GSEA dot plot

`dds_GSE206648.rds` and `vsd_GSE206648.rds` let you reload the analysis with `readRDS()` without re-running it.
