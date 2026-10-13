###############################################################################
#  GSE206648 - ILC2 in Dengue Fever : DESeq2 differential expression pipeline
#
#  Dataset : Hardman / Malavige / Ogg lab (Oxford), Fonseka et al. 2022
#            Blood ILC2 (Lin- CD3- CD45+ CRTH2+ IL-7Ra+), sorted, SmartSeq2,
#            Illumina HiSeq 4000, 75 bp PE
#  Groups  : Healthy (HC, 4 donors) | Dengue fever (DF, 4) | DHF (4)
#            Each donor was sequenced in technical triplicate (e.g. DF1, DF1b, DF1c)
#
#  Inputs (gzipped, NO need to extract - R reads .gz directly):
#     GSE206648_raw_counts_GRCh38.p13_NCBI.tsv.gz   (NCBI raw counts, Entrez IDs)
#     GSE206648_series_matrix.txt.gz                 (sample metadata)
#     Human.GRCh38.p13.annot.tsv.gz                  (gene annotation)
#
#  How to run:
#     RStudio : open this file, set DATA_DIR below, then "Source" (Ctrl+Shift+S)
#     Terminal: Rscript GSE206648_DESeq2_pipeline.R "C:/path/to/Dengue_virus_analysis"
#
#  Output  : <DATA_DIR>/results/
#              main/              core results (start here)
#              qc/                quality-control plots
#              gsea/              Hallmark pathway analysis
#              paper_comparison/  optional check against the published paper
###############################################################################

## ---------------------------------------------------------------------------
## 0. SETTINGS - edit these
## ---------------------------------------------------------------------------
DATA_DIR <- "" # put your own directory where you downloaded the data
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) DATA_DIR <- args[1]

COLLAPSE_TECH_REPS <- TRUE   # sum technical triplicates per donor (recommended)
PADJ_CUTOFF        <- 0.05   # FDR threshold
LFC_CUTOFF         <- 1      # |log2 fold change| threshold (1 = 2-fold)
MIN_COUNT          <- 10     # pre-filter: >= MIN_COUNT reads ...
SET.SEED           <- 42

## ---------------------------------------------------------------------------
## 1. PACKAGES (installed automatically on first run)
## ---------------------------------------------------------------------------
cran_pkgs <- c("ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "msigdbr")
bioc_pkgs <- c("DESeq2", "apeglm", "fgsea")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE))
  BiocManager::install(p, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})
HAS_APEGLM <- requireNamespace("apeglm", quietly = TRUE)
HAS_GSEA   <- requireNamespace("fgsea", quietly = TRUE) &&
  requireNamespace("msigdbr", quietly = TRUE)
set.seed(SET.SEED)

# Output folders
#   results/main/              -> the core analysis: DE tables, donor PCA,
#                                 volcano plots, top-gene heatmaps.
#                                 Start here.
#   results/qc/                -> quality-control plots and tables
#                                 (library sizes, replicate agreement,
#                                 dispersion, MA plots, ILC2 marker check).
#                                 Use these to judge whether the data is
#                                 trustworthy; they are not results.
#   results/gsea/              -> Hallmark pathway analysis (hypothesis-free,
#                                 tests all 50 Hallmark gene sets).
#   results/paper_comparison/  -> type I interferon gene panel. This is a
#                                 hypothesis-driven check taken from
#                                 Fonseka et al. 2022. It is NOT part of the
#                                 independent analysis; use it only to
#                                 compare with the published findings.
OUT          <- file.path(DATA_DIR, "results")
OUT_TAB      <- file.path(OUT, "main", "tables")
OUT_PLOT     <- file.path(OUT, "main", "plots")
OUT_QC       <- file.path(OUT, "qc")
OUT_GSEA     <- file.path(OUT, "gsea")
OUT_PAPER    <- file.path(OUT, "paper_comparison")
for (d in c(OUT_TAB, OUT_PLOT, OUT_QC, OUT_GSEA, OUT_PAPER))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

GROUP_COLS <- c(Healthy = "#4C72B0", DF = "#DD8452", DHF = "#C44E52")

## ---------------------------------------------------------------------------
## 2. LOAD DATA (straight from the .gz files)
## ---------------------------------------------------------------------------
msg("Loading data from ", DATA_DIR)
f_counts <- file.path(DATA_DIR, "GSE206648_raw_counts_GRCh38.p13_NCBI.tsv.gz")
f_meta   <- file.path(DATA_DIR, "GSE206648_series_matrix.txt.gz")
f_annot  <- file.path(DATA_DIR, "Human.GRCh38.p13.annot.tsv.gz")
stopifnot(file.exists(f_counts), file.exists(f_meta), file.exists(f_annot))

counts <- read.delim(gzfile(f_counts), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

annot <- read.delim(gzfile(f_annot), quote = "", stringsAsFactors = FALSE)
annot <- annot[, c("GeneID", "Symbol", "Description", "GeneType", "EnsemblGeneID")]
annot$GeneID <- as.character(annot$GeneID)

# --- parse sample metadata from the series matrix (no GEOquery needed) ---
sm <- readLines(gzfile(f_meta))
get_field <- function(tag) {
  lines <- grep(paste0("^", tag, "\t"), sm, value = TRUE)
  lapply(lines, function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
}
gsm   <- get_field("!Sample_geo_accession")[[1]]
title <- get_field("!Sample_title")[[1]]
chars <- get_field("!Sample_characteristics_ch1")
disease <- sub("^disease state: ", "",
               chars[[grep("disease state", sapply(chars, `[`, 1))]])

meta <- data.frame(
  sample  = gsm,
  title   = title,
  disease = disease,
  stringsAsFactors = FALSE
)
meta$condition <- factor(c("healthy" = "Healthy",
                           "dengue fever" = "DF",
                           "dengue hemorrhagic fever" = "DHF")[meta$disease],
                         levels = c("Healthy", "DF", "DHF"))
meta$donor     <- sub("[bc]$", "", meta$title)          # DF1b -> DF1
meta$replicate <- ifelse(grepl("b$", meta$title), "b",
                         ifelse(grepl("c$", meta$title), "c", "a"))
rownames(meta) <- meta$sample

# keep samples present in BOTH files (NCBI dropped GSM6260025 / DHF2 from counts)
missing <- setdiff(meta$sample, colnames(counts))
if (length(missing))
  msg("NOTE: in metadata but not in count matrix (excluded): ",
      paste(missing, meta[missing, "title"], collapse = ", "))
common <- intersect(colnames(counts), meta$sample)
counts <- counts[, common]
meta   <- meta[common, ]
stopifnot(identical(colnames(counts), rownames(meta)))

write.csv(meta, file.path(OUT_TAB, "00_sample_metadata.csv"), row.names = FALSE)
msg("Count matrix: ", nrow(counts), " genes x ", ncol(counts), " libraries")
print(table(meta$condition))

## ---------------------------------------------------------------------------
## 3. RAW-DATA QC (per library, before collapsing)          -> results/qc/
##    Checks data quality only; nothing here is a biological result.
## ---------------------------------------------------------------------------
msg("Library-level QC")
qc <- data.frame(meta, lib_size = colSums(counts),
                 genes_detected = colSums(counts > 0))
write.csv(qc, file.path(OUT_QC, "library_QC.csv"), row.names = FALSE)

p <- ggplot(qc, aes(x = reorder(title, as.integer(condition)),
                    y = lib_size / 1e6, fill = condition)) +
  geom_col() + scale_fill_manual(values = GROUP_COLS) +
  labs(x = NULL, y = "Assigned reads (millions)", title = "Library sizes") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave(file.path(OUT_QC, "library_sizes.png"), p, width = 10, height = 4.5, dpi = 150)

# technical-replicate concordance (log counts correlation)
lc <- log2(counts + 1)
cor_mat <- cor(lc[rowSums(counts >= MIN_COUNT) >= 3, ], method = "pearson")
ann <- data.frame(condition = meta$condition, donor = meta$donor,
                  row.names = meta$title)
dimnames(cor_mat) <- list(meta$title, meta$title)
pheatmap(cor_mat, annotation_col = ann["condition"],
         annotation_colors = list(condition = GROUP_COLS),
         main = "Library correlation (log2 counts) - technical replicates should cluster",
         filename = file.path(OUT_QC, "library_correlation.png"),
         width = 11, height = 10)

# PCA on all libraries, before collapsing (shows technical replicate spread)
dds_lib <- DESeqDataSetFromMatrix(counts, meta, design = ~ condition)
dds_lib <- dds_lib[rowSums(counts(dds_lib) >= MIN_COUNT) >= 3, ]
vsd_lib <- vst(dds_lib, blind = TRUE)
pc <- plotPCA(vsd_lib, intgroup = c("condition", "donor"), returnData = TRUE)
pv <- round(100 * attr(pc, "percentVar"))
p <- ggplot(pc, aes(PC1, PC2, colour = condition, label = donor)) +
  geom_point(size = 3) + geom_text_repel(size = 3, show.legend = FALSE) +
  scale_colour_manual(values = GROUP_COLS) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       title = "PCA - all libraries (technical replicates labelled by donor)") +
  theme_bw()
ggsave(file.path(OUT_QC, "PCA_all_libraries.png"), p, width = 7.5, height = 6, dpi = 150)

## ---------------------------------------------------------------------------
## 4. BUILD DESeq2 OBJECT (collapse technical replicates -> biological n)
## ---------------------------------------------------------------------------
# Technical replicates are not independent biological samples. Treating them as
# such inflates n and gives far too many false positives. DESeq2's
# collapseReplicates() sums their counts so each donor becomes one sample.
dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ condition)
if (COLLAPSE_TECH_REPS) {
  msg("Collapsing technical replicates by donor")
  dds <- collapseReplicates(dds, groupby = dds$donor, run = dds$title)
  colnames(dds) <- as.character(dds$donor)
}
msg("Samples entering DE: ", ncol(dds))
print(table(dds$condition))

# pre-filter: gene must have >= MIN_COUNT reads in at least the smallest group
min_n <- min(table(dds$condition))
keep  <- rowSums(counts(dds) >= MIN_COUNT) >= min_n
dds   <- dds[keep, ]
msg("Genes kept after filtering: ", nrow(dds))

## ---------------------------------------------------------------------------
## 5. RUN DESeq2
## ---------------------------------------------------------------------------
dds$condition <- relevel(dds$condition, ref = "Healthy")
dds <- DESeq(dds)

png(file.path(OUT_QC, "dispersion_estimates.png"), width = 7, height = 6, units = "in", res = 150)
plotDispEsts(dds, main = "DESeq2 dispersion estimates")
dev.off()

norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = FALSE)
add_annot <- function(df, id_col = "GeneID") {
  out <- merge(annot, df, by.x = "GeneID", by.y = id_col, all.y = TRUE, sort = FALSE)
  out$Symbol[is.na(out$Symbol) | out$Symbol == ""] <- out$GeneID[is.na(out$Symbol) | out$Symbol == ""]
  out
}
write.csv(add_annot(data.frame(GeneID = rownames(norm_counts), norm_counts, check.names = FALSE)),
          file.path(OUT_TAB, "01_normalized_counts.csv"), row.names = FALSE)
write.csv(add_annot(data.frame(GeneID = rownames(vsd), assay(vsd), check.names = FALSE)),
          file.path(OUT_TAB, "02_vst_expression.csv"), row.names = FALSE)

sym <- setNames(annot$Symbol, annot$GeneID)
sym_of <- function(ids) { s <- sym[ids]; ifelse(is.na(s) | s == "", ids, s) }

## ---------------------------------------------------------------------------
## 6. SAMPLE-LEVEL PLOTS (donor level)
## ---------------------------------------------------------------------------
pc <- plotPCA(vsd, intgroup = "condition", returnData = TRUE, ntop = 500)
pv <- round(100 * attr(pc, "percentVar"))
p <- ggplot(pc, aes(PC1, PC2, colour = condition, label = name)) +
  geom_point(size = 4) + geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = GROUP_COLS) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"),
       title = "PCA - donors (top 500 variable genes)") +
  theme_bw()
ggsave(file.path(OUT_PLOT, "01_PCA_donors.png"), p, width = 7, height = 5.5, dpi = 150)

# QC (-> results/qc/): how similar each donor is to every other donor
sd_mat <- as.matrix(dist(t(assay(vsd))))
ann_d  <- data.frame(condition = dds$condition, row.names = colnames(dds))
pheatmap(sd_mat, annotation_col = ann_d, annotation_colors = list(condition = GROUP_COLS),
         color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         main = "Sample-to-sample Euclidean distance (VST)",
         filename = file.path(OUT_QC, "sample_distance_heatmap.png"), width = 7, height = 6)

# QC (-> results/qc/): ILC2 identity check - canonical ILC2 genes should be
# high in every sample, other lineage markers low. Checks the cell sort, not the disease.
ilc2_markers <- c("GATA3", "IL7R", "PTGDR2", "KLRB1", "IL1RL1", "IL17RB",
                  "HPGDS", "KIT", "RORA", "IL2RA", "CD3E", "CD19", "CD14")
idx <- match(ilc2_markers, sym_of(rownames(vsd)))
if (any(!is.na(idx))) {
  m <- assay(vsd)[na.omit(idx), , drop = FALSE]
  rownames(m) <- sym_of(rownames(m))
  pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = ann_d, annotation_colors = list(condition = GROUP_COLS),
           main = "ILC2 identity check (CD3E/CD19/CD14 should be low)",
           filename = file.path(OUT_QC, "ILC2_marker_check.png"), width = 7, height = 5)
}

## ---------------------------------------------------------------------------
## 7. DIFFERENTIAL EXPRESSION - three contrasts
## ---------------------------------------------------------------------------
contrasts <- list(
  DF_vs_Healthy  = c("condition", "DF",  "Healthy"),
  DHF_vs_Healthy = c("condition", "DHF", "Healthy"),
  DHF_vs_DF      = c("condition", "DHF", "DF")
)

# shrunken LFC (apeglm) needs the comparison as a model coefficient, so for
# DHF vs DF we refit with DF as the reference level
get_shrunk <- function(dds, ctr) {
  ref <- ctr[3]; tst <- ctr[2]
  d <- dds
  if (levels(d$condition)[1] != ref) {
    d$condition <- relevel(d$condition, ref = ref)
    d <- nbinomWaldTest(d)
  }
  coef <- paste0("condition_", tst, "_vs_", ref)
  if (HAS_APEGLM) lfcShrink(d, coef = coef, type = "apeglm", quiet = TRUE)
  else            lfcShrink(d, coef = coef, type = "normal", quiet = TRUE)
}

volcano <- function(df, title, file) {
  df <- df[!is.na(df$padj), ]
  df$status <- "NS"
  df$status[df$padj < PADJ_CUTOFF & df$log2FoldChange >=  LFC_CUTOFF] <- "Up"
  df$status[df$padj < PADJ_CUTOFF & df$log2FoldChange <= -LFC_CUTOFF] <- "Down"
  top <- head(df[df$status != "NS", ][order(df$padj[df$status != "NS"]), ], 20)
  n_up <- sum(df$status == "Up"); n_dn <- sum(df$status == "Down")
  p <- ggplot(df, aes(log2FoldChange, -log10(pvalue), colour = status)) +
    geom_point(size = 1, alpha = 0.6) +
    scale_colour_manual(values = c(Down = "#3B6FB6", NS = "grey75", Up = "#C0392B")) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = 2, colour = "grey40") +
    geom_text_repel(data = top, aes(label = Symbol), size = 3, colour = "black",
                    max.overlaps = 30) +
    labs(title = title,
         subtitle = sprintf("padj < %.2f & |log2FC| >= %g :  %d up, %d down",
                            PADJ_CUTOFF, LFC_CUTOFF, n_up, n_dn),
         x = "log2 fold change (shrunken)", y = "-log10(p-value)") +
    theme_bw()
  ggsave(file, p, width = 7.5, height = 6.5, dpi = 150)
}

all_res   <- list()
summary_tab <- data.frame()
for (nm in names(contrasts)) {
  ctr <- contrasts[[nm]]
  msg("Contrast: ", nm)
  res_raw <- results(dds, contrast = ctr, alpha = PADJ_CUTOFF)
  res_shr <- get_shrunk(dds, ctr)
  
  df <- data.frame(GeneID          = rownames(res_raw),
                   baseMean        = res_raw$baseMean,
                   log2FoldChange  = res_shr$log2FoldChange,   # shrunken (for ranking/plots)
                   lfcSE           = res_shr$lfcSE,
                   log2FC_unshrunk = res_raw$log2FoldChange,
                   stat            = res_raw$stat,
                   pvalue          = res_raw$pvalue,
                   padj            = res_raw$padj)
  df <- add_annot(df)
  df <- df[order(df$padj, df$pvalue), ]
  all_res[[nm]] <- df
  
  sig <- subset(df, !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) >= LFC_CUTOFF)
  write.csv(df,  file.path(OUT_TAB, paste0("DE_", nm, "_all_genes.csv")), row.names = FALSE)
  write.csv(sig, file.path(OUT_TAB, paste0("DE_", nm, "_significant.csv")), row.names = FALSE)
  
  summary_tab <- rbind(summary_tab, data.frame(
    contrast = nm,
    padj_lt_cutoff = sum(df$padj < PADJ_CUTOFF, na.rm = TRUE),
    sig_up   = sum(sig$log2FoldChange > 0),
    sig_down = sum(sig$log2FoldChange < 0)))
  
  volcano(df, paste("Volcano:", gsub("_", " ", nm)),
          file.path(OUT_PLOT, paste0("02_volcano_", nm, ".png")))
  
  # QC (-> results/qc/): MA plot - checks fold-change shrinkage behaved sensibly
  png(file.path(OUT_QC, paste0("MA_", nm, ".png")), width = 7, height = 5.5, units = "in", res = 150)
  plotMA(res_shr, alpha = PADJ_CUTOFF, main = paste("MA plot (shrunken):", nm), ylim = c(-6, 6))
  dev.off()
  
  # heatmap of top 40 DE genes for this contrast (all donors shown)
  top_ids <- head(df$GeneID[!is.na(df$padj)], 40)
  m <- assay(vsd)[top_ids, , drop = FALSE]
  m <- t(scale(t(m)))
  rownames(m) <- sym_of(rownames(m))
  pheatmap(m, annotation_col = ann_d, annotation_colors = list(condition = GROUP_COLS),
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = paste("Top 40 genes by padj:", nm, "(row z-score, VST)"),
           fontsize_row = 7,
           filename = file.path(OUT_PLOT, paste0("03_heatmap_top40_", nm, ".png")),
           width = 7, height = 9)
}
write.csv(summary_tab, file.path(OUT_TAB, "03_DE_summary.csv"), row.names = FALSE)
print(summary_tab)

## ---------------------------------------------------------------------------
## 8. PAPER COMPARISON: type I IFN gene panel             -> results/paper_comparison/
##    Hypothesis-driven: the gene list comes from the published study
##    (Fonseka et al. 2022). Not part of the independent analysis - use it
##    only to compare your results with theirs. Safe to delete this section.
## ---------------------------------------------------------------------------
msg("Type I IFN-stimulated gene panel")
isg <- c("IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
         "ISG15", "ISG20", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2",
         "IFITM1", "IFITM3", "STAT1", "STAT2", "IRF7", "IRF9", "XAF1", "HERC5",
         "EPSTI1", "USP18", "BST2", "CMPK2", "LY6E", "SIGLEC1", "IFNAR1", "IFNAR2")
ids <- rownames(vsd)[sym_of(rownames(vsd)) %in% isg]
if (length(ids) > 2) {
  m <- assay(vsd)[ids, , drop = FALSE]
  rownames(m) <- sym_of(rownames(m))
  m <- m[order(rownames(m)), ]
  pheatmap(t(scale(t(m))), cluster_cols = FALSE,
           annotation_col = ann_d, annotation_colors = list(condition = GROUP_COLS),
           color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
           main = "Type I IFN-stimulated genes (row z-score)",
           filename = file.path(OUT_PAPER, "ISG_heatmap.png"), width = 7, height = 8)
  
  # ISG score per donor = mean z-score of ISGs
  score <- colMeans(t(scale(t(m))))
  sdf <- data.frame(donor = names(score), ISG_score = score,
                    condition = dds$condition[match(names(score), colnames(dds))])
  write.csv(sdf, file.path(OUT_PAPER, "ISG_score_per_donor.csv"), row.names = FALSE)
  p <- ggplot(sdf, aes(condition, ISG_score, colour = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 3) +
    scale_colour_manual(values = GROUP_COLS) +
    labs(y = "Mean ISG z-score", x = NULL, title = "Type I IFN signature score") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(OUT_PAPER, "ISG_score.png"), p, width = 5, height = 5, dpi = 150)
  
  isg_tab <- do.call(rbind, lapply(names(all_res), function(nm) {
    d <- all_res[[nm]][all_res[[nm]]$Symbol %in% isg, c("Symbol", "log2FoldChange", "padj")]
    if (nrow(d)) d$contrast <- nm
    d
  }))
  write.csv(isg_tab, file.path(OUT_PAPER, "ISG_results_all_contrasts.csv"), row.names = FALSE)
}

## ---------------------------------------------------------------------------
## 9. GSEA - MSigDB Hallmark gene sets (fgsea, pre-ranked by Wald stat) -> results/gsea/
##    Hypothesis-free: all 50 Hallmark pathways are tested, none chosen in advance.
## ---------------------------------------------------------------------------
if (HAS_GSEA) {
  msg("GSEA (Hallmark)")
  hm <- tryCatch({
    h <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
    h
  }, error = function(e) tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", category = "H"),
    error = function(e2) NULL))
  if (!is.null(hm)) {
    pathways <- split(as.character(hm$ncbi_gene), hm$gs_name)
    gsea_all <- list()
    for (nm in names(all_res)) {
      d <- all_res[[nm]]
      d <- d[!is.na(d$stat), ]
      ranks <- sort(setNames(d$stat, d$GeneID), decreasing = TRUE)
      fg <- fgsea::fgsea(pathways, ranks, minSize = 15, maxSize = 500)
      fg <- fg[order(fg$padj), ]
      fg$leadingEdge <- vapply(fg$leadingEdge,
                               function(x) paste(sym_of(x), collapse = ";"), "")
      write.csv(as.data.frame(fg), file.path(OUT_GSEA, paste0("GSEA_hallmark_", nm, ".csv")),
                row.names = FALSE)
      fg$contrast <- nm
      gsea_all[[nm]] <- fg
    }
    g <- do.call(rbind, gsea_all)
    keep_p <- unique(g$pathway[g$padj < 0.25])
    if (length(keep_p)) {
      g <- g[g$pathway %in% keep_p, ]
      g$pathway <- gsub("HALLMARK_", "", g$pathway)
      p <- ggplot(g, aes(contrast, reorder(pathway, NES), colour = NES, size = -log10(padj))) +
        geom_point() +
        scale_colour_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B") +
        labs(x = NULL, y = NULL, title = "Hallmark GSEA (pathways with FDR < 0.25 in any contrast)") +
        theme_bw()
      ggsave(file.path(OUT_GSEA, "GSEA_hallmark_dotplot.png"), p,
             width = 8, height = max(4, 0.25 * length(keep_p) + 2), dpi = 150)
    }
  } else msg("msigdbr gene sets unavailable - GSEA skipped")
} else msg("fgsea/msigdbr not installed - GSEA skipped")

## ---------------------------------------------------------------------------
## 10. SAVE
## ---------------------------------------------------------------------------
saveRDS(dds, file.path(OUT, "dds_GSE206648.rds"))
saveRDS(vsd, file.path(OUT, "vsd_GSE206648.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
msg("Done. Results in: ", normalizePath(OUT))