###############################################################################
#  02_GSEA_hallmark.R  -  Hallmark GSEA for GSE206648 (ILC2 in dengue)
#
#  Hypothesis-free pathway analysis: all 50 MSigDB Hallmark gene sets are
#  tested; none are chosen in advance.
#
#  Input : results/main/tables/DE_<contrast>_all_genes.csv  (from the main pipeline)
#  Method: preranked GSEA (fgsea), genes ranked by the DESeq2 Wald statistic
#  Output: results/gsea/
###############################################################################

## ---- 1. Packages ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("fgsea",   quietly = TRUE)) BiocManager::install("fgsea", update = FALSE, ask = FALSE)
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
library(fgsea)
library(msigdbr)
library(ggplot2)

set.seed(42)                                   # fgsea uses random sampling
IN_DIR  <- file.path("results", "main", "tables")
OUT_DIR <- file.path("results", "gsea")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

contrasts <- c("DF_vs_Healthy", "DHF_vs_Healthy", "DHF_vs_DF")

## ---- 2. Hallmark gene sets ---------------------------------------------------
# msigdbr changed its arguments in version 10 ("collection" instead of
# "category"), so try the new form first and fall back to the old one.
hm <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
               error = function(e) msigdbr(species = "Homo sapiens", category = "H"))

# our count matrix uses NCBI (Entrez) Gene IDs; the column name differs by version
id_col <- intersect(c("ncbi_gene", "entrez_gene"), colnames(hm))[1]
pathways <- split(as.character(hm[[id_col]]), hm$gs_name)
cat("Loaded", length(pathways), "Hallmark gene sets\n")

## ---- 3. Run GSEA per contrast -----------------------------------------------
all_gsea <- list()
for (nm in contrasts) {
  de <- read.csv(file.path(IN_DIR, paste0("DE_", nm, "_all_genes.csv")))
  de <- de[!is.na(de$stat), ]
  de$GeneID <- as.character(de$GeneID)
  
  # ranked list: Wald stat, high = up in the first group, low = down
  ranks <- sort(setNames(de$stat, de$GeneID), decreasing = TRUE)
  
  res <- fgsea(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500)
  res <- res[order(res$padj), ]
  
  # leading-edge genes: Entrez IDs -> gene symbols (easier to read)
  id2sym <- setNames(de$Symbol, de$GeneID)
  res$leadingEdge <- vapply(res$leadingEdge,
                            function(x) paste(id2sym[x], collapse = ";"), "")
  
  write.csv(as.data.frame(res),
            file.path(OUT_DIR, paste0("GSEA_hallmark_", nm, ".csv")), row.names = FALSE)
  
  res$contrast <- nm
  all_gsea[[nm]] <- res
  cat(nm, ":", sum(res$padj < 0.25), "pathways with FDR < 0.25\n")
}

## ---- 4. Summary dot plot -----------------------------------------------------
# shows every pathway with FDR < 0.25 in at least one contrast
g <- do.call(rbind, all_gsea)
keep <- unique(g$pathway[g$padj < 0.25])
if (length(keep) > 0) {
  g <- g[g$pathway %in% keep, ]
  g$pathway  <- sub("HALLMARK_", "", g$pathway)
  g$contrast <- factor(g$contrast, levels = contrasts)
  p <- ggplot(g, aes(contrast, reorder(pathway, NES), colour = NES, size = -log10(padj))) +
    geom_point() +
    scale_colour_gradient2(low = "#3B6FB6", mid = "white", high = "#C0392B") +
    labs(x = NULL, y = NULL, size = "-log10(FDR)",
         title = "Hallmark GSEA (FDR < 0.25 in at least one comparison)") +
    theme_bw()
  ggsave(file.path(OUT_DIR, "GSEA_hallmark_dotplot.png"), p,
         width = 8, height = max(4, 0.3 * length(keep) + 2), dpi = 150)
} else {
  cat("No pathway reached FDR < 0.25 - no dot plot made\n")
}