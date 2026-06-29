#!/usr/bin/env Rscript

################################################################################
# Script: heatmap_36_genes.R
#
# Description:
# This script generates a heatmap of 36 differentially expressed genes identified
# after fitting a linear model corrected for PC1 and PC2.
#
# Gene names are converted from RefSeq IDs to HGNC symbols for improved
# biological interpretability.
################################################################################

library(pheatmap)

# -------------------------
# LOAD EXPRESSION MATRIX
# -------------------------
expr <- read.csv("matriz_normalizada_para_modelos_nb.csv", check.names = FALSE)

rownames(expr) <- make.unique(expr[,1])
expr <- expr[,-1]

# convertir a data.frame con IDs
expr_df <- as.data.frame(expr)
expr_df$gene <- rownames(expr_df)

# limpiar IDs
expr_df$gene_clean <- sub("\\..*$", "", expr_df$gene)

# -------------------------
# LOAD ANNOTATION (HGNC)
# -------------------------
annot <- read.csv("Supplementary_Table_DEG_PC1_PC2_HGNC.csv")

annot$gene_clean <- sub("\\..*$", "", annot$gene)

# -------------------------
# MERGE REAL (CLAVE)
# -------------------------
merged <- merge(annot, expr_df, by = "gene_clean")

# -------------------------
# USAR SYMBOLS COMO ROW NAMES
# -------------------------
rownames(merged) <- merged$gene_symbol

# eliminar columnas no expresión
expr_36 <- merged[, !(colnames(merged) %in% c(
  "gene_clean", "gene.x", "gene.y",
  "gene_symbol", "logFC", "pvalue", "padj", "regulation"
))]

# -------------------------
# LOG
# -------------------------
expr_log <- log2(expr_36 + 1)

# -------------------------
# METADATA
# -------------------------
meta <- read.csv("metadata_pca.csv")

annotation_col <- data.frame(group = meta$group)
rownames(annotation_col) <- meta$sample

annotation_colors <- list(
  group = c("Active" = "#1f77b4", "Sedentary" = "#d62728")
)

# -------------------------
# HEATMAP
# -------------------------
pheatmap(
  expr_log,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  fontsize_row = 8,
  main = "Heatmap – 36 genes (PC1+PC2 corrected)",
  filename = "heatmap_36_genes.png",
  width = 6,
  height = 8
)