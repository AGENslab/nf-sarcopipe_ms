#!/usr/bin/env Rscript

############################################################
# Figure 4D – Heatmap of 36 PC1+PC2-corrected genes
#
# Description:
# This script generates a heatmap of 36 differentially expressed
# genes identified after fitting a linear model corrected for PC1
# and PC2. Gene names are converted from RefSeq IDs to HGNC symbols
# for improved biological interpretability.
#
# Inputs:
# - Normalized expression matrix CSV
# - DEG annotation CSV with HGNC symbols
# - Metadata CSV with sample and group columns
# - Output PNG file
#
# Outputs:
# - Heatmap PNG
#
# Usage:
# Rscript Figure4D_heatmap_36_genes.R \
#   matriz_normalizada_para_modelos_nb.csv \
#   Supplementary_Table_DEG_PC1_PC2_HGNC.csv \
#   metadata_pca.csv \
#   heatmap_36_genes.png
#
############################################################

suppressPackageStartupMessages({
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript Figure4D_heatmap_36_genes.R ",
    "<expression_matrix.csv> <deg_annotation.csv> <metadata.csv> <output.png>",
    call. = FALSE
  )
}

expr_file <- args[1]
annot_file <- args[2]
metadata_file <- args[3]
out_png <- args[4]

input_files <- c(expr_file, annot_file, metadata_file)
missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Input file(s) not found: ",
    paste(missing_files, collapse = ", "),
    call. = FALSE
  )
}

out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD EXPRESSION MATRIX
# -------------------------
expr <- read.csv(expr_file, check.names = FALSE)

rownames(expr) <- make.unique(expr[, 1])
expr <- expr[, -1]

expr_df <- as.data.frame(expr)
expr_df$gene <- rownames(expr_df)

expr_df$gene_clean <- sub("\\..*$", "", expr_df$gene)

# -------------------------
# LOAD ANNOTATION HGNC
# -------------------------
annot <- read.csv(annot_file)

required_annot_cols <- c("gene", "gene_symbol")
missing_annot_cols <- setdiff(required_annot_cols, colnames(annot))

if (length(missing_annot_cols) > 0) {
  stop(
    "Annotation file is missing required columns: ",
    paste(missing_annot_cols, collapse = ", "),
    call. = FALSE
  )
}

annot$gene_clean <- sub("\\..*$", "", annot$gene)

# -------------------------
# MERGE EXPRESSION AND ANNOTATION
# -------------------------
merged <- merge(annot, expr_df, by = "gene_clean")

# -------------------------
# USE HGNC SYMBOLS AS ROW NAMES
# -------------------------
rownames(merged) <- merged$gene_symbol

expr_36 <- merged[, !(colnames(merged) %in% c(
  "gene_clean",
  "gene.x",
  "gene.y",
  "gene_symbol",
  "logFC",
  "pvalue",
  "padj",
  "regulation"
))]

# -------------------------
# LOG TRANSFORMATION
# -------------------------
expr_log <- log2(expr_36 + 1)

# -------------------------
# METADATA
# -------------------------
meta <- read.csv(metadata_file)

required_meta_cols <- c("sample", "group")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_meta_cols, collapse = ", "),
    call. = FALSE
  )
}

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
  filename = out_png,
  width = 6,
  height = 8
)

cat("Heatmap written to:", out_png, "\n")