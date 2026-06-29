#!/usr/bin/env Rscript

############################################################
# Figure S3C – PCA of 36 corrected genes
#
# Description:
# Performs Principal Component Analysis (PCA) using a subset of
# 36 differentially expressed genes identified from a linear model
# corrected for PC1 and PC2.
#
# This PCA is a post-hoc exploratory visualization because the genes
# were pre-selected based on statistical significance. It should be
# interpreted as visual confirmation of robust biological signal, not
# independent evidence of group separation.
#
# Inputs:
# 1) Expression matrix CSV:
#    - genes x samples
#    - first column contains gene identifiers
# 2) Gene list TXT:
#    - one selected gene per line
# 3) Metadata CSV:
#    - must include columns:
#      sample
#      group
#
# Outputs:
# - PCA_36_genes.png
#
# Usage:
# Rscript FigureS3C_PCA_36_corrected_genes.R \
#   <matriz_normalizada_para_modelos_nb.csv> \
#   <genes_36.txt> \
#   <metadata_pca.csv> \
#   <PCA_36_genes.png>
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript FigureS3C_PCA_36_corrected_genes.R ",
    "<expression_matrix.csv> <genes_36.txt> <metadata_pca.csv> <out_png>"
  )
}

expr_file <- args[1]
genes_file <- args[2]
metadata_file <- args[3]
out_png <- args[4]

if (!file.exists(expr_file)) {
  stop("Expression matrix file not found: ", expr_file)
}

if (!file.exists(genes_file)) {
  stop("Gene list file not found: ", genes_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# Load expression matrix
# -------------------------
expr <- read.csv(expr_file, check.names = FALSE)

if (ncol(expr) < 2) {
  stop("Expression matrix must contain one gene column and at least one sample column.")
}

rownames(expr) <- make.unique(expr[, 1])
expr <- expr[, -1, drop = FALSE]

# -------------------------
# Load selected genes
# -------------------------
genes <- read.table(genes_file, stringsAsFactors = FALSE)$V1
expr_36 <- expr[rownames(expr) %in% genes, , drop = FALSE]

if (nrow(expr_36) == 0) {
  stop("No selected genes were found in the expression matrix.")
}

# -------------------------
# Log transformation
# -------------------------
expr_log <- log2(expr_36 + 1)

# -------------------------
# PCA
# -------------------------
pca <- prcomp(t(expr_log), center = TRUE, scale. = TRUE)
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# -------------------------
# Metadata
# -------------------------
meta <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_meta_cols <- c("sample", "group")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing required columns in metadata: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = "sample")

if (any(is.na(pca_df$group))) {
  stop("Some PCA samples do not have matching group information in metadata.")
}

pca_df$group <- factor(pca_df$group, levels = c("Sedentary", "Active"))

# -------------------------
# Plot
# -------------------------
p <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  stat_ellipse(aes(group = group), linewidth = 1, level = 0.95) +
  geom_point(size = 4, alpha = 1) +
  scale_color_manual(values = c("Sedentary" = "#d62728", "Active" = "#1f77b4")) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA – 36 genes (post-correction)",
    x = paste0("PC1 (", round(var_exp[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "% variance)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# -------------------------
# Save
# -------------------------
png(
  out_png,
  width = 1800,
  height = 1500,
  res = 300,
  type = "cairo",
  bg = "white"
)

print(p)
dev.off()

cat("PCA plot written to:", out_png, "\n")
