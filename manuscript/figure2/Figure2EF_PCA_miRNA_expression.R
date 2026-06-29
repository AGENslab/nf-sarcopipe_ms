#!/usr/bin/env Rscript

############################################################
# Figure 2E/F – PCA of miRNA expression
#
# Description:
# This script performs PCA analysis of miRNA expression using
# variance-stabilized counts. It visualizes the global expression
# structure across samples and adds 95% confidence ellipses by
# condition.
#
# Inputs:
# - Counts matrix TSV file
# - Sample metadata TSV file with columns: sample, condition
# - Output prefix
#
# Outputs:
# - <out_prefix>_pca_mirna_expression.png
#
# Usage:
# Rscript Figure2EF_PCA_miRNA_expression.R \
#   <counts.tsv> \
#   <sample_info.tsv> \
#   <out_prefix>
#
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript Figure2EF_PCA_miRNA_expression.R ",
    "<counts.tsv> <sample_info.tsv> <out_prefix>",
    call. = FALSE
  )
}

counts_file <- args[1]
metadata_file <- args[2]
out_prefix <- args[3]

if (!file.exists(counts_file)) {
  stop("Input counts file not found: ", counts_file, call. = FALSE)
}

if (!file.exists(metadata_file)) {
  stop("Input metadata file not found: ", metadata_file, call. = FALSE)
}

out_plot <- paste0(out_prefix, "_pca_mirna_expression.png")

out_dir <- dirname(out_plot)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD DATA
# -------------------------
counts <- read_tsv(counts_file, show_col_types = FALSE)
meta <- read_tsv(metadata_file, show_col_types = FALSE)

# -------------------------
# FORMAT COUNTS
# -------------------------
counts_df <- as.data.frame(counts)
rownames(counts_df) <- counts_df[, 1]
counts_df <- counts_df[, -1, drop = FALSE]

counts_mat <- as.matrix(counts_df)
mode(counts_mat) <- "numeric"

# -------------------------
# FORMAT METADATA
# -------------------------
meta_df <- as.data.frame(meta)

required_cols <- c("sample", "condition")
if (!all(required_cols %in% colnames(meta_df))) {
  stop(
    "Metadata file must contain columns: sample and condition",
    call. = FALSE
  )
}

rownames(meta_df) <- meta_df$sample
meta_df$condition <- factor(meta_df$condition, levels = c("active", "sedentary"))

# -------------------------
# MATCH SAMPLES
# -------------------------
common_samples <- intersect(colnames(counts_mat), meta_df$sample)

if (length(common_samples) < 2) {
  stop(
    "Fewer than 2 shared samples between counts matrix and metadata",
    call. = FALSE
  )
}

counts_mat <- counts_mat[, common_samples, drop = FALSE]
meta_df <- meta_df[common_samples, , drop = FALSE]

if (!identical(colnames(counts_mat), rownames(meta_df))) {
  stop(
    "Sample order mismatch between counts matrix and metadata",
    call. = FALSE
  )
}

# -------------------------
# BUILD DESEQ2 OBJECT
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),
  colData = meta_df,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds)) > 10, ]

if (nrow(dds) < 2) {
  stop("Too few miRNAs remaining after filtering", call. = FALSE)
}

# -------------------------
# TRANSFORMATION
# -------------------------
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

# -------------------------
# PCA
# -------------------------
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"), 1)

# -------------------------
# AUTO TITLE BASED ON INPUT
# -------------------------
method_label <- if (grepl("mirdeep", tolower(counts_file))) {
  "miRDeep2"
} else if (grepl("brumir", tolower(counts_file))) {
  "BrumiR"
} else {
  "miRNA"
}

# -------------------------
# PLOT
# -------------------------
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 5, alpha = 0.9) +
  stat_ellipse(
    aes(group = condition),
    type = "norm",
    level = 0.95,
    linewidth = 1.2,
    linetype = "solid"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c(
    "active" = "blue",
    "sedentary" = "red"
  )) +
  labs(
    title = paste0("PCA of ", method_label, " miRNA expression"),
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Condition"
  )

# -------------------------
# SAVE
# -------------------------
png(
  filename = out_plot,
  width = 1800,
  height = 1500,
  res = 300,
  bg = "white",
  type = "cairo"
)

print(p)
dev.off()

cat("PCA plot written to:", out_plot, "\n")
