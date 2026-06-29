#!/usr/bin/env Rscript

############################################################
# Figure S1C – Binary heatmap of miRDeep2 novel miRNAs
#
# Description:
# Generates a binary presence/absence heatmap for miRDeep2 novel
# miRNAs detected at p0.3 across samples.
#
# Inputs:
# 1) counts.tsv:
#    - first column: miRNA identifier
#    - remaining columns: sample counts
# 2) sample_info.tsv:
#    - sample
#    - condition
#
# Outputs:
# - PNG heatmap showing presence (1) / absence (0)
#
# Usage:
# Rscript FigureS1C_heatmap_miRDeep2_novel_binary.R \
#   <counts.tsv> <sample_info.tsv> <out.png>
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript FigureS1C_heatmap_miRDeep2_novel_binary.R ",
    "<counts.tsv> <sample_info.tsv> <out.png>"
  )
}

counts_file <- args[1]
meta_file <- args[2]
out_plot <- args[3]

if (!file.exists(counts_file)) {
  stop("Input counts file not found: ", counts_file)
}

if (!file.exists(meta_file)) {
  stop("Input sample metadata file not found: ", meta_file)
}

out_dir <- dirname(out_plot)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -------------------------
# Load data
# -------------------------
counts <- read.delim(counts_file, check.names = FALSE)

meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  as.data.frame()

required_meta_cols <- c("sample", "condition")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing required columns in sample metadata: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

if (ncol(counts) < 2) {
  stop("Counts file must contain one ID column and at least one sample column.")
}

# -------------------------
# Format counts
# -------------------------
rownames(counts) <- counts[, 1]
counts <- counts[, -1, drop = FALSE]

counts <- as.data.frame(
  lapply(counts, function(x) as.numeric(as.character(x)))
)
counts[is.na(counts)] <- 0

# -------------------------
# Binarize
# -------------------------
counts_bin <- counts
counts_bin[counts_bin > 0] <- 1

# -------------------------
# Match metadata
# -------------------------
meta <- meta %>%
  filter(sample %in% colnames(counts_bin))

if (nrow(meta) == 0) {
  stop("No metadata samples match count matrix columns.")
}

counts_bin <- counts_bin[, meta$sample, drop = FALSE]

meta$condition <- factor(meta$condition, levels = c("active", "sedentary"))
rownames(meta) <- meta$sample

# -------------------------
# Optional filter
# Keeps original logic: retain miRNAs present in >=2 samples
# -------------------------
counts_bin <- counts_bin[rowSums(counts_bin) >= 2, ]

if (nrow(counts_bin) == 0) {
  stop("No miRNAs remain after filtering for presence in >=2 samples.")
}

# -------------------------
# Colors
# -------------------------
annotation_col <- data.frame(
  condition = meta$condition
)
rownames(annotation_col) <- rownames(meta)

ann_colors <- list(
  condition = c(
    active = "blue",
    sedentary = "red"
  )
)

# -------------------------
# Heatmap
# -------------------------
png(
  filename = out_plot,
  width = 2200,
  height = 2600,
  res = 300,
  bg = "white",
  type = "cairo"
)

pheatmap(
  as.matrix(counts_bin),
  color = c("white", "black"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 10,
  main = "Binary presence of miRDeep2 novel miRNAs (p0.3)"
)

dev.off()

cat("Binary heatmap written to:", out_plot, "\n")