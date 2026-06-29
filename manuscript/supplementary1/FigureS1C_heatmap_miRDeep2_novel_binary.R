#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Binary heatmap of miRDeep2 novel miRNAs (p0.3)
# Displays presence (1) / absence (0) across samples
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_heatmap_mirdeep2_novel_binary.R <counts.tsv> <sample_info.tsv> <out.png>")
}

counts_file <- args[1]
meta_file   <- args[2]
out_plot    <- args[3]

# -------------------------
# LOAD DATA (robusto)
# -------------------------
counts <- read.delim(counts_file, check.names = FALSE)

meta <- read_tsv(meta_file, show_col_types = FALSE) %>%
  as.data.frame()

# -------------------------
# FORMAT COUNTS
# -------------------------
rownames(counts) <- counts[,1]
counts <- counts[,-1, drop = FALSE]

# asegurar numérico
counts <- as.data.frame(lapply(counts, function(x) as.numeric(as.character(x))))
counts[is.na(counts)] <- 0

# -------------------------
# BINARIZE
# -------------------------
counts_bin <- counts
counts_bin[counts_bin > 0] <- 1

# -------------------------
# MATCH METADATA
# -------------------------
meta <- meta %>%
  filter(sample %in% colnames(counts_bin))

counts_bin <- counts_bin[, meta$sample, drop = FALSE]

meta$condition <- factor(meta$condition, levels = c("active", "sedentary"))
rownames(meta) <- meta$sample

# -------------------------
# OPTIONAL FILTER
# (solo los que aparecen en >=2 muestras)
# -------------------------
counts_bin <- counts_bin[rowSums(counts_bin) >= 2, ]

# -------------------------
# COLORS
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
# HEATMAP
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
