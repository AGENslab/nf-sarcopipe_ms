#!/usr/bin/env Rscript

############################################################
# Figure 5B – Seed prediction summary
#
# Description:
# This script generates a grouped barplot summarizing seed-based
# target prediction results across algorithms.
#
# Inputs:
# - Summary TSV with source, n_miRNAs, n_genes, n_pairs, n_coherent
# - Output PNG file
#
# Outputs:
# - Seed prediction summary barplot PNG
#
# Usage:
# Rscript Figure5B_seed_prediction_summary.R \
#   fig5a_summary.tsv \
#   Figure5B_summary.png
#
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure5B_seed_prediction_summary.R <summary.tsv> <output.png>",
    call. = FALSE
  )
}

input_tsv <- args[1]
output_png <- args[2]

if (!file.exists(input_tsv)) {
  stop("Input TSV not found: ", input_tsv, call. = FALSE)
}

out_dir <- dirname(output_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

df <- read_tsv(input_tsv, show_col_types = FALSE)

required_cols <- c("source", "n_miRNAs", "n_genes", "n_pairs", "n_coherent")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input TSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

long <- pivot_longer(
  df,
  cols = -source,
  names_to = "metric",
  values_to = "value"
)

long$metric <- factor(
  long$metric,
  levels = c("n_miRNAs", "n_genes", "n_pairs", "n_coherent"),
  labels = c("miRNAs evaluated", "Genes with match", "Raw pairs", "Coherent pairs")
)

colors_algorithms <- c(
  "BrumiR" = "#4D4D4D",
  "miRDeep2" = "#1B9E77"
)

p <- ggplot(long, aes(x = metric, y = value, fill = source)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(
    aes(label = value),
    position = position_dodge(width = 0.75),
    vjust = -0.25,
    size = 3.5
  ) +
  scale_fill_manual(values = colors_algorithms) +
  labs(
    title = "Seed-based target prediction summary across algorithms",
    x = NULL,
    y = "Count",
    fill = "Algorithm"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  output_png,
  p,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Seed prediction summary plot written to:", output_png, "\n")