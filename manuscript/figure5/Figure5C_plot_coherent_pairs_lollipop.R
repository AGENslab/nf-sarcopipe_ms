#!/usr/bin/env Rscript

############################################################
# Figure 5C – Coherent miRNA–mRNA pairs lollipop plot
#
# Description:
# This script generates a lollipop plot of coherent miRNA–mRNA
# target pairs identified by seed-based prediction. It keeps all
# BrumiR coherent pairs and the top 15 miRDeep2 pairs ranked by
# adjusted p-value.
#
# Inputs:
# - TSV file with coherent pairs for lollipop plotting
# - Output PNG file
#
# Outputs:
# - Coherent pairs lollipop plot PNG
#
# Usage:
# Rscript Figure5C_plot_coherent_pairs_lollipop.R \
#   coherent_pairs_for_lollipop.tsv \
#   Figure5C_coherent_lollipop.png
#
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure5C_plot_coherent_pairs_lollipop.R ",
    "<coherent_pairs_for_lollipop.tsv> <output.png>",
    call. = FALSE
  )
}

infile <- args[1]
outfile <- args[2]

if (!file.exists(infile)) {
  stop("Input TSV not found: ", infile, call. = FALSE)
}

out_dir <- dirname(outfile)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

df <- read_tsv(infile, show_col_types = FALSE)

required_cols <- c(
  "source",
  "padj",
  "renamed_miRNA",
  "gene_symbol",
  "direction",
  "abs_log2FC"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input TSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -------------------------
# SELECT TOP PAIRS
# -------------------------
md_top <- df %>%
  filter(source == "miRDeep2") %>%
  arrange(padj) %>%
  slice(1:15)

br_all <- df %>%
  filter(source == "BrumiR")

df_plot <- bind_rows(br_all, md_top)

# -------------------------
# PREPARE VARIABLES
# -------------------------
df_plot <- df_plot %>%
  mutate(
    pair_label = paste0(renamed_miRNA, " \u2192 ", gene_symbol),
    signed_value = ifelse(direction == "Up_in_Active", -abs_log2FC, abs_log2FC)
  )

df_plot <- df_plot %>%
  group_by(source) %>%
  arrange(signed_value, .by_group = TRUE) %>%
  mutate(pair_label = factor(pair_label, levels = unique(pair_label))) %>%
  ungroup()

# -------------------------
# PLOT
# -------------------------
p <- ggplot(df_plot, aes(x = signed_value, y = pair_label, color = direction)) +
  geom_segment(
    aes(x = 0, xend = signed_value, y = pair_label, yend = pair_label),
    linewidth = 0.7
  ) +
  geom_point(size = 3) +
  facet_wrap(~source, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(
    "Up_in_Active" = "#3B6FB6",
    "Up_in_Sedentary" = "#C94C4C"
  )) +
  labs(
    title = "Coherent miRNA–mRNA target pairs identified by seed-based prediction",
    x = "miRNA |log2 fold change| (signed by direction)",
    y = NULL,
    color = "miRNA direction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )

ggsave(outfile, p, width = 10, height = 9, dpi = 300)

cat("Coherent pairs lollipop plot written to:", outfile, "\n")
