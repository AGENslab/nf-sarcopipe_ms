#!/usr/bin/env Rscript

############################################################
# Figure S4B – Prioritized predicted target genes
#
# Description:
# Generates a publication-ready plot of prioritized predicted
# target genes based on coherent miRNA support, number of
# seed-matched sites, and miRNA effect size.
#
# Inputs:
# 1) prioritized_predicted_target_genes_top15.tsv
#
# Outputs:
# - FigureS4B_prioritized_predicted_targets.png
#
# Usage:
# Rscript FigureS4B_plot_prioritized_targets.R \
#   <prioritized_predicted_target_genes_top15.tsv> \
#   <FigureS4B_prioritized_predicted_targets.png>
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript FigureS4B_plot_prioritized_targets.R ",
    "<prioritized_predicted_target_genes_top15.tsv> ",
    "<FigureS4B_prioritized_predicted_targets.png>"
  )
}

infile <- args[1]
outfile <- args[2]

if (!file.exists(infile)) {
  stop("Input prioritized target table not found: ", infile)
}

out_dir <- dirname(outfile)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

df <- read_tsv(infile, show_col_types = FALSE)

required_cols <- c(
  "gene_symbol",
  "category",
  "priority_score",
  "total_n_sites"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in input table: ",
    paste(missing_cols, collapse = ", ")
  )
}

df <- df %>%
  mutate(
    gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)),
    category = factor(category, levels = c("shared", "BrumiR_only", "miRDeep2_only"))
  )

p <- ggplot(
  df,
  aes(
    x = priority_score,
    y = gene_symbol,
    color = category,
    size = total_n_sites
  )
) +
  geom_point(alpha = 0.95) +
  scale_color_manual(values = c(
    "shared" = "#7A4EAB",
    "BrumiR_only" = "#4D4D4D",
    "miRDeep2_only" = "#1B9E77"
  )) +
  labs(
    title = "Prioritized predicted target genes from coherent miRNA–mRNA pairs",
    x = "Priority score",
    y = NULL,
    color = "Support category",
    size = "Total\nsites"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold")
  )

ggsave(
  outfile,
  p,
  width = 9,
  height = 6,
  dpi = 300
)

cat("Prioritized target plot written to:", outfile, "\n")
