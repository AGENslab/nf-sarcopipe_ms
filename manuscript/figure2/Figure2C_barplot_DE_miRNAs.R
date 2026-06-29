#!/usr/bin/env Rscript

############################################################
# Figure 2C – Differentially expressed miRNAs
#
# Description:
# Summarizes differentially expressed miRNAs from DESeq2 results
# generated for BrumiR and miRDeep2, and produces a stacked bar
# plot showing the number of miRNAs upregulated in Active and
# Sedentary groups.
#
# Inputs:
# 1. BrumiR DESeq2 results CSV
# 2. miRDeep2 DESeq2 results CSV
#
# Output:
# 1. PNG bar plot
#
# Usage:
# Rscript Figure2C_barplot_DE_miRNAs.R \
#   brumir_deseq2.csv \
#   mirdeep2_deseq2.csv \
#   de_mirna_barplot.png
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: Rscript Figure2C_barplot_DE_miRNAs.R ",
    "<brumir_deseq2.csv> <mirdeep2_deseq2.csv> <out_plot.png>"
  )
}

brumir_file <- args[1]
mirdeep_file <- args[2]
out_plot <- args[3]

if (!file.exists(brumir_file)) {
  stop("BrumiR input file not found: ", brumir_file)
}

if (!file.exists(mirdeep_file)) {
  stop("miRDeep2 input file not found: ", mirdeep_file)
}

out_dir <- dirname(out_plot)
if (!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

brumir <- read_csv(brumir_file, show_col_types = FALSE)
mirdeep <- read_csv(mirdeep_file, show_col_types = FALSE)

summarize_de <- function(df, label) {
  df %>%
    filter(!is.na(padj), padj < 0.05) %>%
    mutate(direction = case_when(
      log2FoldChange > 0 ~ "Up in Active",
      log2FoldChange < 0 ~ "Up in Sedentary"
    )) %>%
    count(direction) %>%
    mutate(method = label)
}

brumir_sum <- summarize_de(brumir, "BrumiR")
mirdeep_sum <- summarize_de(mirdeep, "miRDeep2")

df_plot <- bind_rows(brumir_sum, mirdeep_sum)

p <- ggplot(df_plot, aes(x = method, y = n, fill = direction)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white")
  ) +
  scale_fill_manual(values = c(
    "Up in Active" = "blue",
    "Up in Sedentary" = "red"
  )) +
  labs(
    title = "Differentially expressed miRNAs",
    x = "Method",
    y = "Count",
    fill = "miRNA regulation"
  )

ggsave(out_plot, p, width = 6, height = 5, bg = "white")
