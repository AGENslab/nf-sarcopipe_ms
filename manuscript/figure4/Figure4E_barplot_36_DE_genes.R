#!/usr/bin/env Rscript

################################################################################
# Script: barplot_36_genes_directional.R
#
# Description:
# This script generates a diverging barplot for the 36 differentially expressed
# genes identified after correction for PC1 and PC2.
#
# Objective:
# To visualize the direction and magnitude of differential expression for each
# gene, highlighting genes upregulated in Active individuals on one side and
# genes upregulated in Sedentary individuals on the opposite side.
#
# Methodological note:
# Positive values indicate genes upregulated in Sedentary individuals, whereas
# negative values indicate genes upregulated in Active individuals. Gene labels
# are shown using HGNC symbols for improved biological interpretability.
#
# Inputs:
# - Supplementary_Table_DEG_PC1_PC2_HGNC.csv
#
# Output:
# - barplot_36_genes_directional.png
#
# Author: Natalia Poblete
################################################################################

library(ggplot2)

# -------------------------
# LOAD DATA
# -------------------------
df <- read.csv("Supplementary_Table_DEG_PC1_PC2_HGNC.csv", stringsAsFactors = FALSE)

# -------------------------
# PREPARE DATA
# -------------------------
# Keep only significant genes
df <- df[order(df$logFC), ]

# Factor for plotting order
df$gene_symbol <- factor(df$gene_symbol, levels = df$gene_symbol)

# -------------------------
# BARPLOT
# -------------------------
p <- ggplot(df, aes(x = gene_symbol, y = logFC, fill = regulation)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Up in Active" = "#1f77b4",
    "Up in Sedentary" = "#d62728"
  )) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  theme_classic(base_size = 14) +
  labs(
    title = "Direction and magnitude of differential expression in the 36 corrected genes",
    x = "",
    y = "Log2 fold change",
    fill = "Regulation"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 8)
  )

# -------------------------
# SAVE
# -------------------------
png("barplot_36_genes_directional.png", width = 2400, height = 2200, res = 300, type = "cairo", bg = "white")
print(p)
dev.off()
