#!/usr/bin/env Rscript

###############################################################
# Script: supplementary_pc_sensitivity_plot.R
# Author: Natalia Poblete
# Description:
# This script generates a supplementary barplot illustrating
# the effect of including additional principal components (PCs)
# in the linear model used for differential expression analysis.
#
# The figure shows how the number of significant genes (FDR < 0.05)
# decreases as more PCs are included, highlighting over-adjustment
# and loss of biological signal.
#
# Input:
#   Hardcoded summary of results from sensitivity analysis:
#   - PC1 + PC2       -> 36 genes
#   - PC1 + PC2 + PC3 -> 9 genes
#   - PC1 + PC2 + PC3 + PC4 -> 0 genes
#
# Output:
#   results/figures/supplementary_pc_sensitivity.png
###############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

# -------------------------
# DATA
# -------------------------
df <- data.frame(
  model = c("PC1+PC2", "PC1+PC2+PC3", "PC1+PC2+PC3+PC4"),
  genes = c(36, 9, 0),
  stringsAsFactors = FALSE
)

# -------------------------
# PLOT
# -------------------------
p <- ggplot(df, aes(x = model, y = genes)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(aes(label = genes), vjust = -0.4, size = 5, fontface = "bold") +
  expand_limits(y = max(df$genes) * 1.2) +   # 🔥 ESTE ES EL FIX
  theme_classic(base_size = 14) +
  labs(
    title = "Effect of including additional principal components",
    x = "Linear model",
    y = "Number of significant genes (FDR < 0.05)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------------------------
# SAVE FIGURE (HPC-safe)
# -------------------------
output_path <- "results/figures/supplementary_pc_sensitivity.png"

png(
  filename = output_path,
  width = 1800,
  height = 1200,
  res = 300,
  type = "cairo",
  bg = "white"
)
print(p)
dev.off()

cat("Figure saved to:", output_path, "\n")
