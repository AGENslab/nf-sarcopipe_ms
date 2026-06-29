#!/usr/bin/env Rscript

############################################################
# Figure S3D – PC sensitivity barplot
#
# Description:
# Generates a supplementary barplot illustrating the effect of
# including additional principal components (PCs) in the linear
# model used for differential expression analysis.
#
# The figure shows how the number of significant genes
# (FDR < 0.05) decreases as more PCs are included, highlighting
# over-adjustment and loss of biological signal.
#
# Inputs:
# Hardcoded summary of sensitivity analysis results:
# - PC1 + PC2            -> 36 genes
# - PC1 + PC2 + PC3      -> 9 genes
# - PC1 + PC2 + PC3 + PC4 -> 0 genes
#
# Outputs:
# - supplementary_pc_sensitivity.png
#
# Usage:
# Rscript FigureS3D_sensitivity_barplot.R \
#   <supplementary_pc_sensitivity.png>
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript FigureS3D_sensitivity_barplot.R ",
    "<supplementary_pc_sensitivity.png>"
  )
}

output_path <- args[1]

out_dir <- dirname(output_path)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# Data
# -------------------------
df <- data.frame(
  model = c("PC1+PC2", "PC1+PC2+PC3", "PC1+PC2+PC3+PC4"),
  genes = c(36, 9, 0),
  stringsAsFactors = FALSE
)

# -------------------------
# Plot
# -------------------------
p <- ggplot(df, aes(x = model, y = genes)) +
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(aes(label = genes), vjust = -0.4, size = 5, fontface = "bold") +
  expand_limits(y = max(df$genes) * 1.2) +
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
# Save figure
# -------------------------
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
