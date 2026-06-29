#!/usr/bin/env Rscript

################################################################################
# Script: gene_reduction_barplot.R
#
# Description:
# This script generates a barplot summarizing the progressive reduction in the
# number of candidate genes throughout the mRNA-seq analysis workflow.
#
# Objective:
# To illustrate how the initial set of expressed genes was reduced after
# differential expression analysis before and after correction for PC1 and PC2.
#
# Methodological note:
# A log10 scale is used on the y-axis to highlight the final reduced set of
# candidate genes while preserving the large differences in magnitude across
# analysis steps.
#
# Output:
# - gene_reduction_barplot.png
#
# Author: Natalia Poblete
################################################################################

library(ggplot2)

# -------------------------
# DATA
# -------------------------
df <- data.frame(
  step = c(
    "Initial genes\n(post-filtering)",
    "Significant genes\n(before PC1+PC2 correction)",
    "Significant genes\n(after PC1+PC2 correction)"
  ),
  genes = c(44234, 24081, 36)
)

df$step <- factor(df$step, levels = df$step)

# custom label positions so the 36 is more visible
df$label_y <- c(50000, 27000, 55)
df$label_size <- c(5, 5, 6)

# -------------------------
# BARPLOT
# -------------------------
p <- ggplot(df, aes(x = step, y = genes, fill = step)) +
  geom_bar(stat = "identity", width = 0.65) +
  
  geom_text(
    aes(y = label_y, label = format(genes, big.mark = ",")),
    size = df$label_size,
    fontface = "bold"
  ) +
  
  scale_fill_manual(values = c("#9ecae1", "#6baed6", "#de2d26")) +
  
  scale_y_log10(
    breaks = c(10, 100, 1000, 10000, 50000),
    labels = c("10", "100", "1,000", "10,000", "50,000")
  ) +
  
  theme_classic(base_size = 14) +
  labs(
    title = "Progressive reduction of candidate genes after confounding correction",
    subtitle = "Log-scale representation highlighting the final robust gene set",
    x = "",
    y = "Number of genes (log10 scale)"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# -------------------------
# SAVE
# -------------------------
png("gene_reduction_barplot.png", width = 2200, height = 1400, res = 300, type = "cairo", bg = "white")
print(p)
dev.off()
