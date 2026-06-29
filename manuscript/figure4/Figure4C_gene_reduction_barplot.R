#!/usr/bin/env Rscript

############################################################
# Figure 4C – Gene reduction barplot
#
# Description:
# This script generates a barplot summarizing the progressive
# reduction in the number of candidate genes throughout the
# mRNA-seq analysis workflow.
#
# Inputs:
# - Output PNG file
#
# Outputs:
# - gene_reduction_barplot.png or user-defined output PNG
#
# Usage:
# Rscript Figure4C_gene_reduction_barplot.R \
#   gene_reduction_barplot.png
#
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(
    "Usage: Rscript Figure4C_gene_reduction_barplot.R <output.png>",
    call. = FALSE
  )
}

out_png <- args[1]

out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

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
png(
  out_png,
  width = 2200,
  height = 1400,
  res = 300,
  type = "cairo",
  bg = "white"
)

print(p)
dev.off()

cat("Gene reduction barplot written to:", out_png, "\n")