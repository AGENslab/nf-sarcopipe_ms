#!/usr/bin/env Rscript

############################################################
# Figure 4E – Barplot of 36 differentially expressed genes
#
# Description:
# This script generates a diverging barplot for the 36
# differentially expressed genes identified after correction for
# PC1 and PC2. Positive values indicate genes upregulated in
# Sedentary individuals, whereas negative values indicate genes
# upregulated in Active individuals.
#
# Inputs:
# - DEG annotation CSV with HGNC symbols, logFC, and regulation
# - Output PNG file
#
# Outputs:
# - barplot_36_genes_directional.png or user-defined output PNG
#
# Usage:
# Rscript Figure4E_barplot_36_DE_genes.R \
#   Supplementary_Table_DEG_PC1_PC2_HGNC.csv \
#   barplot_36_genes_directional.png
#
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure4E_barplot_36_DE_genes.R ",
    "<Supplementary_Table_DEG_PC1_PC2_HGNC.csv> <output.png>",
    call. = FALSE
  )
}

in_csv <- args[1]
out_png <- args[2]

if (!file.exists(in_csv)) {
  stop("Input CSV not found: ", in_csv, call. = FALSE)
}

out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD DATA
# -------------------------
df <- read.csv(in_csv, stringsAsFactors = FALSE)

required_cols <- c("gene_symbol", "logFC", "regulation")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input CSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -------------------------
# PREPARE DATA
# -------------------------
df <- df[order(df$logFC), ]

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
png(
  out_png,
  width = 2400,
  height = 2200,
  res = 300,
  type = "cairo",
  bg = "white"
)

print(p)
dev.off()

cat("Barplot written to:", out_png, "\n")
