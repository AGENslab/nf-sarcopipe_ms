#!/usr/bin/env Rscript

############################################################
# Figure 4B – Volcano plot corrected by PC1 and PC2
#
# Description:
# This script generates a volcano plot from the gene-wise linear
# model corrected by PC1 and PC2. The top 10 most significant
# genes are labeled using HGNC symbols retrieved with biomaRt.
#
# Inputs:
# - DEG CSV file from Figure4B_linear_model_PC1_PC2.R
# - Output PNG file
#
# Outputs:
# - Volcano plot PNG
#
# Usage:
# Rscript Figure4B_volcano_PC1_PC2.R \
#   DEG_PC1_PC2_corrected.csv \
#   volcano_PC1_PC2_labeled.png
#
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(biomaRt)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure4B_volcano_PC1_PC2.R ",
    "<DEG_PC1_PC2_corrected.csv> <output.png>",
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
df <- read.csv(in_csv)

required_cols <- c("gene", "logFC", "padj")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input CSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -------------------------
# CALCULATE VALUES
# -------------------------
df$neglog10_padj <- -log10(df$padj)

df$status <- "Not significant"
df$status[df$padj < 0.05 & df$logFC > 0] <- "Up in Sedentary"
df$status[df$padj < 0.05 & df$logFC < 0] <- "Up in Active"

# -------------------------
# CONVERT REFSEQ → HGNC
# -------------------------
refseq_clean <- sub("\\..*$", "", df$gene)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(
  attributes = c("refseq_mrna", "hgnc_symbol"),
  filters = "refseq_mrna",
  values = unique(refseq_clean),
  mart = ensembl
)

annot <- annot[annot$hgnc_symbol != "", ]
annot <- annot[!duplicated(annot$refseq_mrna), ]

symbol_map <- annot$hgnc_symbol
names(symbol_map) <- annot$refseq_mrna

gene_symbols <- symbol_map[refseq_clean]
gene_symbols[is.na(gene_symbols)] <- df$gene[is.na(gene_symbols)]

df$gene_symbol <- gene_symbols

# -------------------------
# TOP 10 GENES
# -------------------------
df_sig <- df[df$padj < 0.05, ]
df_sig <- df_sig[order(df_sig$padj), ]

top10 <- df_sig[1:10, ]

# -------------------------
# PLOT
# -------------------------
p <- ggplot(df, aes(logFC, neglog10_padj, color = status)) +
  geom_point(size = 1.2) +
  geom_text_repel(
    data = top10,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = c(
    "Up in Sedentary" = "#d62728",
    "Up in Active" = "#1f77b4",
    "Not significant" = "grey70"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano plot (linear model corrected by PC1 + PC2)",
    x = "log2 Fold Change",
    y = "-log10 FDR"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# -------------------------
# SAVE
# -------------------------
png(
  out_png,
  width = 2100,
  height = 1500,
  res = 300,
  type = "cairo",
  bg = "white"
)

print(p)
dev.off()

cat("Volcano plot written to:", out_png, "\n")
