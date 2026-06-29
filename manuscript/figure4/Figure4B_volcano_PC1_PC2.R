#!/usr/bin/env Rscript

################################################################################
# Script: volcano_PC1_PC2.R
#
# Description:
# Volcano plot with top 10 most significant genes labeled using HGNC symbols.
################################################################################

library(ggplot2)
library(ggrepel)
library(biomaRt)

# -------------------------
# LOAD DATA
# -------------------------
df <- read.csv("DEG_PC1_PC2_corrected.csv")

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
  filters    = "refseq_mrna",
  values     = unique(refseq_clean),
  mart       = ensembl
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
png("volcano_PC1_PC2_labeled.png", width = 2100, height = 1500, res = 300, type = "cairo", bg = "white")
print(p)
dev.off()
