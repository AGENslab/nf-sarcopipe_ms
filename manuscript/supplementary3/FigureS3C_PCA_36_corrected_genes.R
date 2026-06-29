#!/usr/bin/env Rscript

################################################################################
# Script: PCA_36_genes.R
#
# Description:
# This script performs Principal Component Analysis (PCA) using a subset of
# 36 differentially expressed genes identified from a linear model corrected
# for PC1 and PC2.
#
# Objective:
# To visualize the separation between Active and Sedentary groups based on
# genes that remain significant after correcting for tissue effects (PC1)
# and intra-group variability (PC2).
#
# Methodological note:
# This PCA represents a post-hoc exploratory visualization, as the genes were
# pre-selected based on statistical significance. Therefore, it should not be
# interpreted as independent evidence of group separation, but rather as a
# visual confirmation of robust biological signal.
#
# Inputs:
# - matriz_normalizada_para_modelos_nb.csv : gene expression matrix (genes x samples)
# - genes_36.txt : list of significant genes (FDR < 0.05)
# - metadata_pca.csv : sample metadata (must include "group")
#
# Output:
# - PCA_36_genes.png : PCA plot (PC1 vs PC2)
#
# Author: Natalia Poblete
################################################################################

library(ggplot2)
library(dplyr)

# -------------------------
# LOAD EXPRESSION MATRIX
# -------------------------
expr <- read.csv("matriz_normalizada_para_modelos_nb.csv", check.names = FALSE)
rownames(expr) <- make.unique(expr[,1])
expr <- expr[,-1]

# -------------------------
# LOAD SELECTED GENES
# -------------------------
genes <- read.table("genes_36.txt", stringsAsFactors = FALSE)$V1
expr_36 <- expr[rownames(expr) %in% genes, , drop = FALSE]

# -------------------------
# LOG TRANSFORMATION
# -------------------------
expr_log <- log2(expr_36 + 1)

# -------------------------
# PCA
# -------------------------
pca <- prcomp(t(expr_log), center = TRUE, scale. = TRUE)
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

# -------------------------
# METADATA
# -------------------------
meta <- read.csv("metadata_pca.csv", stringsAsFactors = FALSE)

pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, meta, by = "sample")
pca_df$group <- factor(pca_df$group, levels = c("Sedentary", "Active"))

# -------------------------
# PLOT
# -------------------------
p <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  stat_ellipse(aes(group = group), linewidth = 1, level = 0.95) +
  geom_point(size = 4, alpha = 1) +
  scale_color_manual(values = c("Sedentary" = "#d62728", "Active" = "#1f77b4")) +
  theme_classic(base_size = 14) +
  labs(
    title = "PCA – 36 genes (post-correction)",
    x = paste0("PC1 (", round(var_exp[1],1), "% variance)"),
    y = paste0("PC2 (", round(var_exp[2],1), "% variance)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# -------------------------
# SAVE (HPC-safe)
# -------------------------
png("PCA_36_genes.png", width = 1800, height = 1500, res = 300, type = "cairo", bg = "white")
print(p)
dev.off()
