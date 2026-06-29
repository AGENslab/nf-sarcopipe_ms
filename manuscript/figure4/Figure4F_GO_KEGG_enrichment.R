#!/usr/bin/env Rscript

############################################################
# Figure 4F – GO enrichment by regulation direction
#
# Description:
# This script performs exploratory GO biological process
# enrichment separately for genes upregulated in Active and
# Sedentary individuals among the 36 PC1+PC2-corrected genes.
#
# Inputs:
# - CSV file containing columns: direction, ENTREZID
# - Output prefix
#
# Outputs:
# - <out_prefix>_active_GO_dotplot.png
# - <out_prefix>_sedentary_GO_dotplot.png
#
# Usage:
# Rscript Figure4F_GO_KEGG_enrichment.R \
#   DEG_36_for_networks.csv \
#   GO_enrichment_36_genes
#
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure4F_GO_KEGG_enrichment.R ",
    "<DEG_36_for_networks.csv> <out_prefix>",
    call. = FALSE
  )
}

input_csv <- args[1]
out_prefix <- args[2]

if (!file.exists(input_csv)) {
  stop("Input CSV not found: ", input_csv, call. = FALSE)
}

out_active <- paste0(out_prefix, "_active_GO_dotplot.png")
out_sedentary <- paste0(out_prefix, "_sedentary_GO_dotplot.png")

out_dir <- dirname(out_active)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD DATA
# -------------------------
df <- read.csv(input_csv)

required_cols <- c("direction", "ENTREZID")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input CSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -------------------------
# SPLIT GROUPS
# -------------------------
active <- df %>%
  filter(direction == "Up_in_Active")

sedentary <- df %>%
  filter(direction == "Up_in_Sedentary")

active_genes <- unique(active$ENTREZID[!is.na(active$ENTREZID)])
sedentary_genes <- unique(sedentary$ENTREZID[!is.na(sedentary$ENTREZID)])

# -------------------------
# GO ACTIVE
# -------------------------
ego_active <- enrichGO(
  gene = active_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.1,
  readable = TRUE
)

# -------------------------
# GO SEDENTARY
# -------------------------
ego_sedentary <- enrichGO(
  gene = sedentary_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.1,
  readable = TRUE
)

# -------------------------
# PLOTS
# -------------------------
p_active <- dotplot(
  ego_active,
  showCategory = 6,
  title = "Active (immune signature)"
)

p_sedentary <- dotplot(
  ego_sedentary,
  showCategory = 6,
  title = "Sedentary (metabolic signature)"
)

png(
  filename = out_active,
  width = 2200,
  height = 1600,
  res = 300,
  type = "cairo",
  bg = "white"
)
print(p_active)
dev.off()

png(
  filename = out_sedentary,
  width = 2200,
  height = 1600,
  res = 300,
  type = "cairo",
  bg = "white"
)
print(p_sedentary)
dev.off()

cat("Active GO dotplot written to:", out_active, "\n")
cat("Sedentary GO dotplot written to:", out_sedentary, "\n")