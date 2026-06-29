#!/usr/bin/env Rscript

############################################################
# Figure 5E – Retrieve 3'UTR sequences for sarcopenia gene set
#
# Description:
# This script downloads annotated 3'UTR sequences for a curated
# sarcopenia/exercise-related gene set using Ensembl BioMart.
# For genes with multiple 3'UTR entries, the longest available
# sequence is retained.
#
# Inputs:
# - Text file containing gene symbols, one per line
# - Output TSV file
#
# Outputs:
# - TSV file with columns: gene_symbol, utr_3
#
# Usage:
# Rscript Figure5E_get_sarcopenia_3UTR.R \
#   input_sarcopenia_gene_set.txt \
#   sarcopenia_gene_set_3UTR.tsv
#
############################################################

suppressPackageStartupMessages({
  library(biomaRt)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure5E_get_sarcopenia_3UTR.R ",
    "<input_sarcopenia_gene_set.txt> <sarcopenia_gene_set_3UTR.tsv>",
    call. = FALSE
  )
}

genes_file <- args[1]
out_file <- args[2]

if (!file.exists(genes_file)) {
  stop("Input gene set file not found: ", genes_file, call. = FALSE)
}

out_dir <- dirname(out_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD GENES
# -------------------------
genes <- read_lines(genes_file)
genes <- unique(genes[genes != ""])

# -------------------------
# QUERY BIOMART
# -------------------------
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)

utr_df <- getSequence(
  id = genes,
  type = "hgnc_symbol",
  seqType = "3utr",
  mart = mart
)

colnames(utr_df) <- tolower(colnames(utr_df))

symbol_col <- intersect(
  colnames(utr_df),
  c("hgnc_symbol", "external_gene_name", "gene_symbol")
)[1]

utr_col <- "3utr"

if (is.na(symbol_col) || !utr_col %in% colnames(utr_df)) {
  stop(
    "BioMart output does not contain expected gene symbol or 3'UTR columns.",
    call. = FALSE
  )
}

# -------------------------
# FORMAT OUTPUT
# -------------------------
utr_df <- utr_df %>%
  transmute(
    gene_symbol = .data[[symbol_col]],
    utr_3 = .data[[utr_col]]
  ) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    !is.na(utr_3),
    utr_3 != "",
    utr_3 != "Sequence unavailable"
  ) %>%
  distinct() %>%
  mutate(utr_len = nchar(utr_3)) %>%
  arrange(gene_symbol, desc(utr_len)) %>%
  group_by(gene_symbol) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene_symbol, utr_3)

# -------------------------
# SAVE
# -------------------------
write_tsv(utr_df, out_file)

cat("Genes requested:", length(genes), "\n")
cat("Genes with 3'UTR retrieved:", nrow(utr_df), "\n")
cat("Output:", out_file, "\n")
