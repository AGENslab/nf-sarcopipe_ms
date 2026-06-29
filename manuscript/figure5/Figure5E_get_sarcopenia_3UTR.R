# ============================================================
# get_sarcopenia_gene_set_3utr.R
# Description:
# Downloads annotated 3'UTR sequences for the curated
# sarcopenia/exercise gene set using Ensembl BioMart.
#
# Input:
#   ../input/input_sarcopenia_gene_set.txt
#
# Output:
#   ../output/sarcopenia_gene_set_3UTR.tsv
# ============================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(readr)
  library(dplyr)
})

genes_file <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/input/input_sarcopenia_gene_set.txt"
out_file   <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/sarcopenia_gene_set_3UTR.tsv"

genes <- read_lines(genes_file)
genes <- unique(genes[genes != ""])

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

symbol_col <- intersect(colnames(utr_df), c("hgnc_symbol", "external_gene_name", "gene_symbol"))[1]
utr_col <- "3utr"

utr_df <- utr_df %>%
  transmute(
    gene_symbol = .data[[symbol_col]],
    utr_3 = .data[[utr_col]]
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "", !is.na(utr_3), utr_3 != "", utr_3 != "Sequence unavailable") %>%
  distinct() %>%
  mutate(utr_len = nchar(utr_3)) %>%
  arrange(gene_symbol, desc(utr_len)) %>%
  group_by(gene_symbol) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene_symbol, utr_3)

write_tsv(utr_df, out_file)

cat("Genes requested:", length(genes), "\n")
cat("Genes with 3'UTR retrieved:", nrow(utr_df), "\n")
cat("Output:", out_file, "\n")
