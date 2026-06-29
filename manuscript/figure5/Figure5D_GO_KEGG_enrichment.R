#!/usr/bin/env Rscript

############################################################
# Figure 5D – GO and KEGG enrichment of coherent target genes
#
# Description:
# This script performs GO and KEGG enrichment on the union of
# coherent miRNA–mRNA target genes predicted from BrumiR and
# miRDeep2. Gene symbols are mapped to Entrez IDs before
# enrichment analysis.
#
# Inputs:
# - Text file with coherent target gene symbols, one per line
# - Output GO TSV file
# - Output KEGG TSV file
# - Output gene mapping TSV file
#
# Outputs:
# - union_coherent_GO.tsv
# - union_coherent_KEGG.tsv
# - union_coherent_gene_mapping.tsv
#
# Usage:
# Rscript Figure5D_GO_KEGG_enrichment.R \
#   coherent_target_genes_union.txt \
#   union_coherent_GO.tsv \
#   union_coherent_KEGG.tsv \
#   union_coherent_gene_mapping.tsv
#
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript Figure5D_GO_KEGG_enrichment.R ",
    "<coherent_target_genes_union.txt> ",
    "<union_coherent_GO.tsv> ",
    "<union_coherent_KEGG.tsv> ",
    "<union_coherent_gene_mapping.tsv>",
    call. = FALSE
  )
}

input_file <- args[1]
out_go <- args[2]
out_kegg <- args[3]
out_mapping <- args[4]

if (!file.exists(input_file)) {
  stop("Input gene list not found: ", input_file, call. = FALSE)
}

output_files <- c(out_go, out_kegg, out_mapping)
for (output_file in output_files) {
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

genes <- read_lines(input_file)
genes <- unique(genes[genes != ""])

mapping <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Hs.eg.db
)

write_tsv(mapping, out_mapping)

entrez <- unique(mapping$ENTREZID)

ego <- enrichGO(
  gene = entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

ekegg <- enrichKEGG(
  gene = entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)

if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  write_tsv(as.data.frame(ego), out_go)
} else {
  write_tsv(tibble(), out_go)
}

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  write_tsv(as.data.frame(ekegg), out_kegg)
} else {
  write_tsv(tibble(), out_kegg)
}

cat("Genes input:", length(genes), "\n")
cat("Mapped genes:", nrow(mapping), "\n")
cat("GO written to:", out_go, "\n")
cat("KEGG written to:", out_kegg, "\n")
cat("Mapping written to:", out_mapping, "\n")