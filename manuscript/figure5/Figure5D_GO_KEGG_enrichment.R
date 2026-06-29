# ============================================================
# run_union_target_enrichment.R
# Description:
# Performs GO and KEGG enrichment on the UNION of coherent
# miRNA–mRNA target genes predicted from BrumiR and miRDeep2.
# Input:
#   ../output/coherent_target_genes_union.txt
# Outputs:
#   ../output/union_coherent_GO.tsv
#   ../output/union_coherent_KEGG.tsv
#   ../output/union_coherent_gene_mapping.tsv
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(tibble)
})

input_file <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/coherent_target_genes_union.txt"
out_go <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/union_coherent_GO.tsv"
out_kegg <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/union_coherent_KEGG.tsv"
out_mapping <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/union_coherent_gene_mapping.tsv"

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
