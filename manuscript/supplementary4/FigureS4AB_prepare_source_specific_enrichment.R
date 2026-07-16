#!/usr/bin/env Rscript

############################################################
# Figure 5D – GO, KEGG and Reactome enrichment
#
# Description
# -----------
# Performs functional enrichment using the unique target genes
# contained in the final validated miRNA–mRNA regulatory network.
#
# Only interactions supported by:
# - miRanda for BrumiR-derived de novo miRNAs
# - multiMiR for miRDeep2 annotated miRNAs
#
# are retained.
#
# Analyses
# --------
# - Gene Ontology Biological Process
# - KEGG pathways
# - Reactome pathways
#
# No genes, interaction counts, or pathway names are hardcoded.
#
# Usage
# -----
# Rscript Figure5D_GO_KEGG_Reactome_enrichment.R \
#   <cytoscape_edges.tsv> \
#   <out_target_genes.tsv> \
#   <out_gene_mapping.tsv> \
#   <out_go_bp.tsv> \
#   <out_kegg.tsv> \
#   <out_reactome.tsv> \
#   <out_summary.tsv>
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript Figure5D_GO_KEGG_Reactome_enrichment.R ",
      "<cytoscape_edges.tsv> ",
      "<out_target_genes.tsv> ",
      "<out_gene_mapping.tsv> ",
      "<out_go_bp.tsv> ",
      "<out_kegg.tsv> ",
      "<out_reactome.tsv> ",
      "<out_summary.tsv>\n"
    ),
    call. = FALSE
  )
}

input_edges <- args[1]
out_target_genes <- args[2]
out_mapping <- args[3]
out_go <- args[4]
out_kegg <- args[5]
out_reactome <- args[6]
out_summary <- args[7]

if (!file.exists(input_edges)) {
  stop(
    "Input Cytoscape edge table not found: ",
    input_edges,
    call. = FALSE
  )
}

for (path in c(
  out_target_genes,
  out_mapping,
  out_go,
  out_kegg,
  out_reactome,
  out_summary
)) {
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


require_columns <- function(data, required, label) {
  missing_columns <- setdiff(required, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}


result_to_table <- function(result) {
  if (is.null(result)) {
    return(tibble())
  }

  result_df <- as.data.frame(result)

  if (nrow(result_df) == 0) {
    return(tibble())
  }

  as_tibble(result_df)
}


write_enrichment <- function(result, output_path) {
  result_df <- result_to_table(result)
  write_tsv(result_df, output_path)
  invisible(result_df)
}


count_terms <- function(result_df) {
  if (nrow(result_df) == 0) {
    return(0L)
  }

  nrow(result_df)
}


count_significant <- function(result_df, cutoff = 0.05) {
  if (
    nrow(result_df) == 0 ||
    !"p.adjust" %in% colnames(result_df)
  ) {
    return(0L)
  }

  sum(
    !is.na(result_df$p.adjust) &
      result_df$p.adjust < cutoff
  )
}


# ----------------------------------------------------------
# Read final validated Cytoscape edges
# ----------------------------------------------------------

edges <- read_tsv(
  input_edges,
  show_col_types = FALSE,
  progress = FALSE
)

require_columns(
  edges,
  c(
    "target_node",
    "validated_by"
  ),
  "Cytoscape edge table"
)

validated_edges <- edges %>%
  filter(
    validated_by %in% c(
      "miRanda",
      "multiMiR"
    )
  )

if (nrow(validated_edges) == 0) {
  stop(
    "No miRanda- or multiMiR-validated interactions were found.",
    call. = FALSE
  )
}

target_genes <- validated_edges %>%
  transmute(
    SYMBOL = trimws(as.character(target_node))
  ) %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  distinct() %>%
  arrange(SYMBOL)

if (nrow(target_genes) == 0) {
  stop(
    "No target genes were recovered from the validated network.",
    call. = FALSE
  )
}

write_tsv(
  target_genes,
  out_target_genes
)


# ----------------------------------------------------------
# Map SYMBOL to ENTREZID
# ----------------------------------------------------------

mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = target_genes$SYMBOL,
  keytype = "SYMBOL",
  columns = c(
    "SYMBOL",
    "ENTREZID"
  )
) %>%
  as_tibble() %>%
  filter(
    !is.na(SYMBOL),
    !is.na(ENTREZID),
    SYMBOL != "",
    ENTREZID != ""
  ) %>%
  distinct(
    SYMBOL,
    ENTREZID
  ) %>%
  arrange(
    SYMBOL,
    ENTREZID
  )

if (nrow(mapping) == 0) {
  stop(
    "None of the validated target genes mapped to Entrez IDs.",
    call. = FALSE
  )
}

write_tsv(
  mapping,
  out_mapping
)

mapped_symbols <- unique(mapping$SYMBOL)

unmapped_symbols <- setdiff(
  target_genes$SYMBOL,
  mapped_symbols
)

entrez_ids <- unique(
  as.character(mapping$ENTREZID)
)


# ----------------------------------------------------------
# GO Biological Process
# ----------------------------------------------------------

go_result <- enrichGO(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)


# ----------------------------------------------------------
# KEGG
# ----------------------------------------------------------

kegg_result <- tryCatch(
  {
    result <- enrichKEGG(
      gene = entrez_ids,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )

    if (
      !is.null(result) &&
      nrow(as.data.frame(result)) > 0
    ) {
      result <- setReadable(
        result,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID"
      )
    }

    result
  },
  error = function(error) {
    warning(
      "KEGG enrichment failed: ",
      conditionMessage(error)
    )
    NULL
  }
)


# ----------------------------------------------------------
# Reactome
# ----------------------------------------------------------

reactome_result <- tryCatch(
  {
    enrichPathway(
      gene = entrez_ids,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    )
  },
  error = function(error) {
    warning(
      "Reactome enrichment failed: ",
      conditionMessage(error)
    )
    NULL
  }
)


# ----------------------------------------------------------
# Export complete enrichment tables
# ----------------------------------------------------------

go_table <- write_enrichment(
  go_result,
  out_go
)

kegg_table <- write_enrichment(
  kegg_result,
  out_kegg
)

reactome_table <- write_enrichment(
  reactome_result,
  out_reactome
)


# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------

summary_table <- tibble(
  metric = c(
    "validated_edges",
    "unique_target_genes",
    "mapped_target_genes",
    "unmapped_target_genes",
    "GO_BP_terms_total",
    "GO_BP_terms_padj_lt_0.05",
    "KEGG_terms_total",
    "KEGG_terms_padj_lt_0.05",
    "Reactome_terms_total",
    "Reactome_terms_padj_lt_0.05"
  ),
  value = c(
    nrow(validated_edges),
    nrow(target_genes),
    length(mapped_symbols),
    length(unmapped_symbols),
    count_terms(go_table),
    count_significant(go_table),
    count_terms(kegg_table),
    count_significant(kegg_table),
    count_terms(reactome_table),
    count_significant(reactome_table)
  )
)

write_tsv(
  summary_table,
  out_summary
)


# ----------------------------------------------------------
# Console report
# ----------------------------------------------------------

cat("===== FUNCTIONAL ENRICHMENT =====\n")
cat("Validated edges:", nrow(validated_edges), "\n")
cat("Unique target genes:", nrow(target_genes), "\n")
cat("Mapped target genes:", length(mapped_symbols), "\n")
cat("Unmapped target genes:", length(unmapped_symbols), "\n")

if (length(unmapped_symbols) > 0) {
  cat(
    "Unmapped symbols:",
    paste(unmapped_symbols, collapse = ", "),
    "\n"
  )
}

cat(
  "GO BP terms:",
  count_terms(go_table),
  "| BH < 0.05:",
  count_significant(go_table),
  "\n"
)

cat(
  "KEGG terms:",
  count_terms(kegg_table),
  "| BH < 0.05:",
  count_significant(kegg_table),
  "\n"
)

cat(
  "Reactome terms:",
  count_terms(reactome_table),
  "| BH < 0.05:",
  count_significant(reactome_table),
  "\n"
)

cat("Written:", out_target_genes, "\n")
cat("Written:", out_mapping, "\n")
cat("Written:", out_go, "\n")
cat("Written:", out_kegg, "\n")
cat("Written:", out_reactome, "\n")
cat("Written:", out_summary, "\n")
