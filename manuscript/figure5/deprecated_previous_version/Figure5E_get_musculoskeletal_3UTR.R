#!/usr/bin/env Rscript

############################################################
# Figure 5E – Retrieve 3'UTR sequences for a musculoskeletal
# gene set
#
# Description
# -----------
# Retrieves annotated human 3'UTR sequences from Ensembl
# BioMart for a user-provided gene list.
#
# For genes with multiple available 3'UTR sequences, the
# longest sequence is retained.
#
# Inputs
# ------
# 1) Text file with one HGNC gene symbol per line
# 2) Output TSV with retrieved 3'UTR sequences
# 3) Optional output TSV listing genes without a retrieved 3'UTR
#
# Outputs
# -------
# - gene_symbol
# - utr_3
# - utr_length
#
# Usage
# -----
# Rscript Figure5E_get_musculoskeletal_3UTR.R \
#   Figure5E_musculoskeletal_gene_symbols.txt \
#   Figure5E_musculoskeletal_3UTRs.tsv \
#   Figure5E_musculoskeletal_missing_3UTRs.tsv
############################################################

suppressPackageStartupMessages({
  library(biomaRt)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 3) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript Figure5E_get_musculoskeletal_3UTR.R ",
      "<gene_symbols.txt> <out_3UTR.tsv> [out_missing.tsv]\n"
    ),
    call. = FALSE
  )
}

genes_file <- args[1]
out_file <- args[2]

out_missing <- if (length(args) == 3) {
  args[3]
} else {
  file.path(
    dirname(out_file),
    "Figure5E_musculoskeletal_missing_3UTRs.tsv"
  )
}

if (!file.exists(genes_file)) {
  stop(
    "Input gene-symbol file not found: ",
    genes_file,
    call. = FALSE
  )
}

dir.create(
  dirname(out_file),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(out_missing),
  recursive = TRUE,
  showWarnings = FALSE
)

genes <- read_lines(
  genes_file,
  progress = FALSE
) %>%
  trimws() %>%
  toupper()

genes <- unique(
  genes[
    !is.na(genes) &
      genes != ""
  ]
)

if (length(genes) == 0) {
  stop(
    "No valid gene symbols were found in: ",
    genes_file,
    call. = FALSE
  )
}

cat("Genes requested:", length(genes), "\n")

connect_ensembl <- function() {
  mirrors <- c(
    "useast",
    "www",
    "asia"
  )

  last_error <- NULL

  for (mirror_name in mirrors) {
    cat(
      "Trying Ensembl mirror:",
      mirror_name,
      "\n"
    )

    mart <- tryCatch(
      {
        useEnsembl(
          biomart = "genes",
          dataset = "hsapiens_gene_ensembl",
          mirror = mirror_name
        )
      },
      error = function(error) {
        last_error <<- error
        NULL
      }
    )

    if (!is.null(mart)) {
      return(mart)
    }
  }

  stop(
    "Could not connect to any Ensembl BioMart mirror. Last error: ",
    conditionMessage(last_error),
    call. = FALSE
  )
}

mart <- connect_ensembl()

utr_raw <- tryCatch(
  {
    getSequence(
      id = genes,
      type = "hgnc_symbol",
      seqType = "3utr",
      mart = mart
    )
  },
  error = function(error) {
    stop(
      "BioMart 3'UTR query failed: ",
      conditionMessage(error),
      call. = FALSE
    )
  }
)

if (nrow(utr_raw) == 0) {
  stop(
    "BioMart returned no 3'UTR sequences for the requested genes.",
    call. = FALSE
  )
}

colnames(utr_raw) <- tolower(
  colnames(utr_raw)
)

symbol_candidates <- c(
  "hgnc_symbol",
  "external_gene_name",
  "gene_symbol"
)

symbol_col <- symbol_candidates[
  symbol_candidates %in% colnames(utr_raw)
][1]

utr_candidates <- c(
  "3utr",
  "3_utr",
  "utr3"
)

utr_col <- utr_candidates[
  utr_candidates %in% colnames(utr_raw)
][1]

if (
  is.na(symbol_col) ||
  is.na(utr_col)
) {
  stop(
    paste0(
      "BioMart output does not contain the expected gene-symbol ",
      "or 3'UTR columns. Columns returned: ",
      paste(
        colnames(utr_raw),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

utr_table <- utr_raw %>%
  transmute(
    gene_symbol = toupper(
      trimws(
        as.character(.data[[symbol_col]])
      )
    ),
    utr_3 = toupper(
      gsub(
        "\\s+",
        "",
        as.character(.data[[utr_col]])
      )
    )
  ) %>%
  mutate(
    utr_3 = gsub(
      "T",
      "U",
      utr_3,
      fixed = TRUE
    )
  ) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != "",
    gene_symbol %in% genes,
    !is.na(utr_3),
    utr_3 != "",
    utr_3 != "SEQUENCEUNAVAILABLE"
  ) %>%
  mutate(
    utr_length = nchar(utr_3)
  ) %>%
  filter(
    utr_length > 0
  ) %>%
  distinct(
    gene_symbol,
    utr_3,
    .keep_all = TRUE
  ) %>%
  arrange(
    gene_symbol,
    desc(utr_length)
  ) %>%
  group_by(
    gene_symbol
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  arrange(
    gene_symbol
  )

retrieved_genes <- unique(
  utr_table$gene_symbol
)

missing_genes <- setdiff(
  genes,
  retrieved_genes
)

missing_table <- tibble(
  gene_symbol = missing_genes,
  status = "3UTR_not_retrieved_from_Ensembl_BioMart"
)

write_tsv(
  utr_table,
  out_file
)

write_tsv(
  missing_table,
  out_missing
)

cat(
  "Genes with retrieved 3'UTR:",
  nrow(utr_table),
  "\n"
)

cat(
  "Genes without retrieved 3'UTR:",
  length(missing_genes),
  "\n"
)

if (length(missing_genes) > 0) {
  cat(
    "Missing genes:",
    paste(
      missing_genes,
      collapse = ", "
    ),
    "\n"
  )
}

cat("Written:", out_file, "\n")
cat("Written:", out_missing, "\n")
