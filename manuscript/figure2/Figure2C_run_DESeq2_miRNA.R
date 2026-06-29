#!/usr/bin/env Rscript

# ============================================================
# run_deseq2_mirna.R
# ============================================================
#
# Description:
# This script performs differential expression analysis for miRNA count matrices
# using DESeq2. It is designed to run two analyses in one execution:
#   1) BrumiR count matrix
#   2) miRDeep2 count matrix
#
# Inputs:
#   --brumir_counts   TSV file with BrumiR counts
#   --mirdeep2_counts TSV file with miRDeep2 counts
#   --sample_info     TSV file with sample metadata
#   --outdir          Output directory
#
# Expected sample_info columns:
#   - sample
#   - condition
#
# Notes:
# - The first level in the condition factor is used as reference.
# - The second level is used as the comparison group.
# - The script writes full DESeq2 results, significant results,
#   significant results with abs(log2FC) >= 1, and VST matrices.
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list(
    brumir_counts = NULL,
    mirdeep2_counts = NULL,
    sample_info = NULL,
    outdir = "."
  )

  i <- 1
  while (i <= length(x)) {
    key <- x[i]

    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }

    if (i == length(x)) {
      stop("Missing value for argument: ", key)
    }

    value <- x[i + 1]
    name <- sub("^--", "", key)

    if (!name %in% names(out)) {
      stop("Unknown argument: ", key)
    }

    out[[name]] <- value
    i <- i + 2
  }

  out
}

opt <- parse_args(args)

if (is.null(opt$brumir_counts) || is.null(opt$mirdeep2_counts) || is.null(opt$sample_info)) {
  stop(
    paste(
      "Usage:",
      "Rscript run_deseq2_mirna.R",
      "--brumir_counts brumir.counts.tsv",
      "--mirdeep2_counts miRDeep2.counts.tsv",
      "--sample_info sample_info.tsv",
      "--outdir results_dir",
      sep = "\n"
    )
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
})

run_one_deseq <- function(counts_file, sample_file, out_prefix) {

  counts <- read.delim(counts_file, check.names = FALSE)
  samples <- read.delim(sample_file, sep = "\t", stringsAsFactors = FALSE)

  required_cols <- c("sample", "condition")
  missing_cols <- setdiff(required_cols, colnames(samples))
  if (length(missing_cols) > 0) {
    stop(
      "Sample info file is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # first col = feature id
  rownames(counts) <- counts[[1]]
  counts <- counts[, -1, drop = FALSE]

  # keep only shared samples, same order
  common <- intersect(colnames(counts), samples$sample)
  if (length(common) < 2) {
    stop("Fewer than 2 shared samples between counts and sample_info for: ", counts_file)
  }

  counts <- counts[, common, drop = FALSE]
  samples <- samples[match(common, samples$sample), , drop = FALSE]

  samples$condition <- factor(samples$condition)

  if (length(levels(samples$condition)) != 2) {
    stop("This script expects exactly 2 conditions in sample_info for: ", counts_file)
  }

  ref_cond <- levels(samples$condition)[1]
  test_cond <- levels(samples$condition)[2]
  samples$condition <- relevel(samples$condition, ref = ref_cond)

  # numeric matrix
  count_mat <- as.matrix(counts)
  mode(count_mat) <- "integer"

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = samples,
    design    = ~ condition
  )

  # minimal filtering
  dds <- dds[rowSums(counts(dds)) >= 10, ]

  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", test_cond, ref_cond))
  res <- as.data.frame(res)
  res$feature <- rownames(res)
  res <- res[, c("feature", setdiff(colnames(res), "feature"))]

  write.csv(res, paste0(out_prefix, "_results.csv"), row.names = FALSE)

  sig <- subset(res, !is.na(padj) & padj < 0.05)
  write.csv(sig, paste0(out_prefix, "_sig.csv"), row.names = FALSE)

  sig_lfc <- subset(res, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1)
  write.csv(sig_lfc, paste0(out_prefix, "_sig_absLFC1.csv"), row.names = FALSE)

  vsd <- vst(dds, blind = TRUE)
  vsd_mat <- assay(vsd)
  vsd_mat <- cbind(feature = rownames(vsd_mat), as.data.frame(vsd_mat))
  write.csv(vsd_mat, paste0(out_prefix, "_vst.csv"), row.names = FALSE)

  cat("Done:", out_prefix, "\n")
  cat("Reference condition:", ref_cond, "\n")
  cat("Comparison condition:", test_cond, "\n")
  cat("Samples used:", ncol(count_mat), "\n")
  cat("Features tested:", nrow(dds), "\n")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

run_one_deseq(
  counts_file = opt$brumir_counts,
  sample_file = opt$sample_info,
  out_prefix  = file.path(opt$outdir, "brumir_deseq2")
)

run_one_deseq(
  counts_file = opt$mirdeep2_counts,
  sample_file = opt$sample_info,
  out_prefix  = file.path(opt$outdir, "mirdeep2_deseq2")
)
