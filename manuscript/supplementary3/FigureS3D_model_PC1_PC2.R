#!/usr/bin/env Rscript

############################################################
# Figure S3D – Linear model corrected by PC1 and PC2
#
# Description:
# Runs a gene-wise linear model on variance-stabilized mRNA-seq
# expression values while correcting for PC1 and PC2.
#
# The script:
# - reads Salmon quant.genes.sf files
# - builds a count matrix
# - filters low-expression genes
# - runs DESeq2 and VST
# - computes PCA using the top 5000 most variable genes
# - fits one linear model per gene:
#     expression ~ group + PC1 + PC2
# - adjusts p-values using Benjamini-Hochberg FDR
#
# Inputs:
# --base_dir
#     Directory containing one folder per sample with quant.genes.sf files.
#
# --samplesheet
#     CSV/TSV file with columns:
#     - sample
#     - condition
#
# --out_csv
#     Output CSV path.
#
# Outputs:
# - DEG_PC1_PC2_corrected.csv
#
# Usage:
# Rscript FigureS3D_model_PC1_PC2.R \
#   --base_dir /path/to/star_salmon \
#   --samplesheet metadata.tsv \
#   --out_csv DEG_PC1_PC2_corrected.csv
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(matrixStats)
})

# -------------------------
# Argument parsing
# -------------------------
parse_args <- function(x) {
  args <- list(
    base_dir = NULL,
    samplesheet = NULL,
    out_csv = "DEG_PC1_PC2_corrected.csv"
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

    if (!name %in% names(args)) {
      stop("Unknown argument: ", key)
    }

    args[[name]] <- value
    i <- i + 2
  }

  args
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(args$base_dir) || is.null(args$samplesheet)) {
  stop(
    paste(
      "Missing required arguments.",
      "Usage:",
      "Rscript FigureS3D_model_PC1_PC2.R",
      "--base_dir /path/to/star_salmon",
      "--samplesheet metadata.tsv",
      "--out_csv DEG_PC1_PC2_corrected.csv",
      sep = "\n"
    )
  )
}

if (!dir.exists(args$base_dir)) {
  stop("Base directory not found: ", args$base_dir)
}

if (!file.exists(args$samplesheet)) {
  stop("Samplesheet not found: ", args$samplesheet)
}

out_dir <- dirname(args$out_csv)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# Read samplesheet
# -------------------------
read_samplesheet <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "csv") {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  required_cols <- c("sample", "condition")
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Samplesheet is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  df$sample <- trimws(df$sample)
  df$condition <- trimws(df$condition)

  if (any(df$sample == "")) {
    stop("Samplesheet contains empty sample names.")
  }

  if (any(df$condition == "")) {
    stop("Samplesheet contains empty condition values.")
  }

  if (length(unique(df$condition)) != 2) {
    stop("This script expects exactly 2 conditions.")
  }

  df
}

samples_df <- read_samplesheet(args$samplesheet)

samples <- samples_df$sample
files <- file.path(args$base_dir, samples, "quant.genes.sf")
names(files) <- samples

condition <- samples_df$condition
coldata <- data.frame(
  row.names = samples,
  condition = factor(condition)
)

# -------------------------
# Check input files
# -------------------------
missing_files <- files[!file.exists(files)]

if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

# -------------------------
# Read counts from Salmon
# -------------------------
counts_list <- lapply(files, function(f) {
  x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  required_cols <- c("Name", "NumReads")
  missing_cols <- setdiff(required_cols, colnames(x))

  if (length(missing_cols) > 0) {
    stop(
      "File is missing required columns Name and NumReads: ",
      f
    )
  }

  x <- x[, c("Name", "NumReads")]
  x
})

counts <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), counts_list)
rownames(counts) <- counts$Name
counts <- counts[, -1, drop = FALSE]
colnames(counts) <- samples
counts[is.na(counts)] <- 0
counts <- round(counts)

# -------------------------
# Filtering
# -------------------------
keep <- rowSums(counts >= 10) >= 2
counts_filt <- counts[keep, , drop = FALSE]

# -------------------------
# DESeq2 + VST
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# -------------------------
# PCA with top 5000 most variable genes
# -------------------------
gene_variance <- rowVars(vsd_mat)
top5000_idx <- order(gene_variance, decreasing = TRUE)[
  1:min(5000, length(gene_variance))
]
vsd_top5000 <- vsd_mat[top5000_idx, , drop = FALSE]

pca <- prcomp(t(vsd_top5000), scale. = FALSE)

meta <- data.frame(
  sample = colnames(vsd_mat),
  group = coldata$condition,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

# -------------------------
# Gene-wise linear model
# -------------------------
pvals <- numeric(nrow(vsd_mat))
logFC <- numeric(nrow(vsd_mat))

group <- factor(meta$group, levels = c("active", "sedentary"))

for (i in seq_len(nrow(vsd_mat))) {
  y <- as.numeric(vsd_mat[i, ])

  model <- lm(y ~ group + PC1 + PC2, data = meta)
  coef_summary <- summary(model)$coefficients

  if ("groupsedentary" %in% rownames(coef_summary)) {
    pvals[i] <- coef_summary["groupsedentary", "Pr(>|t|)"]
    logFC[i] <- coef_summary["groupsedentary", "Estimate"]
  } else {
    pvals[i] <- NA
    logFC[i] <- NA
  }
}

padj <- p.adjust(pvals, method = "BH")

results <- data.frame(
  gene = rownames(vsd_mat),
  logFC = logFC,
  pvalue = pvals,
  padj = padj
)

write.csv(results, args$out_csv, row.names = FALSE)

cat("Modelo completado\n")
cat("Genes analizados:", nrow(results), "\n")
cat("Significativos FDR < 0.05:", sum(results$padj < 0.05, na.rm = TRUE), "\n")
cat("Written:", args$out_csv, "\n")
