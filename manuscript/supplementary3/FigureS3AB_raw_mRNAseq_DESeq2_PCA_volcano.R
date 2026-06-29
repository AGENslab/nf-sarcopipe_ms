#!/usr/bin/env Rscript

############################
# nf-core mRNA-seq downstream analysis
# DESeq2 workflow
# Author: NP
############################

# DESCRIPTION:
# This script performs downstream differential expression analysis starting from
# nf-core/rnaseq gene-level quantification files (quant.genes.sf).
#
# REQUIRED INPUTS:
#   --base_dir     Directory containing per-sample nf-core/rnaseq results
#   --samplesheet  CSV/TSV with at least:
#                    sample,condition
#                  Optional column:
#                    quant_file
#   --outdir       Output directory
#
# OPTIONAL INPUT:
#   --quant_pattern  Relative path from base_dir to each quantification file.
#                    Use "{sample}" as placeholder.
#                    Default: "{sample}/quant.genes.sf"
#
# EXAMPLE:
#   Rscript run_nfcore_deseq2_downstream.R \
#     --base_dir /path/to/nfcore/rnaseq/results/star_salmon \
#     --samplesheet /path/to/samplesheet.tsv \
#     --outdir paper_raw_mrna
#
# NOTES:
# - The samplesheet must contain a binary comparison in the column "condition".
# - The first condition encountered is used as reference level.
# - If "quant_file" is provided in the samplesheet, it is used directly.
# - Otherwise, quant files are inferred as:
#       file.path(base_dir, "{sample}/quant.genes.sf")

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(matrixStats)
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})

############################
# ARGUMENT PARSING
############################

parse_args <- function(x) {
  args <- list(
    base_dir = NULL,
    samplesheet = NULL,
    outdir = "paper_raw_mrna",
    quant_pattern = "{sample}/quant.genes.sf"
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
      "Rscript run_nfcore_deseq2_downstream.R",
      "--base_dir /path/to/star_salmon",
      "--samplesheet /path/to/samplesheet.tsv",
      "--outdir paper_raw_mrna",
      sep = "\n"
    )
  )
}

############################
# OUTPUT DIRECTORIES
############################

outdir <- args$outdir
figdir <- file.path(outdir, "figures")
tabledir <- file.path(outdir, "tables")
objdir <- file.path(outdir, "objects")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabledir, showWarnings = FALSE, recursive = TRUE)
dir.create(objdir, showWarnings = FALSE, recursive = TRUE)

############################
# READ SAMPLE SHEET
############################

read_samplesheet <- function(path) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "csv") {
    df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  required <- c("sample", "condition")
  missing_cols <- setdiff(required, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Samplesheet is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (nrow(df) < 2) {
    stop("Samplesheet must contain at least 2 samples.")
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
    stop("This script currently expects exactly 2 conditions in the samplesheet.")
  }

  df
}

samples_df <- read_samplesheet(args$samplesheet)

############################
# INPUT FILES
############################

build_quant_path <- function(base_dir, sample_name, quant_pattern) {
  rel <- gsub("\\{sample\\}", sample_name, args$quant_pattern)
  file.path(base_dir, rel)
}

if ("quant_file" %in% colnames(samples_df)) {
  files <- samples_df$quant_file
  files <- ifelse(grepl("^/", files), files, file.path(args$base_dir, files))
} else {
  files <- vapply(
    samples_df$sample,
    function(s) build_quant_path(args$base_dir, s, args$quant_pattern),
    character(1)
  )
}

samples <- samples_df$sample
names(files) <- samples

condition_levels <- unique(samples_df$condition)
condition <- factor(samples_df$condition, levels = condition_levels)

coldata <- data.frame(
  row.names = samples,
  condition = condition
)

############################
# CHECK INPUTS
############################

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

############################
# BUILD COUNT MATRIX
############################

cat("Building count matrix...\n")

counts_list <- lapply(files, function(f) {
  x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  needed <- c("Name", "NumReads")
  if (!all(needed %in% colnames(x))) {
    stop("File does not contain required columns Name and NumReads: ", f)
  }

  x[, c("Name", "NumReads")]
})

counts <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), counts_list)
rownames(counts) <- counts$Name
counts <- counts[, -1, drop = FALSE]
colnames(counts) <- samples
counts[is.na(counts)] <- 0
counts <- round(as.matrix(counts))

write.csv(
  counts,
  file.path(tabledir, "counts_raw_sin_filtrar_para_transparencia.csv")
)

############################
# FILTER LOW EXPRESSION
############################

cat("Filtering low expression genes...\n")

keep <- rowSums(counts >= 10) >= 2
counts_filt <- counts[keep, , drop = FALSE]

write.csv(
  counts_filt,
  file.path(tabledir, "counts_filtrados_para_DESeq2.csv")
)

############################
# DESEQ2 ANALYSIS
############################

cat("Running DESeq2...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)
res_refseq <- results(dds)

saveRDS(dds, file.path(objdir, "objeto_dds_DESeq2.rds"))
saveRDS(res_refseq, file.path(objdir, "objeto_resultados_DESeq2_refseq.rds"))

write.csv(
  as.data.frame(res_refseq),
  file.path(tabledir, "resultados_DESeq2_refseq.csv")
)

############################
# REFSEQ TO GENE SYMBOL
############################

cat("Annotating genes with biomaRt...\n")

refseq_full <- rownames(dds)
refseq_clean <- sub("\\..*$", "", refseq_full)

ensembl <- NULL
for (m in c("www", "uswest", "asia")) {
  try({
    ensembl <- biomaRt::useEnsembl(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      mirror = m
    )
  }, silent = TRUE)
  if (!is.null(ensembl)) break
}

if (is.null(ensembl)) {
  stop("Ensembl is not available from www/uswest/asia mirrors.")
}

annot <- getBM(
  attributes = c("refseq_mrna", "hgnc_symbol"),
  filters = "refseq_mrna",
  values = unique(refseq_clean),
  mart = ensembl
)

annot <- annot[annot$hgnc_symbol != "", , drop = FALSE]
annot <- annot[!duplicated(annot$refseq_mrna), , drop = FALSE]

symbol_map <- annot$hgnc_symbol
names(symbol_map) <- annot$refseq_mrna

gene_symbols <- symbol_map[refseq_clean]
gene_symbols[is.na(gene_symbols)] <- refseq_full[is.na(gene_symbols)]

rownames(dds) <- gene_symbols
res <- results(dds)

saveRDS(res, file.path(objdir, "objeto_resultados_DESeq2_gene_symbol.rds"))

write.csv(
  as.data.frame(res),
  file.path(tabledir, "resultados_DESeq2_gene_symbol.csv")
)

############################
# NORMALIZED MATRICES
############################

cat("Exporting normalized matrices...\n")

matriz_normalizada_para_modelos_nb <- counts(dds, normalized = TRUE)
write.csv(
  matriz_normalizada_para_modelos_nb,
  file.path(tabledir, "matriz_normalizada_para_modelos_nb.csv")
)

vsd <- vst(dds, blind = FALSE)
para_correlacion_mirnas_matrix <- assay(vsd)

write.csv(
  para_correlacion_mirnas_matrix,
  file.path(tabledir, "para_correlacion_mirnas_matrix.csv")
)

############################
# SIGNIFICANT GENES
############################

cat("Exporting significant genes...\n")

res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)

res_sig <- res_df %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write.csv(
  res_sig,
  file.path(tabledir, "genes_significativos_FDR0.05_absLFC1.csv"),
  row.names = FALSE
)

############################
# SUMMARY BARPLOT
############################

cat("Generating summary barplot...\n")

condition_ref <- condition_levels[1]
condition_test <- condition_levels[2]

genes_test_up <- res_sig %>%
  filter(log2FoldChange > 1) %>%
  pull(gene_symbol)

genes_ref_up <- res_sig %>%
  filter(log2FoldChange < -1) %>%
  pull(gene_symbol)

n_total_DE <- nrow(res_sig)

summary_df <- data.frame(
  Group = c("Total DE genes", paste0("Up in ", condition_test), paste0("Up in ", condition_ref)),
  Count = c(n_total_DE, length(genes_test_up), length(genes_ref_up))
)

p_resumen_DE <- ggplot(summary_df, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Summary of Differentially Expressed Genes",
    y = "Number of genes",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figdir, "resumen_genes_DE.png"),
  p_resumen_DE,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(figdir, "resumen_genes_DE.pdf"),
  p_resumen_DE,
  width = 8,
  height = 5
)

############################
# VOLCANO PLOT
############################

cat("Generating volcano plot...\n")

res_df <- as.data.frame(res)
res_df$gene_symbol <- rownames(res_df)
res_df$neglog10_padj <- -log10(res_df$padj)

res_df$status <- "Not significant"
res_df$status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- paste0("Up in ", condition_test)
res_df$status[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- paste0("Up in ", condition_ref)

res_df$status <- factor(
  res_df$status,
  levels = c(paste0("Up in ", condition_test), paste0("Up in ", condition_ref), "Not significant")
)

top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(desc(neglog10_padj)) %>%
  head(10)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj, color = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = "grey50"
  ) +
  coord_cartesian(ylim = c(0, 500)) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Volcano plot – ", condition_test, " vs ", condition_ref),
    x = paste0("log2 Fold Change (", condition_test, " / ", condition_ref, ")"),
    y = "-log10 adjusted p-value",
    color = "Gene status"
  )

ggsave(
  file.path(figdir, "volcano_plot_condition_comparison.png"),
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figdir, "volcano_plot_condition_comparison.pdf"),
  p_volcano,
  width = 8,
  height = 6
)

############################
# TOP 20 SIGNIFICANT GENES
############################

cat("Exporting top 20 significant genes...\n")

top_by_padj <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(desc(neglog10_padj)) %>%
  select(gene_symbol, log2FoldChange, padj, neglog10_padj) %>%
  head(20)

write.csv(
  top_by_padj,
  file.path(tabledir, "top20_genes_mas_significativos.csv"),
  row.names = FALSE
)

############################
# TOP 10 GENE BARPLOT
############################

cat("Generating top 10 gene barplot...\n")

top10_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(desc(neglog10_padj)) %>%
  pull(gene_symbol) %>%
  unique() %>%
  head(10)

norm_counts <- counts(dds, normalized = TRUE)
top10_genes <- intersect(top10_genes, rownames(norm_counts))

plot_df <- as.data.frame(norm_counts[top10_genes, , drop = FALSE])
plot_df$gene <- rownames(plot_df)

plot_df <- plot_df %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(
    data.frame(sample = rownames(coldata), condition = coldata$condition),
    by = "sample"
  )

plot_df$condition <- factor(plot_df$condition, levels = condition_levels)

p_top10 <- ggplot(plot_df, aes(x = gene, y = expression, fill = condition)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Normalized expression of top 10 differentially expressed genes",
    x = "Gene",
    y = "Normalized expression (DESeq2)",
    fill = "Condition"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figdir, "barplot_top10_genes.png"),
  p_top10,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figdir, "barplot_top10_genes.pdf"),
  p_top10,
  width = 10,
  height = 6
)

############################
# PCA
############################

cat("Generating PCA...\n")

vsd_mat <- assay(vsd)
gene_variance <- matrixStats::rowVars(vsd_mat)

top_n <- min(5000, nrow(vsd_mat))
top_idx <- order(gene_variance, decreasing = TRUE)[seq_len(top_n)]
vsd_top <- vsd[top_idx, ]

pca <- prcomp(t(assay(vsd_top)), scale. = FALSE)
var_exp <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = colData(vsd_top)$condition
)

pca_df$condition <- factor(pca_df$condition, levels = condition_levels)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = condition), linewidth = 1, level = 0.95) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("PCA of Top ", top_n, " Most Variable Genes"),
    subtitle = "Variance-stabilized expression (DESeq2)",
    x = paste0("PC1 (", var_exp[1], "% variance)"),
    y = paste0("PC2 (", var_exp[2], "% variance)"),
    color = "Condition"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(figdir, "PCA_top_variable_genes.png"),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figdir, "PCA_top_variable_genes.pdf"),
  p_pca,
  width = 8,
  height = 6
)

############################
# HEATMAP
############################

cat("Generating heatmap...\n")

mat_top <- assay(vsd_top)

annotation_col <- data.frame(condition = colData(vsd_top)$condition)
rownames(annotation_col) <- colnames(mat_top)
annotation_col$condition <- factor(annotation_col$condition, levels = condition_levels)

stopifnot(all(rownames(annotation_col) == colnames(mat_top)))

pheatmap(
  mat_top,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_col = annotation_col,
  fontsize_col = 10,
  main = paste0("Heatmap – Top ", top_n, " most variable genes"),
  filename = file.path(figdir, "heatmap_top_variable_genes.png"),
  width = 2000 / 300,
  height = 1800 / 300
)

pheatmap(
  mat_top,
  scale = "row",
  show_rownames = FALSE,
  show_colnames = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  annotation_col = annotation_col,
  fontsize_col = 10,
  main = paste0("Heatmap – Top ", top_n, " most variable genes"),
  filename = file.path(figdir, "heatmap_top_variable_genes.pdf"),
  width = 9,
  height = 8
)

cat("Analysis completed successfully.\n")
