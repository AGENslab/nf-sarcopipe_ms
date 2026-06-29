#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Combined heatmap of miRNA expression from miRDeep2 and BrumiR-RF.
# Selects the top differentially expressed miRNAs from each method,
# applies variance-stabilizing transformation to count matrices,
# scales rows as z-scores, and draws a combined annotated heatmap.
# ------------------------------------------------------------

.libPaths(c("~/Rlibs", .libPaths()))

suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript plot_heatmap_mirna_expression.R <brumir_counts.tsv> <mirdeep_counts.tsv> <sample_info.tsv> <brumir_deseq2.csv> <mirdeep_deseq2.csv> <out_prefix>")
}

brumir_counts_file <- args[1]
mirdeep_counts_file <- args[2]
sample_info_file <- args[3]
brumir_deseq_file <- args[4]
mirdeep_deseq_file <- args[5]
out_prefix <- args[6]

out_plot <- paste0(out_prefix, "_heatmap_mirna_expression.png")
out_tsv  <- paste0(out_prefix, "_heatmap_mirna_expression.tsv")

# -------------------------
# HELPERS
# -------------------------
read_counts_matrix <- function(path) {
  x <- read_tsv(path, show_col_types = FALSE)
  x <- as.data.frame(x)
  rownames(x) <- x[,1]
  x <- x[,-1, drop = FALSE]
  x <- as.matrix(x)
  mode(x) <- "numeric"
  x
}

vst_matrix <- function(counts_mat, meta_df) {
  common_samples <- intersect(colnames(counts_mat), meta_df$sample)
  counts_mat <- counts_mat[, common_samples, drop = FALSE]
  meta_sub <- meta_df[match(common_samples, meta_df$sample), , drop = FALSE]
  rownames(meta_sub) <- meta_sub$sample
  meta_sub$condition <- factor(meta_sub$condition, levels = c("active", "sedentary"))

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_mat),
    colData = meta_sub,
    design = ~ condition
  )

  dds <- dds[rowSums(counts(dds)) > 10, ]
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  assay(vsd)
}

row_zscore <- function(mat) {
  t(scale(t(mat)))
}

# -------------------------
# LOAD INPUTS
# -------------------------
meta <- read_tsv(sample_info_file, show_col_types = FALSE) %>%
  as.data.frame()

br_counts <- read_counts_matrix(brumir_counts_file)
md_counts <- read_counts_matrix(mirdeep_counts_file)

br_deseq <- read_csv(brumir_deseq_file, show_col_types = FALSE) %>%
  mutate(method = "BrumiR-RF")

md_deseq <- read_csv(mirdeep_deseq_file, show_col_types = FALSE) %>%
  mutate(method = "miRDeep2")

# -------------------------
# SELECT TOP FEATURES
# -------------------------
top_br <- br_deseq %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 25) %>%
  pull(feature)

top_md <- md_deseq %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 25) %>%
  pull(feature)

# -------------------------
# VST MATRICES
# -------------------------
br_vst <- vst_matrix(br_counts, meta)
md_vst <- vst_matrix(md_counts, meta)

top_br <- intersect(top_br, rownames(br_vst))
top_md <- intersect(top_md, rownames(md_vst))

br_sub <- br_vst[top_br, , drop = FALSE]
md_sub <- md_vst[top_md, , drop = FALSE]

# Prefix rownames so the origin is explicit
rownames(br_sub) <- paste0("BR_", rownames(br_sub))
rownames(md_sub) <- paste0("MD_", rownames(md_sub))

# -------------------------
# COMBINE
# -------------------------
common_samples <- intersect(colnames(br_sub), colnames(md_sub))
br_sub <- br_sub[, common_samples, drop = FALSE]
md_sub <- md_sub[, common_samples, drop = FALSE]

combined <- rbind(md_sub, br_sub)
combined_z <- row_zscore(combined)
combined_z[is.na(combined_z)] <- 0

# Save numeric matrix
write_tsv(
  data.frame(feature = rownames(combined_z), combined_z, check.names = FALSE),
  out_tsv
)

# -------------------------
# ANNOTATIONS
# -------------------------
annotation_col <- meta %>%
  filter(sample %in% common_samples) %>%
  select(sample, condition) %>%
  as.data.frame()
rownames(annotation_col) <- annotation_col$sample
annotation_col <- annotation_col[common_samples, "condition", drop = FALSE]

annotation_row <- data.frame(
  Method = c(rep("miRDeep2", nrow(md_sub)), rep("BrumiR-RF", nrow(br_sub)))
)
rownames(annotation_row) <- rownames(combined_z)

ann_colors <- list(
  condition = c(active = "blue", sedentary = "red"),
  Method = c(`miRDeep2` = "#4D4D4D", `BrumiR-RF` = "#B8860B")
)

# -------------------------
# PLOT
# -------------------------
png(
  filename = out_plot,
  width = 2400,
  height = 2600,
  res = 300,
  bg = "white",
  type = "cairo"
)

pheatmap(
  combined_z,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  scale = "none",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 9,
  main = "Combined heatmap of miRNA expression"
)

dev.off()

cat("Heatmap written to:", out_plot, "\n")
cat("Heatmap matrix written to:", out_tsv, "\n")
