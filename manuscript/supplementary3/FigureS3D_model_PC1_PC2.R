#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(matrixStats)
})

# -------------------------
# PATHS Y MUESTRAS
# -------------------------
base_dir <- "/mnt/beegfs/home/npoblete/AGENsLab/Sarcopenia/analysis/mrnaobj1/nfcore/rnaseq/results/star_salmon"

samples <- c(
  "SRR1424731","SRR1424738","SRR1424739","SRR1424741","SRR1424745","SRR1424754",
  "SRR13442895","SRR13442897","SRR13442899","SRR13442901","SRR13442903","SRR13442905","SRR13442907"
)

files <- file.path(base_dir, samples, "quant.genes.sf")
names(files) <- samples

condition <- c(rep("sedentary", 6), rep("active", 7))
coldata <- data.frame(row.names = samples, condition = factor(condition))

# -------------------------
# LEER COUNTS DESDE SALMON
# -------------------------
counts_list <- lapply(files, function(f) {
  x <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
# FILTRADO
# -------------------------
keep <- rowSums(counts >= 10) >= 2
counts_filt <- counts[keep, , drop = FALSE]

# -------------------------
# DESEQ2 + VST
# -------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# -------------------------
# PCA con top 5000 genes más variables
# -------------------------
gene_variance <- rowVars(vsd_mat)
top5000_idx <- order(gene_variance, decreasing = TRUE)[1:min(5000, length(gene_variance))]
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
# MODELO LINEAL POR GEN
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

write.csv(results, "DEG_PC1_PC2_corrected.csv", row.names = FALSE)

cat("Modelo completado\n")
cat("Genes analizados:", nrow(results), "\n")
cat("Significativos FDR < 0.05:", sum(results$padj < 0.05, na.rm = TRUE), "\n")
