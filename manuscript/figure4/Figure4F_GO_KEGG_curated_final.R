#!/usr/bin/env Rscript

############################################################
# Figure 4F – Curated GO and KEGG enrichment
#
# Description:
# This script performs exploratory GO and KEGG enrichment
# analysis using the 36 corrected mRNAs and generates a curated
# dotplot. Irrelevant terms, such as viral/infectious disease
# pathways and other categories not directly related to exercise,
# aging, muscle remodeling, or sarcopenia, are removed using
# predefined filtering patterns.
#
# Inputs:
# - CSV file containing the 36 corrected genes with ENTREZID column
# - Output PNG file
#
# Outputs:
# - Curated GO + KEGG enrichment dotplot PNG
#
# Usage:
# Rscript Figure4F_GO_KEGG_curated_final.R \
#   DEG_36_for_networks.csv \
#   GO_KEGG_combined_curated.png
#
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure4F_GO_KEGG_curated_final.R ",
    "<DEG_36_for_networks.csv> <output.png>",
    call. = FALSE
  )
}

input_file <- args[1]
output_file <- args[2]

if (!file.exists(input_file)) {
  stop("Input CSV not found: ", input_file, call. = FALSE)
}

out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD GENES
# -------------------------
df <- read.csv(input_file, stringsAsFactors = FALSE)

if (!"ENTREZID" %in% colnames(df)) {
  stop("Input CSV must contain column: ENTREZID", call. = FALSE)
}

genes_entrez <- unique(na.omit(df$ENTREZID))
genes_entrez <- as.character(genes_entrez)

cat("Genes with Entrez ID:", length(genes_entrez), "\n")

# -------------------------
# GO ENRICHMENT
# -------------------------
ego <- enrichGO(
  gene = genes_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# -------------------------
# KEGG ENRICHMENT
# -------------------------
ekegg <- enrichKEGG(
  gene = genes_entrez,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# -------------------------
# CONVERT TO DATA FRAMES
# -------------------------
go_df <- as.data.frame(ego)
kegg_df <- as.data.frame(ekegg)

if (nrow(go_df) == 0 && nrow(kegg_df) == 0) {
  stop("No enrichment terms found.", call. = FALSE)
}

# -------------------------
# CURATE TERMS
# -------------------------
remove_pattern <- paste(
  c(
    "viral",
    "virus",
    "chagas",
    "hepatitis",
    "influenza",
    "herpes",
    "epstein",
    "cytomegalo",
    "tuberculosis",
    "malaria",
    "measles",
    "pertussis",
    "salmonella",
    "toxoplasmosis",
    "leishmania",
    "amoebiasis",
    "covid",
    "hiv",
    "aids",
    "papilloma",
    "carcinoma",
    "cancer",
    "tumor",
    "leukemia",
    "glioma",
    "melanoma",
    "infection"
  ),
  collapse = "|"
)

if (nrow(go_df) > 0) {
  go_df <- go_df %>%
    filter(!str_detect(tolower(Description), remove_pattern))
}

if (nrow(kegg_df) > 0) {
  kegg_df <- kegg_df %>%
    filter(!str_detect(tolower(Description), remove_pattern))
}

# -------------------------
# OPTIONAL BIOLOGICAL RELEVANCE FILTER
# -------------------------
keep_pattern <- paste(
  c(
    "muscle",
    "calcium",
    "oxidative",
    "stress",
    "migration",
    "chemotaxis",
    "lymphocyte",
    "natural killer",
    "immune",
    "cytokine",
    "inflammatory",
    "homeostasis",
    "metabolic",
    "granzyme",
    "cell death",
    "toll-like",
    "dna-sensing",
    "th1",
    "th2"
  ),
  collapse = "|"
)

if (nrow(go_df) > 6) {
  go_df <- go_df %>%
    filter(str_detect(tolower(Description), keep_pattern))
}

if (nrow(kegg_df) > 6) {
  kegg_df <- kegg_df %>%
    filter(str_detect(tolower(Description), keep_pattern))
}

# -------------------------
# SELECT TOP TERMS
# -------------------------
go_plot <- go_df %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 8) %>%
  mutate(Type = "GO")

kegg_plot <- kegg_df %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 5) %>%
  mutate(Type = "KEGG")

plot_df <- bind_rows(go_plot, kegg_plot)

if (nrow(plot_df) == 0) {
  stop("No curated terms left after filtering.", call. = FALSE)
}

plot_df$GeneRatio_num <- sapply(plot_df$GeneRatio, function(x) {
  ratio <- strsplit(x, "/")[[1]]
  as.numeric(ratio[1]) / as.numeric(ratio[2])
})

plot_df$Description <- factor(
  plot_df$Description,
  levels = rev(plot_df$Description[order(plot_df$GeneRatio_num)])
)

# -------------------------
# PLOT
# -------------------------
plot_df$Description <- stringr::str_wrap(plot_df$Description, width = 40)

p <- ggplot(
  plot_df,
  aes(x = GeneRatio_num, y = Description, color = Type, size = Count)
) +
  geom_point() +
  scale_color_manual(values = c("GO" = "#1f77b4", "KEGG" = "#d62728")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Exploratory GO + KEGG enrichment",
    subtitle = "36 corrected genes (curated terms)",
    x = "Gene ratio",
    y = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

# -------------------------
# SAVE
# -------------------------
png(
  filename = output_file,
  width = 2600,
  height = 1800,
  res = 300,
  type = "cairo",
  bg = "white"
)

print(p)
dev.off()

cat("Curated GO + KEGG plot written to:", output_file, "\n")
