#!/usr/bin/env Rscript

###############################################################
# Script: go_kegg_36_mrnas_curated.R
# Author: Natalia Poblete
# Description:
# This script performs exploratory GO and KEGG enrichment
# analysis using the 36 corrected mRNAs and generates a dotplot
# similar to the original figure, but removing irrelevant terms
# such as viral/infectious disease pathways and other categories
# not directly related to exercise, aging, muscle remodeling,
# or sarcopenia.
#
# Input:
#   ../results/tables/DEG_36_for_networks.csv
#
# Output:
#   ../results/figures/GO_KEGG_combined_curated.png
###############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# -------------------------
# PATHS
# -------------------------
base_dir <- "/mnt/beegfs/home/npoblete/AGENsLab/Sarcopenia/analysis/mrnaobj1/nfcore/rnaseq_confounding_analysis"
input_file <- file.path(base_dir, "results/tables/DEG_36_for_networks.csv")
output_file <- file.path(base_dir, "results/figures/GO_KEGG_combined_curated.png")

# -------------------------
# LOAD GENES
# -------------------------
df <- read.csv(input_file, stringsAsFactors = FALSE)

# use Entrez IDs when available
genes_entrez <- unique(na.omit(df$ENTREZID))
genes_entrez <- as.character(genes_entrez)

cat("Genes with Entrez ID:", length(genes_entrez), "\n")

# -------------------------
# GO ENRICHMENT
# -------------------------
ego <- enrichGO(
  gene          = genes_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.2,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# -------------------------
# KEGG ENRICHMENT
# -------------------------
ekegg <- enrichKEGG(
  gene          = genes_entrez,
  organism      = "hsa",
  keyType       = "ncbi-geneid",
  pvalueCutoff  = 0.2,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.2
)

# -------------------------
# CONVERT TO DATA FRAMES
# -------------------------
go_df <- as.data.frame(ego)
kegg_df <- as.data.frame(ekegg)

if (nrow(go_df) == 0 && nrow(kegg_df) == 0) {
  stop("No enrichment terms found.")
}

# -------------------------
# CURATE TERMS
# remove pathways not relevant to aging/exercise/sarcopenia
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
# OPTIONAL:
# keep only biologically relevant concepts if enough remain
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
  stop("No curated terms left after filtering.")
}

# GeneRatio as numeric
plot_df$GeneRatio_num <- sapply(plot_df$GeneRatio, function(x) {
  x <- strsplit(x, "/")[[1]]
  as.numeric(x[1]) / as.numeric(x[2])
})

# order for plotting
plot_df$Description <- factor(
  plot_df$Description,
  levels = rev(plot_df$Description[order(plot_df$GeneRatio_num)])
)

# -------------------------
# PLOT
# -------------------------

plot_df$Description <- stringr::str_wrap(plot_df$Description, width = 40)

p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description, color = Type, size = Count)) +
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
    #  TÍTULOS
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    
    #  EJE Y
    axis.text.y = element_text(size = 13, color = "black"), 
    axis.title.y = element_text(size = 14, face = "bold"),
    
    #  EJE X
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    
    #  LEYENDA
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )
    

# -------------------------
# SAVE (HPC-safe)
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
