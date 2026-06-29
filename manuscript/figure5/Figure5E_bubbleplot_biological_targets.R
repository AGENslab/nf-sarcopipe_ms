#!/usr/bin/env Rscript

############################################################
# Figure 5E – Biological prioritized target bubbleplot
#
# Description:
# This script generates a publication-ready dot plot of
# biologically prioritized BrumiR target genes across
# sarcopenia/exercise-related pathways. Category names are
# abbreviated in facet strips for readability, while full
# biological meanings are preserved through the color legend
# and figure text.
#
# Inputs:
# - TSV file with biologically prioritized BrumiR targets
# - Output PNG file
#
# Outputs:
# - Biological prioritized target bubbleplot PNG
#
# Usage:
# Rscript Figure5E_bubbleplot_biological_targets.R \
#   brumir_biological_targets_top.tsv \
#   Figure5E_biological_prioritized_targets.png
#
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure5E_bubbleplot_biological_targets.R ",
    "<brumir_biological_targets_top.tsv> <output.png>",
    call. = FALSE
  )
}

infile <- args[1]
outfile <- args[2]

if (!file.exists(infile)) {
  stop("Input TSV not found: ", infile, call. = FALSE)
}

out_dir <- dirname(outfile)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD DATA
# -------------------------
df <- read_tsv(infile, show_col_types = FALSE)

required_cols <- c(
  "renamed_miRNA",
  "gene_symbol",
  "category",
  "priority_score",
  "n_sites"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input TSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -------------------------
# PREPARE DATA
# -------------------------
candidate_levels <- c(
  "hsa-miR-660-5p",
  "known_seed_A",
  "hsa-miR-novel_A",
  "hsa-miR-novel_B",
  "hsa-miR-novel_C",
  "hsa-miR-novel_D"
)

df <- df %>%
  mutate(
    renamed_miRNA = factor(renamed_miRNA, levels = candidate_levels),
    gene_label = paste0(gene_symbol, "  "),
    category = factor(
      category,
      levels = c(
        "Inflammation / immunity",
        "Muscle atrophy / growth",
        "Fibrosis / ECM remodeling",
        "Senescence / damage",
        "Autophagy / degradation",
        "Myogenesis / muscle structure"
      ),
      labels = c(
        "Inflamm./Imm.",
        "Atrophy/Growth",
        "Fibrosis/ECM",
        "Senescence",
        "Autophagy",
        "Myogenesis"
      )
    )
  ) %>%
  arrange(category, priority_score) %>%
  mutate(gene_label = factor(gene_label, levels = unique(gene_label)))

category_colors <- c(
  "Inflamm./Imm." = "#c94c4c",
  "Atrophy/Growth" = "#d17c2f",
  "Fibrosis/ECM" = "#7a4eab",
  "Senescence" = "#8c564b",
  "Autophagy" = "#4c9a8a",
  "Myogenesis" = "#3b6fb6"
)

# -------------------------
# PLOT
# -------------------------
p <- ggplot(
  df,
  aes(x = renamed_miRNA, y = gene_label, size = n_sites, color = category)
) +
  geom_point(alpha = 0.95) +
  scale_color_manual(values = category_colors) +
  scale_size_continuous(range = c(3, 10)) +
  labs(
    title = "Biologically prioritized predicted targets of BrumiR candidates",
    x = "BrumiR candidate miRNAs",
    y = NULL,
    color = "Biological category",
    size = "Seed-matched\nsites"
  ) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text.y = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# -------------------------
# SAVE
# -------------------------
ggsave(outfile, p, width = 9, height = 10, dpi = 300)

cat("Biological target bubbleplot written to:", outfile, "\n")
