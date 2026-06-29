#!/usr/bin/env Rscript

############################################################
# Figure S2F – Candidate seed heatmap
#
# Description:
# Plots the 7-mer seed-space analysis of selected candidate
# miRNAs, highlighting canonical and non-canonical matches
# to known human miRNA seeds.
#
# Inputs:
# 1) candidate_7mer_long.tsv
#
# Outputs:
# - Supplementary_Figure_S2_7mer_heatmap.png
# - Supplementary_Figure_S2_7mer_heatmap.pdf
#
# Usage:
# Rscript FigureS2F_candidate_seed_heatmap.R \
#   <candidate_7mer_long.tsv> \
#   <out_png> \
#   <out_pdf>
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    "Usage: Rscript FigureS2F_candidate_seed_heatmap.R ",
    "<candidate_7mer_long.tsv> <out_png> <out_pdf>"
  )
}

input_tsv <- args[1]
out_png <- args[2]
out_pdf <- args[3]

if (!file.exists(input_tsv)) {
  stop("Input candidate 7-mer table not found: ", input_tsv)
}

out_png_dir <- dirname(out_png)
out_pdf_dir <- dirname(out_pdf)

if (!dir.exists(out_png_dir)) {
  dir.create(out_png_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(out_pdf_dir)) {
  dir.create(out_pdf_dir, recursive = TRUE, showWarnings = FALSE)
}

df <- read_tsv(input_tsv, show_col_types = FALSE)

required_cols <- c(
  "candidate",
  "kmer_start_pos_1based",
  "is_canonical_seed_2_8",
  "kmer_status"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in input table: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Order candidates
candidate_order <- c(
  "hsa_miR-660-5p",
  "known_seed_A",
  "hsa-miR-novel_A",
  "hsa-miR-novel_B",
  "hsa-miR-novel_C",
  "hsa-miR-novel_D"
)

df <- df %>%
  mutate(
    candidate = factor(candidate, levels = rev(candidate_order)),
    fill_group = case_when(
      is_canonical_seed_2_8 == "yes" & kmer_status == "known_seed_match" ~ "canonical_known",
      is_canonical_seed_2_8 == "yes" & kmer_status == "novel_seed" ~ "canonical_novel",
      is_canonical_seed_2_8 == "no" & kmer_status == "known_seed_match" ~ "noncanonical_known",
      TRUE ~ "noncanonical_novel"
    )
  )

p <- ggplot(df, aes(x = kmer_start_pos_1based, y = candidate, fill = fill_group)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c(
      "canonical_known" = "#D73027",
      "canonical_novel" = "#FC8D59",
      "noncanonical_known" = "#4575B4",
      "noncanonical_novel" = "#D9D9D9"
    ),
    labels = c(
      "Canonical known seed",
      "Canonical novel seed",
      "Non-canonical known seed",
      "No known seed match"
    ),
    name = ""
  ) +
  labs(
    title = "7-mer seed-space analysis of BrumiR candidate miRNAs",
    subtitle = "All possible 7-mers compared against known human canonical miRNA seeds",
    x = "7-mer start position within mature miRNA",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  out_png,
  p,
  width = 9,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_pdf,
  p,
  width = 9,
  height = 4.8,
  bg = "white"
)

cat("Candidate seed heatmap PNG written to:", out_png, "\n")
cat("Candidate seed heatmap PDF written to:", out_pdf, "\n")