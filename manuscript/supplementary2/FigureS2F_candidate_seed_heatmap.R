#!/usr/bin/env Rscript

############################################################
# Figure S2F – Final de novo candidate 7-mer seed heatmap
#
# Description
# -----------
# Visualizes every possible 7-mer in the final BrumiR de novo
# candidate sequences and indicates whether each 7-mer matches
# a known human canonical miRNA seed.
#
# Candidate names and their order are read directly from the
# input table. No candidate IDs or candidate counts are
# hardcoded.
#
# Usage
# -----
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

if (length(args) != 3) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript FigureS2F_candidate_seed_heatmap.R ",
      "<candidate_7mer_long.tsv> <out_png> <out_pdf>\n"
    )
  )
}

input_tsv <- args[1]
out_png <- args[2]
out_pdf <- args[3]

if (!file.exists(input_tsv)) {
  stop(
    "Input candidate 7-mer table not found: ",
    input_tsv
  )
}

dir.create(
  dirname(out_png),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(out_pdf),
  recursive = TRUE,
  showWarnings = FALSE
)

df <- read_tsv(
  input_tsv,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "provisional_miRNA_ID",
  "kmer_start_pos_1based",
  "is_canonical_seed_2_8",
  "kmer_status"
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns in input table: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (nrow(df) == 0) {
  stop("Input candidate 7-mer table contains no rows.")
}

candidate_order <- df %>%
  distinct(provisional_miRNA_ID) %>%
  pull(provisional_miRNA_ID)

df <- df %>%
  mutate(
    provisional_miRNA_ID = factor(
      provisional_miRNA_ID,
      levels = rev(candidate_order)
    ),
    fill_group = case_when(
      is_canonical_seed_2_8 == "yes" &
        kmer_status == "known_seed_match" ~
        "canonical_known",

      is_canonical_seed_2_8 == "yes" &
        kmer_status == "novel_seed" ~
        "canonical_novel",

      is_canonical_seed_2_8 == "no" &
        kmer_status == "known_seed_match" ~
        "noncanonical_known",

      TRUE ~ "noncanonical_novel"
    )
  )

n_candidates <- length(candidate_order)

plot_height <- max(
  4.8,
  1.0 + 0.55 * n_candidates
)

p <- ggplot(
  df,
  aes(
    x = kmer_start_pos_1based,
    y = provisional_miRNA_ID,
    fill = fill_group
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_x_continuous(
    breaks = sort(
      unique(df$kmer_start_pos_1based)
    )
  ) +
  scale_fill_manual(
    values = c(
      "canonical_known" = "#D73027",
      "canonical_novel" = "#FC8D59",
      "noncanonical_known" = "#4575B4",
      "noncanonical_novel" = "#D9D9D9"
    ),
    breaks = c(
      "canonical_known",
      "canonical_novel",
      "noncanonical_known",
      "noncanonical_novel"
    ),
    labels = c(
      "Canonical seed matching a known miRNA seed",
      "Canonical seed without a known miRNA seed match",
      "Non-canonical 7-mer matching a known miRNA seed",
      "No known miRNA seed match"
    ),
    name = NULL
  ) +
  labs(
    title = paste0(
      "7-mer seed-space analysis of final de novo miRNA candidates"
    ),
    subtitle = paste0(
      "All possible 7-mers were compared with known human ",
      "canonical miRNA seeds"
    ),
    x = "7-mer start position within the mature miRNA",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(
      size = 10
    )
  )

ggsave(
  filename = out_png,
  plot = p,
  width = 10.5,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 10.5,
  height = plot_height,
  bg = "white"
)

cat(
  "Candidates plotted:",
  n_candidates,
  "\n"
)

cat(
  "Candidate seed heatmap PNG written to:",
  out_png,
  "\n"
)

cat(
  "Candidate seed heatmap PDF written to:",
  out_pdf,
  "\n"
)