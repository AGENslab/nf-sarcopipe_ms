suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# plot_candidate_7mer_heatmap.R
# Description:
# Plots the 7-mer seed-space analysis of selected candidate
# miRNAs, highlighting canonical and non-canonical matches
# to known human miRNA seeds.
#
# Input:
#   kmer_seed_analysis/candidate_7mer_long.tsv
#
# Output:
#   kmer_seed_analysis/Supplementary_Figure_S2_7mer_heatmap.png
# ============================================================

df <- read_tsv("kmer_seed_analysis/candidate_7mer_long.tsv", show_col_types = FALSE)

# order candidates
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
  "kmer_seed_analysis/Supplementary_Figure_S2_7mer_heatmap.png",
  p, width = 9, height = 4.8, dpi = 300, bg = "white"
)

ggsave(
  "kmer_seed_analysis/Supplementary_Figure_S2_7mer_heatmap.pdf",
  p, width = 9, height = 4.8, bg = "white"
)
