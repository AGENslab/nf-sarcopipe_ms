# ============================================================
# plot_biological_prioritized_targets.R
# Description:
# Generates a publication-ready dot plot of biologically
# prioritized BrumiR target genes across sarcopenia/exercise-
# related pathways for Figure 5d. Category names are abbreviated
# in facet strips for better readability, while full biological
# meanings are preserved through the color legend and figure text.
#
# Input:
#   ../output/brumir_biological_targets_top.tsv
#
# Output:
#   ../plots/Figure5d_biological_prioritized_targets.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

infile <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/brumir_biological_targets_top.tsv"
outfile <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5d_biological_prioritized_targets.png"

df <- read_tsv(infile, show_col_types = FALSE)

# order of BrumiR candidates
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

p <- ggplot(df, aes(x = renamed_miRNA, y = gene_label, size = n_sites, color = category)) +
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

ggsave(outfile, p, width = 9, height = 10, dpi = 300)
