# ============================================================
# plot_prioritized_predicted_target_genes.R
# Description:
# Generates a publication-ready plot of prioritized predicted
# target genes based on coherent miRNA support, number of
# seed-matched sites, and miRNA effect size. Intended for
# Figure 5d.
#
# Input:
#   ../output/prioritized_predicted_target_genes_top15.tsv
#
# Output:
#   ../plots/Figure5d_prioritized_predicted_targets.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

infile <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/prioritized_predicted_target_genes_top15.tsv"
outfile <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5d_prioritized_predicted_targets.png"

df <- read_tsv(infile, show_col_types = FALSE)

df <- df %>%
  mutate(
    gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)),
    category = factor(category, levels = c("shared", "BrumiR_only", "miRDeep2_only"))
  )

p <- ggplot(df, aes(x = priority_score, y = gene_symbol, color = category, size = total_n_sites)) +
  geom_point(alpha = 0.95) +
  scale_color_manual(values = c(
    "shared" = "#7A4EAB",
    "BrumiR_only" = "#4D4D4D",
    "miRDeep2_only" = "#1B9E77"
  )) +
  labs(
    title = "Prioritized predicted target genes from coherent miRNA–mRNA pairs",
    x = "Priority score",
    y = NULL,
    color = "Support category",
    size = "Total\nsites"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.y = element_text(face = "bold", size = 10),
    legend.title = element_text(face = "bold")
  )

ggsave(outfile, p, width = 9, height = 6, dpi = 300)
