# ============================================================
# plot_union_target_enrichment.R
# Description:
# Generates GO and KEGG bubble plots for the UNION of coherent
# miRNA–mRNA target genes predicted from BrumiR and miRDeep2.
# This script is intended for Figure 5c–d.
# Inputs:
#   ../output/union_coherent_GO.tsv
#   ../output/union_coherent_KEGG.tsv
# Outputs:
#   ../plots/Figure5c_union_GO.png
#   ../plots/Figure5d_union_KEGG.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

go_file <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/union_coherent_GO.tsv"
kegg_file <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/union_coherent_KEGG.tsv"

go_plot <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5c_union_GO.png"
kegg_plot <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5d_union_KEGG.png"

plot_bubble <- function(df, title, outfile) {
  if (nrow(df) == 0) {
    message("No enrichment results for: ", outfile)
    return(NULL)
  }

  top_n <- min(10, nrow(df))

  df2 <- df %>%
    slice_min(order_by = p.adjust, n = top_n) %>%
    mutate(
      GeneRatio_num = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
    ) %>%
    arrange(GeneRatio_num)

  df2$Description <- factor(df2$Description, levels = df2$Description)

  p <- ggplot(df2, aes(x = GeneRatio_num, y = Description, size = Count, color = p.adjust)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#C94C4C",   # más significativo
      high = "#3B6FB6",  # menos significativo
      trans = "reverse"
    ) +
    labs(
      title = title,
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10),
      legend.title = element_text(face = "bold")
    )

  ggsave(outfile, p, width = 8, height = 5, dpi = 300)
}

go <- read_tsv(go_file, show_col_types = FALSE)
kegg <- read_tsv(kegg_file, show_col_types = FALSE)

plot_bubble(
  go,
  "Union of coherent target genes: GO enrichment",
  go_plot
)

plot_bubble(
  kegg,
  "Union of coherent target genes: KEGG enrichment",
  kegg_plot
)
