# ============================================================
# plot_filtered_union_GO.R
# Description:
# Generates a filtered and publication-ready GO bubble plot
# for the union of coherent miRNA–mRNA predicted target genes.
# Only biologically interpretable GO terms relevant to the
# study narrative are retained for visualization.
# Input:
#   ../output/union_coherent_GO.tsv
# Output:
#   ../plots/Figure5c_filtered_union_GO.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

go_file <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/union_coherent_GO.tsv"
go_plot <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5c_filtered_union_GO.png"

go <- read_tsv(go_file, show_col_types = FALSE)

# términos a conservar
keep_terms <- c(
  "regulation of lymphocyte migration",
  "lymphocyte migration",
  "regulation of lymphocyte chemotaxis",
  "lymphocyte chemotaxis",
  "natural killer cell chemotaxis",
  "stimulatory C-type lectin receptor signaling pathway",
  "cellular response to lectin",
  "response to lectin"
)

go2 <- go %>%
  filter(Description %in% keep_terms) %>%
  mutate(
    GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    Description_wrapped = str_wrap(Description, width = 38)
  ) %>%
  arrange(GeneRatio_num)

# ordenar eje Y
go2$Description_wrapped <- factor(go2$Description_wrapped, levels = go2$Description_wrapped)

p <- ggplot(go2, aes(x = GeneRatio_num, y = Description_wrapped, size = Count, color = p.adjust)) +
  geom_point(alpha = 0.95) +
  scale_color_gradient(
    low = "#C94C4C",
    high = "#3B6FB6",
    trans = "reverse"
  ) +
  coord_cartesian(xlim = c(0.075, 0.116)) +
  labs(
    title = "Filtered GO enrichment of coherent predicted target genes",
    x = "Gene ratio",
    y = NULL,
    color = "Adjusted p-value",
    size = "Count"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(go_plot, p, width = 11, height = 5.5, dpi = 300)
