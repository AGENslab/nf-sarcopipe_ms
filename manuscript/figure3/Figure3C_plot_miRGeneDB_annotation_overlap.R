# ============================================================
# Figure: miRNA annotation and overlap between BrumiR and miRDeep2
# ============================================================

library(ggplot2)
library(dplyr)

df <- data.frame(
  algorithm = c("BrumiR", "BrumiR",
                "miRDeep2", "miRDeep2", "miRDeep2"),
  category = c("Shared annotated", "Putative novel",
               "Shared annotated", "Unique annotated", "Putative novel"),
  count = c(141, 79,
            141, 238, 562)
)

df$algorithm <- factor(df$algorithm, levels = c("miRDeep2", "BrumiR"))
df$category <- factor(
  df$category,
  levels = c("Shared annotated", "Unique annotated", "Putative novel")
)

p <- ggplot(df, aes(x = algorithm, y = count, fill = category)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.25) +
  coord_flip() +
  scale_fill_manual(
    breaks = c("Putative novel", "Unique annotated", "Shared annotated"),
    values = c(
      "Putative novel" = "#BFBFBF",
      "Unique annotated" = "#DD8452",
      "Shared annotated" = "#4C72B0"
    )
  ) +
  labs(
    title = "miRNA annotation against miRGeneDB",
    subtitle = "≥98% identity, ≥75% query coverage",
    x = "",
    y = "Number of miRNAs",
    fill = ""
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black")
  )

print(p)

ggsave(
  "Figure_miRGeneDB_annotation_overlap_clean.png",
  p,
  width = 8,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Figure_miRGeneDB_annotation_overlap_clean.pdf",
  p,
  width = 8,
  height = 4.8,
  bg = "white"
)
