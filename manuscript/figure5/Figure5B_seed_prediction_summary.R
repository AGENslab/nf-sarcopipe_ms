library(ggplot2)
library(readr)
library(tidyr)

# input
df <- read_tsv(
  "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/fig5a_summary.tsv",
  show_col_types = FALSE
)

# formato largo
long <- pivot_longer(df, cols = -source, names_to = "metric", values_to = "value")

# etiquetas bonitas
long$metric <- factor(
  long$metric,
  levels = c("n_miRNAs", "n_genes", "n_pairs", "n_coherent"),
  labels = c("miRNAs evaluated", "Genes with match", "Raw pairs", "Coherent pairs")
)

colors_algorithms <- c(
  "BrumiR" = "#4D4D4D",   # gris elegante
  "miRDeep2" = "#1B9E77"  # azul petróleo/verde profesional
)

# plot
p <- ggplot(long, aes(x = metric, y = value, fill = source)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = value),
            position = position_dodge(width = 0.75),
            vjust = -0.25, size = 3.5) +
  scale_fill_manual(values = colors_algorithms) +
  labs(
    title = "Seed-based target prediction summary across algorithms",
    x = NULL,
    y = "Count",
    fill = "Algorithm"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# guardar
ggsave(
  "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5a_summary.png",
  p,
  width = 9,
  height = 5,
  dpi = 300
)