library(ggplot2)
library(readr)
library(dplyr)

infile <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/output/coherent_pairs_for_lollipop.tsv"
outfile <- "/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/plots/Figure5b_coherent_lollipop.png"

df <- read_tsv(infile, show_col_types = FALSE)

# 🔥 seleccionar top miRDeep2 (15 mejores por padj)
md_top <- df %>%
  filter(source == "miRDeep2") %>%
  arrange(padj) %>%
  slice(1:15)

# mantener todos los BrumiR
br_all <- df %>%
  filter(source == "BrumiR")

df_plot <- bind_rows(br_all, md_top)

# construir variables
df_plot <- df_plot %>%
  mutate(
    pair_label = paste0(renamed_miRNA, " \u2192 ", gene_symbol),
    signed_value = ifelse(direction == "Up_in_Active", -abs_log2FC, abs_log2FC)
  )

# ordenar por valor
df_plot <- df_plot %>%
  group_by(source) %>%
  arrange(signed_value, .by_group = TRUE) %>%
  mutate(pair_label = factor(pair_label, levels = unique(pair_label))) %>%
  ungroup()

# 🎨 plot mejorado
p <- ggplot(df_plot, aes(x = signed_value, y = pair_label, color = direction)) +
  geom_segment(aes(x = 0, xend = signed_value, y = pair_label, yend = pair_label),
               linewidth = 0.7) +
  geom_point(size = 3) +
  facet_wrap(~source, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(
    "Up_in_Active" = "#3B6FB6",     # azul
    "Up_in_Sedentary" = "#C94C4C"   # rojo
  )) +
  labs(
    title = "Coherent miRNA–mRNA target pairs identified by seed-based prediction",
    x = "miRNA |log2 fold change| (signed by direction)",
    y = NULL,
    color = "miRNA direction"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 10, face = "bold"),  # 🔥 mejora clave
    axis.text.x = element_text(size = 10),
    legend.position = "right"
  )

ggsave(outfile, p, width = 10, height = 9, dpi = 300)
