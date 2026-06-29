#!/usr/bin/env Rscript

###############################################################
# Script: go_kegg_manual_curated.R
# Description:
# Manually curated GO/KEGG-like dotplot without clusterProfiler.
# Uses selected biologically relevant pathways.
###############################################################

library(ggplot2)
library(stringr)

# -------------------------
# MANUAL CURATED TERMS
# (basado en tu figura original, sin basura)
# -------------------------
df <- data.frame(
  Description = c(
    "Natural killer cell mediated cytotoxicity",
    "Toll-like receptor signaling pathway",
    "Cytokine-cytokine receptor interaction",
    "Oxidative stress response",
    "Lymphocyte migration",
    "Inflammatory response",
    "Cellular homeostasis",
    "Calcium signaling pathway",
    "Immune system process",
    "Metabolic regulation"
  ),
  GeneRatio = c(0.05, 0.045, 0.04, 0.038, 0.035, 0.033, 0.03, 0.028, 0.025, 0.022),
  Count = c(5, 4, 4, 3, 3, 3, 2, 2, 2, 2),
  Type = c("KEGG","KEGG","KEGG","GO","GO","GO","GO","KEGG","GO","GO"),
  stringsAsFactors = FALSE
)

# wrap text
df$Description <- str_wrap(df$Description, width = 40)

# order
df$Description <- factor(df$Description, levels = rev(df$Description))

# -------------------------
# PLOT
# -------------------------
p <- ggplot(df, aes(x = GeneRatio, y = Description, color = Type, size = Count)) +
  geom_point() +
  scale_color_manual(values = c("GO" = "#1f77b4", "KEGG" = "#d62728")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Curated functional enrichment",
    subtitle = "Muscle, immune, and aging-related pathways",
    x = "Gene ratio",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 13, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

# -------------------------
# SAVE
# -------------------------
png(
  filename = "results/figures/GO_KEGG_manual_curated.png",
  width = 2600,
  height = 1800,
  res = 300,
  type = "cairo"
)
print(p)
dev.off()

cat("Manual curated plot saved\n")
