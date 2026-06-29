#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Combined volcano plot for miRNA differential expression.
# This script merges DESeq2 results from miRDeep2 and BrumiR
# and visualizes log2 fold change versus adjusted p-value.
# Points are colored by regulation direction and shaped by method.
# The top 20 most significant miRNAs are labeled using ggrepel.
# ------------------------------------------------------------

.libPaths(c("~/Rlibs", .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript plot_volcano_mirna_expression.R <brumir_deseq2.csv> <mirdeep2_deseq2.csv> <out_prefix>")
}

brumir_file  <- args[1]
mirdeep_file <- args[2]
out_prefix   <- args[3]

out_plot <- paste0(out_prefix, "_volcano_mirna_expression.png")
out_tsv  <- paste0(out_prefix, "_volcano_mirna_expression.tsv")

# -------------------------
# LOAD DATA
# -------------------------
brumir <- read_csv(brumir_file, show_col_types = FALSE) %>%
  mutate(method = "BrumiR")

mirdeep <- read_csv(mirdeep_file, show_col_types = FALSE) %>%
  mutate(method = "miRDeep2")

# -------------------------
# MERGE AND CLASSIFY
# -------------------------
df <- bind_rows(brumir, mirdeep) %>%
  mutate(
    padj_plot = ifelse(is.na(padj) | padj <= 0, 1, padj),
    neglog10_padj = -log10(padj_plot),
    regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Up in Active",
      !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Up in Sedentary",
      TRUE ~ "Not significant"
    )
  )

write_tsv(df, out_tsv)

# -------------------------
# TOP 10 PER METHOD
# -------------------------
df_top <- df %>%
  filter(!is.na(padj)) %>%
  group_by(method) %>%
  arrange(padj, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

# wider x-axis
x_limit <- max(abs(df$log2FoldChange), na.rm = TRUE) * 1.35

# -------------------------
# PLOT
# -------------------------
p <- ggplot(df, aes(x = log2FoldChange, y = neglog10_padj)) +
  geom_point(aes(color = regulation, shape = method), size = 2.8, alpha = 0.85) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "solid",
    linewidth = 0.5,
    color = "grey60"
  ) +
  geom_text_repel(
    data = df_top,
    aes(label = feature),
    size = 4.2,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 123
  ) +
  coord_cartesian(xlim = c(-x_limit, x_limit)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.grid.major  = element_line(color = "grey85"),
    panel.grid.minor  = element_blank()
  ) +
  scale_color_manual(values = c(
    "Up in Active" = "blue",
    "Up in Sedentary" = "red",
    "Not significant" = "grey70"
  )) +
  scale_shape_manual(values = c(
    "BrumiR" = 17,
    "miRDeep2" = 16
  )) +
  labs(
    title = "Combined volcano plot of miRNA differential expression",
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    color = "miRNA regulation",
    shape = "Method"
  )

# -------------------------
# SAVE
# -------------------------
out_dir <- dirname(out_plot)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

png(
  filename = out_plot,
  width = 2600,
  height = 1900,
  res = 300,
  bg = "white",
  type = "cairo"
)
print(p)
dev.off()

cat("Combined volcano plot written to:", out_plot, "\n")
cat("Combined volcano table written to:", out_tsv, "\n")
