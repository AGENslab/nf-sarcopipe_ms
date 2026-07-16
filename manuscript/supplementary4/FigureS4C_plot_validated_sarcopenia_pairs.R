#!/usr/bin/env Rscript

############################################################
# Figure 5E – Validated musculoskeletal miRNA-target pairs
#
# Description
# -----------
# Generates a bubbleplot comparing validated de novo and
# annotated miRNA interactions with musculoskeletal genes
# targeted by the final de novo miRNA candidates.
#
# Color:
# - de novo miRNAs
# - annotated miRNAs
#
# Point size:
# - number of canonical seed-matched sites
#
# Point shape:
# - validation method
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript Figure5E_plot_validated_musculoskeletal_pairs.R ",
      "<input.tsv> <output.png> <output.pdf>\n"
    ),
    call. = FALSE
  )
}

input_file <- args[1]
out_png <- args[2]
out_pdf <- args[3]

if (!file.exists(input_file)) {
  stop(
    "Input table not found: ",
    input_file,
    call. = FALSE
  )
}

dir.create(
  dirname(out_png),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(out_pdf),
  recursive = TRUE,
  showWarnings = FALSE
)

df <- read_tsv(
  input_file,
  show_col_types = FALSE
)

required_cols <- c(
  "gene_symbol",
  "category",
  "renamed_miRNA",
  "miRNA_type",
  "n_sites",
  "log2FoldChange",
  "validated_by"
)

missing_cols <- setdiff(
  required_cols,
  colnames(df)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

df <- df %>%
  mutate(
    n_sites = as.numeric(n_sites),

    miRNA_type = factor(
      miRNA_type,
      levels = c(
        "de_novo",
        "annotated"
      ),
      labels = c(
        "De novo",
        "Annotated"
      )
    ),

    pair_label = paste0(
      renamed_miRNA,
      " \u2192 ",
      gene_symbol
    ),

    pair_label = str_wrap(
      pair_label,
      width = 34
    ),

    category_label = str_wrap(
      category,
      width = 24
    )
  ) %>%
  arrange(
    category,
    gene_symbol,
    miRNA_type,
    padj
  ) %>%
  mutate(
    pair_label = factor(
      pair_label,
      levels = rev(unique(pair_label))
    )
  )

type_colors <- c(
  "De novo" = "#6A3D9A",
  "Annotated" = "#E7298A"
)

validation_shapes <- c(
  "miRanda" = 21,
  "multiMiR" = 22
)

p <- ggplot(
  df,
  aes(
    x = abs(log2FoldChange),
    y = pair_label,
    size = n_sites,
    fill = miRNA_type,
    shape = validated_by
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = abs(log2FoldChange),
      y = pair_label,
      yend = pair_label
    ),
    linewidth = 0.55,
    color = "grey70"
  ) +
  geom_point(
    color = "grey20",
    stroke = 0.4,
    alpha = 0.95
  ) +
  facet_grid(
    category_label ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(
    values = type_colors,
    name = "miRNA type"
  ) +
  scale_shape_manual(
    values = validation_shapes,
    name = "Supported by"
  ) +
  scale_size_continuous(
    range = c(3.5, 8),
    breaks = sort(unique(df$n_sites)),
    name = "Seed-matched\nsites"
  ) +
  labs(
    title = "Validated miRNA interactions with musculoskeletal genes",
    subtitle = paste0(
      "De novo interactions supported by miRanda and annotated ",
      "interactions supported by multiMiR"
    ),
    x = "Absolute miRNA log2 fold change",
    y = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    strip.text.y = element_text(
      face = "bold",
      size = 10
    ),
    strip.background = element_rect(
      fill = "grey96",
      color = "grey65"
    ),
    axis.text.y = element_text(
      size = 9
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

plot_height <- max(
  6,
  2.5 + 0.38 * nrow(df)
)

ggsave(
  out_png,
  p,
  width = 10.5,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_pdf,
  p,
  width = 10.5,
  height = plot_height,
  bg = "white"
)

cat("Figure 5E PNG written:", out_png, "\n")
cat("Figure 5E PDF written:", out_pdf, "\n")
