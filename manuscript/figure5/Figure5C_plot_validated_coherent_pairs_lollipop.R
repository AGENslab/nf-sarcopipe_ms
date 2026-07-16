#!/usr/bin/env Rscript

############################################################
# Figure 5C – Top validated miRNA–mRNA interactions
#
# Description
# -----------
# Generates a lollipop plot showing the top validated
# miRNA–mRNA interactions from each discovery branch.
#
# The input must be the final Cytoscape edge table containing
# only coherent interactions independently supported by:
#
# - miRanda for BrumiR de novo candidates
# - multiMiR for miRDeep2 annotated miRNAs
#
# Ranking
# -------
# BrumiR pairs:
#   1. lowest adjusted p value
#   2. highest miRanda score
#   3. largest absolute miRNA log2 fold change
#
# miRDeep2 pairs:
#   1. lowest adjusted p value
#   2. largest number of supporting multiMiR databases
#   3. largest absolute miRNA log2 fold change
#
# No candidate IDs, target genes, pair counts, or numerical
# results are hardcoded.
#
# Usage
# -----
# Rscript Figure5C_plot_validated_pairs_lollipop.R \
#   <validated_edges.tsv> \
#   <top_n_per_direction> \
#   <out_tsv> \
#   <out_png> \
#   <out_pdf>
#
# top_n_per_direction:
# Number of most positive and most negative pairs selected
# independently per discovery branch.
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript Figure5C_plot_validated_pairs_lollipop.R ",
      "<validated_edges.tsv> <top_n_per_source> ",
      "<out_tsv> <out_png> <out_pdf>\n"
    ),
    call. = FALSE
  )
}

input_file <- args[1]
top_n <- suppressWarnings(as.integer(args[2]))
out_tsv <- args[3]
out_png <- args[4]
out_pdf <- args[5]

if (!file.exists(input_file)) {
  stop(
    "Input validated edge table not found: ",
    input_file,
    call. = FALSE
  )
}

if (is.na(top_n) || top_n < 1) {
  stop(
    "top_n_per_source must be a positive integer.",
    call. = FALSE
  )
}

dir.create(
  dirname(out_tsv),
  recursive = TRUE,
  showWarnings = FALSE
)

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
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "source_node",
  "target_node",
  "miRNA_source",
  "miRNA_direction",
  "gene_direction",
  "miRNA_log2FoldChange",
  "miRNA_padj",
  "validated_by",
  "validation_database_list",
  "validation_database_count",
  "miranda_score",
  "miranda_energy"
)

missing_columns <- setdiff(
  required_columns,
  colnames(df)
)

if (length(missing_columns) > 0) {
  stop(
    "Input table is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(df) == 0) {
  stop(
    "Input validated edge table contains no interactions.",
    call. = FALSE
  )
}

# ----------------------------------------------------------
# Normalize numerical fields and remove nonvalidated records
# ----------------------------------------------------------

df <- df %>%
  mutate(
    miRNA_log2FoldChange = as.numeric(miRNA_log2FoldChange),
    miRNA_padj = as.numeric(miRNA_padj),
    validation_database_count = suppressWarnings(
      as.numeric(validation_database_count)
    ),
    miranda_score = suppressWarnings(
      as.numeric(miranda_score)
    ),
    miranda_energy = suppressWarnings(
      as.numeric(miranda_energy)
    ),
    validation_database_count = coalesce(
      validation_database_count,
      0
    ),
    abs_log2FC = abs(miRNA_log2FoldChange)
  ) %>%
  filter(
    validated_by %in% c("miRanda", "multiMiR")
  )

if (nrow(df) == 0) {
  stop(
    "No miRanda- or multiMiR-supported interactions were found.",
    call. = FALSE
  )
}

# ----------------------------------------------------------
# Select the most positive and most negative interactions
# within each discovery branch
# ----------------------------------------------------------

select_extremes <- function(data, n_each) {

  positive <- data %>%
    filter(miRNA_log2FoldChange > 0) %>%
    arrange(
      desc(miRNA_log2FoldChange),
      miRNA_padj,
      desc(validation_database_count),
      desc(miranda_score),
      source_node,
      target_node
    ) %>%
    slice_head(n = n_each) %>%
    mutate(
      regulation_group = "Most positive",
      source_rank = row_number()
    )

  negative <- data %>%
    filter(miRNA_log2FoldChange < 0) %>%
    arrange(
      miRNA_log2FoldChange,
      miRNA_padj,
      desc(validation_database_count),
      desc(miranda_score),
      source_node,
      target_node
    ) %>%
    slice_head(n = n_each) %>%
    mutate(
      regulation_group = "Most negative",
      source_rank = row_number()
    )

  bind_rows(
    negative,
    positive
  )
}

brumir_top <- df %>%
  filter(miRNA_source == "BrumiR") %>%
  select_extremes(top_n)

mirdeep2_top <- df %>%
  filter(miRNA_source == "miRDeep2") %>%
  select_extremes(top_n)

df_plot <- bind_rows(
  brumir_top,
  mirdeep2_top
)

if (nrow(df_plot) == 0) {
  stop(
    "No interactions remained after selecting positive and negative extremes.",
    call. = FALSE
  )
}

# ----------------------------------------------------------
# Plot variables
# ----------------------------------------------------------

df_plot <- df_plot %>%
  mutate(
    pair_label = paste0(
      source_node,
      " \u2192 ",
      target_node
    ),

    signed_log2FC = miRNA_log2FoldChange,

    direction_label = case_when(
      miRNA_direction == "Up_in_Active" ~
        "miRNA up in Active",

      miRNA_direction == "Up_in_Sedentary" ~
        "miRNA up in Sedentary",

      TRUE ~ miRNA_direction
    ),

    validation_label = case_when(
      validated_by == "miRanda" ~
        paste0(
          "miRanda score = ",
          format(
            miranda_score,
            trim = TRUE,
            scientific = FALSE
          ),
          "; energy = ",
          format(
            miranda_energy,
            trim = TRUE,
            scientific = FALSE
          )
        ),

      validated_by == "multiMiR" ~
        paste0(
          validation_database_count,
          " supporting database(s)"
        ),

      TRUE ~ validated_by
    )
  )

# Unique factor levels prevent duplicated pair labels across facets
df_plot <- df_plot %>%
  mutate(
    pair_key = paste(
      miRNA_source,
      pair_label,
      sep = "__"
    )
  ) %>%
  arrange(
    miRNA_source,
    signed_log2FC,
    source_rank
  )

pair_levels <- unique(df_plot$pair_key)

df_plot <- df_plot %>%
  mutate(
    pair_key = factor(
      pair_key,
      levels = pair_levels
    )
  )

pair_labels <- setNames(
  df_plot$pair_label,
  df_plot$pair_key
)

# ----------------------------------------------------------
# Export exact plotting table
# ----------------------------------------------------------

output_table <- df_plot %>%
  select(
    miRNA_source,
    source_rank,
    miRNA = source_node,
    target_gene = target_node,
    miRNA_direction,
    gene_direction,
    miRNA_log2FoldChange,
    miRNA_padj,
    validated_by,
    validation_database_list,
    validation_database_count,
    miranda_score,
    miranda_energy,
    validation_label
  )

write_tsv(
  output_table,
  out_tsv
)

# ----------------------------------------------------------
# Plot
# ----------------------------------------------------------

p <- ggplot(
  df_plot,
  aes(
    x = signed_log2FC,
    y = pair_key,
    color = direction_label
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = signed_log2FC,
      yend = pair_key
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5
  ) +
  facet_grid(
    miRNA_source ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(
    labels = pair_labels
  ) +
  scale_color_manual(
    values = c(
      "miRNA up in Active" = "#3B6FB6",
      "miRNA up in Sedentary" = "#C94C4C"
    ),
    name = "miRNA regulation"
  ) +
  labs(
    title = paste0(
      "Top validated miRNA\u2013mRNA regulatory interactions"
    ),
    subtitle = paste0(
      "Coherent pairs supported by miRanda or multiMiR"
    ),
    x = "miRNA log2 fold change",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(
      size = 9,
      face = "bold"
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

plot_height <- max(
  7,
  2.8 + 0.36 * nrow(df_plot)
)

ggsave(
  filename = out_png,
  plot = p,
  width = 11,
  height = plot_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 11,
  height = plot_height,
  bg = "white"
)

cat("===== FIGURE 5C =====\n")

cat(
  "Pairs requested per direction and source:",
  top_n,
  "\n"
)

count_summary <- df_plot %>%
  count(miRNA_source, name = "n_selected")

print(count_summary)

cat("Written:", out_tsv, "\n")
cat("Written:", out_png, "\n")
cat("Written:", out_pdf, "\n")
