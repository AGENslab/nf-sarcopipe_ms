#!/usr/bin/env Rscript

############################################################
# Figure 5B – Progressive prioritization of miRNA–mRNA pairs
#
# Description
# -----------
# Summarizes the progression from final differentially
# expressed miRNAs to seed-matched, coherent, and independently
# validated miRNA–mRNA interactions.
#
# All values are calculated from nf-Sarcopipe output tables.
# No candidate counts or interaction counts are hardcoded.
#
# Inputs
# ------
# 1) Final BrumiR de novo catalog
# 2) Final unique miRDeep2 known catalog
# 3) BrumiR seed-matched pairs
# 4) miRDeep2 seed-matched pairs
# 5) Coherent miRNA–mRNA pairs
# 6) Final validated Cytoscape edges
# 7) Output TSV
# 8) Output PNG
# 9) Output PDF
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript Figure5B_seed_target_prioritization.R ",
      "<brumir_final.tsv> <mirdeep2_final.tsv> ",
      "<brumir_seed.tsv> <mirdeep2_seed.tsv> ",
      "<coherent.tsv> <validated_edges.tsv> ",
      "<out_tsv> <out_png> <out_pdf>\n"
    )
  )
}

brumir_final_file <- args[1]
mirdeep2_final_file <- args[2]
brumir_seed_file <- args[3]
mirdeep2_seed_file <- args[4]
coherent_file <- args[5]
validated_edges_file <- args[6]
out_tsv <- args[7]
out_png <- args[8]
out_pdf <- args[9]

input_files <- c(
  brumir_final_file,
  mirdeep2_final_file,
  brumir_seed_file,
  mirdeep2_seed_file,
  coherent_file,
  validated_edges_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing input file(s): ",
    paste(missing_files, collapse = ", ")
  )
}

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

read_table <- function(path) {
  read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
}

require_columns <- function(df, required, label) {
  missing <- setdiff(required, colnames(df))

  if (length(missing) > 0) {
    stop(
      label,
      " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

brumir_final <- read_table(brumir_final_file)
mirdeep2_final <- read_table(mirdeep2_final_file)
brumir_seed <- read_table(brumir_seed_file)
mirdeep2_seed <- read_table(mirdeep2_seed_file)
coherent <- read_table(coherent_file)
validated_edges <- read_table(validated_edges_file)

require_columns(
  brumir_final,
  c("candidate"),
  "BrumiR final catalog"
)

require_columns(
  mirdeep2_final,
  c("representative_miRNA"),
  "miRDeep2 final catalog"
)

require_columns(
  brumir_seed,
  c("renamed_miRNA", "gene_symbol"),
  "BrumiR seed-match table"
)

require_columns(
  mirdeep2_seed,
  c("renamed_miRNA", "gene_symbol"),
  "miRDeep2 seed-match table"
)

require_columns(
  coherent,
  c("source", "renamed_miRNA", "gene_symbol"),
  "Coherent-pair table"
)

require_columns(
  validated_edges,
  c("miRNA_source", "source_node", "target_node", "validated_by"),
  "Validated Cytoscape edge table"
)

count_rows <- function(df) {
  nrow(df)
}

count_unique_nonempty <- function(values) {
  values <- trimws(as.character(values))
  length(unique(values[!is.na(values) & values != ""]))
}

count_source_rows <- function(df, source_column, source_value) {
  sum(df[[source_column]] == source_value, na.rm = TRUE)
}

summary_df <- tibble(
  Step = c(
    "Final miRNAs\nevaluated",
    "miRNAs with\nseed match",
    "Genes with\nseed match",
    "Raw seed-matched\npairs",
    "Coherent\npairs",
    "Validated\npairs"
  ),

  BrumiR = c(
    count_unique_nonempty(brumir_final$candidate),
    count_unique_nonempty(brumir_seed$renamed_miRNA),
    count_unique_nonempty(brumir_seed$gene_symbol),
    count_rows(brumir_seed),
    count_source_rows(coherent, "source", "BrumiR"),
    count_source_rows(validated_edges, "miRNA_source", "BrumiR")
  ),

  miRDeep2 = c(
    count_unique_nonempty(mirdeep2_final$representative_miRNA),
    count_unique_nonempty(mirdeep2_seed$renamed_miRNA),
    count_unique_nonempty(mirdeep2_seed$gene_symbol),
    count_rows(mirdeep2_seed),
    count_source_rows(coherent, "source", "miRDeep2"),
    count_source_rows(validated_edges, "miRNA_source", "miRDeep2")
  )
)

step_order <- summary_df$Step

plot_df <- summary_df %>%
  pivot_longer(
    cols = c("BrumiR", "miRDeep2"),
    names_to = "Algorithm",
    values_to = "Count"
  ) %>%
  mutate(
    Step = factor(Step, levels = step_order),
    Algorithm = factor(
      Algorithm,
      levels = c("BrumiR", "miRDeep2")
    )
  )

write_tsv(
  summary_df,
  out_tsv
)

p <- ggplot(
  plot_df,
  aes(
    x = Step,
    y = Count,
    fill = Algorithm
  )
) +
  geom_col(
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = Count),
    position = position_dodge(width = 0.82),
    vjust = -0.35,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_y_log10(
    breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000),
    labels = scales::label_number(big.mark = ","),
    expand = expansion(mult = c(0.03, 0.16))
  ) +
  scale_fill_manual(
    values = c(
      "BrumiR" = "#7A6A8A",
      "miRDeep2" = "#1B9E77"
    )
  ) +
  labs(
    title = "Progressive prioritization of miRNA–mRNA interactions",
    subtitle = paste0(
      "Canonical seed matching, inverse-expression filtering, ",
      "and independent validation"
    ),
    x = NULL,
    y = "Count (log10 scale)",
    fill = "Discovery method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      face = "bold"
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = out_png,
  plot = p,
  width = 10.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 10.5,
  height = 6.5,
  bg = "white"
)

cat("===== FIGURE 5B =====\n")
print(summary_df)
cat("Written:", out_tsv, "\n")
cat("Written:", out_png, "\n")
cat("Written:", out_pdf, "\n")
