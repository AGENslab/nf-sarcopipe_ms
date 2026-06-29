#!/usr/bin/env Rscript

############################################################
# Figure S4A – Union KEGG and GO enrichment
#
# Description:
# Generates GO and KEGG bubble plots for the union of coherent
# miRNA–mRNA target genes predicted from BrumiR and miRDeep2.
#
# Inputs:
# 1) union_coherent_GO.tsv
# 2) union_coherent_KEGG.tsv
#
# Outputs:
# 1) FigureS4A_union_GO.png
# 2) FigureS4A_union_KEGG.png
#
# Usage:
# Rscript FigureS4A_union_KEGG_enrichment.R \
#   <union_coherent_GO.tsv> \
#   <union_coherent_KEGG.tsv> \
#   <FigureS4A_union_GO.png> \
#   <FigureS4A_union_KEGG.png>
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    "Usage: Rscript FigureS4A_union_KEGG_enrichment.R ",
    "<union_GO.tsv> <union_KEGG.tsv> <out_GO.png> <out_KEGG.png>"
  )
}

go_file <- args[1]
kegg_file <- args[2]
go_plot <- args[3]
kegg_plot <- args[4]

if (!file.exists(go_file)) {
  stop("GO enrichment file not found: ", go_file)
}

if (!file.exists(kegg_file)) {
  stop("KEGG enrichment file not found: ", kegg_file)
}

for (outdir in unique(c(dirname(go_plot), dirname(kegg_plot)))) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
}

plot_bubble <- function(df, title, outfile) {

  if (nrow(df) == 0) {
    message("No enrichment results for: ", outfile)
    return(invisible(NULL))
  }

  required_cols <- c(
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )

  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
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

  df2$Description <- factor(
    df2$Description,
    levels = df2$Description
  )

  p <- ggplot(
    df2,
    aes(
      x = GeneRatio_num,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#C94C4C",
      high = "#3B6FB6",
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

  ggsave(
    outfile,
    p,
    width = 8,
    height = 5,
    dpi = 300
  )

  message("Written: ", outfile)
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
