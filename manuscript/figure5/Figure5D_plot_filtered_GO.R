#!/usr/bin/env Rscript

############################################################
# Figure 5D – Filtered GO enrichment plot
#
# Description:
# This script generates a filtered and publication-ready GO
# bubble plot for the union of coherent miRNA–mRNA predicted
# target genes. Only biologically interpretable GO terms relevant
# to the study narrative are retained for visualization.
#
# Inputs:
# - GO enrichment TSV file
# - Output PNG file
#
# Outputs:
# - Filtered GO enrichment bubble plot PNG
#
# Usage:
# Rscript Figure5D_plot_filtered_GO.R \
#   union_coherent_GO.tsv \
#   Figure5D_filtered_union_GO.png
#
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript Figure5D_plot_filtered_GO.R ",
    "<union_coherent_GO.tsv> <output.png>",
    call. = FALSE
  )
}

go_file <- args[1]
go_plot <- args[2]

if (!file.exists(go_file)) {
  stop("Input GO TSV not found: ", go_file, call. = FALSE)
}

out_dir <- dirname(go_plot)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# -------------------------
# LOAD DATA
# -------------------------
go <- read_tsv(go_file, show_col_types = FALSE)

required_cols <- c("Description", "GeneRatio", "Count", "p.adjust")
missing_cols <- setdiff(required_cols, colnames(go))

if (length(missing_cols) > 0) {
  stop(
    "Input GO TSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# -------------------------
# FILTER TERMS
# -------------------------
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
    GeneRatio_num = sapply(
      strsplit(GeneRatio, "/"),
      function(x) as.numeric(x[1]) / as.numeric(x[2])
    ),
    Description_wrapped = str_wrap(Description, width = 38)
  ) %>%
  arrange(GeneRatio_num)

go2$Description_wrapped <- factor(
  go2$Description_wrapped,
  levels = go2$Description_wrapped
)

# -------------------------
# PLOT
# -------------------------
p <- ggplot(
  go2,
  aes(x = GeneRatio_num, y = Description_wrapped, size = Count, color = p.adjust)
) +
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

# -------------------------
# SAVE
# -------------------------
ggsave(go_plot, p, width = 11, height = 5.5, dpi = 300)

cat("Filtered GO plot written to:", go_plot, "\n")
