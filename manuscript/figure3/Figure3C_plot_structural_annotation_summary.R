#!/usr/bin/env Rscript

############################################################
# Figure 3C – miRGeneDB annotation of miRNA candidate sets
#
# Input TSV columns:
#   feature_set
#   category
#   count
#
# Usage:
# Rscript Figure3C_plot_structural_annotation_summary.R \
#   Figure3C_annotation_summary.tsv \
#   Figure3C_annotation_summary
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript Figure3C_plot_structural_annotation_summary.R ",
    "<summary.tsv> <out_prefix>",
    call. = FALSE
  )
}

input_file <- args[1]
out_prefix <- args[2]

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

out_png <- paste0(out_prefix, ".png")
out_pdf <- paste0(out_prefix, ".pdf")

out_dir <- dirname(out_png)

if (!dir.exists(out_dir) && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE)
}

df <- read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required <- c("feature_set", "category", "count")
missing <- setdiff(required, colnames(df))

if (length(missing) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

df$feature_set <- factor(
  df$feature_set,
  levels = c(
    "miRDeep2 core",
    "BrumiR-RF core",
    "BrumiR-RF supported"
  )
)

df$category <- factor(
  df$category,
  levels = c(
    "Not annotated in MirGeneDB",
    "Unique annotated",
    "Shared annotated"
  )
)

category_count <- function(sub, category_name) {

    x <- sub$count[sub$category == category_name]

    if(length(x)==0) return(0)

    x[1]
}

side_labels <- do.call(
    rbind,
    lapply(levels(df$feature_set), function(s){

        sub <- df[df$feature_set==s,]

        shared <- category_count(sub,"Shared annotated")
        unique <- category_count(sub,"Unique annotated")
        notann <- category_count(sub,"Not annotated in MirGeneDB")

        total <- shared + unique + notann

        data.frame(
            feature_set=s,
            total=total,
            label=paste0(
                "n=",total,
                "   shared=",shared,
                "   unique=",unique,
                "   not=",notann
            )
        )

    })
)

side_labels$feature_set <- factor(
    side_labels$feature_set,
    levels=levels(df$feature_set)

)

p <- ggplot(
  df,
  aes(
    x = feature_set,
    y = count,
    fill = category
  )
) +
  geom_col(
    width = 0.65,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = side_labels,
    aes(
        x = feature_set,
        y = total,
        label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    fontface = "bold",
    size = 3.5
) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Not annotated in MirGeneDB" = "#BFBFBF",
      "Unique annotated" = "#DD8452",
      "Shared annotated" = "#4C72B0"
    ),
    breaks = c(
      "Not annotated in MirGeneDB",
      "Unique annotated",
      "Shared annotated"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.55))
  ) +
  labs(
    title = "miRNA annotation against MirGeneDB",
    subtitle = paste(
      "≥98% identity, ≥75% query coverage, 0 mismatches;",
      "BrumiR-RF supported candidates passed BrumiR2Reference"
    ),
    x = NULL,
    y = "Number of miRNA candidates",
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.margin = margin(10, 35, 10, 10)
  )

ggsave(
  out_png,
  p,
  width = 9,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_pdf,
  p,
  width = 9,
  height = 5.5,
  bg = "white"
)

cat("Plot written to:", out_png, "\n")
cat("Plot written to:", out_pdf, "\n")
