#!/usr/bin/env Rscript

############################################################
# Figure 3B – BrumiR2Reference support by group
#
# Description:
# This script plots the number of supported and not-supported
# BrumiR core candidates across biological groups using a
# stacked barplot.
#
# Inputs:
# - TSV file containing columns: group, metric, count
# - Output PNG file
# - Plot title
#
# Outputs:
# - A stacked barplot PNG showing supported vs not-supported
#   candidates by group
#
# Usage:
# Rscript Figure3B_plot_BrumiR2Reference.R \
#   <input.tsv> \
#   <output.png> \
#   <title>
#
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat(
    "Usage:\n",
    "  Rscript Figure3B_plot_BrumiR2Reference.R <input.tsv> <output.png> <title>\n",
    sep = ""
  )
  quit(status = 1)
}

in_tsv <- args[1]
out_png <- args[2]
ttl <- args[3]

if (!file.exists(in_tsv)) {
  stop("Input TSV not found: ", in_tsv, call. = FALSE)
}

out_dir <- dirname(out_png)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

df <- read.table(
  in_tsv,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

required_cols <- c("group", "metric", "count")
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Input TSV is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

sub <- df[df$metric %in% c("supported", "not_supported"), , drop = FALSE]

if (nrow(sub) == 0) {
  stop(
    "No rows found with metric values 'supported' or 'not_supported'.",
    call. = FALSE
  )
}

sub$group <- factor(sub$group, levels = c("athlete", "sedentary", "shared"))
sub$metric <- factor(sub$metric, levels = c("supported", "not_supported"))

png(out_png, width = 1400, height = 900, res = 150)

par(mar = c(11, 5, 4, 2), xpd = NA)

mat <- tapply(sub$count, list(sub$metric, sub$group), sum)
mat[is.na(mat)] <- 0

bp <- barplot(
  mat,
  beside = FALSE,
  col = c("grey40", "grey85"),
  ylab = "Number of BrumiR core candidates",
  main = ttl,
  ylim = c(0, max(colSums(mat)) * 1.15)
)

tot <- colSums(mat)

text(
  bp,
  tot + max(tot) * 0.03,
  labels = tot,
  cex = 1
)

text(
  bp,
  mat[1, ] / 2,
  labels = mat[1, ],
  col = "white",
  cex = 0.9
)

text(
  bp,
  mat[1, ] + mat[2, ] / 2,
  labels = mat[2, ],
  cex = 0.9
)

legend(
  "bottom",
  inset = c(0, -0.35),
  legend = c("supported", "not_supported"),
  fill = c("grey40", "grey85"),
  horiz = TRUE,
  bty = "n"
)

dev.off()