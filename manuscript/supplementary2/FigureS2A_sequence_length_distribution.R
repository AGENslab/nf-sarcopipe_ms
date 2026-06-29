#!/usr/bin/env Rscript

############################################################
# Figure S2A – Sequence length distribution
#
# Description:
# Compares the length distribution of BrumiR core candidate
# sequences and miRDeep2 known core sequences.
#
# The script reads two FASTA files, computes sequence lengths,
# filters them to a user-defined interval, writes a summary TSV,
# and generates a grouped barplot.
#
# Inputs:
# 1) BrumiR core FASTA
# 2) miRDeep2 core FASTA
# 3) output PNG
# 4) output TSV
# 5) plot title
# 6) minimum sequence length
# 7) maximum sequence length
#
# Outputs:
# - A TSV summarizing sequence counts per length
# - A PNG grouped barplot of length distributions
#
# Usage:
# Rscript FigureS2A_sequence_length_distribution.R \
#   <brumir_fasta> <mirdeep2_fasta> <out_png> <out_tsv> \
#   <title> <min_len> <max_len>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  cat(
    "Usage:\n",
    "  Rscript FigureS2A_sequence_length_distribution.R ",
    "<brumir_fasta> <mirdeep2_fasta> <out_png> <out_tsv> ",
    "<title> <min_len> <max_len>\n",
    sep = ""
  )
  quit(status = 1)
}

br_fa <- args[1]
md_fa <- args[2]
out_png <- args[3]
out_tsv <- args[4]
title <- args[5]
min_len <- as.integer(args[6])
max_len <- as.integer(args[7])

if (!file.exists(br_fa)) {
  stop("Input BrumiR FASTA not found: ", br_fa)
}

if (!file.exists(md_fa)) {
  stop("Input miRDeep2 FASTA not found: ", md_fa)
}

if (is.na(min_len) || is.na(max_len)) {
  stop("Minimum and maximum sequence lengths must be integers.")
}

if (min_len > max_len) {
  stop("Minimum sequence length cannot be greater than maximum sequence length.")
}

out_png_dir <- dirname(out_png)
out_tsv_dir <- dirname(out_tsv)

if (!dir.exists(out_png_dir)) {
  dir.create(out_png_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(out_tsv_dir)) {
  dir.create(out_tsv_dir, recursive = TRUE, showWarnings = FALSE)
}

read_fasta_lengths <- function(path) {
  fasta_lines <- readLines(path, warn = FALSE)
  lengths <- integer(0)
  current_sequence <- ""

  for (line in fasta_lines) {
    if (startsWith(line, ">")) {
      if (nchar(current_sequence) > 0) {
        lengths <- c(lengths, nchar(current_sequence))
        current_sequence <- ""
      }
    } else {
      current_sequence <- paste0(
        current_sequence,
        gsub("[^ACGTUacgtu]", "", line)
      )
    }
  }

  if (nchar(current_sequence) > 0) {
    lengths <- c(lengths, nchar(current_sequence))
  }

  lengths
}

br_len <- read_fasta_lengths(br_fa)
md_len <- read_fasta_lengths(md_fa)

br_len <- br_len[br_len >= min_len & br_len <= max_len]
md_len <- md_len[md_len >= min_len & md_len <= max_len]

lengths <- min_len:max_len

count_vec <- function(values, lengths) {
  tab <- table(factor(values, levels = lengths))
  as.integer(tab)
}

df <- data.frame(
  length = lengths,
  BrumiR_candidates = count_vec(br_len, lengths),
  miRDeep2_known = count_vec(md_len, lengths),
  stringsAsFactors = FALSE
)

write.table(
  df,
  file = out_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -------- Plot: legend in a dedicated bottom panel --------

png(out_png, width = 2400, height = 1800, res = 250, bg = "white")

layout(matrix(c(1, 2), nrow = 2), heights = c(4.5, 1.5))

col_brumir <- "#5E3C99"
col_mirdeep <- "#1B9E77"

# Panel 1: grouped bars
par(mar = c(5.5, 6, 4, 2) + 0.1)

mat <- rbind(df$BrumiR_candidates, df$miRDeep2_known)
colnames(mat) <- as.character(df$length)

barplot(
  mat,
  beside = TRUE,
  names.arg = df$length,
  col = c(col_brumir, col_mirdeep),
  border = "black",
  las = 1,
  main = title,
  xlab = "Sequence length (nt)",
  ylab = "Number of core sequences"
)

# Panel 2: legend only
par(mar = c(0.5, 2, 0.5, 2))
plot.new()

legend(
  "center",
  legend = c(
    "BrumiR-RF high-confidence candidate miRNAs (core set; CD-HIT clustered)",
    "miRDeep2 known miRNAs (miRBase-annotated; core set)"
  ),
  fill = c(col_brumir, col_mirdeep),
  border = "black",
  bty = "n",
  cex = 1.1
)

dev.off()

cat("Length distribution plot written to:", out_png, "\n")
cat("Length distribution table written to:", out_tsv, "\n")