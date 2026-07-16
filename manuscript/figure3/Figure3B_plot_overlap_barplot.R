#!/usr/bin/env Rscript

# ============================================================
# plot_figure2B_overlap_barplot.R
# ============================================================
#
# Description:
# This script generates the Figure 2B barplot summarizing overlap between
# BrumiR core sequences and miRDeep2 known core sequences.
#
# It combines:
# - exact overlap counts from overlap_counts.tsv
# - near-match overlap counts inferred from a CD-HIT-2D .clstr file
#
# The script also writes a summary TSV containing:
# - total BrumiR sequences
# - total miRDeep2 sequences
# - exact shared sequences
# - near-match shared sequences
# - BrumiR-only exact sequences
# - miRDeep2-only exact sequences
# - sequences added by near-match clustering
#
# Inputs:
#   1. exact_counts.tsv
#   2. brumir_fasta         (kept for interface consistency; currently not used)
#   3. mirdeep2_fasta       (kept for interface consistency; currently not used)
#   4. overlap_similar.clstr
#   5. out_png
#   6. plot title
#   7. optional CD-HIT identity (default: 0.95)
#
# Outputs:
# - figure2B overlap barplot PNG
# - overlap_figure2B_summary.tsv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  cat(
    "Usage:\n",
    "  Rscript plot_figure2B_overlap_barplot.R <exact_counts.tsv> <brumir_fasta> <mirdeep2_fasta> <nearmatch.clstr> <out_png> <title> [cdhit_identity]\n\n",
    "Example:\n",
    "  Rscript plot_figure2B_overlap_barplot.R overlap_counts.tsv brumir.fa mirdeep.fa overlap_similar_0.95.clstr out.png \"Title\" 0.95\n",
    sep = ""
  )
  quit(status = 1)
}

exact_counts   <- args[1]
br_fa          <- args[2]
md_fa          <- args[3]
clstr_path     <- args[4]
out_png        <- args[5]
plot_title     <- args[6]
cdhit_identity <- ifelse(length(args) >= 7, as.numeric(args[7]), 0.95)

# ----------------------------
# Helpers
# ----------------------------

read_tsv1 <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
}

# Parse CD-HIT .clstr and count clusters that contain both:
# - BrumiR headers: cluster_*
# - miRDeep2/miRBase headers: typically hsa-* or anything not cluster_*
count_shared_from_clstr <- function(clstr_path) {
  x <- readLines(clstr_path, warn = FALSE)

  shared <- 0
  has_br <- FALSE
  has_md <- FALSE

  flush_cluster <- function() {
    if (has_br && has_md) {
      shared <<- shared + 1
    }
  }

  for (line in x) {
    if (startsWith(line, ">Cluster")) {
      flush_cluster()
      has_br <- FALSE
      has_md <- FALSE
      next
    }

    m <- regexpr(">[^ .]+", line)
    if (m[1] != -1) {
      name <- substr(line, m[1] + 1, m[1] + attr(m, "match.length") - 1)
      if (startsWith(name, "cluster_")) {
        has_br <- TRUE
      } else {
        has_md <- TRUE
      }
    }
  }

  flush_cluster()
  shared
}

# ----------------------------
# 1) Read exact counts
# ----------------------------

ex <- read_tsv1(exact_counts)

if (nrow(ex) < 1) {
  stop("exact_counts.tsv has no rows: ", exact_counts)
}

brN_exact    <- ex$brumir_n[1]
mdN_exact    <- ex$mirdeep2_n[1]
shared_exact <- ex$shared_n[1]

# ----------------------------
# 2) Read near-match overlap from .clstr
# ----------------------------

if (!file.exists(clstr_path)) {
  stop("Expected .clstr not found: ", clstr_path)
}

shared_95 <- count_shared_from_clstr(clstr_path)

br_only_exact <- brN_exact - shared_exact
md_only_exact <- mdN_exact - shared_exact
add_gain      <- max(shared_95 - shared_exact, 0)

# Summary TSV
sumdf <- data.frame(
  brumir_n = brN_exact,
  mirdeep2_n = mdN_exact,
  shared_exact = shared_exact,
  shared_nearmatch = shared_95,
  brumir_only_exact = br_only_exact,
  mirdeep2_only_exact = md_only_exact,
  added_by_nearmatch = add_gain,
  stringsAsFactors = FALSE
)

write.table(
  sumdf,
  "overlap_figure2B_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------
# 3) Plot
# ----------------------------

col_brumir  <- "#6A3D9A"
col_overlap <- "#BDBDBD"
col_added   <- "#1B9E77"
bar_cols    <- c(col_brumir, col_overlap, col_added)

png(out_png, width = 1800, height = 1200, res = 250, bg = "white")

vals <- c(br_only_exact, shared_exact, add_gain)
labels <- c(
  "BrumiR-only\n(exact)",
  "Overlap\n(exact 100%)",
  sprintf("Added by\nnear-match %.2f", cdhit_identity)
)

par(mar = c(7, 5, 4, 2) + 0.1)

bp <- barplot(
  vals,
  names.arg = labels,
  col = bar_cols,
  border = "black",
  las = 1,
  main = plot_title,
  ylab = "Number of sequences",
  cex.names = 0.9
)

text(x = bp, y = vals, labels = vals, pos = 3, cex = 0.9)

mtext(
  sprintf(
    "BrumiR total=%d | miRDeep2 total=%d | shared_exact=%d | shared_%.2f=%d",
    brN_exact, mdN_exact, shared_exact, cdhit_identity, shared_95
  ),
  side = 1,
  line = 4.5,
  cex = 0.9
)

dev.off()
