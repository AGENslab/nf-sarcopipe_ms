#!/usr/bin/env bash

# Filter top BrumiR2Reference-supported candidates for Supplementary Figure S2B.
# Criteria used for manual table assembly:
#   - de novo exact match = 1
#   - no known exact match
#   - ranked by BrumiR2Reference/RNAfold MFE and adjusted p-value when available

INPUT="brumir2reference_candidates.tsv"
OUTPUT="FigureS2B_top4_BrumiR2Reference_candidates.tsv"

awk 'BEGIN{FS=OFS="\t"} NR==1 || ($0 ~ /cluster_384|cluster_389|cluster_428|cluster_487/)' "$INPUT" > "$OUTPUT"
