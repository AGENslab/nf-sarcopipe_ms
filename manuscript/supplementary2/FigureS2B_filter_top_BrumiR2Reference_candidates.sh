#!/usr/bin/env bash

###############################################################################
# Figure S2B – Filter top BrumiR2Reference-supported candidates
#
# Description:
# Selects the manually curated BrumiR2Reference-supported candidate miRNAs
# used in Supplementary Figure S2B.
#
# Selection criteria (unchanged):
#   - de novo exact match = 1
#   - no known exact match
#   - ranked by BrumiR2Reference/RNAfold MFE and adjusted p-value
#   - retain clusters:
#       cluster_384
#       cluster_389
#       cluster_428
#       cluster_487
#
# Inputs:
#   1) BrumiR2Reference candidate table (TSV)
#
# Outputs:
#   - FigureS2B_top4_BrumiR2Reference_candidates.tsv
#
# Usage:
#   bash FigureS2B_filter_top_BrumiR2Reference_candidates.sh \
#       input.tsv output.tsv
###############################################################################

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "  bash FigureS2B_filter_top_BrumiR2Reference_candidates.sh <input.tsv> <output.tsv>"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: Input file not found: $INPUT"
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT")"

awk '
BEGIN { FS=OFS="\t" }
NR==1 ||
$0 ~ /cluster_384|cluster_389|cluster_428|cluster_487/
' "$INPUT" > "$OUTPUT"

echo "Filtered candidate table written to: $OUTPUT"
