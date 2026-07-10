#!/usr/bin/env bash

# ============================================================
# suppl_table_2_count_unique_miRDeep2_entities.sh
# ============================================================
#
# Description
# -----------
# Quantifies the non-redundant annotated and novel miRNA
# entities identified by miRDeep2 across multiple sequencing
# libraries.
#
# Annotated miRNAs are counted using unique miRBase identifiers
# extracted from the first column of the combined known-miRNA
# count matrix.
#
# Novel miRNAs are counted using non-redundant consensus mature
# sequences extracted from the per-library
# novel_miRNAs_*.tsv files. This avoids counting repeated
# predictions of the same mature sequence across libraries.
#
# The script does not modify the original miRDeep2 results.
#
# Definitions
# -----------
# Unique annotated miRNA:
#   A distinct miRBase miRNA identifier represented in the
#   combined miRDeep2 known-miRNA count matrix.
#
# Unique novel mature sequence:
#   A distinct nucleotide sequence reported in the
#   "consensus mature sequence" column across all per-library
#   novel-miRNA result tables.
#
# Usage
# -----
# bash Figure2A_count_unique_miRDeep2_entities.sh \
#     --known-dir /path/to/miRDeep2_known \
#     --novel-dir /path/to/miRDeep2_novel \
#     --out /path/to/miRDeep2_unique_counts.tsv
#
# Required arguments
# ------------------
# --known-dir DIR
#     Directory containing:
#       miRDeep2.known.counts.tsv
#
# --novel-dir DIR
#     Directory containing:
#       novel_miRNAs_*.tsv
#
# Optional arguments
# ------------------
# --out FILE
#     Output summary table.
#     Default: miRDeep2_unique_counts.tsv in the current
#     working directory.
#
# --help
#     Display the help message.
#
# Output
# ------
# A tab-separated table containing:
#
# category    uniqueness_definition              unique_count
# known       unique_miRBase_ID                   N
# novel       unique_consensus_mature_sequence    N
#
# Manuscript use
# --------------
# The resulting values summarize the non-redundant miRDeep2
# entities reported in the Results section associated with
# Figure 2.
#
# Requirements
# ------------
# bash, awk, sort, wc, find
#
# ============================================================

set -euo pipefail

PROGRAM_NAME="$(basename "$0")"

KNOWN_DIR=""
NOVEL_DIR=""
OUTFILE="miRDeep2_unique_counts.tsv"


usage() {
    cat <<EOF
Usage:
  bash ${PROGRAM_NAME} \\
      --known-dir DIR \\
      --novel-dir DIR \\
      [--out FILE]

Required:
  --known-dir DIR   Directory containing miRDeep2.known.counts.tsv
  --novel-dir DIR   Directory containing novel_miRNAs_*.tsv

Optional:
  --out FILE        Output TSV file
                    Default: miRDeep2_unique_counts.tsv
  --help            Show this help message
EOF
}


error_exit() {
    echo "ERROR: $*" >&2
    exit 1
}


while [[ $# -gt 0 ]]; do
    case "$1" in
        --known-dir)
            [[ $# -ge 2 ]] || error_exit "--known-dir requires a directory."
            KNOWN_DIR="$2"
            shift 2
            ;;

        --novel-dir)
            [[ $# -ge 2 ]] || error_exit "--novel-dir requires a directory."
            NOVEL_DIR="$2"
            shift 2
            ;;

        --out)
            [[ $# -ge 2 ]] || error_exit "--out requires a file path."
            OUTFILE="$2"
            shift 2
            ;;

        --help|-h)
            usage
            exit 0
            ;;

        *)
            error_exit "Unknown argument: $1. Use --help for usage."
            ;;
    esac
done


# ------------------------------------------------------------
# Validate arguments and input files
# ------------------------------------------------------------

[[ -n "$KNOWN_DIR" ]] ||
    error_exit "Missing required argument: --known-dir"

[[ -n "$NOVEL_DIR" ]] ||
    error_exit "Missing required argument: --novel-dir"

[[ -d "$KNOWN_DIR" ]] ||
    error_exit "Known-miRNA directory not found: $KNOWN_DIR"

[[ -d "$NOVEL_DIR" ]] ||
    error_exit "Novel-miRNA directory not found: $NOVEL_DIR"

KNOWN_COUNTS="${KNOWN_DIR%/}/miRDeep2.known.counts.tsv"

[[ -s "$KNOWN_COUNTS" ]] ||
    error_exit "Known-miRNA count matrix not found or empty: $KNOWN_COUNTS"


mapfile -d '' NOVEL_FILES < <(
    find "$NOVEL_DIR" \
        -maxdepth 1 \
        -type f \
        -name 'novel_miRNAs_*.tsv' \
        -print0 |
    sort -z
)

[[ ${#NOVEL_FILES[@]} -gt 0 ]] ||
    error_exit \
        "No files matching novel_miRNAs_*.tsv were found in: $NOVEL_DIR"


# Create the output directory when the output includes a path.
OUTDIR="$(dirname "$OUTFILE")"
mkdir -p "$OUTDIR"


# Temporary workspace removed automatically when the script exits.
TMPDIR_SCRIPT="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_SCRIPT"' EXIT

KNOWN_IDS="${TMPDIR_SCRIPT}/known_unique_ids.txt"
NOVEL_SEQUENCES="${TMPDIR_SCRIPT}/novel_unique_sequences.txt"


# ------------------------------------------------------------
# Count unique annotated miRNAs
# ------------------------------------------------------------
#
# The first column of miRDeep2.known.counts.tsv contains the
# annotated miRBase identifier. The header is excluded, blank
# values are ignored, and identifiers are deduplicated.
# ------------------------------------------------------------

awk -F '\t' '
    NR == 1 {
        next
    }

    {
        id = $1
        sub(/\r$/, "", id)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", id)

        if (id != "") {
            print id
        }
    }
' "$KNOWN_COUNTS" |
LC_ALL=C sort -u \
> "$KNOWN_IDS"

KNOWN_UNIQUE_COUNT="$(wc -l < "$KNOWN_IDS" | tr -d '[:space:]')"


# ------------------------------------------------------------
# Count unique novel mature sequences
# ------------------------------------------------------------
#
# For every input file, the column is located by its header
# name rather than by a fixed column number. Sequences are:
#
#   - converted to uppercase;
#   - stripped of whitespace and carriage returns;
#   - filtered to retain nucleotide sequences only;
#   - deduplicated across all sequencing libraries.
#
# The column lookup is repeated for every file because FNR
# returns to 1 when awk begins reading a new file.
# ------------------------------------------------------------

awk -F '\t' '
    FNR == 1 {
        mature_col = 0

        for (i = 1; i <= NF; i++) {
            header = $i
            sub(/\r$/, "", header)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", header)

            if (tolower(header) == "consensus mature sequence") {
                mature_col = i
                break
            }
        }

        if (mature_col == 0) {
            print "ERROR: Column \"consensus mature sequence\" not found in " \
                  FILENAME > "/dev/stderr"
            exit 2
        }

        next
    }

    {
        sequence = $mature_col
        sub(/\r$/, "", sequence)
        gsub(/[[:space:]]/, "", sequence)
        sequence = toupper(sequence)

        # Accept RNA or DNA notation and exclude missing values,
        # punctuation, or malformed entries.
        if (sequence ~ /^[ACGTUN]+$/) {
            print sequence
        }
    }
' "${NOVEL_FILES[@]}" |
LC_ALL=C sort -u \
> "$NOVEL_SEQUENCES"

NOVEL_UNIQUE_COUNT="$(
    wc -l < "$NOVEL_SEQUENCES" |
    tr -d '[:space:]'
)"


# ------------------------------------------------------------
# Validate calculated counts
# ------------------------------------------------------------

[[ "$KNOWN_UNIQUE_COUNT" -gt 0 ]] ||
    error_exit "The annotated-miRNA unique count is zero."

[[ "$NOVEL_UNIQUE_COUNT" -gt 0 ]] ||
    error_exit "The novel mature-sequence unique count is zero."


# ------------------------------------------------------------
# Write summary table
# ------------------------------------------------------------

{
    printf "category\tuniqueness_definition\tunique_count\n"
    printf "known\tunique_miRBase_ID\t%s\n" "$KNOWN_UNIQUE_COUNT"
    printf "novel\tunique_consensus_mature_sequence\t%s\n" \
        "$NOVEL_UNIQUE_COUNT"
} > "$OUTFILE"


# ------------------------------------------------------------
# Print execution summary
# ------------------------------------------------------------

printf "\n"
printf "============================================================\n"
printf "miRDeep2 non-redundant entity summary\n"
printf "============================================================\n"
printf "Known-miRNA matrix:        %s\n" "$KNOWN_COUNTS"
printf "Novel input directory:    %s\n" "$NOVEL_DIR"
printf "Novel files processed:    %s\n" "${#NOVEL_FILES[@]}"
printf "Unique annotated miRNAs:  %s\n" "$KNOWN_UNIQUE_COUNT"
printf "Unique novel sequences:   %s\n" "$NOVEL_UNIQUE_COUNT"
printf "Output table:             %s\n" "$OUTFILE"
printf "============================================================\n"
printf "\n"
