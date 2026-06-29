#!/usr/bin/env python3

"""
Figure 3C – Shared miRGeneDB table

Description
-----------
This script builds a table of shared miRNAs between BrumiR and
miRDeep2 based on miRGeneDB annotation. It loads BLAST-like
alignment tables, filters hits present in a shared miRGeneDB ID
list, keeps the best hit per subject based on bitscore, and merges
BrumiR and miRDeep2 annotations.

Inputs
------
--brumir
    BrumiR vs miRGeneDB TSV file.

--mirdeep
    miRDeep2 vs miRGeneDB TSV file.

--shared
    Text file containing shared miRGeneDB IDs, one per line.

--out
    Output TSV file.

Outputs
-------
Shared miRNA table TSV with BrumiR and miRDeep2 miRGeneDB matches.

Usage
-----
python3 Figure3C_build_shared_table.py \\
  --brumir brumir_core220_vs_miRGeneDB.tsv \\
  --mirdeep md_core_vs_miRGeneDB.tsv \\
  --shared shared_miRGeneDB_hits.txt \\
  --out Supplementary_Table_S4_shared_miRNAs.tsv
"""

import argparse
from pathlib import Path

import pandas as pd


BLAST_COLUMNS = [
    "query",
    "subject",
    "identity",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build shared miRGeneDB annotation table for BrumiR and miRDeep2."
    )

    parser.add_argument(
        "--brumir",
        required=True,
        help="BrumiR vs miRGeneDB TSV file.",
    )
    parser.add_argument(
        "--mirdeep",
        required=True,
        help="miRDeep2 vs miRGeneDB TSV file.",
    )
    parser.add_argument(
        "--shared",
        required=True,
        help="Text file with shared miRGeneDB IDs, one per line.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file.",
    )

    return parser.parse_args()


def main() -> None:
    """Build shared BrumiR/miRDeep2 miRGeneDB table."""
    args = parse_args()

    brumir_file = validate_file(args.brumir, "BrumiR miRGeneDB TSV")
    mirdeep_file = validate_file(args.mirdeep, "miRDeep2 miRGeneDB TSV")
    shared_file = validate_file(args.shared, "Shared miRGeneDB ID file")
    out_file = Path(args.out)

    out_file.parent.mkdir(parents=True, exist_ok=True)

    shared_ids = set(shared_file.read_text().splitlines())

    brumir_hits = pd.read_csv(brumir_file, sep="\t", names=BLAST_COLUMNS)
    mirdeep_hits = pd.read_csv(mirdeep_file, sep="\t", names=BLAST_COLUMNS)

    brumir_shared = brumir_hits[brumir_hits["subject"].isin(shared_ids)].copy()
    mirdeep_shared = mirdeep_hits[mirdeep_hits["subject"].isin(shared_ids)].copy()

    brumir_shared = (
        brumir_shared
        .sort_values("bitscore", ascending=False)
        .drop_duplicates("subject")
    )

    mirdeep_shared = (
        mirdeep_shared
        .sort_values("bitscore", ascending=False)
        .drop_duplicates("subject")
    )

    merged = pd.merge(
        brumir_shared,
        mirdeep_shared,
        on="subject",
        suffixes=("_BrumiR", "_miRDeep2"),
    )

    final = merged[[
        "subject",
        "query_BrumiR",
        "identity_BrumiR",
        "length_BrumiR",
        "query_miRDeep2",
        "identity_miRDeep2",
        "length_miRDeep2",
    ]]

    final.columns = [
        "miRGeneDB_ID",
        "BrumiR_cluster",
        "BrumiR_identity",
        "BrumiR_alignment_length",
        "miRDeep2_miRNA",
        "miRDeep2_identity",
        "miRDeep2_alignment_length",
    ]

    final.to_csv(out_file, sep="\t", index=False)

    print("Written:", out_file)
    print("Total shared miRNAs:", len(final))


if __name__ == "__main__":
    main()
