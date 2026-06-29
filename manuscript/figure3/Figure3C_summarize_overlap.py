#!/usr/bin/env python3

"""
Figure 3C – Summarize miRGeneDB overlap

Description
-----------
This script builds a compact summary table for the exact-match
overlap against human miRGeneDB for BrumiR core220 sequences and
miRDeep2 core941 sequences.

It counts:
- total input sequences from FASTA headers
- matched sequences from unique query IDs in BLAST-like TSV files
- unmatched putative novel sequences as total minus matched

Inputs
------
--brumir_fasta
    FASTA file containing BrumiR core220 sequences.

--mirdeep_fasta
    FASTA file containing miRDeep2 core941 sequences.

--brumir_blast
    BrumiR vs miRGeneDB BLAST-like TSV file.

--mirdeep_blast
    miRDeep2 vs miRGeneDB BLAST-like TSV file.

--out
    Output summary TSV file.

Outputs
-------
TSV file with columns:
algorithm, total, matched_100pct, unmatched_putative_novel

Usage
-----
python3 Figure3C_summarize_overlap.py \\
  --brumir_fasta brumir_core220_sequences.fa \\
  --mirdeep_fasta md_p08_unique_sequences.fa \\
  --brumir_blast brumir_core220_vs_miRGeneDB.tsv \\
  --mirdeep_blast md_core_vs_miRGeneDB.tsv \\
  --out overlap_miRGeneDB_summary.tsv
"""

import argparse
from pathlib import Path


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def count_fasta_headers(path: Path) -> int:
    """Count FASTA records by counting header lines."""
    return sum(1 for line in path.open() if line.startswith(">"))


def count_blast_queries(path: Path) -> int:
    """Count unique query IDs in a BLAST-like TSV file."""
    hits = set()

    with path.open() as blast_handle:
        for line in blast_handle:
            line = line.strip()

            if line:
                hits.add(line.split("\t")[0])

    return len(hits)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize miRGeneDB overlap for BrumiR and miRDeep2 core sequences."
    )

    parser.add_argument(
        "--brumir_fasta",
        required=True,
        help="FASTA file containing BrumiR core220 sequences.",
    )
    parser.add_argument(
        "--mirdeep_fasta",
        required=True,
        help="FASTA file containing miRDeep2 core941 sequences.",
    )
    parser.add_argument(
        "--brumir_blast",
        required=True,
        help="BrumiR vs miRGeneDB BLAST-like TSV file.",
    )
    parser.add_argument(
        "--mirdeep_blast",
        required=True,
        help="miRDeep2 vs miRGeneDB BLAST-like TSV file.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output summary TSV file.",
    )

    return parser.parse_args()


def main() -> None:
    """Build miRGeneDB overlap summary table."""
    args = parse_args()

    brumir_fasta = validate_file(args.brumir_fasta, "BrumiR FASTA")
    mirdeep_fasta = validate_file(args.mirdeep_fasta, "miRDeep2 FASTA")
    brumir_blast = validate_file(args.brumir_blast, "BrumiR BLAST TSV")
    mirdeep_blast = validate_file(args.mirdeep_blast, "miRDeep2 BLAST TSV")
    out_file = Path(args.out)

    out_file.parent.mkdir(parents=True, exist_ok=True)

    brumir_total = count_fasta_headers(brumir_fasta)
    mirdeep_total = count_fasta_headers(mirdeep_fasta)

    brumir_matched = count_blast_queries(brumir_blast)
    mirdeep_matched = count_blast_queries(mirdeep_blast)

    with out_file.open("w") as out_handle:
        out_handle.write("algorithm\ttotal\tmatched_100pct\tunmatched_putative_novel\n")
        out_handle.write(
            f"BrumiR_core220\t{brumir_total}\t{brumir_matched}\t"
            f"{brumir_total - brumir_matched}\n"
        )
        out_handle.write(
            f"miRDeep2_core941\t{mirdeep_total}\t{mirdeep_matched}\t"
            f"{mirdeep_total - mirdeep_matched}\n"
        )

    print("Written:", out_file)
    print(
        f"BrumiR_core220: total={brumir_total}, "
        f"matched_100pct={brumir_matched}, "
        f"unmatched={brumir_total - brumir_matched}"
    )
    print(
        f"miRDeep2_core941: total={mirdeep_total}, "
        f"matched_100pct={mirdeep_matched}, "
        f"unmatched={mirdeep_total - mirdeep_matched}"
    )


if __name__ == "__main__":
    main()
