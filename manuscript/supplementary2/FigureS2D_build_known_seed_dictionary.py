#!/usr/bin/env python3
"""
Figure S2D – Build known miRNA seed dictionary

Description
-----------
Extract canonical 7-mer seeds from known human miRNAs and build a seed
reference dictionary.

The canonical seed is extracted from positions 2–8 of the mature miRNA
sequence using Python slicing s[1:8].

Inputs
------
--fasta
    FASTA file containing known mature human miRNAs.

Outputs
-------
--out
    TSV file with columns:
    - miRNA_id
    - family
    - mature_seq
    - canonical_seed

Usage
-----
python FigureS2D_build_known_seed_dictionary.py \\
  --fasta path/to/hsa_mature.mirdeep2.OK.fa \\
  --out path/to/known_seed_dictionary.tsv
"""

import argparse
import sys
from pathlib import Path


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build known miRNA canonical seed dictionary."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA with known mature human miRNAs.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV path.",
    )

    return parser.parse_args()


def validate_input_file(path: Path, label: str) -> None:
    """Validate that an input file exists."""
    if not path.exists():
        sys.exit(f"ERROR: {label} not found: {path}")
    if not path.is_file():
        sys.exit(f"ERROR: {label} is not a file: {path}")


def family_from_id(mirna_id: str) -> str:
    """Return a simple family approximation for summary."""
    family = mirna_id.replace("hsa-", "")
    return family


def add_record(rows: list[tuple[str, str, str, str]], header: str, seq_parts: list[str]) -> None:
    """Add one FASTA record to the output rows when sequence length is >= 8."""
    mature_seq = "".join(seq_parts).upper().replace("T", "U")

    if len(mature_seq) >= 8:
        seed = mature_seq[1:8]  # canonical positions 2–8
        rows.append(
            (
                header,
                family_from_id(header),
                mature_seq,
                seed,
            )
        )


def read_fasta_seed_rows(fasta_path: Path):
    """Read FASTA and extract known miRNA seed rows."""
    rows = []
    header = None
    seq_parts = []

    with fasta_path.open("r") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    add_record(rows, header, seq_parts)

                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

        if header is not None:
            add_record(rows, header, seq_parts)

    return rows


def write_seed_dictionary(rows, output_path: Path) -> None:
    """Write seed dictionary TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w") as handle:
        handle.write("miRNA_id\tfamily\tmature_seq\tcanonical_seed\n")

        for row in rows:
            handle.write("\t".join(row) + "\n")


def main():
    args = parse_args()

    fasta_path = Path(args.fasta)
    output_path = Path(args.out)

    validate_input_file(fasta_path, "Input FASTA")

    rows = read_fasta_seed_rows(fasta_path)
    write_seed_dictionary(rows, output_path)

    print("Written:", output_path)
    print("Total known miRNAs:", len(rows))


if __name__ == "__main__":
    main()
