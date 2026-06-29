#!/usr/bin/env python3

"""
Figure 3C – Extract miRDeep2 core sequences from miRGeneDB

Description
-----------
This script extracts mature miRNA sequences from a human miRGeneDB
FASTA file for the unique miRDeep2 p=0.8 core miRNA IDs.

Inputs
------
--ids
    Text file containing unique miRDeep2 miRNA IDs, one per line.

--ref_fasta
    Human miRGeneDB mature FASTA file.

--out_fasta
    Output FASTA file.

Outputs
-------
FASTA file containing mature miRGeneDB sequences for the requested
miRDeep2 core miRNA IDs.

Usage
-----
python3 Figure3C_extract_miRDeep2_miRGeneDB.py \\
  --ids md_p08_unique_ids.txt \\
  --ref_fasta mirgenedb_hsa_mature.fa \\
  --out_fasta md_p08_unique_sequences.fa
"""

import argparse
from pathlib import Path


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract miRDeep2 core mature sequences from miRGeneDB FASTA."
    )

    parser.add_argument(
        "--ids",
        required=True,
        help="Text file containing unique miRDeep2 miRNA IDs.",
    )
    parser.add_argument(
        "--ref_fasta",
        required=True,
        help="Human miRGeneDB mature FASTA file.",
    )
    parser.add_argument(
        "--out_fasta",
        required=True,
        help="Output FASTA file.",
    )

    return parser.parse_args()


def main() -> None:
    """Extract mature miRNA sequences from miRGeneDB FASTA."""
    args = parse_args()

    ids_file = validate_file(args.ids, "miRNA ID file")
    ref_fasta = validate_file(args.ref_fasta, "miRGeneDB mature FASTA")
    out_fasta = Path(args.out_fasta)

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    wanted = {
        line.strip()
        for line in ids_file.read_text().splitlines()
        if line.strip()
    }

    header = None
    sequence_parts = []
    keep = False
    written = 0

    with ref_fasta.open() as fasta_handle, out_fasta.open("w") as out_handle:
        for raw_line in fasta_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None and keep:
                    out_handle.write(f">{header}\n{''.join(sequence_parts)}\n")
                    written += 1

                header = line[1:].split()[0]
                keep = header in wanted
                sequence_parts = []

            else:
                sequence_parts.append(line)

        if header is not None and keep:
            out_handle.write(f">{header}\n{''.join(sequence_parts)}\n")
            written += 1

    print("Requested IDs:", len(wanted))
    print("Written sequences:", written)
    print("Output:", out_fasta)


if __name__ == "__main__":
    main()
