#!/usr/bin/env python3

"""
Figure 3C – Extract BrumiR core FASTA sequences

Description
-----------
This script extracts representative centroid sequences corresponding
to the BrumiR-RF 0.95 core union cluster list and writes them to a
FASTA file for downstream annotation against miRGeneDB.

Inputs
------
--core_ids
    Text file containing BrumiR core union cluster IDs, one per line.

--cluster_fasta
    FASTA file containing named BrumiR clustered centroid sequences.

--out_fasta
    Output FASTA file.

Outputs
-------
FASTA file containing only BrumiR core union sequences.

Usage
-----
python3 Figure3C_extract_BrumiR_core.py \\
  --core_ids brumir.0.95.core_union.txt \\
  --cluster_fasta all.candidates_clustered.cluster_named.fasta \\
  --out_fasta brumir_core220_sequences.fa
"""

import argparse
from pathlib import Path
from typing import TextIO


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract BrumiR core union sequences from clustered FASTA."
    )

    parser.add_argument(
        "--core_ids",
        required=True,
        help="Text file with BrumiR core union cluster IDs, one per line.",
    )
    parser.add_argument(
        "--cluster_fasta",
        required=True,
        help="FASTA file with named BrumiR clustered centroid sequences.",
    )
    parser.add_argument(
        "--out_fasta",
        required=True,
        help="Output FASTA file.",
    )

    return parser.parse_args()


def flush_record(header: str, sequence: str, wanted_ids: set[str], out_handle: TextIO) -> None:
    """Write a FASTA record if its first header token is in the wanted ID set."""
    if header and header.split()[0] in wanted_ids:
        out_handle.write(f">{header}\n{sequence}\n")


def main() -> None:
    """Extract selected FASTA records."""
    args = parse_args()

    core_ids_file = validate_file(args.core_ids, "Core ID file")
    cluster_fasta = validate_file(args.cluster_fasta, "Cluster FASTA")
    out_fasta = Path(args.out_fasta)

    out_fasta.parent.mkdir(parents=True, exist_ok=True)

    wanted_ids = {
        line.strip()
        for line in core_ids_file.read_text().splitlines()
        if line.strip()
    }

    header = None
    sequence_parts = []

    with cluster_fasta.open() as fasta_handle, out_fasta.open("w") as out_handle:
        for raw_line in fasta_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    flush_record(header, "".join(sequence_parts), wanted_ids, out_handle)

                header = line[1:]
                sequence_parts = []

            else:
                sequence_parts.append(line)

        if header is not None:
            flush_record(header, "".join(sequence_parts), wanted_ids, out_handle)

    print("Written:", out_fasta)


if __name__ == "__main__":
    main()
