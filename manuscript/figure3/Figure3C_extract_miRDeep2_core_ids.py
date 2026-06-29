#!/usr/bin/env python3

"""
Figure 3C – Extract unique miRDeep2 core IDs

Description
-----------
This script extracts a unique list of miRDeep2 known miRNA IDs
from the p=0.8 core sets table, removing duplicated entries and
writing the identifiers to a text file for downstream annotation.

Inputs
------
--input
    miRDeep2 p0.8 core sets TSV file.

--output
    Output text file containing one unique miRNA ID per line.

Outputs
-------
Unique miRDeep2 ID list.

Usage
-----
python3 Figure3C_extract_miRDeep2_core_ids.py \
    --input miRDeep2.known.p08.core_sets.tsv \
    --output md_p08_unique_ids.txt
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
        description="Extract unique miRDeep2 IDs from a p0.8 core set table."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="miRDeep2 p0.8 core sets TSV file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output text file containing unique miRNA IDs.",
    )

    return parser.parse_args()


def main() -> None:
    """Extract unique miRDeep2 IDs."""
    args = parse_args()

    input_file = validate_file(args.input, "Core set TSV")
    output_file = Path(args.output)

    output_file.parent.mkdir(parents=True, exist_ok=True)

    ids = set()

    with input_file.open() as infile:
        next(infile)

        for line in infile:
            fields = line.rstrip("\n").split("\t")

            if fields and fields[0].strip():
                ids.add(fields[0].strip())

    with output_file.open("w") as outfile:
        for mirna_id in sorted(ids):
            outfile.write(f"{mirna_id}\n")

    print("Unique IDs:", len(ids))
    print("Written:", output_file)


if __name__ == "__main__":
    main()
