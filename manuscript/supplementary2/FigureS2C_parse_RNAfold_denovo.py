#!/usr/bin/env python3
"""
Figure S2C – Parse RNAfold output for de novo candidates

Description
-----------
Parse raw RNAfold output for de novo BrumiR-supported candidates, merge it
with the corresponding candidate catalog, and compute a simple ranking score
to prioritize candidates for downstream inspection.

The merge is performed by row order, assuming that:
1. the catalog rows
2. the precursor FASTA entries
3. the RNAfold output blocks

all preserve the same order.

Inputs
------
--catalog_tsv
    Input catalog TSV with de novo candidate metadata.

--rnafold_raw
    Raw RNAfold output file.

--out_all
    Output TSV with all merged candidates and RNAfold annotations.

--out_top
    Output TSV with the top-ranked candidates.

--top_n
    Number of top candidates to export. Default: 4.

Outputs
-------
1. A full TSV with catalog + RNAfold-derived structure and score.
2. A top-N TSV filtered to candidates with parsed structure.

Usage
-----
python FigureS2C_parse_RNAfold_denovo.py \\
  --catalog_tsv path/to/catalog.tsv \\
  --rnafold_raw path/to/rnafold.raw.txt \\
  --out_all path/to/all_candidates_rnafold.tsv \\
  --out_top path/to/top_candidates_rnafold.tsv \\
  --top_n 4

Notes
-----
- The script expects RNAfold blocks in the form:
    >header
    sequence
    structure (MFE)

- Ranking is unchanged and based on:
    - mature length between 21 and 24 nt
    - precursor length between 60 and 120 nt
    - more negative RNAfold MFE
"""

import argparse
import csv
import re
import sys
from pathlib import Path


MFE_RE = re.compile(r"\(([-]?\d+(?:\.\d+)?)\)\s*$")


def validate_input_file(path: Path, label: str) -> None:
    """Validate that an input file exists."""
    if not path.exists():
        sys.exit(f"ERROR: {label} not found: {path}")
    if not path.is_file():
        sys.exit(f"ERROR: {label} is not a file: {path}")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Parse RNAfold output for de novo candidates and rank top structures."
        )
    )
    parser.add_argument(
        "--catalog_tsv",
        required=True,
        help="Input candidate catalog TSV.",
    )
    parser.add_argument(
        "--rnafold_raw",
        required=True,
        help="Raw RNAfold output.",
    )
    parser.add_argument(
        "--out_all",
        required=True,
        help="Output TSV with all parsed candidates.",
    )
    parser.add_argument(
        "--out_top",
        required=True,
        help="Output TSV with top-ranked candidates.",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=4,
        help="Number of top candidates to export. Default: 4.",
    )

    return parser.parse_args()


def read_catalog(path: Path):
    """Read candidate catalog TSV."""
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames

    if fieldnames is None:
        sys.exit(f"ERROR: Catalog TSV appears to be empty or lacks a header: {path}")

    return rows, fieldnames


def read_rnafold_blocks(path: Path):
    """
    Parse RNAfold raw output and return ordered blocks.

    Returns
    -------
    list of dict
        Each dictionary contains:
        - header_raw
        - rnafold_seq
        - structure
        - rnafold_mfe
    """
    blocks = []

    with path.open("r") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip() != ""]

    line_index = 0

    while line_index < len(lines):
        if not lines[line_index].startswith(">"):
            line_index += 1
            continue

        header = lines[line_index][1:].strip()

        if line_index + 2 >= len(lines):
            break

        sequence = lines[line_index + 1].strip()
        structure_line = lines[line_index + 2].rstrip()

        mfe_match = MFE_RE.search(structure_line)
        mfe = ""
        structure = structure_line

        if mfe_match:
            mfe = mfe_match.group(1)
            structure = structure_line[:mfe_match.start()].rstrip()

        blocks.append(
            {
                "header_raw": header,
                "rnafold_seq": sequence,
                "structure": structure,
                "rnafold_mfe": mfe,
            }
        )

        line_index += 3

    return blocks


def to_int(value):
    """Convert a value to int when possible."""
    try:
        return int(float(value))
    except Exception:
        return None


def to_float(value):
    """Convert a value to float when possible."""
    try:
        return float(value)
    except Exception:
        return None


def score_row(row):
    """
    Compute ranking score.

    Original scoring logic preserved:
    - +10 if mature length is in [21, 24]
    - +10 if precursor length is in [60, 120]
    - add capped contribution from negative MFE
    """
    mature_len = to_int(row.get("mature_len"))
    precursor_len = to_int(row.get("precursor_len"))
    mfe = to_float(row.get("rnafold_mfe"))

    mature_ok = mature_len is not None and 21 <= mature_len <= 24
    precursor_ok = precursor_len is not None and 60 <= precursor_len <= 120

    score = 0.0
    score += 10.0 if mature_ok else 0.0
    score += 10.0 if precursor_ok else 0.0

    if mfe is not None:
        score += max(0.0, min(40.0, -mfe))

    return mature_ok, precursor_ok, score


def write_tsv(path: Path, rows, fieldnames) -> None:
    """Write rows to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    catalog_path = Path(args.catalog_tsv)
    rnafold_path = Path(args.rnafold_raw)
    out_all_path = Path(args.out_all)
    out_top_path = Path(args.out_top)

    validate_input_file(catalog_path, "Catalog TSV")
    validate_input_file(rnafold_path, "RNAfold raw output")

    if args.top_n < 1:
        sys.exit("ERROR: --top_n must be greater than or equal to 1.")

    catalog_rows, catalog_fields = read_catalog(catalog_path)
    rnafold_blocks = read_rnafold_blocks(rnafold_path)

    if len(catalog_rows) == 0:
        sys.exit(f"ERROR: Catalog TSV contains no rows: {catalog_path}")

    if len(rnafold_blocks) == 0:
        sys.exit(f"ERROR: No RNAfold blocks were parsed from: {rnafold_path}")

    # Merge by order, assuming RNAfold blocks correspond to catalog rows.
    n_merged = min(len(catalog_rows), len(rnafold_blocks))
    merged_rows = []

    for index in range(n_merged):
        row = catalog_rows[index]
        block = rnafold_blocks[index]

        row["rnafold_seq"] = block["rnafold_seq"]
        row["structure"] = block["structure"]
        row["rnafold_mfe"] = block["rnafold_mfe"]
        row["rnafold_header_raw"] = block["header_raw"]

        mature_ok, precursor_ok, score = score_row(row)
        row["mature_len_ok"] = str(mature_ok)
        row["prec_len_ok"] = str(precursor_ok)
        row["score"] = f"{score:.2f}"

        merged_rows.append(row)

    # If lengths mismatch, still write the subset that could be merged.
    fieldnames = catalog_fields[:] if catalog_fields else []
    extra_fields = [
        "rnafold_header_raw",
        "rnafold_seq",
        "structure",
        "rnafold_mfe",
        "mature_len_ok",
        "prec_len_ok",
        "score",
    ]

    for column in extra_fields:
        if column not in fieldnames:
            fieldnames.append(column)

    write_tsv(out_all_path, merged_rows, fieldnames)

    # Top N candidates with parsed structure.
    parsed_rows = [
        row for row in merged_rows
        if (row.get("structure", "") or "").strip() != ""
    ]

    parsed_rows.sort(
        key=lambda row: float(row.get("score", "0") or 0.0),
        reverse=True,
    )

    top_rows = parsed_rows[:args.top_n]
    write_tsv(out_top_path, top_rows, fieldnames)

    if len(top_rows) == 0:
        raise SystemExit(
            "ERROR: Top TSV has 0 rows. No structure was parsed. "
            "Check RNAfold raw output and FASTA generation."
        )

    if out_top_path.stat().st_size == 0:
        raise SystemExit("ERROR: out_top file is empty.")

    if len(catalog_rows) != len(rnafold_blocks):
        print(
            "WARNING: Catalog rows and RNAfold blocks have different lengths. "
            f"Merged {n_merged} rows by order.",
            file=sys.stderr,
        )

    print(f"RNAfold full table written to: {out_all_path}")
    print(f"RNAfold top table written to: {out_top_path}")


if __name__ == "__main__":
    main()
