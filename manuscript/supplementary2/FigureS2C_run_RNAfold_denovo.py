#!/usr/bin/env python3
"""
Figure S2C – Run RNAfold on de novo candidates

Description
-----------
Run RNAfold on de novo BrumiR-supported precursor sequences and compute
simple structural metrics to rank candidates for downstream inspection.

For each de novo candidate, it calculates:
- RNAfold secondary structure
- minimum free energy (MFE)
- MFE normalized by precursor length
- fraction of paired bases within the mature miRNA window
- a simple potency score based on structural stability and mature-region pairing

Inputs
------
--catalog_tsv
    Input TSV file containing de novo candidate information.
    Expected columns:
    - is_denovo_exact
    - precursor_seq
    - mature_start
    - mature_end

--out_all
    Output TSV with all scored de novo candidates.
    Default: denovo_rnafold_all.tsv

--out_top
    Output TSV with the top-ranked candidates.
    Default: denovo_rnafold_top4.tsv

--top_n
    Number of top candidates to export.
    Default: 4

Outputs
-------
1. A TSV containing all de novo candidates with RNAfold-derived metrics.
2. A TSV containing only the top N ranked candidates.

Usage
-----
python FigureS2C_run_RNAfold_denovo.py \\
  --catalog_tsv path/to/denovo_passfilter_catalog.tsv \\
  --out_all path/to/denovo_rnafold_all.tsv \\
  --out_top path/to/denovo_rnafold_top4.tsv \\
  --top_n 4

Requirements
------------
- RNAfold must be available in the execution environment.
- pandas must be installed.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def validate_input_file(path: Path, label: str) -> None:
    """Validate that an input file exists."""
    if not path.exists():
        sys.exit(f"ERROR: {label} not found: {path}")
    if not path.is_file():
        sys.exit(f"ERROR: {label} is not a file: {path}")


def validate_columns(df: pd.DataFrame, required_columns: list[str], path: Path) -> None:
    """Validate required columns in the input catalog."""
    missing_columns = [col for col in required_columns if col not in df.columns]

    if missing_columns:
        sys.exit(
            "ERROR: Missing required columns in "
            f"{path}: {', '.join(missing_columns)}"
        )


def run_rnafold(seq: str) -> tuple[str, float]:
    """
    Run RNAfold on a single precursor sequence.

    RNAfold output format:
      SEQUENCE
      DOTBRACKET ( -12.30 )

    Returns
    -------
    struct : str
        Dot-bracket secondary structure.
    mfe : float
        Minimum free energy.
    """
    process = subprocess.run(
        ["RNAfold", "--noPS"],
        input=(seq.strip() + "\n").encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )

    output_lines = process.stdout.decode().strip().splitlines()
    last_line = output_lines[-1]

    # Example: "(((...))) (-12.30)"
    structure = last_line.split()[0]
    mfe_value = last_line.split()[-1].strip("()")

    return structure, float(mfe_value)


def paired_fraction_in_window(struct: str, start1: int, end1: int) -> float:
    """
    Compute the fraction of paired positions in a 1-based inclusive window
    over the dot-bracket structure string.
    """
    window = struct[start1 - 1:end1]

    if len(window) == 0:
        return 0.0

    paired = sum(1 for char in window if char in "()")
    return paired / len(window)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Run RNAfold on de novo candidates and rank them using simple "
            "structural metrics."
        )
    )
    parser.add_argument(
        "--catalog_tsv",
        required=True,
        help="Input denovo_passfilter_catalog.tsv.",
    )
    parser.add_argument(
        "--out_all",
        default="denovo_rnafold_all.tsv",
        help="Output TSV with all candidates. Default: denovo_rnafold_all.tsv.",
    )
    parser.add_argument(
        "--out_top",
        default="denovo_rnafold_top4.tsv",
        help="Output TSV with top-ranked candidates. Default: denovo_rnafold_top4.tsv.",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=4,
        help="Number of top-ranked candidates to export. Default: 4.",
    )

    return parser.parse_args()


def main():
    args = parse_args()

    catalog_path = Path(args.catalog_tsv)
    out_all_path = Path(args.out_all)
    out_top_path = Path(args.out_top)

    validate_input_file(catalog_path, "Catalog TSV")

    if shutil.which("RNAfold") is None:
        sys.exit("ERROR: RNAfold was not found in PATH.")

    if args.top_n < 1:
        sys.exit("ERROR: --top_n must be greater than or equal to 1.")

    out_all_path.parent.mkdir(parents=True, exist_ok=True)
    out_top_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(catalog_path, sep="\t")

    required_columns = [
        "is_denovo_exact",
        "precursor_seq",
        "mature_start",
        "mature_end",
    ]
    validate_columns(df, required_columns, catalog_path)

    df = df[df["is_denovo_exact"] == 1].copy()

    if df.empty:
        print(
            "ERROR: No de novo rows found in catalog (is_denovo_exact == 1).",
            file=sys.stderr,
        )
        sys.exit(1)

    structures = []
    mfes = []
    mfe_per_nt = []
    paired_mature = []

    for _, row in df.iterrows():
        sequence = row["precursor_seq"]
        structure, mfe = run_rnafold(sequence)

        structures.append(structure)
        mfes.append(mfe)
        mfe_per_nt.append(mfe / len(sequence))

        paired_fraction = paired_fraction_in_window(
            structure,
            int(row["mature_start"]),
            int(row["mature_end"]),
        )
        paired_mature.append(paired_fraction)

    df["rnafold_structure"] = structures
    df["rnafold_mfe"] = mfes
    df["rnafold_mfe_per_nt"] = mfe_per_nt
    df["paired_frac_mature"] = paired_mature

    # Original ranking logic preserved:
    # more negative MFE/nt and more mature-region pairing -> higher score.
    df["potency_score"] = (-df["rnafold_mfe_per_nt"]) + df["paired_frac_mature"]

    df.sort_values(
        ["potency_score", "rnafold_mfe_per_nt"],
        ascending=[False, True],
        inplace=True,
    )

    df.to_csv(out_all_path, sep="\t", index=False)

    top = df.head(args.top_n).copy()
    top.to_csv(out_top_path, sep="\t", index=False)

    print(f"[OK] RNAfold computed for {len(df)} de novo candidates.", file=sys.stderr)
    print(f"[OK] Wrote: {out_all_path}", file=sys.stderr)
    print(f"[OK] Wrote TOP{args.top_n}: {out_top_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
