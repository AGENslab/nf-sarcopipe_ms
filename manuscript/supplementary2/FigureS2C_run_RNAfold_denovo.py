#!/usr/bin/env python3
"""
run_rnafold_denovo.py

Description
-----------
This script runs RNAfold on de novo BrumiR-supported precursor sequences and computes
simple structural metrics to rank candidates for downstream inspection.

For each de novo candidate, it calculates:
- RNAfold secondary structure (dot-bracket notation)
- minimum free energy (MFE)
- MFE normalized by precursor length
- fraction of paired bases within the mature miRNA window
- a simple potency score based on structural stability and mature-region pairing

Inputs
------
--catalog_tsv   Input TSV file containing de novo candidate information.
                Expected columns include:
                - is_denovo_exact
                - precursor_seq
                - mature_start
                - mature_end

--out_all       Output TSV with all scored de novo candidates
--out_top       Output TSV with the top-ranked candidates
--top_n         Number of top candidates to export (default: 4)

Outputs
-------
1. A TSV containing all de novo candidates with RNAfold-derived metrics
2. A TSV containing only the top N ranked candidates

Requirements
------------
- RNAfold must be available in the execution environment
- pandas must be installed
"""

import argparse
import subprocess
import sys

import pandas as pd


def run_rnafold(seq: str) -> tuple[str, float]:
    """
    Run RNAfold on a single precursor sequence.

    RNAfold output format:
      SEQUENCE
      DOTBRACKET ( -12.30 )

    Returns
    -------
    struct : str
        Dot-bracket secondary structure
    mfe : float
        Minimum free energy
    """
    p = subprocess.run(
        ["RNAfold", "--noPS"],
        input=(seq.strip() + "\n").encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )

    out = p.stdout.decode().strip().splitlines()
    last = out[-1]

    # Example: "(((...))) (-12.30)"
    struct = last.split()[0]
    mfe_str = last.split()[-1].strip("()")

    return struct, float(mfe_str)


def paired_fraction_in_window(struct: str, start1: int, end1: int) -> float:
    """
    Compute the fraction of paired positions in a 1-based inclusive window
    over the dot-bracket structure string.
    """
    window = struct[start1 - 1:end1]

    if len(window) == 0:
        return 0.0

    paired = sum(1 for c in window if c in "()")
    return paired / len(window)


def main():
    ap = argparse.ArgumentParser(
        description="Run RNAfold on de novo candidates and rank them using simple structural metrics."
    )
    ap.add_argument("--catalog_tsv", required=True, help="Input denovo_passfilter_catalog.tsv")
    ap.add_argument("--out_all", default="denovo_rnafold_all.tsv", help="Output TSV with all candidates")
    ap.add_argument("--out_top", default="denovo_rnafold_top4.tsv", help="Output TSV with top-ranked candidates")
    ap.add_argument("--top_n", type=int, default=4, help="Number of top-ranked candidates to export")
    args = ap.parse_args()

    df = pd.read_csv(args.catalog_tsv, sep="\t")
    df = df[df["is_denovo_exact"] == 1].copy()

    if df.empty:
        print("[ERROR] No de novo rows found in catalog (is_denovo_exact == 1).", file=sys.stderr)
        sys.exit(1)

    structs = []
    mfes = []
    mfe_per_nt = []
    paired_mature = []

    for _, row in df.iterrows():
        seq = row["precursor_seq"]
        struct, mfe = run_rnafold(seq)

        structs.append(struct)
        mfes.append(mfe)
        mfe_per_nt.append(mfe / len(seq))

        pm = paired_fraction_in_window(
            struct,
            int(row["mature_start"]),
            int(row["mature_end"])
        )
        paired_mature.append(pm)

    df["rnafold_structure"] = structs
    df["rnafold_mfe"] = mfes
    df["rnafold_mfe_per_nt"] = mfe_per_nt
    df["paired_frac_mature"] = paired_mature

    # Simple ranking score:
    # more negative MFE/nt and more mature-region pairing -> higher score
    df["potency_score"] = (-df["rnafold_mfe_per_nt"]) + df["paired_frac_mature"]

    df.sort_values(
        ["potency_score", "rnafold_mfe_per_nt"],
        ascending=[False, True],
        inplace=True
    )

    df.to_csv(args.out_all, sep="\t", index=False)

    top = df.head(args.top_n).copy()
    top.to_csv(args.out_top, sep="\t", index=False)

    print(f"[OK] RNAfold computed for {len(df)} de novo candidates.", file=sys.stderr)
    print(f"[OK] Wrote: {args.out_all}", file=sys.stderr)
    print(f"[OK] Wrote TOP{args.top_n}: {args.out_top}", file=sys.stderr)


if __name__ == "__main__":
    main()
