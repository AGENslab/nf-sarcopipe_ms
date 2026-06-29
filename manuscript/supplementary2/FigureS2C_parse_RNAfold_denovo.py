#!/usr/bin/env python3
"""
parse_rnafold_denovo.py

Description
-----------
This script parses raw RNAfold output for de novo BrumiR-supported candidates,
merges it with the corresponding candidate catalog, and computes a simple ranking
score to prioritize candidates for downstream inspection.

The merge is performed by row order, assuming that:
1. the catalog rows
2. the precursor FASTA entries
3. the RNAfold output blocks

all preserve the same order.

Inputs
------
--catalog_tsv   Input catalog TSV with de novo candidate metadata
--rnafold_raw   Raw RNAfold output file
--out_all       Output TSV with all merged candidates and RNAfold annotations
--out_top       Output TSV with the top-ranked candidates
--top_n         Number of top candidates to export (default: 4)

Outputs
-------
1. A full TSV with catalog + RNAfold-derived structure and score
2. A top-N TSV filtered to candidates with parsed structure

Notes
-----
- The script expects RNAfold blocks in the form:
    >header
    sequence
    structure (MFE)
- Ranking is based on:
    - mature length between 21 and 24 nt
    - precursor length between 60 and 120 nt
    - more negative RNAfold MFE
"""

import argparse
import csv
import re
from pathlib import Path

mfe_re = re.compile(r"\(([-]?\d+(?:\.\d+)?)\)\s*$")


def parse_args():
    p = argparse.ArgumentParser(
        description="Parse RNAfold output for de novo candidates and rank top structures."
    )
    p.add_argument("--catalog_tsv", required=True, help="Input candidate catalog TSV")
    p.add_argument("--rnafold_raw", required=True, help="Raw RNAfold output")
    p.add_argument("--out_all", required=True, help="Output TSV with all parsed candidates")
    p.add_argument("--out_top", required=True, help="Output TSV with top-ranked candidates")
    p.add_argument("--top_n", type=int, default=4, help="Number of top candidates to export")
    return p.parse_args()


def read_catalog(path):
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        return list(r), r.fieldnames


def read_rnafold_blocks(path):
    """
    Parse RNAfold raw output and return ordered blocks:
      {header_raw, rnafold_seq, structure, rnafold_mfe}
    """
    blocks = []

    with open(path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip() != ""]

    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue

        header = lines[i][1:].strip()

        if i + 2 >= len(lines):
            break

        seq = lines[i + 1].strip()
        struct_line = lines[i + 2].rstrip()

        m = mfe_re.search(struct_line)
        mfe = ""
        struct = struct_line

        if m:
            mfe = m.group(1)
            struct = struct_line[:m.start()].rstrip()

        blocks.append({
            "header_raw": header,
            "rnafold_seq": seq,
            "structure": struct,
            "rnafold_mfe": mfe
        })

        i += 3

    return blocks


def to_int(x):
    try:
        return int(float(x))
    except Exception:
        return None


def to_float(x):
    try:
        return float(x)
    except Exception:
        return None


def score_row(row):
    """
    Compute a simple ranking score based on:
    - mature length in [21, 24]
    - precursor length in [60, 120]
    - more negative MFE
    """
    mlen = to_int(row.get("mature_len"))
    plen = to_int(row.get("precursor_len"))
    mfe = to_float(row.get("rnafold_mfe"))

    mature_ok = (mlen is not None and 21 <= mlen <= 24)
    prec_ok = (plen is not None and 60 <= plen <= 120)

    score = 0.0
    score += 10.0 if mature_ok else 0.0
    score += 10.0 if prec_ok else 0.0

    if mfe is not None:
        score += max(0.0, min(40.0, -mfe))

    return mature_ok, prec_ok, score


def main():
    args = parse_args()

    cat_rows, cat_fields = read_catalog(args.catalog_tsv)
    blocks = read_rnafold_blocks(args.rnafold_raw)

    # Merge by order, assuming RNAfold blocks correspond to catalog rows
    n = min(len(cat_rows), len(blocks))
    merged = []

    for idx in range(n):
        row = cat_rows[idx]
        blk = blocks[idx]

        row["rnafold_seq"] = blk["rnafold_seq"]
        row["structure"] = blk["structure"]
        row["rnafold_mfe"] = blk["rnafold_mfe"]
        row["rnafold_header_raw"] = blk["header_raw"]

        mature_ok, prec_ok, sc = score_row(row)
        row["mature_len_ok"] = str(mature_ok)
        row["prec_len_ok"] = str(prec_ok)
        row["score"] = f"{sc:.2f}"

        merged.append(row)

    # If lengths mismatch, still write the subset that could be merged
    fieldnames = cat_fields[:] if cat_fields else []
    extra = [
        "rnafold_header_raw",
        "rnafold_seq",
        "structure",
        "rnafold_mfe",
        "mature_len_ok",
        "prec_len_ok",
        "score"
    ]

    for col in extra:
        if col not in fieldnames:
            fieldnames.append(col)

    with open(args.out_all, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for row in merged:
            w.writerow(row)

    # Top N candidates with parsed structure
    ok = [r for r in merged if (r.get("structure", "") or "").strip() != ""]
    ok.sort(key=lambda r: float(r.get("score", "0") or 0.0), reverse=True)
    top = ok[:args.top_n]

    with open(args.out_top, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for row in top:
            w.writerow(row)

    # Fail here if nothing usable was parsed
    if len(top) == 0:
        raise SystemExit(
            "ERROR: Top TSV has 0 rows (no structure parsed). Check RNAfold raw output and FASTA generation."
        )

    if Path(args.out_top).stat().st_size == 0:
        raise SystemExit("ERROR: out_top file empty")


if __name__ == "__main__":
    main()
