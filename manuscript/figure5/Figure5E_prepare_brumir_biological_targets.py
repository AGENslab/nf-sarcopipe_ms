#!/usr/bin/env python3

"""
Figure 5E – Prepare BrumiR biological targets

Description
-----------
This script filters BrumiR-derived predicted targets from the
curated sarcopenia/exercise gene set, assigns biological categories,
and prioritizes top targets per BrumiR candidate for Figure 5E.

Inputs
------
--input
    TSV file with sarcopenia-related seed matches.

--out_all
    Output TSV with all BrumiR biological targets.

--out_top
    Output TSV with top BrumiR biological targets.

Outputs
-------
brumir_biological_targets.tsv
brumir_biological_targets_top.tsv

Usage
-----
python3 Figure5E_prepare_brumir_biological_targets.py \\
  --input seed_matches_sarcopenia.tsv \\
  --out_all brumir_biological_targets.tsv \\
  --out_top brumir_biological_targets_top.tsv
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


CATEGORY_MAP = {
    "IL6": "Inflammation / immunity",
    "TNF": "Inflammation / immunity",
    "NFKB1": "Inflammation / immunity",
    "RELA": "Inflammation / immunity",
    "IL1B": "Inflammation / immunity",
    "CXCL8": "Inflammation / immunity",
    "CCL2": "Inflammation / immunity",
    "FOXO3": "Muscle atrophy / growth",
    "FOXO1": "Muscle atrophy / growth",
    "TRIM63": "Muscle atrophy / growth",
    "FBXO32": "Muscle atrophy / growth",
    "MSTN": "Muscle atrophy / growth",
    "IGF1": "Muscle atrophy / growth",
    "AKT1": "Muscle atrophy / growth",
    "MTOR": "Muscle atrophy / growth",
    "RPS6KB1": "Muscle atrophy / growth",
    "EIF4EBP1": "Muscle atrophy / growth",
    "COL1A1": "Fibrosis / ECM remodeling",
    "COL3A1": "Fibrosis / ECM remodeling",
    "COL5A1": "Fibrosis / ECM remodeling",
    "TGFB1": "Fibrosis / ECM remodeling",
    "SMAD2": "Fibrosis / ECM remodeling",
    "SMAD3": "Fibrosis / ECM remodeling",
    "SMAD4": "Fibrosis / ECM remodeling",
    "MMP2": "Fibrosis / ECM remodeling",
    "MMP9": "Fibrosis / ECM remodeling",
    "TIMP1": "Fibrosis / ECM remodeling",
    "TIMP2": "Fibrosis / ECM remodeling",
    "SIRT1": "Senescence / damage",
    "SIRT3": "Senescence / damage",
    "CDKN1A": "Senescence / damage",
    "CDKN2A": "Senescence / damage",
    "TP53": "Senescence / damage",
    "LMNB1": "Senescence / damage",
    "GDF11": "Senescence / damage",
    "BECN1": "Autophagy / degradation",
    "MAP1LC3B": "Autophagy / degradation",
    "ATG5": "Autophagy / degradation",
    "ATG7": "Autophagy / degradation",
    "SQSTM1": "Autophagy / degradation",
    "ULK1": "Autophagy / degradation",
    "MYOD1": "Myogenesis / muscle structure",
    "MYOG": "Myogenesis / muscle structure",
    "PAX7": "Myogenesis / muscle structure",
    "MYH1": "Myogenesis / muscle structure",
    "MYH2": "Myogenesis / muscle structure",
    "MYH7": "Myogenesis / muscle structure",
    "ACTA1": "Myogenesis / muscle structure",
    "DES": "Myogenesis / muscle structure",
}


FIELDNAMES = [
    "renamed_miRNA",
    "gene_symbol",
    "category",
    "n_sites",
    "abs_log2FC",
    "priority_score",
    "direction",
]


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def prepare_output(path: str) -> Path:
    """Create parent directory for an output file if needed."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare prioritized BrumiR biological targets for Figure 5E."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="TSV file with sarcopenia-related seed matches.",
    )
    parser.add_argument(
        "--out_all",
        required=True,
        help="Output TSV with all BrumiR biological targets.",
    )
    parser.add_argument(
        "--out_top",
        required=True,
        help="Output TSV with top BrumiR biological targets.",
    )

    return parser.parse_args()


def main() -> None:
    """Prepare BrumiR biological target tables."""
    args = parse_args()

    input_file = validate_file(args.input, "Sarcopenia seed-match table")
    out_all = prepare_output(args.out_all)
    out_top = prepare_output(args.out_top)

    rows = []

    with input_file.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_cols = {
            "source",
            "gene_symbol",
            "renamed_miRNA",
            "n_sites",
            "log2FoldChange",
            "direction",
        }
        missing_cols = required_cols - set(reader.fieldnames or [])

        if missing_cols:
            raise SystemExit(
                "ERROR: Input table missing required columns: "
                f"{', '.join(sorted(missing_cols))}"
            )

        for row in reader:
            if row["source"].strip() != "BrumiR":
                continue

            gene = row["gene_symbol"].strip()
            mirna = row["renamed_miRNA"].strip()

            if CATEGORY_MAP.get(gene, "Other") == "Other":
                continue

            try:
                n_sites = int(row["n_sites"])
            except ValueError:
                n_sites = 0

            try:
                abs_fc = abs(float(row["log2FoldChange"]))
            except ValueError:
                abs_fc = 0.0

            priority_score = n_sites + abs_fc

            rows.append({
                "renamed_miRNA": mirna,
                "gene_symbol": gene,
                "category": CATEGORY_MAP[gene],
                "n_sites": n_sites,
                "abs_log2FC": round(abs_fc, 3),
                "priority_score": round(priority_score, 3),
                "direction": row["direction"].strip(),
            })

    rows = sorted(
        rows,
        key=lambda row: (
            row["renamed_miRNA"],
            row["category"],
            -row["priority_score"],
            -row["n_sites"],
            row["gene_symbol"],
        ),
    )

    with out_all.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    top_rows = []
    seen = defaultdict(int)

    for row in rows:
        key = (row["renamed_miRNA"], row["category"])

        if seen[key] < 2:
            top_rows.append(row)
            seen[key] += 1

    with out_top.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(top_rows)

    print("All BrumiR biological targets:", len(rows))
    print("Top BrumiR biological targets:", len(top_rows))
    print("Output:", out_top)


if __name__ == "__main__":
    main()
