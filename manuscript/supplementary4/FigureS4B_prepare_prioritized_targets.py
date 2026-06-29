#!/usr/bin/env python3
"""
Figure S4B – Prepare prioritized predicted target genes

Description
-----------
Rank coherent predicted target genes based on convergent miRNA support,
number of seed-matched sites, and miRNA effect size.

This script prepares the input table for Figure S4B.

Inputs
------
--infile
    Coherent seed match table. Expected columns:
    - gene_symbol
    - renamed_miRNA
    - source
    - n_sites
    - log2FoldChange

Outputs
-------
--out_all
    Full prioritized target gene table.

--out_top
    Top-ranked prioritized target gene table.

--top_n
    Number of top genes to export. Default: 15.

Usage
-----
python FigureS4B_prepare_prioritized_targets.py \\
  --infile path/to/seed_matches_coherent.tsv \\
  --out_all path/to/prioritized_predicted_target_genes.tsv \\
  --out_top path/to/prioritized_predicted_target_genes_top15.tsv \\
  --top_n 15
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


FIELDNAMES = [
    "gene_symbol",
    "category",
    "n_supporting_miRNAs",
    "total_n_sites",
    "mean_abs_log2FC",
    "priority_score",
    "supporting_pairs",
]


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare prioritized predicted target genes for Figure S4B."
    )
    parser.add_argument(
        "--infile",
        required=True,
        help="Input coherent seed match TSV.",
    )
    parser.add_argument(
        "--out_all",
        required=True,
        help="Output full prioritized target gene TSV.",
    )
    parser.add_argument(
        "--out_top",
        required=True,
        help="Output top prioritized target gene TSV.",
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=15,
        help="Number of top genes to export. Default: 15.",
    )

    return parser.parse_args()


def validate_input_file(path: Path, label: str) -> None:
    """Validate that an input file exists."""
    if not path.exists():
        sys.exit(f"ERROR: {label} not found: {path}")
    if not path.is_file():
        sys.exit(f"ERROR: {label} is not a file: {path}")


def validate_columns(reader: csv.DictReader, path: Path) -> None:
    """Validate required columns in input table."""
    required_cols = {
        "gene_symbol",
        "renamed_miRNA",
        "source",
        "n_sites",
        "log2FoldChange",
    }

    missing_cols = required_cols - set(reader.fieldnames or [])

    if missing_cols:
        sys.exit(
            "ERROR: Missing required columns in "
            f"{path}: {', '.join(sorted(missing_cols))}"
        )


def read_and_rank_targets(infile: Path):
    """Read coherent seed matches and compute prioritization metrics."""
    genes = defaultdict(
        lambda: {
            "supporting_miRNAs": set(),
            "sources": set(),
            "total_sites": 0,
            "abs_log2fc_values": [],
            "pairs": [],
        }
    )

    with infile.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        validate_columns(reader, infile)

        for row in reader:
            gene = row["gene_symbol"].strip()
            mirna = row["renamed_miRNA"].strip()
            source = row["source"].strip()

            try:
                n_sites = int(row["n_sites"])
            except Exception:
                n_sites = 0

            try:
                abs_fc = abs(float(row["log2FoldChange"]))
            except Exception:
                abs_fc = 0.0

            genes[gene]["supporting_miRNAs"].add(mirna)
            genes[gene]["sources"].add(source)
            genes[gene]["total_sites"] += n_sites
            genes[gene]["abs_log2fc_values"].append(abs_fc)
            genes[gene]["pairs"].append(f"{mirna} ({source})")

    output_rows = []

    for gene, info in genes.items():
        n_support = len(info["supporting_miRNAs"])
        total_sites = info["total_sites"]

        mean_abs_fc = (
            sum(info["abs_log2fc_values"]) / len(info["abs_log2fc_values"])
            if info["abs_log2fc_values"]
            else 0.0
        )

        if info["sources"] == {"BrumiR"}:
            category = "BrumiR_only"
        elif info["sources"] == {"miRDeep2"}:
            category = "miRDeep2_only"
        else:
            category = "shared"

        priority_score = n_support + total_sites + mean_abs_fc

        output_rows.append(
            {
                "gene_symbol": gene,
                "category": category,
                "n_supporting_miRNAs": n_support,
                "total_n_sites": total_sites,
                "mean_abs_log2FC": round(mean_abs_fc, 3),
                "priority_score": round(priority_score, 3),
                "supporting_pairs": "; ".join(sorted(info["pairs"])),
            }
        )

    output_rows = sorted(
        output_rows,
        key=lambda row: (
            row["priority_score"],
            row["n_supporting_miRNAs"],
            row["total_n_sites"],
        ),
        reverse=True,
    )

    return output_rows


def write_tsv(path: Path, rows) -> None:
    """Write rows to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    infile = Path(args.infile)
    out_all = Path(args.out_all)
    out_top = Path(args.out_top)

    validate_input_file(infile, "Input coherent seed match TSV")

    if args.top_n < 1:
        sys.exit("ERROR: --top_n must be greater than or equal to 1.")

    output_rows = read_and_rank_targets(infile)

    if not output_rows:
        sys.exit("ERROR: No prioritized target genes were generated.")

    write_tsv(out_all, output_rows)

    top_n = min(args.top_n, len(output_rows))
    write_tsv(out_top, output_rows[:top_n])

    print("Total prioritized genes:", len(output_rows))
    print("Top genes written to:", out_top)


if __name__ == "__main__":
    main()
