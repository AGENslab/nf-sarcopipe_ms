#!/usr/bin/env python3

"""
Figure 5C – Prepare coherent pairs for lollipop plot

Description
-----------
This script prepares coherent miRNA–mRNA seed-match pairs for the
Figure 5C lollipop plot. It keeps all BrumiR coherent pairs and
selects the top miRDeep2 miRNAs by adjusted p-value separately for
Up_in_Active and Up_in_Sedentary directions.

Inputs
------
--mirna_table
    TSV file with miRNA seed summary information.

--coherent_table
    TSV file with coherent miRNA–mRNA pairs.

--out_lollipop
    Output TSV file prepared for lollipop plotting.

Outputs
-------
coherent_pairs_for_lollipop.tsv

Usage
-----
python3 Figure5C_prepare_coherent_pairs.py \
  --mirna_table all_miRNA_seed_summary.tsv \
  --coherent_table seed_matches_coherent.tsv \
  --out_lollipop coherent_pairs_for_lollipop.tsv
"""

import argparse
import csv
from pathlib import Path


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


def read_mirna_stats(mirna_table: Path) -> dict:
    """Read miRNA statistics indexed by source and renamed miRNA ID."""
    stats = {}

    with mirna_table.open(encoding="utf-8", errors="ignore") as handle:
        next(handle)

        for line in handle:
            fields = line.rstrip("\n\r").split("\t")

            if len(fields) < 7:
                continue

            source, original_id, renamed_miRNA, seed, log2fc, padj, direction = fields[:7]

            stats[(source.strip(), renamed_miRNA.strip())] = {
                "log2FoldChange": float(log2fc),
                "padj": float(padj),
                "direction": direction.strip(),
            }

    return stats


def read_coherent(coherent_table: Path) -> list:
    """Read coherent miRNA–mRNA pairs."""
    rows = []

    with coherent_table.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            rows.append({
                key.strip(): value.strip()
                for key, value in row.items()
            })

    return rows


def select_top_mirdeep2(stats: dict) -> set:
    """Select top miRDeep2 miRNAs by padj for each regulation direction."""
    active = []
    sedentary = []

    for (source, mirna), values in stats.items():
        if source != "miRDeep2":
            continue

        record = (
            mirna,
            values["padj"],
            values["log2FoldChange"],
            values["direction"],
        )

        if values["direction"] == "Up_in_Active":
            active.append(record)
        elif values["direction"] == "Up_in_Sedentary":
            sedentary.append(record)

    active = sorted(active, key=lambda x: x[1])[:10]
    sedentary = sorted(sedentary, key=lambda x: x[1])[:10]

    keep = set([record[0] for record in active] + [record[0] for record in sedentary])

    return keep


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare coherent miRNA–mRNA pairs for Figure 5C lollipop plot."
    )

    parser.add_argument(
        "--mirna_table",
        required=True,
        help="TSV file with miRNA seed summary information.",
    )
    parser.add_argument(
        "--coherent_table",
        required=True,
        help="TSV file with coherent miRNA–mRNA pairs.",
    )
    parser.add_argument(
        "--out_lollipop",
        required=True,
        help="Output TSV file for lollipop plotting.",
    )

    return parser.parse_args()


def main() -> None:
    """Prepare selected coherent pairs for lollipop plotting."""
    args = parse_args()

    mirna_table = validate_file(args.mirna_table, "miRNA seed summary table")
    coherent_table = validate_file(args.coherent_table, "Coherent pairs table")
    out_lollipop = prepare_output(args.out_lollipop)

    stats = read_mirna_stats(mirna_table)
    coherent = read_coherent(coherent_table)

    keep_mirdeep2 = select_top_mirdeep2(stats)

    selected = []

    for row in coherent:
        source = row["source"]
        mirna = row["renamed_miRNA"]

        if source == "BrumiR":
            selected.append(row)
        elif source == "miRDeep2" and mirna in keep_mirdeep2:
            selected.append(row)

    for row in selected:
        key = (row["source"], row["renamed_miRNA"])
        row["abs_log2FC"] = abs(float(stats[key]["log2FoldChange"]))

    selected.sort(
        key=lambda row: (
            row["source"],
            row["direction"],
            float(row["padj"]),
        )
    )

    with out_lollipop.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(selected[0].keys()) if selected else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    print("Selected coherent pairs for lollipop:", len(selected))
    print("Output:", out_lollipop)


if __name__ == "__main__":
    main()