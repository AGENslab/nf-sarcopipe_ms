#!/usr/bin/env python3

"""
Figure 5E – Annotate sarcopenia-related target genes

Description
-----------
This script annotates and prioritizes predicted miRNA target genes
from a curated sarcopenia/exercise-related gene set. Genes are grouped
into biologically meaningful categories and ranked based on:
1. number of supporting miRNAs
2. total number of seed-matched sites
3. mean absolute miRNA log2 fold change

Inputs
------
--input
    TSV file with sarcopenia-related seed matches.

--out_all
    Output TSV with all annotated sarcopenia targets.

--out_top
    Output TSV with top 20 prioritized sarcopenia targets.

Outputs
-------
annotated_sarcopenia_targets.tsv
annotated_sarcopenia_targets_top20.tsv

Usage
-----
python3 Figure5E_annotate_sarcopenia_targets.py \\
  --input seed_matches_sarcopenia.tsv \\
  --out_all annotated_sarcopenia_targets.tsv \\
  --out_top annotated_sarcopenia_targets_top20.tsv
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
        description="Annotate and prioritize sarcopenia-related miRNA target genes."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="TSV file with sarcopenia-related seed matches.",
    )
    parser.add_argument(
        "--out_all",
        required=True,
        help="Output TSV with all annotated sarcopenia targets.",
    )
    parser.add_argument(
        "--out_top",
        required=True,
        help="Output TSV with top 20 prioritized sarcopenia targets.",
    )

    return parser.parse_args()


def main() -> None:
    """Annotate and prioritize sarcopenia-related target genes."""
    args = parse_args()

    input_file = validate_file(args.input, "Sarcopenia seed-match table")
    out_all = prepare_output(args.out_all)
    out_top = prepare_output(args.out_top)

    genes = defaultdict(lambda: {
        "category": "Other",
        "supporting_miRNAs": set(),
        "sources": set(),
        "total_sites": 0,
        "abs_log2fc_values": [],
        "pairs": [],
    })

    with input_file.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_cols = {
            "gene_symbol",
            "renamed_miRNA",
            "source",
            "n_sites",
            "log2FoldChange",
        }
        missing_cols = required_cols - set(reader.fieldnames or [])

        if missing_cols:
            raise SystemExit(
                "ERROR: Input table missing required columns: "
                f"{', '.join(sorted(missing_cols))}"
            )

        for row in reader:
            gene = row["gene_symbol"].strip()
            mirna = row["renamed_miRNA"].strip()
            source = row["source"].strip()

            try:
                n_sites = int(row["n_sites"])
            except ValueError:
                n_sites = 0

            try:
                abs_fc = abs(float(row["log2FoldChange"]))
            except ValueError:
                abs_fc = 0.0

            genes[gene]["category"] = CATEGORY_MAP.get(gene, "Other")
            genes[gene]["supporting_miRNAs"].add(mirna)
            genes[gene]["sources"].add(source)
            genes[gene]["total_sites"] += n_sites
            genes[gene]["abs_log2fc_values"].append(abs_fc)
            genes[gene]["pairs"].append(f"{mirna} ({source})")

    out_rows = []

    for gene, info in genes.items():
        n_support = len(info["supporting_miRNAs"])
        total_sites = info["total_sites"]

        mean_abs_fc = (
            sum(info["abs_log2fc_values"]) / len(info["abs_log2fc_values"])
            if info["abs_log2fc_values"] else 0.0
        )

        if info["sources"] == {"BrumiR"}:
            support_class = "BrumiR_only"
        elif info["sources"] == {"miRDeep2"}:
            support_class = "miRDeep2_only"
        else:
            support_class = "shared"

        priority_score = n_support + total_sites + mean_abs_fc

        out_rows.append({
            "gene_symbol": gene,
            "category": info["category"],
            "support_class": support_class,
            "n_supporting_miRNAs": n_support,
            "total_n_sites": total_sites,
            "mean_abs_log2FC": round(mean_abs_fc, 3),
            "priority_score": round(priority_score, 3),
            "supporting_pairs": "; ".join(sorted(info["pairs"])),
        })

    out_rows = sorted(
        out_rows,
        key=lambda row: (
            row["category"] == "Other",
            -row["priority_score"],
            -row["n_supporting_miRNAs"],
            -row["total_n_sites"],
        ),
    )

    fieldnames = [
        "gene_symbol",
        "category",
        "support_class",
        "n_supporting_miRNAs",
        "total_n_sites",
        "mean_abs_log2FC",
        "priority_score",
        "supporting_pairs",
    ]

    with out_all.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    prioritized = [row for row in out_rows if row["category"] != "Other"]

    if len(prioritized) < 20:
        prioritized = out_rows[:20]
    else:
        prioritized = prioritized[:20]

    with out_top.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(prioritized)

    print("Total annotated genes:", len(out_rows))
    print("Top prioritized genes written to:", out_top)


if __name__ == "__main__":
    main()
