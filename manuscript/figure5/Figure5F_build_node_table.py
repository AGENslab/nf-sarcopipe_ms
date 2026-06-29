#!/usr/bin/env python3

"""
Figure 5F – Build Cytoscape node table

Description
-----------
This script builds a node annotation table for Cytoscape from the
integrated miRNA–mRNA network table.

Each node is classified as:
- miRNA or gene
- source
- biological category for genes
- direction for miRNAs

Inputs
------
--network
    TSV file containing the integrated Cytoscape network edge table.

--out
    Output TSV node annotation table.

Outputs
-------
cytoscape_node_table.tsv

Usage
-----
python3 Figure5F_build_node_table.py \\
  --network cytoscape_network_table.tsv \\
  --out cytoscape_node_table.tsv
"""

import argparse
import csv
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
    "node_id",
    "node_type",
    "source",
    "category",
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
        description="Build Cytoscape node annotation table from network edge table."
    )

    parser.add_argument(
        "--network",
        required=True,
        help="TSV file containing Cytoscape network edges.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV node annotation table.",
    )

    return parser.parse_args()


def main() -> None:
    """Build Cytoscape node annotation table."""
    args = parse_args()

    network_file = validate_file(args.network, "Cytoscape network table")
    out_file = prepare_output(args.out)

    nodes = {}

    with network_file.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_cols = {
            "miRNA",
            "target_gene",
            "source",
            "dataset_type",
            "miRNA_direction",
        }
        missing_cols = required_cols - set(reader.fieldnames or [])

        if missing_cols:
            raise SystemExit(
                "ERROR: Network table missing required columns: "
                f"{', '.join(sorted(missing_cols))}"
            )

        for row in reader:
            mirna = row["miRNA"].strip()
            gene = row["target_gene"].strip()
            source = row["source"].strip()
            dataset_type = row["dataset_type"].strip()
            direction = row["miRNA_direction"].strip()

            if mirna not in nodes:
                nodes[mirna] = {
                    "node_id": mirna,
                    "node_type": "miRNA",
                    "source": source,
                    "category": "miRNA",
                    "direction": direction,
                }

            if gene not in nodes:
                nodes[gene] = {
                    "node_id": gene,
                    "node_type": "gene",
                    "source": dataset_type,
                    "category": CATEGORY_MAP.get(gene, "Other"),
                    "direction": "",
                }

    with out_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(nodes.values())

    print("Nodes written:", len(nodes))
    print("Output:", out_file)


if __name__ == "__main__":
    main()
