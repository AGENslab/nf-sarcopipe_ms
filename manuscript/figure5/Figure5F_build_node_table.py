#!/usr/bin/env python3
# ============================================================
# build_cytoscape_node_table.py
# Description:
# Builds a node annotation table for Cytoscape from the
# integrated miRNA–mRNA network table.
#
# Each node is classified as:
#   - miRNA or gene
#   - source (BrumiR / miRDeep2 / predicted / dataset)
#   - biological category (for genes)
#   - direction (for miRNAs)
#
# Input:
#   ../output/cytoscape_network_table.tsv
#
# Output:
#   ../output/cytoscape_node_table.tsv
# ============================================================

from pathlib import Path
import csv
from collections import defaultdict

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
OUTPUT = BASE / "output"

NETWORK = OUTPUT / "cytoscape_network_table.tsv"
OUT = OUTPUT / "cytoscape_node_table.tsv"

# ------------------------------------------------------------
# Category map (same as before, IMPORTANT consistency)
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    nodes = {}

    with open(NETWORK, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for r in reader:
            mirna = r["miRNA"].strip()
            gene = r["target_gene"].strip()
            source = r["source"].strip()
            dataset_type = r["dataset_type"].strip()
            direction = r["miRNA_direction"].strip()

            # -------- miRNA node --------
            if mirna not in nodes:
                nodes[mirna] = {
                    "node_id": mirna,
                    "node_type": "miRNA",
                    "source": source,
                    "category": "miRNA",
                    "direction": direction
                }

            # -------- gene node --------
            if gene not in nodes:
                nodes[gene] = {
                    "node_id": gene,
                    "node_type": "gene",
                    "source": dataset_type,
                    "category": CATEGORY_MAP.get(gene, "Other"),
                    "direction": ""
                }

    # write
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["node_id", "node_type", "source", "category", "direction"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(nodes.values())

    print("Nodes written:", len(nodes))
    print("Output:", OUT)

if __name__ == "__main__":
    main()
