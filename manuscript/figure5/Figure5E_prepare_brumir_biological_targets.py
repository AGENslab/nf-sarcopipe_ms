#!/usr/bin/env python3
# ============================================================
# prepare_brumir_biological_targets.py
# Description:
# Filters BrumiR-derived predicted targets from the curated
# sarcopenia/exercise gene set, assigns biological categories,
# and prioritizes top targets per BrumiR candidate for Figure 5d.
#
# Input:
#   ../output/seed_matches_sarcopenia.tsv
#
# Outputs:
#   ../output/brumir_biological_targets.tsv
#   ../output/brumir_biological_targets_top.tsv
# ============================================================

from pathlib import Path
import csv
from collections import defaultdict

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
OUTPUT = BASE / "output"

INFILE = OUTPUT / "seed_matches_sarcopenia.tsv"
OUT_ALL = OUTPUT / "brumir_biological_targets.tsv"
OUT_TOP = OUTPUT / "brumir_biological_targets_top.tsv"

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

def main():
    rows = []

    with open(INFILE, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["source"].strip() != "BrumiR":
                continue

            gene = r["gene_symbol"].strip()
            mirna = r["renamed_miRNA"].strip()

            if CATEGORY_MAP.get(gene, "Other") == "Other":
                continue

            try:
                n_sites = int(r["n_sites"])
            except:
                n_sites = 0

            try:
                abs_fc = abs(float(r["log2FoldChange"]))
            except:
                abs_fc = 0.0

            priority_score = n_sites + abs_fc

            rows.append({
                "renamed_miRNA": mirna,
                "gene_symbol": gene,
                "category": CATEGORY_MAP[gene],
                "n_sites": n_sites,
                "abs_log2FC": round(abs_fc, 3),
                "priority_score": round(priority_score, 3),
                "direction": r["direction"].strip()
            })

    rows = sorted(
        rows,
        key=lambda x: (x["renamed_miRNA"], x["category"], -x["priority_score"], -x["n_sites"], x["gene_symbol"])
    )

    with open(OUT_ALL, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["renamed_miRNA", "gene_symbol", "category", "n_sites", "abs_log2FC", "priority_score", "direction"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    # top 2 genes per category per miRNA
    top_rows = []
    seen = defaultdict(int)

    for r in rows:
        key = (r["renamed_miRNA"], r["category"])
        if seen[key] < 2:
            top_rows.append(r)
            seen[key] += 1

    with open(OUT_TOP, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["renamed_miRNA", "gene_symbol", "category", "n_sites", "abs_log2FC", "priority_score", "direction"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(top_rows)

    print("All BrumiR biological targets:", len(rows))
    print("Top BrumiR biological targets:", len(top_rows))
    print("Output:", OUT_TOP)

if __name__ == "__main__":
    main()
