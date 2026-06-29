#!/usr/bin/env python3
# ============================================================
# annotate_sarcopenia_targets.py
# Description:
# Annotates and prioritizes predicted miRNA target genes from a
# curated sarcopenia/exercise-related gene set. Genes are grouped
# into biologically meaningful categories and ranked based on:
#   1) number of supporting miRNAs
#   2) total number of seed-matched sites
#   3) mean absolute miRNA log2 fold change
#
# Input:
#   ../output/seed_matches_sarcopenia.tsv
#
# Outputs:
#   ../output/annotated_sarcopenia_targets.tsv
#   ../output/annotated_sarcopenia_targets_top20.tsv
# ============================================================

from pathlib import Path
import csv
from collections import defaultdict

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
OUTPUT = BASE / "output"

INFILE = OUTPUT / "seed_matches_sarcopenia.tsv"
OUT_ALL = OUTPUT / "annotated_sarcopenia_targets.tsv"
OUT_TOP = OUTPUT / "annotated_sarcopenia_targets_top20.tsv"

# ------------------------------------------------------------
# Biological categories
# ------------------------------------------------------------
CATEGORY_MAP = {
    # Inflammation / immunity
    "IL6": "Inflammation / immunity",
    "TNF": "Inflammation / immunity",
    "NFKB1": "Inflammation / immunity",
    "RELA": "Inflammation / immunity",
    "IL1B": "Inflammation / immunity",
    "CXCL8": "Inflammation / immunity",
    "CCL2": "Inflammation / immunity",

    # Muscle atrophy / growth
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

    # Fibrosis / ECM remodeling
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

    # Senescence / damage
    "SIRT1": "Senescence / damage",
    "SIRT3": "Senescence / damage",
    "CDKN1A": "Senescence / damage",
    "CDKN2A": "Senescence / damage",
    "TP53": "Senescence / damage",
    "LMNB1": "Senescence / damage",
    "GDF11": "Senescence / damage",

    # Autophagy / degradation
    "BECN1": "Autophagy / degradation",
    "MAP1LC3B": "Autophagy / degradation",
    "ATG5": "Autophagy / degradation",
    "ATG7": "Autophagy / degradation",
    "SQSTM1": "Autophagy / degradation",
    "ULK1": "Autophagy / degradation",

    # Myogenesis / muscle structure
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
    genes = defaultdict(lambda: {
        "category": "Other",
        "supporting_miRNAs": set(),
        "sources": set(),
        "total_sites": 0,
        "abs_log2fc_values": [],
        "pairs": []
    })

    with open(INFILE, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            gene = r["gene_symbol"].strip()
            mirna = r["renamed_miRNA"].strip()
            source = r["source"].strip()

            try:
                n_sites = int(r["n_sites"])
            except:
                n_sites = 0

            try:
                abs_fc = abs(float(r["log2FoldChange"]))
            except:
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

        # simple, defendible priority score
        priority_score = n_support + total_sites + mean_abs_fc

        out_rows.append({
            "gene_symbol": gene,
            "category": info["category"],
            "support_class": support_class,
            "n_supporting_miRNAs": n_support,
            "total_n_sites": total_sites,
            "mean_abs_log2FC": round(mean_abs_fc, 3),
            "priority_score": round(priority_score, 3),
            "supporting_pairs": "; ".join(sorted(info["pairs"]))
        })

    # prioritize category genes first, then score
    out_rows = sorted(
        out_rows,
        key=lambda x: (
            x["category"] == "Other",
            -x["priority_score"],
            -x["n_supporting_miRNAs"],
            -x["total_n_sites"]
        )
    )

    fieldnames = [
        "gene_symbol",
        "category",
        "support_class",
        "n_supporting_miRNAs",
        "total_n_sites",
        "mean_abs_log2FC",
        "priority_score",
        "supporting_pairs"
    ]

    with open(OUT_ALL, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    # top 20 prioritized genes, excluding "Other" if possible
    prioritized = [r for r in out_rows if r["category"] != "Other"]
    if len(prioritized) < 20:
        prioritized = out_rows[:20]
    else:
        prioritized = prioritized[:20]

    with open(OUT_TOP, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(prioritized)

    print("Total annotated genes:", len(out_rows))
    print("Top prioritized genes written to:", OUT_TOP)

if __name__ == "__main__":
    main()
