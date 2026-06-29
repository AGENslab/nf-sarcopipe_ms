#!/usr/bin/env python3
# ============================================================
# prepare_prioritized_predicted_target_genes.py
# Description:
# Ranks coherent predicted target genes based on convergent
# miRNA support, number of seed-matched sites, and miRNA
# effect size. This script prepares the input table for
# Figure 5d.
#
# Input:
#   ../output/seed_matches_coherent.tsv
#
# Outputs:
#   ../output/prioritized_predicted_target_genes.tsv
#   ../output/prioritized_predicted_target_genes_top15.tsv
# ============================================================

from pathlib import Path
import csv
from collections import defaultdict

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
OUTPUT = BASE / "output"

INFILE = OUTPUT / "seed_matches_coherent.tsv"
OUT_ALL = OUTPUT / "prioritized_predicted_target_genes.tsv"
OUT_TOP = OUTPUT / "prioritized_predicted_target_genes_top15.tsv"

def main():
    genes = defaultdict(lambda: {
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

            genes[gene]["supporting_miRNAs"].add(mirna)
            genes[gene]["sources"].add(source)
            genes[gene]["total_sites"] += n_sites
            genes[gene]["abs_log2fc_values"].append(abs_fc)
            genes[gene]["pairs"].append(f"{mirna} ({source})")

    out_rows = []
    for gene, info in genes.items():
        n_support = len(info["supporting_miRNAs"])
        total_sites = info["total_sites"]
        mean_abs_fc = sum(info["abs_log2fc_values"]) / len(info["abs_log2fc_values"]) if info["abs_log2fc_values"] else 0.0

        if info["sources"] == {"BrumiR"}:
            category = "BrumiR_only"
        elif info["sources"] == {"miRDeep2"}:
            category = "miRDeep2_only"
        else:
            category = "shared"

        priority_score = n_support + total_sites + mean_abs_fc

        out_rows.append({
            "gene_symbol": gene,
            "category": category,
            "n_supporting_miRNAs": n_support,
            "total_n_sites": total_sites,
            "mean_abs_log2FC": round(mean_abs_fc, 3),
            "priority_score": round(priority_score, 3),
            "supporting_pairs": "; ".join(sorted(info["pairs"]))
        })

    out_rows = sorted(
        out_rows,
        key=lambda x: (x["priority_score"], x["n_supporting_miRNAs"], x["total_n_sites"]),
        reverse=True
    )

    fieldnames = [
        "gene_symbol",
        "category",
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

    top_n = min(15, len(out_rows))
    with open(OUT_TOP, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows[:top_n])

    print("Total prioritized genes:", len(out_rows))
    print("Top genes written to:", OUT_TOP)

if __name__ == "__main__":
    main()
