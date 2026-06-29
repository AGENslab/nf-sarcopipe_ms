#!/usr/bin/env python3
# ============================================================
# prepare_union_coherent_target_genes.py
# Description:
# Extracts the UNION of coherent target genes from the
# seed-based miRNA–mRNA prediction table.
# Input:
#   ../output/seed_matches_coherent.tsv
# Outputs:
#   ../output/coherent_target_genes_union.txt
#   ../output/coherent_target_genes_BrumiR.txt
#   ../output/coherent_target_genes_miRDeep2.txt
# ============================================================

from pathlib import Path
import csv

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
OUTPUT = BASE / "output"

COHERENT = OUTPUT / "seed_matches_coherent.tsv"
OUT_UNION = OUTPUT / "coherent_target_genes_union.txt"
OUT_BRUMIR = OUTPUT / "coherent_target_genes_BrumiR.txt"
OUT_MIRDEEP2 = OUTPUT / "coherent_target_genes_miRDeep2.txt"

def main():
    br = set()
    md = set()
    union = set()

    with open(COHERENT, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            source = r["source"].strip()
            gene = r["gene_symbol"].strip()
            union.add(gene)
            if source == "BrumiR":
                br.add(gene)
            elif source == "miRDeep2":
                md.add(gene)

    with open(OUT_BRUMIR, "w") as f:
        for g in sorted(br):
            f.write(g + "\n")

    with open(OUT_MIRDEEP2, "w") as f:
        for g in sorted(md):
            f.write(g + "\n")

    with open(OUT_UNION, "w") as f:
        for g in sorted(union):
            f.write(g + "\n")

    print("BrumiR coherent genes:", len(br))
    print("miRDeep2 coherent genes:", len(md))
    print("Union coherent genes:", len(union))
    print("Written:", OUT_UNION)

if __name__ == "__main__":
    main()
