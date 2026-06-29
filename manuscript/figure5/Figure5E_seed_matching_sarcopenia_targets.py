#!/usr/bin/env python3
# ============================================================
# run_seed_matching_sarcopenia_targets.py
# Description:
# Performs seed-based miRNA target prediction against a curated
# sarcopenia/exercise-related gene set using available 3'UTR
# sequences. Applies the same methodology to BrumiR and miRDeep2.
#
# Inputs:
#   ../input/input_sarcopenia_gene_set.txt
#   ../input/input_mRNA_36/all_miRNA_seed_summary.tsv
#   ../output/DEG_36_3UTR_final.tsv
#
# Output:
#   ../output/seed_matches_sarcopenia.tsv
# ============================================================

from pathlib import Path
import csv

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
INPUT = BASE / "input"
OUTPUT = BASE / "output"

GENES = INPUT / "input_sarcopenia_gene_set.txt"
MIRNA = INPUT / "input_mRNA_36" / "all_miRNA_seed_summary.tsv"
UTR = OUTPUT / "sarcopenia_gene_set_3UTR.tsv"

OUT = OUTPUT / "seed_matches_sarcopenia.tsv"

def revcomp_rna(seed: str) -> str:
    comp = str.maketrans({"A": "U", "U": "A", "G": "C", "C": "G", "T": "A"})
    s = seed.upper().replace("T", "U")
    return s.translate(comp)[::-1]

def load_mirnas():
    mirnas = []
    with open(MIRNA, encoding="utf-8", errors="ignore") as f:
        first = True
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if first:
                first = False
                continue  # saltar header contaminado
            parts = line.split("\t")
            if len(parts) < 7:
                continue

            source, original_id, renamed_miRNA, seed, log2fc, padj, direction = parts[:7]
            seed = seed.strip().upper().replace("T", "U")
            if not seed:
                continue

            try:
                log2fc_val = float(log2fc)
            except:
                log2fc_val = 0.0

            mirnas.append({
                "source": source.strip(),
                "original_id": original_id.strip(),
                "renamed_miRNA": renamed_miRNA.strip(),
                "seed": seed,
                "seed_rc": revcomp_rna(seed),
                "log2FoldChange": log2fc_val,
                "padj": padj.strip(),
                "direction": direction.strip()
            })
    return mirnas

def load_utrs():
    utrs = {}
    with open(UTR, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            gene = r["gene_symbol"].strip()
            seq = r["utr_3"].strip().upper().replace("T", "U")
            if gene and seq:
                utrs[gene] = seq
    return utrs

def load_genes():
    genes = set()
    with open(GENES, encoding="utf-8", errors="ignore") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.add(g)
    return genes

def find_all_positions(seq: str, motif: str):
    positions = []
    start = 0
    while True:
        i = seq.find(motif, start)
        if i == -1:
            break
        positions.append(i + 1)  # 1-based
        start = i + 1
    return positions

def main():
    mirnas = load_mirnas()
    utrs = load_utrs()
    genes = load_genes()

    out = []

    for gene in sorted(genes):
        if gene not in utrs:
            continue
        seq = utrs[gene]

        for m in mirnas:
            positions = find_all_positions(seq, m["seed_rc"])
            if positions:
                out.append({
                    "gene_symbol": gene,
                    "renamed_miRNA": m["renamed_miRNA"],
                    "source": m["source"],
                    "seed": m["seed"],
                    "seed_rc": m["seed_rc"],
                    "n_sites": len(positions),
                    "site_positions": ",".join(map(str, positions)),
                    "log2FoldChange": m["log2FoldChange"],
                    "direction": m["direction"]
                })

    with open(OUT, "w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "gene_symbol", "renamed_miRNA", "source", "seed", "seed_rc",
            "n_sites", "site_positions", "log2FoldChange", "direction"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(out)

    print("miRNAs loaded:", len(mirnas))
    print("genes in sarcopenia set:", len(genes))
    print("genes with available 3'UTR:", len(utrs))
    print("matches found:", len(out))
    print("output:", OUT)

if __name__ == "__main__":
    main()