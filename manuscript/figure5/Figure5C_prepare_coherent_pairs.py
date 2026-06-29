#!/usr/bin/env python3
from pathlib import Path
import csv

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
INPUT = BASE / "input" / "input_mRNA_36"
OUTPUT = BASE / "output"

MIRNA_TABLE = INPUT / "all_miRNA_seed_summary.tsv"
COHERENT_TABLE = OUTPUT / "seed_matches_coherent.tsv"

OUT_LOLLIPOP = OUTPUT / "coherent_pairs_for_lollipop.tsv"

def read_mirna_stats():
    stats = {}
    with open(MIRNA_TABLE, encoding="utf-8", errors="ignore") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 7:
                continue
            source, original_id, renamed_miRNA, seed, log2fc, padj, direction = parts[:7]
            stats[(source.strip(), renamed_miRNA.strip())] = {
                "log2FoldChange": float(log2fc),
                "padj": float(padj),
                "direction": direction.strip()
            }
    return stats

def read_coherent():
    rows = []
    with open(COHERENT_TABLE, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append({k.strip(): v.strip() for k, v in r.items()})
    return rows

def select_top_mirdeep2(stats):
    active = []
    sedentary = []

    for (source, mirna), vals in stats.items():
        if source != "miRDeep2":
            continue
        rec = (mirna, vals["padj"], vals["log2FoldChange"], vals["direction"])
        if vals["direction"] == "Up_in_Active":
            active.append(rec)
        elif vals["direction"] == "Up_in_Sedentary":
            sedentary.append(rec)

    active = sorted(active, key=lambda x: x[1])[:10]
    sedentary = sorted(sedentary, key=lambda x: x[1])[:10]

    keep = set([x[0] for x in active] + [x[0] for x in sedentary])
    return keep

def main():
    stats = read_mirna_stats()
    coherent = read_coherent()

    keep_mirdeep2 = select_top_mirdeep2(stats)

    selected = []
    for r in coherent:
        source = r["source"]
        mirna = r["renamed_miRNA"]

        if source == "BrumiR":
            selected.append(r)
        elif source == "miRDeep2" and mirna in keep_mirdeep2:
            selected.append(r)

    # add abs logFC for plotting
    for r in selected:
        key = (r["source"], r["renamed_miRNA"])
        r["abs_log2FC"] = abs(float(stats[key]["log2FoldChange"]))

    # sort
    selected.sort(key=lambda x: (x["source"], x["direction"], float(x["padj"])))

    with open(OUT_LOLLIPOP, "w", encoding="utf-8", newline="") as f:
        fieldnames = list(selected[0].keys()) if selected else []
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected)

    print("Selected coherent pairs for lollipop:", len(selected))
    print("Output:", OUT_LOLLIPOP)

if __name__ == "__main__":
    main()
