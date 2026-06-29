#!/usr/bin/env python3
# ============================================================
# build_cytoscape_network_table_final.py
# Description:
# Builds a final integrated miRNA–mRNA interaction table for
# Cytoscape, combining:
#
# 1. BrumiR candidates (all retained BrumiR miRNAs)
#    - top 5 targets from the dataset-derived target set
#    - top 5 targets from the curated sarcopenia gene set
#
# 2. miRDeep2 candidates
#    - top 5 most Up_in_Sedentary miRNAs
#    - top 5 most Up_in_Active miRNAs
#    - for each selected miRNA:
#         top 5 targets from the dataset-derived target set
#         top 5 targets from the curated sarcopenia gene set
#
# Direction (Up_in_Active / Up_in_Sedentary) is always taken
# from the global miRNA summary table, so it is propagated
# consistently across all predicted targets, including the
# sarcopenia-specific predictions.
#
# Inputs:
#   ../output/seed_matches_36genes_collapsed.tsv
#   ../output/brumir_biological_targets.tsv
#   ../output/seed_matches_36genes_miRDeep2.tsv
#   ../output/seed_matches_sarcopenia.tsv
#   ../input/input_mRNA_36/all_miRNA_seed_summary.tsv
#
# Output:
#   ../output/cytoscape_network_table.tsv
# ============================================================

from pathlib import Path
import csv
from collections import defaultdict

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
OUTPUT = BASE / "output"
INPUT = BASE / "input" / "input_mRNA_36"

BRUMIR_DATASET = OUTPUT / "seed_matches_36genes_collapsed.tsv"
BRUMIR_SARC = OUTPUT / "brumir_biological_targets.tsv"

MIRDEEP_DATASET = OUTPUT / "seed_matches_36genes_miRDeep2.tsv"
MIRDEEP_SARC = OUTPUT / "seed_matches_sarcopenia.tsv"

MIRNA_INFO = INPUT / "all_miRNA_seed_summary.tsv"

OUT = OUTPUT / "cytoscape_network_table.tsv"

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def read_table(path):
    with open(path, encoding="utf-8", errors="ignore") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def safe_int(x, default=0):
    try:
        return int(x)
    except:
        try:
            return int(float(x))
        except:
            return default

def safe_float(x, default=0.0):
    try:
        return float(x)
    except:
        return default

def top_n_by_score(rows, n=5):
    """
    Rank edges by:
    1. number of seed-matched sites
    2. absolute miRNA log2FC
    """
    return sorted(
        rows,
        key=lambda r: (
            safe_int(r.get("n_sites", 0)),
            abs(safe_float(r.get("log2FoldChange", 0)))
        ),
        reverse=True
    )[:n]

# ------------------------------------------------------------
# Global miRNA metadata
# ------------------------------------------------------------

def load_mirna_metadata():
    """
    Reads the global miRNA table and returns:
    - direction by miRNA
    - log2FC by miRNA
    - source by miRNA
    """
    direction = {}
    log2fc = {}
    source_map = {}

    with open(MIRNA_INFO, encoding="utf-8", errors="ignore") as f:
        first = True
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if first:
                first = False
                continue  # skip contaminated header if present

            parts = line.split("\t")
            if len(parts) < 7:
                continue

            source, _, mirna, _, fc, _, dirn = parts[:7]

            mirna = mirna.strip()
            source = source.strip()
            dirn = dirn.strip()

            direction[mirna] = dirn
            log2fc[mirna] = safe_float(fc, 0.0)
            source_map[mirna] = source

    return direction, log2fc, source_map

def get_mirdeep_top(direction_map, log2fc_map, source_map):
    """
    Select top 5 miRDeep2 Up_in_Sedentary and top 5 miRDeep2 Up_in_Active.
    """
    up_sedentary = []
    up_active = []

    for mirna, source in source_map.items():
        if source != "miRDeep2":
            continue

        d = direction_map.get(mirna, "")
        fc = log2fc_map.get(mirna, 0.0)

        if d == "Up_in_Sedentary":
            up_sedentary.append((mirna, fc))
        elif d == "Up_in_Active":
            up_active.append((mirna, fc))

    up_sedentary = sorted(up_sedentary, key=lambda x: x[1], reverse=True)[:5]
    up_active = sorted(up_active, key=lambda x: x[1])[:5]

    selected = set([m for m, _ in up_sedentary + up_active])
    return selected

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    direction_map, log2fc_map, source_map = load_mirna_metadata()

    edges = []
    seen_edges = set()

    # -------------------------
    # BrumiR: all candidates
    # -------------------------
    br_dataset = read_table(BRUMIR_DATASET)
    br_sarc = read_table(BRUMIR_SARC)

    brumir_groups = defaultdict(lambda: {"dataset26": [], "sarcopenia": []})

    for r in br_dataset:
        mirna = r["renamed_miRNA"].strip()
        if source_map.get(mirna, "") == "BrumiR":
            brumir_groups[mirna]["dataset26"].append(r)

    for r in br_sarc:
        mirna = r["renamed_miRNA"].strip()
        if source_map.get(mirna, "") == "BrumiR":
            brumir_groups[mirna]["sarcopenia"].append(r)

    for mirna, subsets in brumir_groups.items():
        for dataset_type in ["dataset26", "sarcopenia"]:
            top = top_n_by_score(subsets[dataset_type], n=5)

            for r in top:
                gene = r["gene_symbol"].strip()
                edge_key = (mirna, gene, "BrumiR", dataset_type)
                if edge_key in seen_edges:
                    continue
                seen_edges.add(edge_key)

                edges.append({
                    "miRNA": mirna,
                    "target_gene": gene,
                    "source": "BrumiR",
                    "dataset_type": dataset_type,
                    "n_sites": r.get("n_sites", ""),
                    "log2FC": log2fc_map.get(mirna, ""),
                    "miRNA_direction": direction_map.get(mirna, "")
                })

    # -------------------------
    # miRDeep2: top 5 up + top 5 down
    # -------------------------
    selected_mirdeep = get_mirdeep_top(direction_map, log2fc_map, source_map)

    md_dataset = read_table(MIRDEEP_DATASET)
    md_sarc = read_table(MIRDEEP_SARC)

    mirdeep_groups = defaultdict(lambda: {"dataset26": [], "sarcopenia": []})

    for r in md_dataset:
        mirna = r["renamed_miRNA"].strip()
        if mirna in selected_mirdeep:
            mirdeep_groups[mirna]["dataset26"].append(r)

    for r in md_sarc:
        mirna = r["renamed_miRNA"].strip()
        if mirna in selected_mirdeep and r["source"].strip() == "miRDeep2":
            mirdeep_groups[mirna]["sarcopenia"].append(r)

    for mirna, subsets in mirdeep_groups.items():
        for dataset_type in ["dataset26", "sarcopenia"]:
            top = top_n_by_score(subsets[dataset_type], n=5)

            for r in top:
                gene = r["gene_symbol"].strip()
                edge_key = (mirna, gene, "miRDeep2", dataset_type)
                if edge_key in seen_edges:
                    continue
                seen_edges.add(edge_key)

                edges.append({
                    "miRNA": mirna,
                    "target_gene": gene,
                    "source": "miRDeep2",
                    "dataset_type": dataset_type,
                    "n_sites": r.get("n_sites", ""),
                    "log2FC": log2fc_map.get(mirna, ""),
                    "miRNA_direction": direction_map.get(mirna, "")
                })

    # -------------------------
    # Write output
    # -------------------------
    edges = sorted(
        edges,
        key=lambda x: (
            x["source"],
            x["miRNA"],
            x["dataset_type"],
            -safe_int(x["n_sites"], 0),
            x["target_gene"]
        )
    )

    with open(OUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "miRNA",
                "target_gene",
                "source",
                "dataset_type",
                "n_sites",
                "log2FC",
                "miRNA_direction"
            ],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(edges)

    print("Edges written:", len(edges))
    print("Output:", OUT)

if __name__ == "__main__":
    main()
