#!/usr/bin/env python3

"""
Figure 5F – Build Cytoscape network table

Description
-----------
This script builds a final integrated miRNA–mRNA interaction table
for Cytoscape by combining dataset-derived and sarcopenia-specific
target predictions from BrumiR and miRDeep2.

Inputs
------
--brumir_dataset
    TSV file with BrumiR dataset-derived collapsed seed matches.

--brumir_sarcopenia
    TSV file with BrumiR sarcopenia/exercise biological targets.

--mirdeep_dataset
    TSV file with miRDeep2 dataset-derived seed matches.

--mirdeep_sarcopenia
    TSV file with miRDeep2 sarcopenia/exercise seed matches.

--mirna_info
    Global miRNA seed summary TSV.

--out
    Output Cytoscape network TSV file.

Outputs
-------
cytoscape_network_table.tsv

Usage
-----
python3 Figure5F_build_network_table.py \
  --brumir_dataset seed_matches_36genes_collapsed.tsv \
  --brumir_sarcopenia brumir_biological_targets.tsv \
  --mirdeep_dataset seed_matches_36genes_miRDeep2.tsv \
  --mirdeep_sarcopenia seed_matches_sarcopenia.tsv \
  --mirna_info all_miRNA_seed_summary.tsv \
  --out cytoscape_network_table.tsv
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


FIELDNAMES = [
    "miRNA",
    "target_gene",
    "source",
    "dataset_type",
    "n_sites",
    "log2FC",
    "miRNA_direction",
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


def read_table(path: Path) -> list:
    """Read a TSV file as a list of dictionaries."""
    with path.open(encoding="utf-8", errors="ignore") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def safe_int(value, default: int = 0) -> int:
    """Safely convert a value to integer."""
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default


def safe_float(value, default: float = 0.0) -> float:
    """Safely convert a value to float."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def top_n_by_score(rows: list, n: int = 5) -> list:
    """
    Rank edges by:
    1. number of seed-matched sites
    2. absolute miRNA log2FC
    """
    return sorted(
        rows,
        key=lambda row: (
            safe_int(row.get("n_sites", 0)),
            abs(safe_float(row.get("log2FoldChange", 0))),
        ),
        reverse=True,
    )[:n]


def load_mirna_metadata(mirna_info: Path) -> tuple:
    """
    Read global miRNA metadata.

    Returns:
    - direction by miRNA
    - log2FC by miRNA
    - source by miRNA
    """
    direction = {}
    log2fc = {}
    source_map = {}

    with mirna_info.open(encoding="utf-8", errors="ignore") as handle:
        first = True

        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if not line:
                continue

            if first:
                first = False
                continue

            fields = line.split("\t")

            if len(fields) < 7:
                continue

            source, _, mirna, _, fc, _, direction_value = fields[:7]

            mirna = mirna.strip()
            source = source.strip()
            direction_value = direction_value.strip()

            direction[mirna] = direction_value
            log2fc[mirna] = safe_float(fc, 0.0)
            source_map[mirna] = source

    return direction, log2fc, source_map


def get_mirdeep_top(direction_map: dict, log2fc_map: dict, source_map: dict) -> set:
    """Select top 5 miRDeep2 Up_in_Sedentary and top 5 Up_in_Active miRNAs."""
    up_sedentary = []
    up_active = []

    for mirna, source in source_map.items():
        if source != "miRDeep2":
            continue

        direction = direction_map.get(mirna, "")
        log2fc = log2fc_map.get(mirna, 0.0)

        if direction == "Up_in_Sedentary":
            up_sedentary.append((mirna, log2fc))
        elif direction == "Up_in_Active":
            up_active.append((mirna, log2fc))

    up_sedentary = sorted(up_sedentary, key=lambda item: item[1], reverse=True)[:5]
    up_active = sorted(up_active, key=lambda item: item[1])[:5]

    return set([mirna for mirna, _ in up_sedentary + up_active])


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build final integrated miRNA-mRNA network table for Cytoscape."
    )

    parser.add_argument("--brumir_dataset", required=True)
    parser.add_argument("--brumir_sarcopenia", required=True)
    parser.add_argument("--mirdeep_dataset", required=True)
    parser.add_argument("--mirdeep_sarcopenia", required=True)
    parser.add_argument("--mirna_info", required=True)
    parser.add_argument("--out", required=True)

    return parser.parse_args()


def main() -> None:
    """Build Cytoscape network table."""
    args = parse_args()

    brumir_dataset = validate_file(args.brumir_dataset, "BrumiR dataset table")
    brumir_sarcopenia = validate_file(args.brumir_sarcopenia, "BrumiR sarcopenia table")
    mirdeep_dataset = validate_file(args.mirdeep_dataset, "miRDeep2 dataset table")
    mirdeep_sarcopenia = validate_file(args.mirdeep_sarcopenia, "miRDeep2 sarcopenia table")
    mirna_info = validate_file(args.mirna_info, "Global miRNA summary table")
    out_file = prepare_output(args.out)

    direction_map, log2fc_map, source_map = load_mirna_metadata(mirna_info)

    edges = []
    seen_edges = set()

    # -------------------------
    # BrumiR candidates
    # -------------------------
    brumir_dataset_rows = read_table(brumir_dataset)
    brumir_sarcopenia_rows = read_table(brumir_sarcopenia)

    brumir_groups = defaultdict(lambda: {"dataset26": [], "sarcopenia": []})

    for row in brumir_dataset_rows:
        mirna = row["renamed_miRNA"].strip()

        if source_map.get(mirna, "") == "BrumiR":
            brumir_groups[mirna]["dataset26"].append(row)

    for row in brumir_sarcopenia_rows:
        mirna = row["renamed_miRNA"].strip()

        if source_map.get(mirna, "") == "BrumiR":
            brumir_groups[mirna]["sarcopenia"].append(row)

    for mirna, subsets in brumir_groups.items():
        for dataset_type in ["dataset26", "sarcopenia"]:
            top_rows = top_n_by_score(subsets[dataset_type], n=5)

            for row in top_rows:
                gene = row["gene_symbol"].strip()
                edge_key = (mirna, gene, "BrumiR", dataset_type)

                if edge_key in seen_edges:
                    continue

                seen_edges.add(edge_key)

                edges.append({
                    "miRNA": mirna,
                    "target_gene": gene,
                    "source": "BrumiR",
                    "dataset_type": dataset_type,
                    "n_sites": row.get("n_sites", ""),
                    "log2FC": log2fc_map.get(mirna, ""),
                    "miRNA_direction": direction_map.get(mirna, ""),
                })

    # -------------------------
    # miRDeep2 candidates
    # -------------------------
    selected_mirdeep = get_mirdeep_top(direction_map, log2fc_map, source_map)

    mirdeep_dataset_rows = read_table(mirdeep_dataset)
    mirdeep_sarcopenia_rows = read_table(mirdeep_sarcopenia)

    mirdeep_groups = defaultdict(lambda: {"dataset26": [], "sarcopenia": []})

    for row in mirdeep_dataset_rows:
        mirna = row["renamed_miRNA"].strip()

        if mirna in selected_mirdeep:
            mirdeep_groups[mirna]["dataset26"].append(row)

    for row in mirdeep_sarcopenia_rows:
        mirna = row["renamed_miRNA"].strip()

        if mirna in selected_mirdeep and row["source"].strip() == "miRDeep2":
            mirdeep_groups[mirna]["sarcopenia"].append(row)

    for mirna, subsets in mirdeep_groups.items():
        for dataset_type in ["dataset26", "sarcopenia"]:
            top_rows = top_n_by_score(subsets[dataset_type], n=5)

            for row in top_rows:
                gene = row["gene_symbol"].strip()
                edge_key = (mirna, gene, "miRDeep2", dataset_type)

                if edge_key in seen_edges:
                    continue

                seen_edges.add(edge_key)

                edges.append({
                    "miRNA": mirna,
                    "target_gene": gene,
                    "source": "miRDeep2",
                    "dataset_type": dataset_type,
                    "n_sites": row.get("n_sites", ""),
                    "log2FC": log2fc_map.get(mirna, ""),
                    "miRNA_direction": direction_map.get(mirna, ""),
                })

    # -------------------------
    # SAVE
    # -------------------------
    edges = sorted(
        edges,
        key=lambda row: (
            row["source"],
            row["miRNA"],
            row["dataset_type"],
            -safe_int(row["n_sites"], 0),
            row["target_gene"],
        ),
    )

    with out_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(edges)

    print("Edges written:", len(edges))
    print("Output:", out_file)


if __name__ == "__main__":
    main()