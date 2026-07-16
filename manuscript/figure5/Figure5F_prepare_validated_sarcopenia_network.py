#!/usr/bin/env python3
"""
Figure 5E and Supplementary Figure 4C
Validated musculoskeletal miRNA-target network

Description
-----------
Combines:

1. BrumiR-derived de novo miRNA-target interactions supported by miRanda.
2. miRDeep2 annotated miRNA-target interactions supported by multiMiR.

The analysis is centered on the genes targeted by the validated de novo
miRNAs. All de novo interactions are retained. Annotated interactions
are ranked independently within each gene, and the requested number of
top annotated pairs per gene is selected.

No miRNA IDs, genes, interaction counts, or results are hardcoded.

Ranking of annotated interactions
---------------------------------
1. lower adjusted p-value
2. higher number of supporting multiMiR databases
3. higher absolute miRNA log2 fold change
4. higher number of seed-matched sites
5. miRNA identifier

Outputs
-------
--out_pairs
    Table used to generate Figure 5E.

--out_edges
    Cytoscape edge table for Supplementary Figure 4C.

--out_nodes
    Cytoscape node table for Supplementary Figure 4C.

--out_summary
    Summary of interactions, miRNAs, and genes.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


PAIR_FIELDS = [
    "gene_symbol",
    "category",
    "renamed_miRNA",
    "original_id",
    "source",
    "miRNA_type",
    "seed",
    "seed_rc",
    "n_sites",
    "site_positions",
    "log2FoldChange",
    "padj",
    "direction",
    "validated_by",
    "validation_database_list",
    "validation_database_count",
    "miranda_score",
    "miranda_energy",
]

EDGE_FIELDS = [
    "source_node",
    "target_node",
    "interaction",
    "miRNA_source",
    "miRNA_type",
    "miRNA_direction",
    "gene_category",
    "seed",
    "seed_rc",
    "n_sites",
    "site_positions",
    "miRNA_log2FoldChange",
    "miRNA_padj",
    "validated_by",
    "validation_database_list",
    "validation_database_count",
    "miranda_score",
    "miranda_energy",
]

NODE_FIELDS = [
    "node_id",
    "node_type",
    "miRNA_source",
    "miRNA_type",
    "direction",
    "gene_category",
    "log2FoldChange",
    "padj",
    "degree",
    "validated_edge_count",
    "validation_methods",
]

SUMMARY_FIELDS = [
    "metric",
    "value",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare validated musculoskeletal interaction tables for "
            "Figure 5E and Supplementary Figure 4C."
        )
    )
    parser.add_argument(
        "--denovo_validated",
        required=True,
        help="miRanda-supported de novo interaction TSV.",
    )
    parser.add_argument(
        "--known_validated",
        required=True,
        help="multiMiR-supported annotated interaction TSV.",
    )
    parser.add_argument(
        "--top_known_per_gene",
        required=True,
        type=int,
        help="Number of annotated interactions retained per de novo target gene.",
    )
    parser.add_argument("--out_pairs", required=True)
    parser.add_argument("--out_edges", required=True)
    parser.add_argument("--out_nodes", required=True)
    parser.add_argument("--out_summary", required=True)
    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")

    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def as_float(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_denovo_rows(
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    normalized = []

    for row in rows:
        mirna = row.get("renamed_miRNA", "").strip()
        gene = row.get("gene_symbol", "").strip()

        if not mirna or not gene:
            continue

        normalized.append({
            "gene_symbol": gene,
            "category": row.get("category", "").strip(),
            "renamed_miRNA": mirna,
            "original_id": row.get("original_id", "").strip(),
            "source": "BrumiR",
            "miRNA_type": "de_novo",
            "seed": row.get("seed", "").strip(),
            "seed_rc": row.get("seed_rc", "").strip(),
            "n_sites": row.get("n_sites", "").strip(),
            "site_positions": row.get("site_positions", "").strip(),
            "log2FoldChange": row.get("log2FoldChange", "").strip(),
            "padj": row.get("padj", "").strip(),
            "direction": row.get("direction", "").strip(),
            "validated_by": "miRanda",
            "validation_database_list": "miRanda",
            "validation_database_count": "1",
            "miranda_score": row.get("miranda_score", "").strip(),
            "miranda_energy": row.get("miranda_energy", "").strip(),
        })

    return normalized


def normalize_known_rows(
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    normalized = []

    for row in rows:
        mirna = row.get("renamed_miRNA", "").strip()
        gene = row.get("gene_symbol", "").strip()

        if not mirna or not gene:
            continue

        normalized.append({
            "gene_symbol": gene,
            "category": row.get("category", "").strip(),
            "renamed_miRNA": mirna,
            "original_id": row.get("original_id", "").strip(),
            "source": "miRDeep2",
            "miRNA_type": "annotated",
            "seed": row.get("seed", "").strip(),
            "seed_rc": row.get("seed_rc", "").strip(),
            "n_sites": row.get("n_sites", "").strip(),
            "site_positions": row.get("site_positions", "").strip(),
            "log2FoldChange": row.get("log2FoldChange", "").strip(),
            "padj": row.get("padj", "").strip(),
            "direction": row.get("direction", "").strip(),
            "validated_by": "multiMiR",
            "validation_database_list": (
                row.get("validation_database_list", "").strip()
            ),
            "validation_database_count": (
                row.get("validation_database_count", "").strip()
            ),
            "miranda_score": "",
            "miranda_energy": "",
        })

    return normalized


def known_rank_key(row: dict[str, str]) -> tuple:
    padj = as_float(
        row.get("padj", ""),
        float("inf"),
    )

    database_count = as_int(
        row.get("validation_database_count", ""),
        0,
    )

    abs_log2fc = abs(
        as_float(
            row.get("log2FoldChange", ""),
            0.0,
        )
    )

    n_sites = as_int(
        row.get("n_sites", ""),
        0,
    )

    return (
        padj,
        -database_count,
        -abs_log2fc,
        -n_sites,
        row.get("renamed_miRNA", ""),
    )


def select_known_per_gene(
    rows: list[dict[str, str]],
    target_genes: set[str],
    top_n: int,
) -> list[dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)

    for row in rows:
        gene = row["gene_symbol"]

        if gene in target_genes:
            grouped[gene].append(row)

    selected = []

    for gene in sorted(grouped):
        ranked = sorted(
            grouped[gene],
            key=known_rank_key,
        )

        selected.extend(
            ranked[:top_n]
        )

    return selected


def build_edges(
    pairs: list[dict[str, str]],
) -> list[dict[str, str]]:
    edges = []

    for row in pairs:
        edges.append({
            "source_node": row["renamed_miRNA"],
            "target_node": row["gene_symbol"],
            "interaction": "targets",
            "miRNA_source": row["source"],
            "miRNA_type": row["miRNA_type"],
            "miRNA_direction": row["direction"],
            "gene_category": row["category"],
            "seed": row["seed"],
            "seed_rc": row["seed_rc"],
            "n_sites": row["n_sites"],
            "site_positions": row["site_positions"],
            "miRNA_log2FoldChange": row["log2FoldChange"],
            "miRNA_padj": row["padj"],
            "validated_by": row["validated_by"],
            "validation_database_list": (
                row["validation_database_list"]
            ),
            "validation_database_count": (
                row["validation_database_count"]
            ),
            "miranda_score": row["miranda_score"],
            "miranda_energy": row["miranda_energy"],
        })

    return edges


def build_nodes(
    edges: list[dict[str, str]],
) -> list[dict[str, object]]:
    node_data: dict[str, dict[str, object]] = {}

    for edge in edges:
        mirna = edge["source_node"]
        gene = edge["target_node"]
        method = edge["validated_by"]

        if mirna not in node_data:
            node_data[mirna] = {
                "node_id": mirna,
                "node_type": "miRNA",
                "miRNA_source": edge["miRNA_source"],
                "miRNA_type": edge["miRNA_type"],
                "direction": edge["miRNA_direction"],
                "gene_category": "",
                "log2FoldChange": edge["miRNA_log2FoldChange"],
                "padj": edge["miRNA_padj"],
                "degree": 0,
                "validated_edge_count": 0,
                "validation_methods": set(),
            }

        if gene not in node_data:
            node_data[gene] = {
                "node_id": gene,
                "node_type": "gene",
                "miRNA_source": "",
                "miRNA_type": "",
                "direction": "",
                "gene_category": edge["gene_category"],
                "log2FoldChange": "",
                "padj": "",
                "degree": 0,
                "validated_edge_count": 0,
                "validation_methods": set(),
            }

        node_data[mirna]["degree"] += 1
        node_data[gene]["degree"] += 1

        node_data[mirna]["validated_edge_count"] += 1
        node_data[gene]["validated_edge_count"] += 1

        node_data[mirna]["validation_methods"].add(method)
        node_data[gene]["validation_methods"].add(method)

    nodes = []

    for node_id in sorted(node_data):
        row = node_data[node_id]

        row["validation_methods"] = ";".join(
            sorted(row["validation_methods"])
        )

        nodes.append(row)

    return nodes


def main() -> None:
    args = parse_args()

    if args.top_known_per_gene < 0:
        raise SystemExit(
            "ERROR: --top_known_per_gene cannot be negative"
        )

    denovo_path = Path(args.denovo_validated)
    known_path = Path(args.known_validated)

    validate_file(
        denovo_path,
        "miRanda-supported de novo table",
    )
    validate_file(
        known_path,
        "multiMiR-supported annotated table",
    )

    denovo_rows = normalize_denovo_rows(
        read_tsv(denovo_path)
    )

    known_rows = normalize_known_rows(
        read_tsv(known_path)
    )

    if not denovo_rows:
        raise SystemExit(
            "ERROR: no validated de novo interactions were loaded"
        )

    denovo_target_genes = {
        row["gene_symbol"]
        for row in denovo_rows
    }

    selected_known = select_known_per_gene(
        known_rows,
        denovo_target_genes,
        args.top_known_per_gene,
    )

    combined_pairs = denovo_rows + selected_known

    combined_pairs.sort(
        key=lambda row: (
            row["category"],
            row["gene_symbol"],
            row["miRNA_type"] != "de_novo",
            known_rank_key(row),
        )
    )

    edges = build_edges(
        combined_pairs
    )

    nodes = build_nodes(
        edges
    )

    denovo_mirnas = {
        row["renamed_miRNA"]
        for row in denovo_rows
    }

    known_mirnas = {
        row["renamed_miRNA"]
        for row in selected_known
    }

    summary = [
        {
            "metric": "de_novo_validated_pairs",
            "value": len(denovo_rows),
        },
        {
            "metric": "de_novo_miRNAs",
            "value": len(denovo_mirnas),
        },
        {
            "metric": "de_novo_target_genes",
            "value": len(denovo_target_genes),
        },
        {
            "metric": "selected_annotated_pairs",
            "value": len(selected_known),
        },
        {
            "metric": "selected_annotated_miRNAs",
            "value": len(known_mirnas),
        },
        {
            "metric": "combined_edges",
            "value": len(edges),
        },
        {
            "metric": "combined_nodes",
            "value": len(nodes),
        },
    ]

    write_tsv(
        Path(args.out_pairs),
        combined_pairs,
        PAIR_FIELDS,
    )

    write_tsv(
        Path(args.out_edges),
        edges,
        EDGE_FIELDS,
    )

    write_tsv(
        Path(args.out_nodes),
        nodes,
        NODE_FIELDS,
    )

    write_tsv(
        Path(args.out_summary),
        summary,
        SUMMARY_FIELDS,
    )

    print("===== VALIDATED MUSCULOSKELETAL NETWORK =====")
    print("Validated de novo pairs:", len(denovo_rows))
    print("De novo miRNAs:", len(denovo_mirnas))
    print("De novo target genes:", len(denovo_target_genes))
    print("Selected annotated pairs:", len(selected_known))
    print("Selected annotated miRNAs:", len(known_mirnas))
    print("Combined network edges:", len(edges))
    print("Combined network nodes:", len(nodes))
    print("Written:", args.out_pairs)
    print("Written:", args.out_edges)
    print("Written:", args.out_nodes)
    print("Written:", args.out_summary)


if __name__ == "__main__":
    main()
