#!/usr/bin/env python3
"""
Figure 5E – Seed matching against musculoskeletal target genes

Performs canonical 7-mer seed matching for the final nf-Sarcopipe
miRNA catalog against annotated 3'UTR sequences from a user-provided
musculoskeletal gene set.

The script identifies input columns by header name. No miRNA IDs,
gene symbols, source names, interaction counts, or results are
hardcoded.

Inputs
------
--mirna
    Final miRNA seed-summary TSV. Required columns:
      source
      original_id
      renamed_miRNA
      seed
      log2FoldChange
      padj
      direction

--utr
    Musculoskeletal 3'UTR TSV. Required columns:
      gene_symbol
      utr_3

--categories
    Gene-set TSV. Required columns:
      gene_symbol
      category

--out
    Output TSV containing one row per miRNA-gene pair with at least
    one exact reverse-complement canonical seed match.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


OUTPUT_FIELDS = [
    "gene_symbol",
    "category",
    "renamed_miRNA",
    "original_id",
    "source",
    "seed",
    "seed_rc",
    "n_sites",
    "site_positions",
    "log2FoldChange",
    "padj",
    "direction",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Scan musculoskeletal gene 3'UTRs for exact reverse-complement "
            "matches to final nf-Sarcopipe canonical miRNA seeds."
        )
    )
    parser.add_argument(
        "--mirna",
        required=True,
        help="Final miRNA seed-summary TSV.",
    )
    parser.add_argument(
        "--utr",
        required=True,
        help="Musculoskeletal 3'UTR TSV.",
    )
    parser.add_argument(
        "--categories",
        required=True,
        help="Gene-category TSV.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV containing seed-matched pairs.",
    )
    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def normalize_rna(sequence: str) -> str:
    return (
        sequence.strip()
        .upper()
        .replace("T", "U")
        .replace(" ", "")
        .replace("\r", "")
        .replace("\n", "")
    )


def reverse_complement_rna(sequence: str) -> str:
    sequence = normalize_rna(sequence)

    complement = str.maketrans(
        {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "N": "N",
        }
    )

    return sequence.translate(complement)[::-1]


def read_categories(path: Path) -> dict[str, str]:
    categories: dict[str, str] = {}

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"gene_symbol", "category"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: category table missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            gene = row["gene_symbol"].strip().upper()
            category = row["category"].strip()

            if gene:
                categories[gene] = category

    if not categories:
        raise SystemExit("ERROR: no gene categories were loaded")

    return categories


def read_utrs(path: Path) -> dict[str, str]:
    utrs: dict[str, str] = {}

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"gene_symbol", "utr_3"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: 3'UTR table missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            gene = row["gene_symbol"].strip().upper()
            sequence = normalize_rna(row["utr_3"])

            if not gene or not sequence:
                continue

            if gene not in utrs or len(sequence) > len(utrs[gene]):
                utrs[gene] = sequence

    if not utrs:
        raise SystemExit("ERROR: no valid 3'UTR sequences were loaded")

    return utrs


def read_mirnas(path: Path) -> list[dict[str, str]]:
    mirnas: list[dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "source",
            "original_id",
            "renamed_miRNA",
            "seed",
            "log2FoldChange",
            "padj",
            "direction",
        }
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: miRNA table missing columns: "
                + ", ".join(sorted(missing))
            )

        seen: set[tuple[str, str, str]] = set()

        for row in reader:
            source = row["source"].strip()
            original_id = row["original_id"].strip()
            renamed_mirna = row["renamed_miRNA"].strip()
            seed = normalize_rna(row["seed"])

            if not source or not renamed_mirna or not seed:
                continue

            if len(seed) != 7:
                raise SystemExit(
                    "ERROR: canonical seed must contain exactly 7 nt: "
                    f"{renamed_mirna} -> {seed}"
                )

            key = (source, original_id, renamed_mirna)

            if key in seen:
                continue

            seen.add(key)

            mirnas.append(
                {
                    "source": source,
                    "original_id": original_id,
                    "renamed_miRNA": renamed_mirna,
                    "seed": seed,
                    "seed_rc": reverse_complement_rna(seed),
                    "log2FoldChange": row["log2FoldChange"].strip(),
                    "padj": row["padj"].strip(),
                    "direction": row["direction"].strip(),
                }
            )

    if not mirnas:
        raise SystemExit("ERROR: no valid miRNAs were loaded")

    return mirnas


def find_all_positions(sequence: str, motif: str) -> list[int]:
    positions: list[int] = []
    start = 0

    while True:
        index = sequence.find(motif, start)

        if index == -1:
            break

        positions.append(index + 1)
        start = index + 1

    return positions


def numeric_or_inf(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def main() -> None:
    args = parse_args()

    mirna_path = Path(args.mirna)
    utr_path = Path(args.utr)
    category_path = Path(args.categories)
    out_path = Path(args.out)

    validate_file(mirna_path, "final miRNA seed-summary table")
    validate_file(utr_path, "musculoskeletal 3'UTR table")
    validate_file(category_path, "gene-category table")

    categories = read_categories(category_path)
    utrs = read_utrs(utr_path)
    mirnas = read_mirnas(mirna_path)

    genes_to_scan = sorted(
        set(categories).intersection(utrs)
    )

    if not genes_to_scan:
        raise SystemExit(
            "ERROR: no gene symbols were shared between category and 3'UTR tables"
        )

    rows: list[dict[str, str | int]] = []
    seen_pairs: set[tuple[str, str, str]] = set()

    for gene in genes_to_scan:
        sequence = utrs[gene]

        for mirna in mirnas:
            positions = find_all_positions(
                sequence,
                mirna["seed_rc"],
            )

            if not positions:
                continue

            pair_key = (
                gene,
                mirna["source"],
                mirna["renamed_miRNA"],
            )

            if pair_key in seen_pairs:
                continue

            seen_pairs.add(pair_key)

            rows.append(
                {
                    "gene_symbol": gene,
                    "category": categories[gene],
                    "renamed_miRNA": mirna["renamed_miRNA"],
                    "original_id": mirna["original_id"],
                    "source": mirna["source"],
                    "seed": mirna["seed"],
                    "seed_rc": mirna["seed_rc"],
                    "n_sites": len(positions),
                    "site_positions": ",".join(map(str, positions)),
                    "log2FoldChange": mirna["log2FoldChange"],
                    "padj": mirna["padj"],
                    "direction": mirna["direction"],
                }
            )

    rows.sort(
        key=lambda row: (
            str(row["source"]),
            numeric_or_inf(str(row["padj"])),
            -abs(
                numeric_or_inf(
                    str(row["log2FoldChange"])
                )
            ),
            str(row["renamed_miRNA"]),
            str(row["gene_symbol"]),
        )
    )

    out_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with out_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    unique_genes = {
        str(row["gene_symbol"])
        for row in rows
    }

    unique_mirnas = {
        (
            str(row["source"]),
            str(row["renamed_miRNA"]),
        )
        for row in rows
    }

    source_counts: dict[str, int] = {}

    for row in rows:
        source = str(row["source"])
        source_counts[source] = source_counts.get(source, 0) + 1

    print("===== MUSCULOSKELETAL SEED MATCHING =====")
    print("miRNAs loaded:", len(mirnas))
    print("genes scanned:", len(genes_to_scan))
    print("miRNA-gene pairs found:", len(rows))
    print("genes with at least one match:", len(unique_genes))
    print("miRNAs with at least one match:", len(unique_mirnas))

    for source in sorted(source_counts):
        print(f"{source} pairs:", source_counts[source])

    print("Written:", out_path)


if __name__ == "__main__":
    main()
