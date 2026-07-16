#!/usr/bin/env python3

"""
Figure 5E – Validate musculoskeletal annotated miRNA–gene pairs with multiMiR

Intersects miRDeep2 seed-matched musculoskeletal miRNA–gene pairs with
the existing nf-Sarcopipe multiMiR target table.

No pipeline files are modified.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


OUT_FIELDS = [
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
    "validation_database_list",
    "validation_database_count",
    "validated_by",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed_matches", required=True)
    parser.add_argument("--multimir", required=True)
    parser.add_argument("--out_validated", required=True)
    parser.add_argument("--out_not_validated", required=True)
    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def read_multimir(path: Path) -> dict[tuple[str, str], set[str]]:
    support: dict[tuple[str, str], set[str]] = defaultdict(set)

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "mature_mirna_id",
            "target_symbol",
            "database",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: multiMiR table missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            mirna = row["mature_mirna_id"].strip()
            gene = row["target_symbol"].strip()
            database = row["database"].strip()

            if mirna and gene:
                support[(mirna, gene)].add(database)

    return support


def read_seed_pairs(path: Path) -> list[dict[str, str]]:
    rows = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
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
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: seed-match table missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            if row["source"].strip() != "miRDeep2":
                continue

            rows.append({
                key: row.get(key, "").strip()
                for key in required
            })

    return rows


def write_tsv(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

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


def main() -> None:
    args = parse_args()

    seed_path = Path(args.seed_matches)
    multimir_path = Path(args.multimir)
    validated_path = Path(args.out_validated)
    rejected_path = Path(args.out_not_validated)

    validate_file(seed_path, "seed-match table")
    validate_file(multimir_path, "multiMiR table")

    multimir = read_multimir(multimir_path)
    seed_pairs = read_seed_pairs(seed_path)

    validated = []
    rejected = []

    for row in seed_pairs:
        key = (
            row["renamed_miRNA"],
            row["gene_symbol"],
        )

        databases = sorted(
            db for db in multimir.get(key, set())
            if db
        )

        if not databases:
            rejected.append(row)
            continue

        validated.append({
            **row,
            "validation_database_list": ";".join(databases),
            "validation_database_count": str(len(databases)),
            "validated_by": "multiMiR",
        })

    validated.sort(
        key=lambda row: (
            row["renamed_miRNA"],
            row["gene_symbol"],
        )
    )

    rejected.sort(
        key=lambda row: (
            row["renamed_miRNA"],
            row["gene_symbol"],
        )
    )

    write_tsv(
        validated_path,
        validated,
        OUT_FIELDS,
    )

    write_tsv(
        rejected_path,
        rejected,
        OUT_FIELDS[:-3],
    )

    unique_mirnas = {
        row["renamed_miRNA"]
        for row in validated
    }

    unique_genes = {
        row["gene_symbol"]
        for row in validated
    }

    print("===== FIGURE 5E MULTIMIR VALIDATION =====")
    print("Seed-predicted annotated pairs:", len(seed_pairs))
    print("Pairs validated by multiMiR:", len(validated))
    print("Pairs not validated:", len(rejected))
    print("Validated annotated miRNAs:", len(unique_mirnas))
    print("Validated target genes:", len(unique_genes))
    print("Written:", validated_path)
    print("Written:", rejected_path)


if __name__ == "__main__":
    main()
