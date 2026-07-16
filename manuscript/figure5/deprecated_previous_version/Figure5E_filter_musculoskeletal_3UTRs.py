#!/usr/bin/env python3
"""
Figure 5E – Filter musculoskeletal gene 3'UTR sequences

Extract the 3'UTR sequences of a user-provided gene set from the same
3'UTR table used by nf-Sarcopipe. This preserves consistency with the
main target-prediction workflow and avoids external BioMart queries.

Inputs
------
--genes
    TSV containing at least:
      gene_symbol
      category

--utr_table
    TSV containing:
      gene_symbol
      utr_3

--out
    Filtered TSV containing the requested genes with available 3'UTRs.

--out_missing
    TSV listing requested genes without an available 3'UTR.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


OUTPUT_FIELDS = [
    "gene_symbol",
    "category",
    "utr_3",
    "utr_length",
]

MISSING_FIELDS = [
    "gene_symbol",
    "category",
    "status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter 3'UTR sequences for a user-provided musculoskeletal "
            "gene set using the nf-Sarcopipe UTR table."
        )
    )
    parser.add_argument(
        "--genes",
        required=True,
        help="Gene-set TSV with gene_symbol and category columns.",
    )
    parser.add_argument(
        "--utr_table",
        required=True,
        help="nf-Sarcopipe UTR TSV with gene_symbol and utr_3 columns.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV with requested genes and available 3'UTRs.",
    )
    parser.add_argument(
        "--out_missing",
        required=True,
        help="Output TSV listing requested genes without 3'UTRs.",
    )
    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def normalize_sequence(sequence: str) -> str:
    return (
        sequence.strip()
        .upper()
        .replace("T", "U")
        .replace(" ", "")
        .replace("\r", "")
        .replace("\n", "")
    )


def read_gene_set(path: Path) -> list[dict[str, str]]:
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
                "ERROR: Gene-set table is missing columns: "
                + ", ".join(sorted(missing))
            )

        rows: list[dict[str, str]] = []
        seen: set[str] = set()

        for row in reader:
            gene = row["gene_symbol"].strip().upper()
            category = row["category"].strip()

            if not gene:
                continue

            if gene in seen:
                raise SystemExit(
                    f"ERROR: duplicated gene in gene-set table: {gene}"
                )

            seen.add(gene)
            rows.append({
                "gene_symbol": gene,
                "category": category,
            })

    if not rows:
        raise SystemExit("ERROR: no genes were loaded from the gene-set table")

    return rows


def read_utrs(path: Path) -> dict[str, str]:
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
                "ERROR: UTR table is missing columns: "
                + ", ".join(sorted(missing))
            )

        utrs: dict[str, str] = {}

        for row in reader:
            gene = row["gene_symbol"].strip().upper()
            sequence = normalize_sequence(row["utr_3"])

            if not gene or not sequence:
                continue

            # If duplicated genes are present, retain the longest sequence.
            if gene not in utrs or len(sequence) > len(utrs[gene]):
                utrs[gene] = sequence

    if not utrs:
        raise SystemExit("ERROR: no valid 3'UTRs were loaded")

    return utrs


def main() -> None:
    args = parse_args()

    genes_path = Path(args.genes)
    utr_path = Path(args.utr_table)
    out_path = Path(args.out)
    missing_path = Path(args.out_missing)

    validate_file(genes_path, "gene-set table")
    validate_file(utr_path, "UTR table")

    genes = read_gene_set(genes_path)
    utrs = read_utrs(utr_path)

    output_rows: list[dict[str, str | int]] = []
    missing_rows: list[dict[str, str]] = []

    for record in genes:
        gene = record["gene_symbol"]
        category = record["category"]

        if gene not in utrs:
            missing_rows.append({
                "gene_symbol": gene,
                "category": category,
                "status": "3UTR_not_available",
            })
            continue

        sequence = utrs[gene]

        output_rows.append({
            "gene_symbol": gene,
            "category": category,
            "utr_3": sequence,
            "utr_length": len(sequence),
        })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    missing_path.parent.mkdir(parents=True, exist_ok=True)

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
        writer.writerows(output_rows)

    with missing_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MISSING_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(missing_rows)

    print("===== MUSCULOSKELETAL 3'UTR FILTERING =====")
    print("Genes requested:", len(genes))
    print("Genes with available 3'UTR:", len(output_rows))
    print("Genes without available 3'UTR:", len(missing_rows))
    print("Written:", out_path)
    print("Written:", missing_path)


if __name__ == "__main__":
    main()
