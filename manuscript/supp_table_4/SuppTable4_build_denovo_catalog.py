#!/usr/bin/env python3
"""
Supplementary Table 4 – Final BrumiR de novo miRNA catalog

Build a manuscript-ready table for the final BrumiR de novo miRNA
candidates by combining:

1. The final BrumiR de novo catalog containing differential-expression
   statistics.
2. The provisional miRNA ID mapping.
3. The 7-mer seed-space summary.

No candidate IDs, sequences, or expected candidate counts are hardcoded.

Output columns
--------------
- provisory_ID
- cluster
- length
- seq
- seed_2_8
- log2FC
- padj
- regulation
- seed_miRBase_match
- miRBase_family
- n_total_7mers
- n_matching_known_7mers
- matched_known_7mers
- matched_families_all_positions
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List


OUTPUT_FIELDS = [
    "provisory_ID",
    "cluster",
    "length",
    "seq",
    "seed_2_8",
    "log2FC",
    "padj",
    "regulation",
    "seed_miRBase_match",
    "miRBase_family",
    "n_total_7mers",
    "n_matching_known_7mers",
    "matched_known_7mers",
    "matched_families_all_positions",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build Supplementary Table 4 for final BrumiR de novo "
            "miRNA candidates."
        )
    )

    parser.add_argument(
        "--denovo_catalog",
        required=True,
        help="Final BrumiR de novo catalog TSV.",
    )
    parser.add_argument(
        "--mapping",
        required=True,
        help="BrumiR provisional-ID mapping TSV.",
    )
    parser.add_argument(
        "--seed_summary",
        required=True,
        help="Final de novo candidate 7-mer summary TSV.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output Supplementary Table 4 TSV.",
    )

    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")

    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if not reader.fieldnames:
            raise SystemExit(f"ERROR: TSV has no header: {path}")

        return list(reader)


def require_columns(
    rows: List[Dict[str, str]],
    required: set[str],
    label: str,
) -> None:
    if not rows:
        raise SystemExit(f"ERROR: {label} contains no rows")

    missing = required - set(rows[0])

    if missing:
        raise SystemExit(
            f"ERROR: {label} is missing columns: "
            + ", ".join(sorted(missing))
        )


def normalize_sequence(sequence: str) -> str:
    return (
        sequence.strip()
        .upper()
        .replace("T", "U")
        .replace(" ", "")
        .replace("-", "")
    )


def index_unique(
    rows: List[Dict[str, str]],
    key_column: str,
    label: str,
) -> Dict[str, Dict[str, str]]:
    indexed: Dict[str, Dict[str, str]] = {}

    for row in rows:
        key = row[key_column].strip()

        if not key:
            raise SystemExit(
                f"ERROR: empty {key_column} in {label}"
            )

        if key in indexed:
            raise SystemExit(
                f"ERROR: duplicated {key_column} '{key}' in {label}"
            )

        indexed[key] = row

    return indexed


def main() -> None:
    args = parse_args()

    catalog_path = Path(args.denovo_catalog)
    mapping_path = Path(args.mapping)
    seed_summary_path = Path(args.seed_summary)
    out_path = Path(args.out)

    validate_file(
        catalog_path,
        "final BrumiR de novo catalog",
    )
    validate_file(
        mapping_path,
        "BrumiR provisional-ID mapping",
    )
    validate_file(
        seed_summary_path,
        "candidate 7-mer summary",
    )

    catalog_rows = read_tsv(catalog_path)
    mapping_rows = read_tsv(mapping_path)
    seed_rows = read_tsv(seed_summary_path)

    require_columns(
        catalog_rows,
        {
            "candidate",
            "mature_seq",
            "log2FoldChange",
            "padj",
            "direction",
        },
        "final BrumiR de novo catalog",
    )

    require_columns(
        mapping_rows,
        {
            "provisional_miRNA_ID",
            "original_cluster",
            "mature_seq",
            "seed_2_8",
        },
        "BrumiR provisional-ID mapping",
    )

    require_columns(
        seed_rows,
        {
            "provisional_miRNA_ID",
            "original_cluster",
            "length",
            "mature_seq",
            "canonical_seed_2_8",
            "seed_miRBase_match",
            "miRBase_family",
            "n_total_7mers",
            "n_matching_known_7mers",
            "matched_known_7mers",
            "matched_families_all_positions",
        },
        "candidate 7-mer summary",
    )

    catalog_by_cluster = index_unique(
        catalog_rows,
        "candidate",
        "final BrumiR de novo catalog",
    )

    seed_by_id = index_unique(
        seed_rows,
        "provisional_miRNA_ID",
        "candidate 7-mer summary",
    )

    output_rows = []

    for mapping in mapping_rows:
        provisional_id = mapping[
            "provisional_miRNA_ID"
        ].strip()

        cluster = mapping[
            "original_cluster"
        ].strip()

        if cluster not in catalog_by_cluster:
            raise SystemExit(
                f"ERROR: cluster '{cluster}' from mapping was not "
                "found in final BrumiR de novo catalog"
            )

        if provisional_id not in seed_by_id:
            raise SystemExit(
                f"ERROR: provisional ID '{provisional_id}' was not "
                "found in candidate 7-mer summary"
            )

        catalog = catalog_by_cluster[cluster]
        seed_row = seed_by_id[provisional_id]

        mapping_seq = normalize_sequence(
            mapping["mature_seq"]
        )
        catalog_seq = normalize_sequence(
            catalog["mature_seq"]
        )
        summary_seq = normalize_sequence(
            seed_row["mature_seq"]
        )

        if not (
            mapping_seq
            == catalog_seq
            == summary_seq
        ):
            raise SystemExit(
                f"ERROR: sequence mismatch for {provisional_id} "
                f"({cluster})"
            )

        mapping_seed = normalize_sequence(
            mapping["seed_2_8"]
        )
        summary_seed = normalize_sequence(
            seed_row["canonical_seed_2_8"]
        )
        computed_seed = mapping_seq[1:8]

        if not (
            mapping_seed
            == summary_seed
            == computed_seed
        ):
            raise SystemExit(
                f"ERROR: seed mismatch for {provisional_id} "
                f"({cluster})"
            )

        if seed_row["original_cluster"].strip() != cluster:
            raise SystemExit(
                f"ERROR: cluster mismatch for {provisional_id}: "
                f"mapping={cluster}; "
                f"seed_summary={seed_row['original_cluster']}"
            )

        output_rows.append({
            "provisory_ID": provisional_id,
            "cluster": cluster,
            "length": str(len(mapping_seq)),
            "seq": mapping_seq,
            "seed_2_8": mapping_seed,
            "log2FC": catalog["log2FoldChange"],
            "padj": catalog["padj"],
            "regulation": catalog["direction"],
            "seed_miRBase_match": (
                seed_row["seed_miRBase_match"]
            ),
            "miRBase_family": (
                seed_row["miRBase_family"]
            ),
            "n_total_7mers": (
                seed_row["n_total_7mers"]
            ),
            "n_matching_known_7mers": (
                seed_row["n_matching_known_7mers"]
            ),
            "matched_known_7mers": (
                seed_row["matched_known_7mers"]
            ),
            "matched_families_all_positions": (
                seed_row[
                    "matched_families_all_positions"
                ]
            ),
        })

    if len(output_rows) != len(mapping_rows):
        raise SystemExit(
            "ERROR: output row count differs from mapping row count"
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
        writer.writerows(output_rows)

    print("===== SUPPLEMENTARY TABLE 4 =====")
    print("Final de novo candidates:", len(output_rows))
    print("Written:", out_path)


if __name__ == "__main__":
    main()
