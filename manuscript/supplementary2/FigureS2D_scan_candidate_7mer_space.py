#!/usr/bin/env python3
"""
Figure S2D – Scan de novo candidate 7-mer space

Generate every possible 7-mer from the final BrumiR de novo candidate
sequences and compare each 7-mer against a dictionary of known human
miRNA canonical seeds.

The script is dataset-independent. Candidate IDs, sequences, groups,
clusters, and candidate counts are read directly from the input TSV.

Expected candidate input columns
--------------------------------
Required:
- provisional_miRNA_ID
- original_cluster
- group
- mature_seq
- seed_2_8

Accepted aliases are also supported for compatibility.

Expected seed dictionary columns
--------------------------------
- miRNA_id
- family
- canonical_seed

Outputs
-------
1. Long-format TSV:
   One row per possible candidate 7-mer.

2. Summary TSV:
   One row per final de novo candidate.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple


LONG_FIELDS = [
    "provisional_miRNA_ID",
    "original_cluster",
    "group",
    "mature_seq",
    "kmer_start_pos_1based",
    "kmer_7mer",
    "is_canonical_seed_2_8",
    "kmer_status",
    "matched_known_miRNAs",
    "matched_families",
]


SUMMARY_FIELDS = [
    "provisional_miRNA_ID",
    "original_cluster",
    "group",
    "length",
    "mature_seq",
    "canonical_seed_2_8",
    "seed_miRBase_match",
    "miRBase_family",
    "n_total_7mers",
    "n_matching_known_7mers",
    "matching_positions",
    "matched_known_7mers",
    "matched_families_all_positions",
]


COLUMN_ALIASES = {
    "provisional_miRNA_ID": [
        "provisional_miRNA_ID",
        "provisory_ID",
        "provisional_id",
        "renamed_miRNA",
        "candidate",
    ],
    "original_cluster": [
        "original_cluster",
        "cluster",
        "original_id",
    ],
    "group": [
        "group",
        "core_membership",
    ],
    "mature_seq": [
        "mature_seq",
        "mature_sequence",
        "sequence",
        "seq",
    ],
    "seed_2_8": [
        "seed_2_8",
        "canonical_seed_2_8",
        "canonical_seed",
        "seed",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate all possible 7-mers from final BrumiR de novo "
            "candidate miRNAs and compare them against known canonical seeds."
        )
    )

    parser.add_argument(
        "--input_tsv",
        required=True,
        help="Final BrumiR provisional-ID mapping or candidate TSV.",
    )
    parser.add_argument(
        "--seed_db",
        required=True,
        help="Known human miRNA seed dictionary TSV.",
    )
    parser.add_argument(
        "--out_long",
        required=True,
        help="Output long-format TSV.",
    )
    parser.add_argument(
        "--out_summary",
        required=True,
        help="Output summary TSV.",
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
        .replace("-", "")
    )


def resolve_column(
    fieldnames: Sequence[str],
    canonical_name: str,
    required: bool = True,
) -> str | None:
    aliases = COLUMN_ALIASES[canonical_name]

    for alias in aliases:
        if alias in fieldnames:
            return alias

    if required:
        raise SystemExit(
            f"ERROR: input candidate TSV lacks a column for "
            f"'{canonical_name}'. Accepted names: {', '.join(aliases)}"
        )

    return None


def load_seed_dictionary(
    seed_db: Path,
) -> Tuple[
    Dict[str, Set[str]],
    Dict[str, Set[str]],
]:
    seed_to_ids: Dict[str, Set[str]] = defaultdict(set)
    seed_to_families: Dict[str, Set[str]] = defaultdict(set)

    with seed_db.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "canonical_seed",
            "miRNA_id",
            "family",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: seed dictionary is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            seed = normalize_sequence(row["canonical_seed"])

            if len(seed) != 7:
                continue

            mirna_id = row["miRNA_id"].strip()
            family = row["family"].strip()

            if mirna_id:
                seed_to_ids[seed].add(mirna_id)

            if family:
                seed_to_families[seed].add(family)

    if not seed_to_ids:
        raise SystemExit(
            f"ERROR: no valid known 7-mer seeds were loaded from {seed_db}"
        )

    return seed_to_ids, seed_to_families


def load_candidates(path: Path) -> List[Dict[str, str]]:
    candidates: List[Dict[str, str]] = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        fieldnames = reader.fieldnames or []

        id_col = resolve_column(
            fieldnames,
            "provisional_miRNA_ID",
        )
        cluster_col = resolve_column(
            fieldnames,
            "original_cluster",
        )
        group_col = resolve_column(
            fieldnames,
            "group",
            required=False,
        )
        sequence_col = resolve_column(
            fieldnames,
            "mature_seq",
        )
        seed_col = resolve_column(
            fieldnames,
            "seed_2_8",
        )

        seen_ids: Set[str] = set()

        for line_number, row in enumerate(reader, start=2):
            provisional_id = row[id_col].strip()
            original_cluster = row[cluster_col].strip()
            group = row[group_col].strip() if group_col else ""
            mature_seq = normalize_sequence(row[sequence_col])
            seed = normalize_sequence(row[seed_col])

            if not provisional_id:
                raise SystemExit(
                    f"ERROR: empty provisional ID at line {line_number}"
                )

            if provisional_id in seen_ids:
                raise SystemExit(
                    f"ERROR: duplicate provisional ID: {provisional_id}"
                )

            if len(mature_seq) < 7:
                raise SystemExit(
                    f"ERROR: mature sequence shorter than 7 nt for "
                    f"{provisional_id}: {mature_seq}"
                )

            computed_seed = mature_seq[1:8]

            if seed != computed_seed:
                raise SystemExit(
                    f"ERROR: seed mismatch for {provisional_id}. "
                    f"Input seed={seed}; sequence positions 2-8={computed_seed}"
                )

            seen_ids.add(provisional_id)

            candidates.append({
                "provisional_miRNA_ID": provisional_id,
                "original_cluster": original_cluster,
                "group": group,
                "mature_seq": mature_seq,
                "seed_2_8": seed,
            })

    if not candidates:
        raise SystemExit(
            f"ERROR: no candidates were loaded from {path}"
        )

    return candidates


def scan_candidate(
    candidate: Dict[str, str],
    seed_to_ids: Dict[str, Set[str]],
    seed_to_families: Dict[str, Set[str]],
) -> Tuple[List[Dict[str, str]], Dict[str, str]]:
    provisional_id = candidate["provisional_miRNA_ID"]
    original_cluster = candidate["original_cluster"]
    group = candidate["group"]
    mature_seq = candidate["mature_seq"]
    canonical_seed = candidate["seed_2_8"]

    long_rows: List[Dict[str, str]] = []

    matched_positions: List[str] = []
    matched_7mers: Set[str] = set()
    matched_families_all: Set[str] = set()

    for index in range(len(mature_seq) - 6):
        kmer = mature_seq[index:index + 7]
        start_position = index + 1
        is_canonical = start_position == 2

        known_ids = sorted(seed_to_ids.get(kmer, set()))
        families = sorted(seed_to_families.get(kmer, set()))

        if known_ids:
            matched_positions.append(str(start_position))
            matched_7mers.add(kmer)
            matched_families_all.update(families)

        long_rows.append({
            "provisional_miRNA_ID": provisional_id,
            "original_cluster": original_cluster,
            "group": group,
            "mature_seq": mature_seq,
            "kmer_start_pos_1based": str(start_position),
            "kmer_7mer": kmer,
            "is_canonical_seed_2_8": (
                "yes" if is_canonical else "no"
            ),
            "kmer_status": (
                "known_seed_match"
                if known_ids
                else "novel_seed"
            ),
            "matched_known_miRNAs": (
                ",".join(known_ids)
                if known_ids
                else "-"
            ),
            "matched_families": (
                ",".join(families)
                if families
                else "-"
            ),
        })

    canonical_families = sorted(
        seed_to_families.get(
            canonical_seed,
            set(),
        )
    )

    summary_row = {
        "provisional_miRNA_ID": provisional_id,
        "original_cluster": original_cluster,
        "group": group,
        "length": str(len(mature_seq)),
        "mature_seq": mature_seq,
        "canonical_seed_2_8": canonical_seed,
        "seed_miRBase_match": (
            "yes"
            if canonical_seed in seed_to_ids
            else "no"
        ),
        "miRBase_family": (
            ",".join(canonical_families)
            if canonical_families
            else "-"
        ),
        "n_total_7mers": str(len(mature_seq) - 6),
        "n_matching_known_7mers": str(len(matched_7mers)),
        "matching_positions": (
            ",".join(matched_positions)
            if matched_positions
            else "-"
        ),
        "matched_known_7mers": (
            ",".join(sorted(matched_7mers))
            if matched_7mers
            else "-"
        ),
        "matched_families_all_positions": (
            ",".join(sorted(matched_families_all))
            if matched_families_all
            else "-"
        ),
    }

    return long_rows, summary_row


def write_tsv(
    path: Path,
    rows: Iterable[Dict[str, str]],
    fieldnames: Sequence[str],
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


def main() -> None:
    args = parse_args()

    input_path = Path(args.input_tsv)
    seed_db_path = Path(args.seed_db)
    out_long_path = Path(args.out_long)
    out_summary_path = Path(args.out_summary)

    validate_file(input_path, "candidate input")
    validate_file(seed_db_path, "known seed dictionary")

    seed_to_ids, seed_to_families = (
        load_seed_dictionary(seed_db_path)
    )

    candidates = load_candidates(input_path)

    long_rows: List[Dict[str, str]] = []
    summary_rows: List[Dict[str, str]] = []

    for candidate in candidates:
        candidate_long, candidate_summary = scan_candidate(
            candidate,
            seed_to_ids,
            seed_to_families,
        )

        long_rows.extend(candidate_long)
        summary_rows.append(candidate_summary)

    write_tsv(
        out_long_path,
        long_rows,
        LONG_FIELDS,
    )

    write_tsv(
        out_summary_path,
        summary_rows,
        SUMMARY_FIELDS,
    )

    print("===== CANDIDATE 7-MER SCAN =====")
    print("Candidates:", len(candidates))
    print("Known canonical seeds:", len(seed_to_ids))
    print("Long-format rows:", len(long_rows))
    print("Written:", out_long_path)
    print("Written:", out_summary_path)


if __name__ == "__main__":
    main()
