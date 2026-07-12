#!/usr/bin/env python3
"""
mirdeep2_unique_mature_sequences.py

Purpose
-------
Count unique mature sequences detected by miRDeep2 across all samples,
separately for known and novel predictions.

IMPORTANT
---------
This script must be run on the published raw per-sample miRDeep2 tables:

  results/miRDeep2_known/known_miRNAs_<sample>.tsv
  results/miRDeep2_novel/novel_miRNAs_<sample>.tsv

Do NOT point it to the Nextflow work directory. The work directory may contain
simplified intermediate tables with only "miRNA" and "count" columns, which do
not contain mature sequences and therefore cannot be used for sequence-level
deduplication.

No downstream filters are applied:
- no p >= 0.8 or p >= 0.3 presence filtering
- no core-set filtering
- no DESeq2 filtering
- no overlap filtering

Sequence normalization
----------------------
Sequences are converted to uppercase RNA alphabet (T -> U). Empty values,
"-", "NA", "N/A", and "none" are ignored.

Outputs
-------
1. <prefix>.summary.tsv
   category, samples_used, raw_rows_with_sequence, unique_mature_sequences

2. <prefix>.known_unique_sequences.tsv
   sequence, n_samples, sample_ids, n_raw_detections

3. <prefix>.novel_unique_sequences.tsv
   sequence, n_samples, sample_ids, n_raw_detections
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple


MISSING = {"", "-", "na", "n/a", "none", "null"}

SEQUENCE_COLUMN_ALIASES = {
    "consensus_mature_sequence",
    "mature_sequence",
    "mature_seq",
    "consensus_mature_seq",
}


def normalize_header(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def normalize_sequence(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None

    cleaned = re.sub(r"\s+", "", value).upper().replace("T", "U")

    if cleaned.lower() in MISSING:
        return None

    if not cleaned:
        return None

    # Mature miRNA sequences should contain nucleotide symbols only.
    if not re.fullmatch(r"[ACGUN]+", cleaned):
        return None

    return cleaned


def infer_sample_id(path: Path, category: str) -> str:
    prefix = "known_miRNAs_" if category == "known" else "novel_miRNAs_"
    name = path.name

    if name.startswith(prefix) and name.endswith(".tsv"):
        return name[len(prefix):-4]

    return path.stem


def split_line(line: str) -> List[str]:
    """Prefer tabs; otherwise tolerate tables separated by 2+ spaces."""
    line = line.rstrip("\n")

    if "\t" in line:
        return line.split("\t")

    return re.split(r"\s{2,}", line.strip())


def find_sequence_column(header: List[str]) -> Optional[int]:
    normalized = [normalize_header(column) for column in header]

    for alias in SEQUENCE_COLUMN_ALIASES:
        if alias in normalized:
            return normalized.index(alias)

    return None


def parse_raw_table(path: Path) -> Tuple[List[str], int]:
    """
    Return normalized mature sequences and total non-empty data-row count.

    Fail explicitly when the table is a simplified miRNA/count table or when no
    mature-sequence column exists.
    """
    with path.open(encoding="utf-8", errors="replace") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]

    if not lines:
        raise ValueError("empty table")

    header = split_line(lines[0])
    normalized_header = [normalize_header(column) for column in header]

    if normalized_header == ["mirna", "count"]:
        raise ValueError(
            'simplified intermediate table with only "miRNA" and "count" columns'
        )

    sequence_column = find_sequence_column(header)

    if sequence_column is None:
        raise ValueError(
            "mature-sequence column not found; columns were: "
            + ", ".join(header)
        )

    sequences: List[str] = []

    for line in lines[1:]:
        cells = split_line(line)

        if sequence_column >= len(cells):
            continue

        sequence = normalize_sequence(cells[sequence_column])

        if sequence is not None:
            sequences.append(sequence)

    return sequences, max(len(lines) - 1, 0)


def collect_tables(directory: Path, category: str) -> List[Path]:
    pattern = (
        "known_miRNAs_*.tsv"
        if category == "known"
        else "novel_miRNAs_*.tsv"
    )

    if not directory.is_dir():
        raise SystemExit(f"ERROR: directory does not exist: {directory}")

    tables = sorted(path for path in directory.glob(pattern) if path.stat().st_size > 0)

    if not tables:
        raise SystemExit(
            f"ERROR: no non-empty {pattern} files found in {directory}"
        )

    return tables


def summarize_category(
    directory: Path,
    category: str,
) -> Tuple[Dict[str, Dict[str, object]], int, int]:
    """
    Return:
      sequence -> {samples: set, detections: int}
      samples successfully used
      raw rows containing a valid sequence
    """
    sequence_stats: Dict[str, Dict[str, object]] = defaultdict(
        lambda: {"samples": set(), "detections": 0}
    )

    samples_used = 0
    rows_with_sequence = 0
    errors: List[str] = []

    for table in collect_tables(directory, category):
        sample_id = infer_sample_id(table, category)

        try:
            sequences, _ = parse_raw_table(table)
        except ValueError as exc:
            errors.append(f"{table}: {exc}")
            continue

        samples_used += 1
        rows_with_sequence += len(sequences)

        for sequence in sequences:
            sample_set = sequence_stats[sequence]["samples"]
            assert isinstance(sample_set, set)
            sample_set.add(sample_id)
            sequence_stats[sequence]["detections"] = (
                int(sequence_stats[sequence]["detections"]) + 1
            )

    if errors:
        print(
            f"WARNING: {len(errors)} {category} table(s) were skipped:",
            file=sys.stderr,
        )
        for error in errors:
            print(f"  - {error}", file=sys.stderr)

    if samples_used == 0:
        raise SystemExit(
            f"ERROR: zero valid raw {category} tables were parsed. "
            "Use results/miRDeep2_known and results/miRDeep2_novel, not work/."
        )

    return dict(sequence_stats), samples_used, rows_with_sequence


def write_unique_table(
    output_path: Path,
    sequence_stats: Dict[str, Dict[str, object]],
) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(
            ["sequence", "n_samples", "sample_ids", "n_raw_detections"]
        )

        for sequence in sorted(sequence_stats):
            samples = sorted(sequence_stats[sequence]["samples"])
            detections = int(sequence_stats[sequence]["detections"])

            writer.writerow(
                [
                    sequence,
                    len(samples),
                    ",".join(samples),
                    detections,
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Count unique known and novel mature miRNA sequences from raw "
            "per-sample miRDeep2 tables."
        )
    )

    parser.add_argument(
        "--known_dir",
        required=True,
        help="Directory containing known_miRNAs_<sample>.tsv raw tables.",
    )
    parser.add_argument(
        "--novel_dir",
        required=True,
        help="Directory containing novel_miRNAs_<sample>.tsv raw tables.",
    )
    parser.add_argument(
        "--out_prefix",
        default="mirdeep2_unique_mature_sequences",
        help="Output prefix, optionally including a directory.",
    )

    args = parser.parse_args()

    known_dir = Path(args.known_dir)
    novel_dir = Path(args.novel_dir)
    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    known_stats, known_samples, known_rows = summarize_category(
        known_dir, "known"
    )
    novel_stats, novel_samples, novel_rows = summarize_category(
        novel_dir, "novel"
    )

    known_output = Path(f"{prefix}.known_unique_sequences.tsv")
    novel_output = Path(f"{prefix}.novel_unique_sequences.tsv")
    summary_output = Path(f"{prefix}.summary.tsv")

    write_unique_table(known_output, known_stats)
    write_unique_table(novel_output, novel_stats)

    with summary_output.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(
            [
                "category",
                "samples_used",
                "raw_rows_with_sequence",
                "unique_mature_sequences",
            ]
        )
        writer.writerow(
            ["known", known_samples, known_rows, len(known_stats)]
        )
        writer.writerow(
            ["novel", novel_samples, novel_rows, len(novel_stats)]
        )

    print(f"Known samples used: {known_samples}")
    print(f"Known raw rows with sequence: {known_rows}")
    print(f"Known unique mature sequences: {len(known_stats)}")
    print(f"Novel samples used: {novel_samples}")
    print(f"Novel raw rows with sequence: {novel_rows}")
    print(f"Novel unique mature sequences: {len(novel_stats)}")
    print(f"[OK] Wrote: {summary_output}")
    print(f"[OK] Wrote: {known_output}")
    print(f"[OK] Wrote: {novel_output}")


if __name__ == "__main__":
    main()
