#!/usr/bin/env python3
"""
FigureS1_prepare_brumir_identity_summary.py

Prepare the BrumiR-RF CD-HIT identity sensitivity summary for
Supplementary Figure 1.

For each identity threshold, combines:
1. CD-HIT cluster sample-presence distribution.
2. BrumiR core-set membership.

Reported values:
- total clusters
- singleton clusters and percentage
- ubiquitous clusters and percentage
- weighted median sample presence
- Active-only core
- Shared core
- Sedentary-only core
- Active core total
- Sedentary core total
- Core union

No manuscript values are hard-coded.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Set, Tuple


def read_cluster_summary(
    path: Path,
    total_samples: int,
) -> Dict[str, int | float]:
    if not path.is_file():
        raise SystemExit(f"[ERROR] File not found: {path}")

    distribution: Dict[int, int] = {}

    with path.open(encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] Empty table: {path}")

        required = {"count_in_n_samples", "n_clusters"}
        if not required.issubset(reader.fieldnames):
            raise SystemExit(
                f"[ERROR] Expected columns {sorted(required)} in {path}. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            n_samples = int(row["count_in_n_samples"])
            n_clusters = int(row["n_clusters"])

            distribution[n_samples] = (
                distribution.get(n_samples, 0) + n_clusters
            )

    total_clusters = sum(distribution.values())
    singletons = distribution.get(1, 0)
    ubiquitous = distribution.get(total_samples, 0)

    midpoint = total_clusters / 2.0
    cumulative = 0
    weighted_median = 0

    for n_samples in sorted(distribution):
        cumulative += distribution[n_samples]

        if cumulative >= midpoint:
            weighted_median = n_samples
            break

    singleton_pct = (
        100.0 * singletons / total_clusters
        if total_clusters
        else 0.0
    )

    ubiquitous_pct = (
        100.0 * ubiquitous / total_clusters
        if total_clusters
        else 0.0
    )

    return {
        "total_clusters": total_clusters,
        "singletons_n": singletons,
        "singletons_pct": singleton_pct,
        "ubiquitous_n": ubiquitous,
        "ubiquitous_pct": ubiquitous_pct,
        "weighted_median_presence": weighted_median,
    }


def read_core_sets(path: Path) -> Tuple[Set[str], Set[str]]:
    if not path.is_file():
        raise SystemExit(f"[ERROR] File not found: {path}")

    active: Set[str] = set()
    sedentary: Set[str] = set()

    with path.open(encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] Empty table: {path}")

        required = {"set", "miRNA"}
        if not required.issubset(reader.fieldnames):
            raise SystemExit(
                f"[ERROR] Expected columns {sorted(required)} in {path}. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            set_name = row["set"].strip()
            candidate = row["miRNA"].strip()

            if not candidate:
                continue

            if set_name == "athlete_core":
                active.add(candidate)
            elif set_name == "sedentary_core":
                sedentary.add(candidate)

    return active, sedentary


def summarize_core_sets(path: Path) -> Dict[str, int]:
    active, sedentary = read_core_sets(path)

    shared = active & sedentary
    active_only = active - sedentary
    sedentary_only = sedentary - active
    union = active | sedentary

    return {
        "active_only": len(active_only),
        "shared_core": len(shared),
        "sedentary_only": len(sedentary_only),
        "active_core_total": len(active),
        "sedentary_core_total": len(sedentary),
        "core_union": len(union),
    }


def parse_row_spec(spec: str) -> Tuple[str, Path, Path]:
    fields = spec.split(",", 2)

    if len(fields) != 3:
        raise SystemExit(
            "[ERROR] Each --row must have the form: "
            "identity,summary_tsv,core_sets_tsv"
        )

    return (
        fields[0].strip(),
        Path(fields[1].strip()),
        Path(fields[2].strip()),
    )


def main() -> None:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--row",
        action="append",
        required=True,
        help="identity,summary_tsv,core_sets_tsv",
    )
    parser.add_argument(
        "--total_samples",
        type=int,
        default=25,
    )
    parser.add_argument("--out", required=True)

    args = parser.parse_args()

    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "cdhit_identity",
        "total_clusters",
        "singletons_n",
        "singletons_pct",
        "ubiquitous_n",
        "ubiquitous_pct",
        "weighted_median_presence",
        "active_only",
        "shared_core",
        "sedentary_only",
        "active_core_total",
        "sedentary_core_total",
        "core_union",
    ]

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
        )
        writer.writeheader()

        for spec in args.row:
            identity, summary_path, core_path = parse_row_spec(spec)

            cluster_stats = read_cluster_summary(
                summary_path,
                args.total_samples,
            )
            core_stats = summarize_core_sets(core_path)

            row = {
                "cdhit_identity": identity,
                **cluster_stats,
                **core_stats,
            }

            row["singletons_pct"] = (
                f"{float(row['singletons_pct']):.1f}"
            )
            row["ubiquitous_pct"] = (
                f"{float(row['ubiquitous_pct']):.1f}"
            )

            writer.writerow(row)

    print(f"[OK] Wrote: {output_path}")


if __name__ == "__main__":
    main()
PY
