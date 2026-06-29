#!/usr/bin/env python3
"""
Figure S1A – Generate CD-HIT supplementary summary table

Description
-----------
Build a compact supplementary summary table from CD-HIT cluster statistics
and BrumiR core-set outputs.

For each CD-HIT identity threshold, this script combines:
1) the cluster-size distribution summary
2) the athlete/sedentary core-set membership table

It reports:
- total number of clusters
- number and percentage of singletons
- number and percentage of ubiquitous clusters
- weighted median number of samples per cluster
- athlete core size
- sedentary core size
- shared core size
- union core size

Inputs
------
1) CD-HIT summary TSV with columns:
   - status
   - identity
   - count_in_n_samples
   - n_clusters

2) Core-sets TSV with columns:
   - set
   - miRNA

Outputs
-------
A single TSV summarizing all requested CD-HIT identity thresholds.

Usage
-----
python FigureS1A_generate_summary_table.py \\
  --rows \\
    0.85,path/to/summary_085.tsv,path/to/core_sets_085.tsv \\
    0.90,path/to/summary_090.tsv,path/to/core_sets_090.tsv \\
    0.95,path/to/summary_095.tsv,path/to/core_sets_095.tsv \\
  --out path/to/FigureS1A_summary_table.tsv
"""

import argparse
import csv
import sys
from pathlib import Path


def validate_input_file(path: Path, label: str) -> None:
    """Validate that an input file exists."""
    if not path.exists():
        sys.exit(f"ERROR: {label} file not found: {path}")
    if not path.is_file():
        sys.exit(f"ERROR: {label} path is not a file: {path}")


def read_summary_tsv(path: Path):
    """
    Read a CD-HIT cluster summary TSV.

    Returns
    -------
    total : int
    singletons : int
    ubiquitous : int
    wmed : int or None
    """
    dist = {}
    total = 0

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"status", "identity", "count_in_n_samples", "n_clusters"}
        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: Missing columns in {path}: {', '.join(sorted(missing))}")

        for row in reader:
            n_samples = int(row["count_in_n_samples"])
            n_clusters = int(row["n_clusters"])
            dist[n_samples] = n_clusters
            total += n_clusters

    half = total / 2.0
    cumulative = 0
    weighted_median = None

    # Original logic preserved: sample-count support is evaluated from 1 to 25.
    for n_samples in range(1, 26):
        cumulative += dist.get(n_samples, 0)
        if weighted_median is None and cumulative >= half:
            weighted_median = n_samples
            break

    singletons = dist.get(1, 0)
    ubiquitous = dist.get(25, 0)

    return total, singletons, ubiquitous, weighted_median


def read_core_sets(path: Path):
    """
    Read a core-set TSV and compute athlete/sedentary overlap statistics.

    Returns
    -------
    athlete_n : int
    sedentary_n : int
    shared_n : int
    union_n : int
    """
    athlete = set()
    sedentary = set()

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"set", "miRNA"}
        missing = required_columns - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"ERROR: Missing columns in {path}: {', '.join(sorted(missing))}")

        for row in reader:
            if row["set"] == "athlete_core":
                athlete.add(row["miRNA"])
            elif row["set"] == "sedentary_core":
                sedentary.add(row["miRNA"])

    shared = athlete & sedentary
    union = athlete | sedentary

    return len(athlete), len(sedentary), len(shared), len(union)


def pct(value: int, denominator: int) -> float:
    """Return percentage, guarding against division by zero."""
    return 100.0 * value / denominator if denominator else 0.0


def parse_row_spec(spec: str):
    """Parse one --rows entry formatted as identity,summary_tsv,core_sets_tsv."""
    fields = [field.strip() for field in spec.split(",")]

    if len(fields) != 3:
        sys.exit(
            "ERROR: Each --rows entry must have exactly three comma-separated fields: "
            "identity,summary_tsv,core_sets_tsv. "
            f"Invalid entry: {spec}"
        )

    identity, summary_path, core_path = fields
    return identity, Path(summary_path), Path(core_path)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build Figure S1A supplementary summary table from CD-HIT cluster "
            "summaries and BrumiR core-set tables."
        )
    )
    parser.add_argument(
        "--rows",
        nargs="+",
        required=True,
        help=(
            "Rows formatted as identity,summary_tsv,core_sets_tsv. "
            "Example: 0.95,results/clustering/0.95/all.0.95.cdhit_summary.tsv,"
            "results/clustering/0.95/brumir.0.95.core_sets.tsv"
        ),
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV path.",
    )

    args = parser.parse_args()

    output_path = Path(args.out)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "cdhit_identity",
        "total_clusters",
        "singletons_n1",
        "singletons_pct",
        "ubiquitous_n25",
        "ubiquitous_pct",
        "weighted_median_samples_per_cluster",
        "athlete_core_n",
        "sedentary_core_n",
        "shared_core_n",
        "core_union_n",
    ]

    output_rows = []

    for row_spec in args.rows:
        identity, summary_path, core_path = parse_row_spec(row_spec)

        validate_input_file(summary_path, "Summary TSV")
        validate_input_file(core_path, "Core-sets TSV")

        total, singletons, ubiquitous, weighted_median = read_summary_tsv(summary_path)
        athlete_n, sedentary_n, shared_n, union_n = read_core_sets(core_path)

        output_rows.append(
            {
                "cdhit_identity": identity,
                "total_clusters": str(total),
                "singletons_n1": str(singletons),
                "singletons_pct": f"{pct(singletons, total):.1f}",
                "ubiquitous_n25": str(ubiquitous),
                "ubiquitous_pct": f"{pct(ubiquitous, total):.1f}",
                "weighted_median_samples_per_cluster": str(weighted_median),
                "athlete_core_n": str(athlete_n),
                "sedentary_core_n": str(sedentary_n),
                "shared_core_n": str(shared_n),
                "core_union_n": str(union_n),
            }
        )

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"OK - wrote {output_path}")


if __name__ == "__main__":
    main()