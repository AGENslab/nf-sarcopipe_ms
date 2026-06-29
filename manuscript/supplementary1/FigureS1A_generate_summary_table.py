#!/usr/bin/env python3
"""
make_cdhit_supp_table.py

Build a compact supplementary summary table from CD-HIT cluster statistics
and BrumiR core-set outputs.

Purpose
-------
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

Expected inputs
---------------
Summary TSV columns:
- status
- identity
- count_in_n_samples
- n_clusters

Core-sets TSV columns:
- set
- miRNA

Output
------
A single TSV summarizing all requested CD-HIT identity thresholds.
"""

import argparse
import csv
from pathlib import Path


def read_summary_tsv(path):
    """
    Read a CD-HIT cluster summary TSV.

    Expected columns:
    - status
    - identity
    - count_in_n_samples
    - n_clusters

    Returns
    -------
    total : int
        Total number of clusters
    singletons : int
        Number of clusters present in exactly 1 sample
    ubiquitous : int
        Number of clusters present in exactly 25 samples
    wmed : int or None
        Weighted median of number of samples per cluster
    """
    dist = {}
    total = 0

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n = int(row["count_in_n_samples"])
            c = int(row["n_clusters"])
            dist[n] = c
            total += c

    half = total / 2.0
    cum = 0
    wmed = None

    # Keeps original logic: assumes sample-count support is evaluated from 1 to 25
    for n in range(1, 26):
        cum += dist.get(n, 0)
        if wmed is None and cum >= half:
            wmed = n
            break

    singletons = dist.get(1, 0)
    ubiquitous = dist.get(25, 0)

    return total, singletons, ubiquitous, wmed


def read_core_sets(path):
    """
    Read a core-set TSV and compute athlete/sedentary overlap statistics.

    Expected columns:
    - set
    - miRNA

    Returns
    -------
    athlete_n : int
    sedentary_n : int
    shared_n : int
    union_n : int
    """
    athlete = set()
    sedentary = set()

    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row["set"] == "athlete_core":
                athlete.add(row["miRNA"])
            elif row["set"] == "sedentary_core":
                sedentary.add(row["miRNA"])

    inter = athlete & sedentary
    union = athlete | sedentary

    return len(athlete), len(sedentary), len(inter), len(union)


def pct(x, denom):
    """Return percentage, guarding against division by zero."""
    return 100.0 * x / denom if denom else 0.0


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a supplementary summary table from CD-HIT cluster summaries "
            "and BrumiR core-set tables."
        )
    )
    parser.add_argument(
        "--rows",
        nargs="+",
        required=True,
        help=(
            "Each row must be provided as: "
            "identity,summary_tsv,core_sets_tsv "
            "Example: "
            "0.95,results/clustering/0.95/all.0.95.cdhit_summary.tsv,"
            "results/clustering/0.95/brumir.0.95.core_sets.tsv"
        ),
    )
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)

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

    rows_out = []

    for spec in args.rows:
        ident, summary_path, core_path = [x.strip() for x in spec.split(",")]

        total, s1, u25, wmed = read_summary_tsv(summary_path)
        a, s, inter, union = read_core_sets(core_path)

        rows_out.append({
            "cdhit_identity": ident,
            "total_clusters": str(total),
            "singletons_n1": str(s1),
            "singletons_pct": f"{pct(s1, total):.1f}",
            "ubiquitous_n25": str(u25),
            "ubiquitous_pct": f"{pct(u25, total):.1f}",
            "weighted_median_samples_per_cluster": str(wmed),
            "athlete_core_n": str(a),
            "sedentary_core_n": str(s),
            "shared_core_n": str(inter),
            "core_union_n": str(union),
        })

    with open(outp, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    print(f"OK - wrote {outp}")


if __name__ == "__main__":
    main()