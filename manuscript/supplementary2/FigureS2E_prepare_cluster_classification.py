#!/usr/bin/env python3

"""
FigureS2E_prepare_cluster_classification.py

Prepare Supplementary Figure S2E by classifying the final
BrumiR-RF candidates according to their original CD-HIT
cluster size.

Inputs
------
--final_candidates
    Figure3E final BrumiR passfilter candidate TSV.

--cluster_support
    cluster_support.tsv containing cluster and n_members.

--out
    Output TSV.

Rules
-----
- Only de novo candidates are retained.
- n_members == 1  -> Singleton
- n_members > 1   -> Cluster
- Final names are assigned to the seven current de novo
  candidates used in the manuscript.

Usage
-----
python3 FigureS2E_prepare_cluster_classification.py \
  --final_candidates Figure3E_final_BrumiR_passfilter_candidates.tsv \
  --cluster_support cluster_support.tsv \
  --out FigureS2E_cluster_classification.tsv
"""

import argparse
import csv
from pathlib import Path


FINAL_NAMES = {
    "cluster_740": "hsa-miR-novel_A",
    "cluster_986": "hsa-miR-novel_B",
    "cluster_42": "hsa-miR-novel_C",
    "cluster_374": "hsa-miR-novel_D",
    "cluster_413": "hsa-miR-novel_E",
    "cluster_465": "hsa-miR-novel_F",
    "cluster_485": "hsa-miR-novel_G",
}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--final_candidates",
        required=True,
        help="Final BrumiR passfilter candidate TSV.",
    )

    parser.add_argument(
        "--cluster_support",
        required=True,
        help="cluster_support.tsv with cluster and n_members.",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV.",
    )

    return parser.parse_args()


def read_cluster_sizes(path):
    sizes = {}

    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"cluster", "n_members"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: cluster support table is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            cluster = row["cluster"].strip()

            try:
                n_members = int(row["n_members"])
            except (TypeError, ValueError):
                continue

            sizes[cluster] = n_members

    return sizes


def main():
    args = parse_args()

    final_candidates_path = Path(args.final_candidates)
    cluster_support_path = Path(args.cluster_support)
    output_path = Path(args.out)

    for path in (final_candidates_path, cluster_support_path):
        if not path.is_file():
            raise SystemExit(f"ERROR: input file not found: {path}")

    cluster_sizes = read_cluster_sizes(cluster_support_path)

    rows = []

    with final_candidates_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "candidate",
            "group",
            "candidate_type",
            "mature_seq",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: final candidate table is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            candidate = row["candidate"].strip()

            if row["candidate_type"].strip() != "de_novo":
                continue

            if candidate not in FINAL_NAMES:
                continue

            if candidate not in cluster_sizes:
                raise SystemExit(
                    f"ERROR: cluster size not found for {candidate}"
                )

            cluster_size = cluster_sizes[candidate]

            rows.append({
                "Final candidate": FINAL_NAMES[candidate],
                "Original cluster": candidate,
                "Group": row["group"],
                "Cluster size": cluster_size,
                "Classification": (
                    "Singleton"
                    if cluster_size == 1
                    else "Cluster"
                ),
                "Mature sequence (5'->3')": row["mature_seq"],
            })

    expected = set(FINAL_NAMES)
    found = {row["Original cluster"] for row in rows}

    missing_candidates = sorted(expected - found)

    if missing_candidates:
        raise SystemExit(
            "ERROR: final de novo candidates not found: "
            + ", ".join(missing_candidates)
        )

    order = {
        candidate: index
        for index, candidate in enumerate(FINAL_NAMES)
    }

    rows.sort(
        key=lambda row: order[row["Original cluster"]]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "Final candidate",
        "Original cluster",
        "Group",
        "Cluster size",
        "Classification",
        "Mature sequence (5'->3')",
    ]

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)

    singleton_count = sum(
        row["Classification"] == "Singleton"
        for row in rows
    )

    cluster_count = sum(
        row["Classification"] == "Cluster"
        for row in rows
    )

    print("Figure S2E cluster classification")
    print(f"De novo candidates: {len(rows)}")
    print(f"Singletons: {singleton_count}")
    print(f"Clusters: {cluster_count}")
    print(f"Written: {output_path}")


if __name__ == "__main__":
    main()
