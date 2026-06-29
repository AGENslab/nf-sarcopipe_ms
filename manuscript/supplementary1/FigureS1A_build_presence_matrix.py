#!/usr/bin/env python3

"""
Figure S1A – Build BrumiR cluster presence matrix

Description
-----------
This script builds a cluster-level BrumiR count matrix from:
1. CD-HIT cluster memberships (.clstr)
2. Per-sample BrumiR-RF high-confidence FASTA headers containing KC counts

It reconstructs cluster abundance after CD-HIT clustering by summing
KC values from BrumiR-RF FASTA headers for all members belonging to
each cluster.

It also generates:
- a counts matrix
- a binary presence/absence matrix
- group-specific core sets based on a minimum within-group presence fraction

Inputs
------
--csv
    Sample metadata CSV with columns: ID, status.

--identity
    Identity label for provenance, for example 0.85.

--clstr
    CD-HIT .clstr file.

--clustered_fasta
    CD-HIT clustered FASTA, kept for provenance.

--rf_fastas
    Per-sample BrumiR-RF high-confidence FASTA files.

--min_presence_frac
    Minimum within-group presence fraction for core sets.

--out_counts
    Output cluster count matrix TSV.

--out_binary
    Output binary presence/absence matrix TSV.

--out_core_sets
    Output core sets TSV.

Outputs
-------
counts TSV
binary TSV
core sets TSV

Usage
-----
python3 FigureS1A_build_presence_matrix.py \\
  --csv files_to_process.csv \\
  --identity 0.95 \\
  --clstr all.candidates_clustered.fasta.clstr \\
  --clustered_fasta all.candidates_clustered.fasta \\
  --rf_fastas sample1.brumir_rf.high_confidence.fasta sample2.brumir_rf.high_confidence.fasta \\
  --min_presence_frac 0.8 \\
  --out_counts brumir_counts.tsv \\
  --out_binary brumir_binary.tsv \\
  --out_core_sets brumir_core_sets.tsv
"""

import argparse
import csv
import math
import re
from collections import OrderedDict, defaultdict
from pathlib import Path


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def prepare_output(path: str) -> Path:
    """Create parent directory for an output file if needed."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def read_samples_csv(csv_path: Path) -> tuple:
    """Read sample metadata CSV."""
    statuses = {}
    groups = defaultdict(list)

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)

        required_cols = {"ID", "status"}
        missing_cols = required_cols - set(reader.fieldnames or [])

        if missing_cols:
            raise SystemExit(
                "ERROR: Sample metadata CSV missing required columns: "
                f"{', '.join(sorted(missing_cols))}"
            )

        for row in reader:
            sample_id = row["ID"]
            status = row["status"]
            statuses[sample_id] = status
            groups[status].append(sample_id)

    return statuses, groups


def parse_rf_fasta_kc(fasta_path: Path, sample_id: str) -> dict:
    """Parse one BrumiR-RF high-confidence FASTA and extract KC counts."""
    kc_map = {}
    kc_re = re.compile(r"\bKC=(\d+)\b")
    seqid_re = re.compile(r"^>(\S+)")

    with fasta_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                match_id = seqid_re.match(line)

                if not match_id:
                    continue

                seqid = match_id.group(1)
                match_kc = kc_re.search(line)

                if not match_kc:
                    kc = 0
                else:
                    kc = int(match_kc.group(1))

                key = f"{sample_id}|{seqid}"
                kc_map[key] = kc

    return kc_map


def parse_rf_fastas_kc(rf_fastas: list[Path]) -> tuple:
    """Parse multiple BrumiR-RF FASTAs and collect KC values."""
    kc_map = {}
    sample_ids_found = []

    for fasta_path in rf_fastas:
        sample_id = fasta_path.name.split(".", 1)[0]
        sample_ids_found.append(sample_id)
        kc_map.update(parse_rf_fasta_kc(fasta_path, sample_id))

    return kc_map, sample_ids_found


def parse_clstr(clstr_path: Path) -> list:
    """Parse a CD-HIT .clstr file."""
    entries = []
    cluster_id = None
    member_pattern = re.compile(r">(.+?)\.\.\.")

    with clstr_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[-1])
                continue

            match = member_pattern.search(line)

            if match and cluster_id is not None:
                member_id = match.group(1).strip()
                is_rep = line.endswith("*")
                entries.append((cluster_id, member_id, is_rep))

    return entries


def write_matrix_tsv(out_path: Path, header_samples: list, feature_to_counts: dict) -> None:
    """Write a matrix-style TSV with features as rows and samples as columns."""
    with out_path.open("w", newline="") as out_handle:
        out_handle.write("miRNA\t" + "\t".join(header_samples) + "\n")

        for feature, counts in feature_to_counts.items():
            row = [feature] + [str(counts.get(sample, 0)) for sample in header_samples]
            out_handle.write("\t".join(row) + "\n")


def write_core_sets(out_path: Path, athlete_core: list, sedentary_core: list) -> None:
    """Write group-specific core feature sets as a two-column TSV."""
    with out_path.open("w", newline="") as out_handle:
        out_handle.write("set\tmiRNA\n")

        for feature in athlete_core:
            out_handle.write(f"athlete_core\t{feature}\n")

        for feature in sedentary_core:
            out_handle.write(f"sedentary_core\t{feature}\n")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build a BrumiR cluster-level count matrix using KC values from "
            "BrumiR-RF FASTA headers and CD-HIT cluster assignments."
        )
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Sample metadata CSV with columns ID and status.",
    )
    parser.add_argument(
        "--identity",
        required=True,
        help="Identity label for provenance, for example 0.85.",
    )
    parser.add_argument(
        "--clstr",
        required=True,
        help="CD-HIT .clstr file.",
    )
    parser.add_argument(
        "--clustered_fasta",
        required=True,
        help="CD-HIT clustered FASTA, kept for provenance.",
    )
    parser.add_argument(
        "--rf_fastas",
        required=True,
        nargs="+",
        help="Per-sample BrumiR-RF high-confidence FASTA files.",
    )
    parser.add_argument(
        "--min_presence_frac",
        type=float,
        default=0.8,
        help="Minimum within-group presence fraction for core sets.",
    )
    parser.add_argument(
        "--out_counts",
        required=True,
        help="Output counts TSV.",
    )
    parser.add_argument(
        "--out_binary",
        required=True,
        help="Output binary TSV.",
    )
    parser.add_argument(
        "--out_core_sets",
        required=True,
        help="Output core sets TSV.",
    )

    return parser.parse_args()


def main() -> None:
    """Build BrumiR count, binary, and core-set matrices."""
    args = parse_args()

    csv_path = validate_file(args.csv, "Sample metadata CSV")
    clstr_path = validate_file(args.clstr, "CD-HIT .clstr file")
    validate_file(args.clustered_fasta, "CD-HIT clustered FASTA")

    rf_fastas = [
        validate_file(path, "BrumiR-RF FASTA")
        for path in args.rf_fastas
    ]

    out_counts = prepare_output(args.out_counts)
    out_binary = prepare_output(args.out_binary)
    out_core_sets = prepare_output(args.out_core_sets)

    statuses, groups = read_samples_csv(csv_path)
    kc_map, rf_samples = parse_rf_fastas_kc(rf_fastas)

    rf_sample_set = set(rf_samples)
    all_samples_ordered = [
        sample_id for sample_id in statuses.keys()
        if sample_id in rf_sample_set
    ]

    athlete_samples = [
        sample_id for sample_id in groups.get("athlete", [])
        if sample_id in rf_sample_set
    ]

    sedentary_samples = [
        sample_id for sample_id in groups.get("sedentary", [])
        if sample_id in rf_sample_set
    ]

    if not all_samples_ordered:
        raise SystemExit(
            "ERROR: No sample IDs from CSV were found among --rf_fastas filenames."
        )

    entries = parse_clstr(clstr_path)

    cluster_to_counts = defaultdict(lambda: defaultdict(int))
    cluster_to_binary = defaultdict(lambda: defaultdict(int))

    for cluster_id, member_id, _is_rep in entries:
        if "|" not in member_id:
            continue

        sample_id = member_id.split("|", 1)[0]

        if sample_id not in rf_sample_set:
            continue

        kc = kc_map.get(member_id, 0)
        feature = f"cluster_{cluster_id}"

        if kc > 0:
            cluster_to_counts[feature][sample_id] += kc
            cluster_to_binary[feature][sample_id] = 1

    features_sorted = sorted(
        cluster_to_counts.keys(),
        key=lambda feature: int(feature.split("_")[1]),
    )

    counts_ordered = OrderedDict(
        (feature, cluster_to_counts[feature])
        for feature in features_sorted
    )

    binary_ordered = OrderedDict(
        (feature, cluster_to_binary[feature])
        for feature in features_sorted
    )

    write_matrix_tsv(out_counts, all_samples_ordered, counts_ordered)
    write_matrix_tsv(out_binary, all_samples_ordered, binary_ordered)

    def core_features(sample_list: list) -> set:
        """Return features present in at least ceil(min_presence_frac * n) samples."""
        if not sample_list:
            return set()

        n_samples = len(sample_list)
        min_n = int(math.ceil(args.min_presence_frac * n_samples))
        core = set()

        for feature, presence in binary_ordered.items():
            present_n = sum(
                1 for sample in sample_list
                if presence.get(sample, 0) == 1
            )

            if present_n >= min_n:
                core.add(feature)

        return core

    athlete_core = sorted(
        core_features(athlete_samples),
        key=lambda feature: int(feature.split("_")[1]),
    )

    sedentary_core = sorted(
        core_features(sedentary_samples),
        key=lambda feature: int(feature.split("_")[1]),
    )

    write_core_sets(out_core_sets, athlete_core, sedentary_core)


if __name__ == "__main__":
    main()
