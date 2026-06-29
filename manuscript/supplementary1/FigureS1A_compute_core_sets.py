#!/usr/bin/env python3

"""
Figure S1A – Compute CD-HIT cluster core sets

Description
-----------
This script generates per-cluster sample presence tables from a
CD-HIT .clstr file and a FASTA containing the corresponding sequence
identifiers.

It summarizes how many distinct samples are represented in each
CD-HIT cluster and exports:
1. a detailed per-sequence table
2. a summary table counting how many clusters appear in exactly N samples

Inputs
------
--status
    Run label, for example athlete, sedentary, or all.

--identity
    Identity label, for example 0.85.

--clstr
    CD-HIT .clstr file.

--fasta
    FASTA containing the same sequence IDs referenced in the .clstr file.

--out_csv
    Output detailed CSV file.

--out_summary
    Output summary TSV file.

Outputs
-------
Detailed CSV with columns:
status, identity, cluster_id, sequence_id, representative, sample, count, sequence

Summary TSV with columns:
status, identity, count_in_n_samples, n_clusters

Usage
-----
python3 FigureS1A_compute_core_sets.py \\
  --status all \\
  --identity 0.95 \\
  --clstr all.candidates_clustered.fasta.clstr \\
  --fasta all.candidates_clustered.fasta \\
  --out_csv cluster_presence.csv \\
  --out_summary cluster_presence_summary.tsv
"""

import argparse
import csv
import re
from collections import defaultdict
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


def read_fasta_as_dict(fasta_path: Path) -> dict:
    """Read a FASTA file into a dictionary: sequence_id -> sequence."""
    sequences = {}
    current_id = None
    current_sequence = []

    with fasta_path.open("r") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if not line:
                continue

            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_sequence)

                current_id = line[1:].split()[0]
                current_sequence = []

            else:
                current_sequence.append(line.strip())

        if current_id is not None:
            sequences[current_id] = "".join(current_sequence)

    return sequences


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
                sequence_id = match.group(1).strip()
                is_rep = line.endswith("*")
                entries.append((cluster_id, sequence_id, is_rep))

    return entries


def split_sample(sequence_id: str) -> str:
    """Extract the sample prefix from a sequence ID formatted as SAMPLE|SEQID."""
    if "|" in sequence_id:
        return sequence_id.split("|", 1)[0]

    return "NA"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate per-cluster sample presence from a CD-HIT .clstr file "
            "and a FASTA with matching sequence identifiers."
        )
    )

    parser.add_argument(
        "--status",
        required=True,
        help="Run label, for example athlete, sedentary, or all.",
    )
    parser.add_argument(
        "--identity",
        required=True,
        help="Identity label, for example 0.85.",
    )
    parser.add_argument(
        "--clstr",
        required=True,
        help="Path to CD-HIT .clstr file.",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="FASTA containing the same sequence IDs referenced in the .clstr file.",
    )
    parser.add_argument(
        "--out_csv",
        required=True,
        help="Output detailed CSV path.",
    )
    parser.add_argument(
        "--out_summary",
        required=True,
        help="Output summary TSV path.",
    )

    return parser.parse_args()


def main() -> None:
    """Generate CD-HIT cluster presence tables."""
    args = parse_args()

    clstr_path = validate_file(args.clstr, "CD-HIT .clstr file")
    fasta_path = validate_file(args.fasta, "Input FASTA")
    out_csv = prepare_output(args.out_csv)
    out_summary = prepare_output(args.out_summary)

    fasta_sequences = read_fasta_as_dict(fasta_path)
    entries = parse_clstr(clstr_path)

    # -------------------------
    # COMPUTE SAMPLE PRESENCE
    # -------------------------
    cluster_samples = defaultdict(set)

    for cluster_id, sequence_id, _is_rep in entries:
        cluster_samples[cluster_id].add(split_sample(sequence_id))

    cluster_count = {
        cluster_id: len(sample_set)
        for cluster_id, sample_set in cluster_samples.items()
    }

    # -------------------------
    # WRITE DETAILED TABLE
    # -------------------------
    with out_csv.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle)
        writer.writerow([
            "status",
            "identity",
            "cluster_id",
            "sequence_id",
            "representative",
            "sample",
            "count",
            "sequence",
        ])

        for cluster_id, sequence_id, is_rep in entries:
            sample = split_sample(sequence_id)
            sequence = fasta_sequences.get(sequence_id, "")

            writer.writerow([
                args.status,
                args.identity,
                cluster_id,
                sequence_id,
                int(is_rep),
                sample,
                cluster_count.get(cluster_id, 0),
                sequence,
            ])

    # -------------------------
    # WRITE SUMMARY TABLE
    # -------------------------
    frequency = defaultdict(int)

    for _cluster_id, count in cluster_count.items():
        frequency[count] += 1

    with out_summary.open("w", newline="") as out_handle:
        out_handle.write("status\tidentity\tcount_in_n_samples\tn_clusters\n")

        for count in sorted(frequency.keys()):
            out_handle.write(
                f"{args.status}\t{args.identity}\t{count}\t{frequency[count]}\n"
            )


if __name__ == "__main__":
    main()
    