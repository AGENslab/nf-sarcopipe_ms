#!/usr/bin/env python3

"""
Figure3C_prepare_annotation_summary.py

Prepare the data supporting Figure 3C.

The script compares three current nf-Sarcopipe feature sets
against the human MirGeneDB mature-miRNA reference:

1. BrumiR-RF core union
2. Structurally supported BrumiR-RF core candidates
3. miRDeep2 known core union

Qualifying annotation criteria:
- nucleotide identity >= 98%
- query coverage >= 75%
- mismatches == 0

A MirGeneDB annotation is classified as shared when the same
MirGeneDB reference ID is recovered by both the complete
BrumiR-RF core and the miRDeep2 core.

No counts or input paths are hardcoded.

Outputs
-------
- BrumiR core FASTA
- BrumiR supported FASTA
- miRDeep2 core FASTA
- raw and filtered BLAST tables
- shared MirGeneDB ID table
- per-query annotation table
- Figure 3C summary TSV
"""

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path


BLAST_COLUMNS = [
    "query",
    "subject",
    "identity",
    "alignment_length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "query_length",
    "subject_length",
]


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--brumir_fasta", required=True)
    parser.add_argument("--brumir_core", required=True)
    parser.add_argument("--brumir_supported_table", required=True)
    parser.add_argument("--mirdeep2_core", required=True)
    parser.add_argument("--mirbase_mature", required=True)
    parser.add_argument("--mirgenedb_mature", required=True)
    parser.add_argument("--blast_bin", required=True)
    parser.add_argument("--outdir", required=True)

    parser.add_argument("--min_identity", type=float, default=98.0)
    parser.add_argument("--min_query_coverage", type=float, default=75.0)
    parser.add_argument("--max_mismatches", type=int, default=0)

    return parser.parse_args()


def normalize_sequence(sequence):
    return sequence.upper().replace("U", "T")


def read_fasta(path):
    records = {}
    name = None
    chunks = []

    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    records[name] = normalize_sequence("".join(chunks))

                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

        if name is not None:
            records[name] = normalize_sequence("".join(chunks))

    return records


def write_fasta(records, path):
    with Path(path).open("w") as handle:
        for identifier in sorted(records, key=natural_key):
            handle.write(f">{identifier}\n{records[identifier]}\n")


def natural_key(value):
    match = re.search(r"(\d+)$", value)

    if match:
        return value[:match.start()], int(match.group(1))

    return value, -1


def read_core_ids(path):
    ids = set()

    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        feature_column = (
            "miRNA"
            if reader.fieldnames and "miRNA" in reader.fieldnames
            else "cluster"
        )

        for row in reader:
            ids.add(row[feature_column].strip())

    return ids


def read_supported_table(path):
    records = {}

    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            identifier = row["Candidate"].strip()
            sequence = normalize_sequence(
                row["Mature sequence (5'->3')"].strip()
            )
            records[identifier] = sequence

    return records


def run_command(command):
    print("+", " ".join(map(str, command)))
    subprocess.run(command, check=True)


def run_blast(query, database, output, blastn):
    run_command([
        str(blastn),
        "-task", "blastn-short",
        "-query", str(query),
        "-db", str(database),
        "-strand", "plus",
        "-dust", "no",
        "-word_size", "7",
        "-evalue", "1000",
        "-outfmt",
        (
            "6 qseqid sseqid pident length mismatch gapopen "
            "qstart qend sstart send evalue bitscore qlen slen"
        ),
        "-out", str(output),
    ])


def filter_hits(
    input_path,
    output_path,
    min_identity,
    min_query_coverage,
    max_mismatches,
):
    rows = []

    with Path(input_path).open() as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
            fieldnames=BLAST_COLUMNS,
        )

        for row in reader:
            identity = float(row["identity"])
            alignment_length = int(row["alignment_length"])
            query_length = int(row["query_length"])
            mismatches = int(row["mismatch"])

            coverage = (
                100.0 * alignment_length / query_length
                if query_length
                else 0.0
            )

            if (
                identity >= min_identity
                and coverage >= min_query_coverage
                and mismatches <= max_mismatches
            ):
                row["query_coverage"] = f"{coverage:.2f}"
                rows.append(row)

    columns = BLAST_COLUMNS + ["query_coverage"]

    with Path(output_path).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    return rows


def subjects_by_query(rows):
    result = defaultdict(set)

    for row in rows:
        result[row["query"]].add(row["subject"])

    return result


def classify_queries(all_queries, query_subjects, shared_subjects):
    classifications = {}

    for query in all_queries:
        subjects = query_subjects.get(query, set())

        if not subjects:
            classifications[query] = "Not annotated in MirGeneDB"
        elif subjects & shared_subjects:
            classifications[query] = "Shared annotated"
        else:
            classifications[query] = "Unique annotated"

    return classifications


def write_detail(path, classifications_by_set, subject_maps):
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "feature_set",
            "query",
            "category",
            "MirGeneDB_IDs",
        ])

        for set_name in classifications_by_set:
            classifications = classifications_by_set[set_name]
            subject_map = subject_maps[set_name]

            for query in sorted(classifications, key=natural_key):
                writer.writerow([
                    set_name,
                    query,
                    classifications[query],
                    ",".join(sorted(subject_map.get(query, set())))
                    or "NA",
                ])


def write_summary(path, classifications_by_set):
    row_order = [
        "BrumiR-RF core",
        "BrumiR-RF supported",
        "miRDeep2 core",
    ]

    category_order = [
        "Not annotated in MirGeneDB",
        "Unique annotated",
        "Shared annotated",
    ]

    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(["feature_set", "category", "count"])

        for set_name in row_order:
            values = classifications_by_set[set_name]

            for category in category_order:
                count = sum(
                    value == category
                    for value in values.values()
                )
                writer.writerow([set_name, category, count])


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    blast_bin = Path(args.blast_bin)
    blastn = blast_bin / "blastn"
    makeblastdb = blast_bin / "makeblastdb"

    for executable in (blastn, makeblastdb):
        if not executable.is_file():
            raise SystemExit(
                f"ERROR: BLAST executable not found: {executable}"
            )

    brumir_all = read_fasta(args.brumir_fasta)
    brumir_core_ids = read_core_ids(args.brumir_core)

    missing_brumir = sorted(brumir_core_ids - brumir_all.keys())

    if missing_brumir:
        raise SystemExit(
            "ERROR: BrumiR core IDs absent from representative FASTA: "
            + ", ".join(missing_brumir[:20])
        )

    brumir_core = {
        identifier: brumir_all[identifier]
        for identifier in brumir_core_ids
    }

    brumir_supported = read_supported_table(
        args.brumir_supported_table
    )

    mirbase = read_fasta(args.mirbase_mature)
    mirdeep2_ids = read_core_ids(args.mirdeep2_core)

    missing_mirdeep2 = sorted(mirdeep2_ids - mirbase.keys())

    if missing_mirdeep2:
        raise SystemExit(
            "ERROR: miRDeep2 core IDs absent from miRBase: "
            + ", ".join(missing_mirdeep2[:20])
        )

    mirdeep2_core = {
        identifier: mirbase[identifier]
        for identifier in mirdeep2_ids
    }

    fasta_paths = {
        "BrumiR-RF core":
            outdir / "Figure3C_BrumiR_RF_core_union.fasta",
        "BrumiR-RF supported":
            outdir / "Figure3C_BrumiR_RF_supported.fasta",
        "miRDeep2 core":
            outdir / "Figure3C_miRDeep2_core_union.fasta",
    }

    record_sets = {
        "BrumiR-RF core": brumir_core,
        "BrumiR-RF supported": brumir_supported,
        "miRDeep2 core": mirdeep2_core,
    }

    for set_name in record_sets:
        write_fasta(record_sets[set_name], fasta_paths[set_name])

    database = outdir / "mirgenedb_hsa_mature_db"

    run_command([
        str(makeblastdb),
        "-in", str(args.mirgenedb_mature),
        "-dbtype", "nucl",
        "-parse_seqids",
        "-out", str(database),
    ])

    raw_paths = {}
    filtered_paths = {}
    filtered_rows = {}

    filename_labels = {
        "BrumiR-RF core": "BrumiR_RF_core",
        "BrumiR-RF supported": "BrumiR_RF_supported",
        "miRDeep2 core": "miRDeep2_core",
    }

    for set_name in record_sets:
        label = filename_labels[set_name]

        raw_paths[set_name] = (
            outdir / f"Figure3C_{label}_vs_MirGeneDB.raw.tsv"
        )

        filtered_paths[set_name] = (
            outdir
            / f"Figure3C_{label}_vs_MirGeneDB.filtered.tsv"
        )

        run_blast(
            fasta_paths[set_name],
            database,
            raw_paths[set_name],
            blastn,
        )

        filtered_rows[set_name] = filter_hits(
            raw_paths[set_name],
            filtered_paths[set_name],
            args.min_identity,
            args.min_query_coverage,
            args.max_mismatches,
        )

    subject_maps = {
        set_name: subjects_by_query(filtered_rows[set_name])
        for set_name in filtered_rows
    }

    brumir_subjects = {
        subject
        for subjects in subject_maps["BrumiR-RF core"].values()
        for subject in subjects
    }

    mirdeep2_subjects = {
        subject
        for subjects in subject_maps["miRDeep2 core"].values()
        for subject in subjects
    }

    shared_subjects = brumir_subjects & mirdeep2_subjects

    shared_path = outdir / "Figure3C_shared_MirGeneDB_IDs.tsv"

    with shared_path.open("w") as handle:
        handle.write("MirGeneDB_ID\n")

        for subject in sorted(shared_subjects):
            handle.write(subject + "\n")

    classifications_by_set = {
        set_name: classify_queries(
            set(record_sets[set_name]),
            subject_maps[set_name],
            shared_subjects,
        )
        for set_name in record_sets
    }

    detail_path = outdir / "Figure3C_annotation_detail.tsv"
    summary_path = outdir / "Figure3C_annotation_summary.tsv"

    write_detail(
        detail_path,
        classifications_by_set,
        subject_maps,
    )

    write_summary(
        summary_path,
        classifications_by_set,
    )

    print()
    print("===== FIGURE 3C SUMMARY =====")

    for set_name in (
        "BrumiR-RF core",
        "BrumiR-RF supported",
        "miRDeep2 core",
    ):
        values = classifications_by_set[set_name]

        print(f"\n{set_name}: {len(values)}")

        for category in (
            "Not annotated in MirGeneDB",
            "Unique annotated",
            "Shared annotated",
        ):
            count = sum(value == category for value in values.values())
            print(f"  {category}: {count}")

    print(f"\nShared MirGeneDB IDs: {len(shared_subjects)}")
    print(f"Summary written: {summary_path}")
    print(f"Detail written: {detail_path}")


if __name__ == "__main__":
    main()
