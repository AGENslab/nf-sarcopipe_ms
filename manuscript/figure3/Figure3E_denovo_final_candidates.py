#!/usr/bin/env python3

"""
Figure 3E – De novo final candidate catalog

Description
-----------
This script builds a catalog of BrumiR core candidates supported by
BrumiR2Reference passfilter results and classifies them as known or
de novo based on exact mature-sequence matching against miRBase and
MirGeneDB.

It produces:
1. A candidate-level catalog of supported mature sequences found inside
   BrumiR2Reference precursor sequences.
2. A summary table of unique supported candidates.
3. A group-level count table for downstream plotting.

Inputs
------
--brumir_core_fasta
    FASTA file with BrumiR core mature sequences.

--passfilter_tsv
    BrumiR2Reference passfilter table.

--core_sets_tsv
    Core sets TSV from Module II.

--mirbase_mature
    miRBase mature FASTA file.

--mirgenedb_mature
    MirGeneDB mature FASTA file.

--out_catalog
    Output TSV with supported candidate catalog.

--out_summary
    Output TSV with summary metrics.

--out_counts
    Output TSV with group-level counts.

Outputs
-------
denovo_passfilter_catalog.tsv
denovo_passfilter_summary.tsv
status_counts_by_group.tsv

Usage
-----
python3 Figure3E_denovo_final_candidates.py \\
  --brumir_core_fasta brumir_core_sequences.fa \\
  --passfilter_tsv brumir2ref.passfilter.tsv \\
  --core_sets_tsv brumir.core_sets.tsv \\
  --mirbase_mature hsa_mature.fa \\
  --mirgenedb_mature mirgenedb_hsa_mature.fa \\
  --out_catalog denovo_passfilter_catalog.tsv \\
  --out_summary denovo_passfilter_summary.tsv \\
  --out_counts status_counts_by_group.tsv
"""

import argparse
import re
import sys
from pathlib import Path


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"[ERROR] {label} not found: {file_path}")

    return file_path


def prepare_output(path: str) -> Path:
    """Create the parent directory for an output file if needed."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path


def read_fasta(path: Path) -> dict:
    """Read FASTA into a dict: header_id -> sequence in RNA alphabet."""
    sequences = {}
    name = None
    buffer = []

    with path.open() as fasta_handle:
        for raw_line in fasta_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(buffer).upper().replace("T", "U")

                name = line[1:].split()[0]
                buffer = []

            else:
                buffer.append(re.sub(r"\s+", "", line))

        if name is not None:
            sequences[name] = "".join(buffer).upper().replace("T", "U")

    return sequences


def fasta_values(path: Path) -> list:
    """Return all sequence values from a FASTA file."""
    return list(read_fasta(path).values())


def load_known(mirbase_path: Path, mirgenedb_path: Path) -> set:
    """Load known mature sequences from miRBase and MirGeneDB into one set."""
    known = set()

    for sequence in fasta_values(mirbase_path):
        if sequence:
            known.add(sequence.upper().replace("T", "U"))

    for sequence in fasta_values(mirgenedb_path):
        if sequence:
            known.add(sequence.upper().replace("T", "U"))

    return known


def parse_core_sets(path: Path) -> dict:
    """
    Read core sets TSV into candidate_id -> set(core_set_labels).

    Expected format:
    set<TAB>miRNA
    """
    mapping = {}

    with path.open() as core_handle:
        _header = core_handle.readline()

        for raw_line in core_handle:
            if not raw_line.strip():
                continue

            fields = raw_line.rstrip("\n").split("\t")

            if len(fields) < 2:
                continue

            core_set, candidate_id = fields[0], fields[1]
            mapping.setdefault(candidate_id, set()).add(core_set)

    return mapping


def core_class(core_sets: set) -> str:
    """Convert core-set labels into a simplified group label."""
    if not core_sets:
        return "NA"

    sets = set(core_sets)
    athlete = "athlete_core" in sets
    sedentary = "sedentary_core" in sets

    if athlete and sedentary:
        return "shared"

    if athlete:
        return "athlete"

    if sedentary:
        return "sedentary"

    if "core_union" in sets:
        return "core_union"

    return "other"


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build a catalog of BrumiR passfilter-supported core candidates "
            "and classify them as known or de novo."
        )
    )

    parser.add_argument(
        "--brumir_core_fasta",
        required=True,
        help="BrumiR core mature FASTA.",
    )
    parser.add_argument(
        "--passfilter_tsv",
        required=True,
        help="BrumiR2Reference passfilter TSV.",
    )
    parser.add_argument(
        "--core_sets_tsv",
        required=True,
        help="BrumiR core sets TSV.",
    )
    parser.add_argument(
        "--mirbase_mature",
        required=True,
        help="miRBase mature FASTA.",
    )
    parser.add_argument(
        "--mirgenedb_mature",
        required=True,
        help="MirGeneDB mature FASTA.",
    )
    parser.add_argument(
        "--out_catalog",
        default="denovo_passfilter_catalog.tsv",
        help="Output catalog TSV.",
    )
    parser.add_argument(
        "--out_summary",
        default="denovo_passfilter_summary.tsv",
        help="Output summary TSV.",
    )
    parser.add_argument(
        "--out_counts",
        default="status_counts_by_group.tsv",
        help="Output count TSV.",
    )

    return parser.parse_args()


def main() -> None:
    """Build de novo passfilter-supported candidate catalog."""
    args = parse_args()

    brumir_core_fasta = validate_file(args.brumir_core_fasta, "BrumiR core FASTA")
    passfilter_tsv = validate_file(args.passfilter_tsv, "Passfilter TSV")
    core_sets_tsv = validate_file(args.core_sets_tsv, "Core sets TSV")
    mirbase_mature = validate_file(args.mirbase_mature, "miRBase mature FASTA")
    mirgenedb_mature = validate_file(args.mirgenedb_mature, "MirGeneDB mature FASTA")

    out_catalog = prepare_output(args.out_catalog)
    out_summary = prepare_output(args.out_summary)
    out_counts_file = prepare_output(args.out_counts)

    core = read_fasta(brumir_core_fasta)
    known = load_known(mirbase_mature, mirgenedb_mature)
    core_sets = parse_core_sets(core_sets_tsv)

    with passfilter_tsv.open() as passfilter_handle:
        header = passfilter_handle.readline().rstrip("\n").split("\t")
        columns = {column: idx for idx, column in enumerate(header)}

        required_columns = ["miRNA", "chr", "start", "stop", "MFE", "Precursor_Seq"]
        for column in required_columns:
            if column not in columns:
                raise SystemExit(
                    f"[ERROR] passfilter missing column '{column}'. Found: {header}"
                )

        pass_rows = []

        for raw_line in passfilter_handle:
            if not raw_line.strip() or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split("\t")

            pass_rows.append({
                "passfilter_miRNA": fields[columns["miRNA"]],
                "chr": fields[columns["chr"]],
                "start": fields[columns["start"]],
                "stop": fields[columns["stop"]],
                "brumir2ref_mfe": fields[columns["MFE"]],
                "precursor_seq": fields[columns["Precursor_Seq"]].upper().replace("T", "U"),
            })

    by_len = {}

    for candidate_id, mature_sequence in core.items():
        by_len.setdefault(len(mature_sequence), []).append(
            (candidate_id, mature_sequence)
        )

    catalog = []
    supported_candidates = set()

    for row in pass_rows:
        precursor = row["precursor_seq"]
        hit = None

        for sequence_length in sorted(by_len.keys()):
            for candidate_id, mature_sequence in by_len[sequence_length]:
                position = precursor.find(mature_sequence)

                if position != -1:
                    hit = (candidate_id, mature_sequence, position)
                    break

            if hit:
                break

        if not hit:
            continue

        candidate_id, mature_sequence, position0 = hit
        supported_candidates.add((candidate_id, mature_sequence))

        sets = core_sets.get(candidate_id, set())
        candidate_group = core_class(sets)
        is_known = "1" if mature_sequence in known else "0"
        is_denovo = "1" if mature_sequence not in known else "0"

        catalog.append({
            "candidate": candidate_id,
            "group": candidate_group,
            "core_sets": ",".join(sorted(sets)) if sets else "NA",
            "passfilter_miRNA": row["passfilter_miRNA"],
            "chr": row["chr"],
            "start": row["start"],
            "stop": row["stop"],
            "brumir2ref_mfe": row["brumir2ref_mfe"],
            "precursor_len": str(len(precursor)),
            "mature_len": str(len(mature_sequence)),
            "mature_start": str(position0 + 1),
            "mature_end": str(position0 + len(mature_sequence)),
            "is_known_exact": is_known,
            "is_denovo_exact": is_denovo,
            "mature_seq": mature_sequence,
            "precursor_seq": precursor,
        })

    def all_core_in_group(group: str) -> set:
        return {
            (candidate_id, mature_sequence)
            for candidate_id, mature_sequence in core.items()
            if core_class(core_sets.get(candidate_id, set())) == group
        }

    groups = ["athlete", "sedentary", "shared"]
    out_counts = []
    supported_set = set(supported_candidates)

    for group in groups:
        core_group = all_core_in_group(group)
        total = len(core_group)
        supported = len(core_group & supported_set)
        not_supported = total - supported

        denovo_core = {
            (candidate_id, mature_sequence)
            for candidate_id, mature_sequence in core_group
            if mature_sequence not in known
        }

        denovo_supported = len(denovo_core & supported_set)

        out_counts.append((group, "core_total", total))
        out_counts.append((group, "supported", supported))
        out_counts.append((group, "not_supported", not_supported))
        out_counts.append((group, "denovo_supported_exact", denovo_supported))

    total_supported = len(supported_set)
    denovo_supported = len({
        candidate for candidate in supported_set
        if candidate[1] not in known
    })
    known_supported = total_supported - denovo_supported

    output_columns = [
        "candidate",
        "group",
        "core_sets",
        "passfilter_miRNA",
        "chr",
        "start",
        "stop",
        "brumir2ref_mfe",
        "precursor_len",
        "mature_len",
        "mature_start",
        "mature_end",
        "is_known_exact",
        "is_denovo_exact",
        "mature_seq",
        "precursor_seq",
    ]

    with out_catalog.open("w") as catalog_handle:
        catalog_handle.write("\t".join(output_columns) + "\n")

        for row in catalog:
            catalog_handle.write(
                "\t".join(row[column] for column in output_columns) + "\n"
            )

    with out_summary.open("w") as summary_handle:
        summary_handle.write("metric\tvalue\n")
        summary_handle.write(f"supported_unique_candidates\t{total_supported}\n")
        summary_handle.write(f"known_exact_supported\t{known_supported}\n")
        summary_handle.write(f"denovo_exact_supported\t{denovo_supported}\n")

    with out_counts_file.open("w") as counts_handle:
        counts_handle.write("group\tmetric\tcount\n")

        for group, metric, count in out_counts:
            counts_handle.write(f"{group}\t{metric}\t{count}\n")

    print("=== SUMMARY (core ∩ passfilter) ===", file=sys.stderr)
    print(f"Supported unique candidates: {total_supported}", file=sys.stderr)
    print(
        f"Known exact supported (miRBase ∪ MirGeneDB): {known_supported}",
        file=sys.stderr,
    )
    print(f"De novo exact supported: {denovo_supported}", file=sys.stderr)
    print(
        f"[OK] Wrote: {out_catalog}, {out_summary}, {out_counts_file}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
