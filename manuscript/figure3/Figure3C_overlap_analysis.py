#!/usr/bin/env python3

"""
Figure 3C – BrumiR-RF and miRDeep2 known overlap analysis

Description
-----------
This script compares all BrumiR-RF candidate sequences before CD-HIT
against miRDeep2 known miRNA sequences to quantify exact sequence
overlap.

Sequences are normalized to uppercase and U is replaced by T before
comparison. The script reports exact shared sequences, BrumiR-only
sequences, miRDeep2-only sequences, and pairwise exact matches.

Inputs
------
--brumir_rf
    BrumiR-RF FASTA file containing all retained RF candidates.

--mirdeep2_known
    miRDeep2 known FASTA file.

--outdir
    Output directory.

Outputs
-------
overlap_exact_matches.tsv
overlap_summary.tsv
overlap_exact_sequences.fasta

Usage
-----
python3 Figure3C_overlap_analysis.py \\
  --brumir_rf path/to/brumir_rf.fasta \\
  --mirdeep2_known path/to/mirdeep2_known.fasta \\
  --outdir path/to/output_dir
"""

import argparse
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, List, Tuple


FastaRecord = Tuple[str, str]


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def read_fasta(path: Path) -> List[FastaRecord]:
    """Read FASTA records and normalize sequences to uppercase DNA alphabet."""
    records: List[FastaRecord] = []
    header = None
    sequence_chunks = []

    with path.open(encoding="utf-8", errors="ignore") as fasta_handle:
        for raw_line in fasta_handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_chunks).upper().replace("U", "T")
                    records.append((header, sequence))

                header = line[1:].strip()
                sequence_chunks = []

            else:
                sequence_chunks.append(line)

        if header is not None:
            sequence = "".join(sequence_chunks).upper().replace("U", "T")
            records.append((header, sequence))

    return records


def write_fasta(records: List[FastaRecord], outpath: Path) -> None:
    """Write FASTA records."""
    with outpath.open("w", encoding="utf-8") as fasta_handle:
        for header, sequence in records:
            fasta_handle.write(f">{header}\n{sequence}\n")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare BrumiR-RF candidates against miRDeep2 known miRNAs "
            "and report exact sequence overlap."
        )
    )

    parser.add_argument(
        "--brumir_rf",
        required=True,
        help="BrumiR-RF FASTA file.",
    )
    parser.add_argument(
        "--mirdeep2_known",
        required=True,
        help="miRDeep2 known FASTA file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )

    return parser.parse_args()


def main() -> None:
    """Run exact overlap analysis between BrumiR-RF and miRDeep2 known miRNAs."""
    args = parse_args()

    brumir_path = validate_file(args.brumir_rf, "BrumiR-RF FASTA")
    mirdeep_path = validate_file(args.mirdeep2_known, "miRDeep2 known FASTA")
    outdir = Path(args.outdir)

    outdir.mkdir(parents=True, exist_ok=True)

    brumir_records = read_fasta(brumir_path)
    mirdeep_records = read_fasta(mirdeep_path)

    print("BrumiR-RF sequences loaded:", len(brumir_records))
    print("miRDeep2 known sequences loaded:", len(mirdeep_records))

    brumir_by_seq: DefaultDict[str, List[str]] = defaultdict(list)
    mirdeep_by_seq: DefaultDict[str, List[str]] = defaultdict(list)

    for header, sequence in brumir_records:
        brumir_by_seq[sequence].append(header)

    for header, sequence in mirdeep_records:
        mirdeep_by_seq[sequence].append(header)

    brumir_set = set(brumir_by_seq.keys())
    mirdeep_set = set(mirdeep_by_seq.keys())

    shared_exact = sorted(brumir_set & mirdeep_set)
    brumir_only = sorted(brumir_set - mirdeep_set)
    mirdeep_only = sorted(mirdeep_set - brumir_set)

    exact_rows = []
    exact_fasta = []

    for sequence in shared_exact:
        for brumir_header in brumir_by_seq[sequence]:
            for mirdeep_header in mirdeep_by_seq[sequence]:
                exact_rows.append({
                    "sequence": sequence,
                    "brumir_rf_id": brumir_header,
                    "mirdeep2_known_id": mirdeep_header,
                    "match_type": "exact",
                })

        exact_fasta.append((f"shared_exact_{len(exact_fasta) + 1}", sequence))

    exact_tsv = outdir / "overlap_exact_matches.tsv"
    with exact_tsv.open("w", encoding="utf-8") as out_handle:
        out_handle.write("sequence\tbrumir_rf_id\tmirdeep2_known_id\tmatch_type\n")

        for row in exact_rows:
            out_handle.write(
                f"{row['sequence']}\t"
                f"{row['brumir_rf_id']}\t"
                f"{row['mirdeep2_known_id']}\t"
                f"{row['match_type']}\n"
            )

    exact_fasta_out = outdir / "overlap_exact_sequences.fasta"
    write_fasta(exact_fasta, exact_fasta_out)

    summary_tsv = outdir / "overlap_summary.tsv"
    with summary_tsv.open("w", encoding="utf-8") as out_handle:
        out_handle.write("metric\tvalue\n")
        out_handle.write(f"BrumiR_RF_total_sequences\t{len(brumir_set)}\n")
        out_handle.write(f"miRDeep2_known_total_sequences\t{len(mirdeep_set)}\n")
        out_handle.write(f"shared_exact_sequences\t{len(shared_exact)}\n")
        out_handle.write(f"BrumiR_RF_only_sequences\t{len(brumir_only)}\n")
        out_handle.write(f"miRDeep2_known_only_sequences\t{len(mirdeep_only)}\n")
        out_handle.write(f"shared_exact_pairwise_matches\t{len(exact_rows)}\n")

    print("Exact shared sequences:", len(shared_exact))
    print("BrumiR-RF only:", len(brumir_only))
    print("miRDeep2 known only:", len(mirdeep_only))
    print("Written:", exact_tsv)
    print("Written:", exact_fasta_out)
    print("Written:", summary_tsv)


if __name__ == "__main__":
    main()
