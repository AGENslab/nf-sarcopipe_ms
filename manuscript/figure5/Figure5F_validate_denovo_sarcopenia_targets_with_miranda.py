#!/usr/bin/env python3
"""
Figure 5E – Validate musculoskeletal de novo miRNA–gene pairs with miRanda

Description
-----------
Validates BrumiR-derived de novo miRNA–gene pairs identified by exact
canonical 7-mer seed matching against a curated musculoskeletal gene set.

The script:

1. Reads the Figure 5E seed-matched interaction table.
2. Retains only BrumiR-derived pairs.
3. Retrieves mature miRNA sequences from a supplied sequence catalog.
4. Writes miRNA and 3'UTR FASTA files.
5. Runs miRanda using the same default thresholds as nf-Sarcopipe:
      score >= 140
      energy <= -20 kcal/mol
6. Intersects miRanda predictions with the original seed-matched pairs.
7. Retains the best miRanda hit per miRNA–gene pair.

No pipeline files, miRNA IDs, genes, interactions, or results are
hardcoded.

Required input columns
----------------------
Seed-match table:
    gene_symbol
    renamed_miRNA
    original_id
    source
    seed
    seed_rc
    n_sites
    site_positions
    log2FoldChange
    padj
    direction

Sequence catalog:
    An identifier column matching original_id or renamed_miRNA
    and one mature-sequence column.

3'UTR table:
    gene_symbol
    utr_3

Outputs
-------
- Figure5E_denovo_miRNAs.fa
- Figure5E_denovo_targets.fa
- Figure5E_miranda.raw.txt
- Figure5E_miranda_all_predictions.tsv
- Figure5E_denovo_seed_pairs_miranda_validated.tsv
- Figure5E_denovo_seed_pairs_not_validated.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path


VALIDATED_FIELDS = [
    "gene_symbol",
    "category",
    "renamed_miRNA",
    "original_id",
    "source",
    "seed",
    "seed_rc",
    "n_sites",
    "site_positions",
    "log2FoldChange",
    "padj",
    "direction",
    "miranda_score",
    "miranda_energy",
    "validated_by",
]

ALL_PREDICTION_FIELDS = [
    "renamed_miRNA",
    "gene_symbol",
    "miranda_score",
    "miranda_energy",
]

ID_COLUMNS = [
    "original_id",
    "candidate",
    "cluster",
    "cluster_id",
    "renamed_miRNA",
    "provisory_id",
    "provisional_id",
    "miRNA",
    "id",
    "ID",
]

SEQUENCE_COLUMNS = [
    "mature_seq",
    "mature_sequence",
    "sequence",
    "seq",
    "Mature sequence (5'->3')",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Figure 5E BrumiR musculoskeletal seed-matched "
            "pairs using miRanda."
        )
    )
    parser.add_argument(
        "--seed_matches",
        required=True,
        help="Figure 5E musculoskeletal seed-match TSV.",
    )
    parser.add_argument(
        "--sequence_catalog",
        required=True,
        help="TSV containing de novo miRNA IDs and mature sequences.",
    )
    parser.add_argument(
        "--utr_table",
        required=True,
        help="TSV containing gene_symbol and utr_3.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )
    parser.add_argument(
        "--score",
        default="140",
        help="Minimum miRanda score. Default: 140.",
    )
    parser.add_argument(
        "--energy",
        default="-20",
        help="Maximum miRanda free energy. Default: -20.",
    )
    parser.add_argument(
        "--miranda_bin",
        default="miranda",
        help="miRanda executable. Default: miranda.",
    )
    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def normalize_rna(sequence: str) -> str:
    return (
        sequence.strip()
        .upper()
        .replace("T", "U")
        .replace(" ", "")
        .replace("-", "")
        .replace("\r", "")
        .replace("\n", "")
    )


def detect_column(fieldnames: list[str], candidates: list[str], label: str) -> str:
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate

    raise SystemExit(
        f"ERROR: could not identify {label}. "
        f"Available columns: {', '.join(fieldnames)}"
    )


def read_seed_pairs(path: Path) -> list[dict[str, str]]:
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "gene_symbol",
            "category",
            "renamed_miRNA",
            "original_id",
            "source",
            "seed",
            "seed_rc",
            "n_sites",
            "site_positions",
            "log2FoldChange",
            "padj",
            "direction",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: seed-match table missing columns: "
                + ", ".join(sorted(missing))
            )

        rows = []

        for row in reader:
            if row["source"].strip() != "BrumiR":
                continue

            gene = row["gene_symbol"].strip()
            mirna = row["renamed_miRNA"].strip()
            original_id = row["original_id"].strip()

            if not gene or not mirna:
                continue

            clean_row = {
                key: row.get(key, "").strip()
                for key in required
            }

            clean_row["gene_symbol"] = gene
            clean_row["renamed_miRNA"] = mirna
            clean_row["original_id"] = original_id
            clean_row["source"] = "BrumiR"

            rows.append(clean_row)

    if not rows:
        raise SystemExit("ERROR: no BrumiR seed-matched pairs were found")

    unique = {}

    for row in rows:
        key = (
            row["renamed_miRNA"],
            row["gene_symbol"],
        )
        unique[key] = row

    return list(unique.values())


def read_sequence_catalog(path: Path) -> dict[str, str]:
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])

        sequence_col = detect_column(
            fieldnames,
            SEQUENCE_COLUMNS,
            "mature-sequence column",
        )

        available_id_columns = [
            column for column in ID_COLUMNS
            if column in fieldnames
        ]

        if not available_id_columns:
            raise SystemExit(
                "ERROR: no usable miRNA identifier column found. "
                f"Available columns: {', '.join(fieldnames)}"
            )

        sequences: dict[str, str] = {}

        for row in reader:
            sequence = normalize_rna(row.get(sequence_col, ""))

            if not sequence:
                continue

            for id_col in available_id_columns:
                identifier = row.get(id_col, "").strip()

                if identifier:
                    sequences[identifier] = sequence

    if not sequences:
        raise SystemExit(
            "ERROR: no mature miRNA sequences were loaded from the catalog"
        )

    return sequences


def read_utrs(path: Path) -> dict[str, str]:
    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"gene_symbol", "utr_3"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: UTR table missing columns: "
                + ", ".join(sorted(missing))
            )

        utrs = {}

        for row in reader:
            gene = row["gene_symbol"].strip()
            sequence = normalize_rna(row["utr_3"])

            if gene and sequence:
                utrs[gene] = sequence

    if not utrs:
        raise SystemExit("ERROR: no 3'UTR sequences were loaded")

    return utrs


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def run_miranda(
    miranda_bin: str,
    mirna_fasta: Path,
    target_fasta: Path,
    raw_output: Path,
    score: str,
    energy: str,
) -> None:
    command = [
        miranda_bin,
        str(mirna_fasta),
        str(target_fasta),
        "-sc",
        str(score),
        "-en",
        str(energy),
        "-out",
        str(raw_output),
    ]

    print("Running:", " ".join(command))
    subprocess.run(command, check=True)


def parse_miranda(path: Path) -> list[dict[str, str]]:
    predictions = []

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            if not line.startswith(">>"):
                continue

            parts = line.replace(">>", "", 1).strip().split()

            if len(parts) < 4:
                continue

            predictions.append(
                {
                    "renamed_miRNA": parts[0],
                    "gene_symbol": parts[1],
                    "miranda_score": parts[2],
                    "miranda_energy": parts[3],
                }
            )

    return predictions


def as_float(value: str, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def select_best_predictions(
    predictions: list[dict[str, str]],
) -> dict[tuple[str, str], dict[str, str]]:
    best = {}

    for prediction in predictions:
        key = (
            prediction["renamed_miRNA"],
            prediction["gene_symbol"],
        )

        if key not in best:
            best[key] = prediction
            continue

        current = best[key]

        new_score = as_float(
            prediction["miranda_score"],
            float("-inf"),
        )
        current_score = as_float(
            current["miranda_score"],
            float("-inf"),
        )

        new_energy = as_float(
            prediction["miranda_energy"],
            float("inf"),
        )
        current_energy = as_float(
            current["miranda_energy"],
            float("inf"),
        )

        if (
            new_score > current_score
            or (
                new_score == current_score
                and new_energy < current_energy
            )
        ):
            best[key] = prediction

    return best


def write_tsv(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()

    seed_path = Path(args.seed_matches)
    catalog_path = Path(args.sequence_catalog)
    utr_path = Path(args.utr_table)
    outdir = Path(args.outdir)

    validate_file(seed_path, "seed-match table")
    validate_file(catalog_path, "de novo sequence catalog")
    validate_file(utr_path, "3'UTR table")

    outdir.mkdir(parents=True, exist_ok=True)

    seed_pairs = read_seed_pairs(seed_path)
    sequence_catalog = read_sequence_catalog(catalog_path)
    utrs = read_utrs(utr_path)

    required_mirnas = sorted(
        {
            row["renamed_miRNA"]
            for row in seed_pairs
        }
    )

    required_genes = sorted(
        {
            row["gene_symbol"]
            for row in seed_pairs
        }
    )

    mirna_sequences = {}
    missing_sequences = []

    for row in seed_pairs:
        renamed = row["renamed_miRNA"]
        original_id = row["original_id"]

        sequence = (
            sequence_catalog.get(original_id)
            or sequence_catalog.get(renamed)
        )

        if sequence:
            mirna_sequences[renamed] = sequence
        else:
            missing_sequences.append(
                f"{renamed} ({original_id})"
            )

    missing_sequences = sorted(set(missing_sequences))

    if missing_sequences:
        raise SystemExit(
            "ERROR: mature sequence not found for: "
            + ", ".join(missing_sequences)
        )

    missing_utrs = [
        gene for gene in required_genes
        if gene not in utrs
    ]

    if missing_utrs:
        raise SystemExit(
            "ERROR: 3'UTR not found for: "
            + ", ".join(missing_utrs)
        )

    mirna_fasta = outdir / "Figure5E_denovo_miRNAs.fa"
    target_fasta = outdir / "Figure5E_denovo_targets.fa"
    raw_output = outdir / "Figure5E_miranda.raw.txt"
    all_output = outdir / "Figure5E_miranda_all_predictions.tsv"
    validated_output = (
        outdir / "Figure5E_denovo_seed_pairs_miranda_validated.tsv"
    )
    rejected_output = (
        outdir / "Figure5E_denovo_seed_pairs_not_validated.tsv"
    )

    write_fasta(
        mirna_fasta,
        [
            (mirna, mirna_sequences[mirna])
            for mirna in required_mirnas
        ],
    )

    write_fasta(
        target_fasta,
        [
            (gene, utrs[gene])
            for gene in required_genes
        ],
    )

    run_miranda(
        miranda_bin=args.miranda_bin,
        mirna_fasta=mirna_fasta,
        target_fasta=target_fasta,
        raw_output=raw_output,
        score=args.score,
        energy=args.energy,
    )

    predictions = parse_miranda(raw_output)
    best_predictions = select_best_predictions(predictions)

    validated_rows = []
    rejected_rows = []

    for seed_pair in seed_pairs:
        key = (
            seed_pair["renamed_miRNA"],
            seed_pair["gene_symbol"],
        )

        prediction = best_predictions.get(key)

        if prediction is None:
            rejected_rows.append(seed_pair)
            continue

        validated_rows.append(
            {
                **seed_pair,
                "miranda_score": prediction["miranda_score"],
                "miranda_energy": prediction["miranda_energy"],
                "validated_by": "miRanda",
            }
        )

    validated_rows.sort(
        key=lambda row: (
            row["renamed_miRNA"],
            row["gene_symbol"],
        )
    )

    rejected_rows.sort(
        key=lambda row: (
            row["renamed_miRNA"],
            row["gene_symbol"],
        )
    )

    write_tsv(
        all_output,
        predictions,
        ALL_PREDICTION_FIELDS,
    )

    write_tsv(
        validated_output,
        validated_rows,
        VALIDATED_FIELDS,
    )

    write_tsv(
        rejected_output,
        rejected_rows,
        VALIDATED_FIELDS[:-3],
    )

    print("===== FIGURE 5E MIRANDA VALIDATION =====")
    print("Seed-predicted BrumiR pairs:", len(seed_pairs))
    print("De novo miRNAs tested:", len(required_mirnas))
    print("Musculoskeletal genes tested:", len(required_genes))
    print("All miRanda predictions:", len(predictions))
    print("Seed pairs validated by miRanda:", len(validated_rows))
    print("Seed pairs not validated:", len(rejected_rows))
    print("Written:", validated_output)
    print("Written:", rejected_output)


if __name__ == "__main__":
    main()
