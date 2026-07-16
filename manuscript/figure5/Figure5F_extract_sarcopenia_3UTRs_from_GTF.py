#!/usr/bin/env python3
"""
Figure 5E – Extract musculoskeletal 3'UTRs from GRCh38

Description
-----------
Extracts annotated 3'UTR sequences for a user-provided gene set
directly from the same GRCh38 genome and Ensembl GTF reference used
by nf-Sarcopipe.

The script:
1. reads gene symbols from a text file;
2. retrieves all annotated three_prime_utr features from the GTF;
3. concatenates 3'UTR segments by transcript, respecting strand;
4. retains the longest complete 3'UTR per gene;
5. reports genes without an annotated 3'UTR.

No genes, coordinates, transcripts, or counts are hardcoded.

Inputs
------
--genes
    Text file containing one HGNC gene symbol per line.

--gtf
    Ensembl GTF annotation containing three_prime_utr features.

--genome
    Reference genome FASTA matching the GTF assembly.

--out
    Output TSV with the longest annotated 3'UTR per gene.

--out_missing
    Output TSV listing genes without a recovered 3'UTR.

Outputs
-------
gene_symbol
transcript_id
chromosome
strand
utr_length
n_utr_segments
utr_3
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


OUTPUT_FIELDS = [
    "gene_symbol",
    "transcript_id",
    "chromosome",
    "strand",
    "utr_length",
    "n_utr_segments",
    "utr_3",
]

MISSING_FIELDS = [
    "gene_symbol",
    "status",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract longest annotated 3'UTR per gene from GRCh38."
    )
    parser.add_argument("--genes", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--out_missing", required=True)
    return parser.parse_args()


def validate_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: {label} is empty: {path}")


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_gene_symbols(path: Path) -> list[str]:
    genes = []
    seen = set()

    with open_text(path) as handle:
        for line in handle:
            gene = line.strip().upper()

            if not gene or gene.startswith("#"):
                continue

            if gene not in seen:
                genes.append(gene)
                seen.add(gene)

    if not genes:
        raise SystemExit("ERROR: no gene symbols were loaded")

    return genes


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current_id = None

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_id = line[1:].split()[0]
                sequences[current_id] = []
            else:
                if current_id is None:
                    raise SystemExit("ERROR: malformed FASTA")
                sequences[current_id].append(line.upper())

    genome = {
        seq_id: "".join(parts)
        for seq_id, parts in sequences.items()
    }

    if not genome:
        raise SystemExit("ERROR: no sequences were loaded from genome FASTA")

    return genome


def parse_attributes(attributes: str) -> dict[str, str]:
    parsed = {}

    for match in re.finditer(r'(\S+)\s+"([^"]*)"', attributes):
        key, value = match.groups()
        parsed[key] = value

    return parsed


def normalize_chromosome(chromosome: str, genome: dict[str, str]) -> str | None:
    candidates = [
        chromosome,
        chromosome.removeprefix("chr"),
        f"chr{chromosome}" if not chromosome.startswith("chr") else chromosome,
    ]

    for candidate in candidates:
        if candidate in genome:
            return candidate

    return None


def reverse_complement(sequence: str) -> str:
    table = str.maketrans(
        {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "N": "N",
        }
    )
    return sequence.translate(table)[::-1]


def load_utr_features(
    gtf_path: Path,
    requested_genes: set[str],
    genome: dict[str, str],
):
    transcripts = defaultdict(
        lambda: {
            "gene_symbol": "",
            "chromosome": "",
            "strand": "",
            "segments": [],
        }
    )

    accepted_features = {
        "three_prime_utr",
        "3UTR",
        "three_prime_UTR",
    }

    with open_text(gtf_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if not raw_line or raw_line.startswith("#"):
                continue

            fields = raw_line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chromosome, _, feature, start, end, _, strand, _, attributes = fields

            if feature not in accepted_features:
                continue

            attrs = parse_attributes(attributes)

            gene_symbol = (
                attrs.get("gene_name")
                or attrs.get("gene_symbol")
                or attrs.get("gene_id")
                or ""
            ).upper()

            transcript_id = attrs.get("transcript_id", "")

            if gene_symbol not in requested_genes or not transcript_id:
                continue

            normalized_chr = normalize_chromosome(chromosome, genome)

            if normalized_chr is None:
                continue

            try:
                start_int = int(start)
                end_int = int(end)
            except ValueError:
                raise SystemExit(
                    f"ERROR: invalid coordinates in GTF line {line_number}"
                )

            record = transcripts[transcript_id]
            record["gene_symbol"] = gene_symbol
            record["chromosome"] = normalized_chr
            record["strand"] = strand
            record["segments"].append((start_int, end_int))

    return transcripts


def assemble_transcript_utrs(
    transcripts,
    genome: dict[str, str],
) -> list[dict[str, str | int]]:
    rows = []

    for transcript_id, record in transcripts.items():
        chromosome = record["chromosome"]
        strand = record["strand"]
        segments = record["segments"]

        if not segments:
            continue

        ordered_segments = sorted(segments)

        sequence_parts = []

        for start, end in ordered_segments:
            # GTF coordinates are 1-based inclusive.
            sequence_parts.append(genome[chromosome][start - 1:end])

        sequence = "".join(sequence_parts)

        if strand == "-":
            sequence = reverse_complement(sequence)

        sequence = sequence.upper().replace("T", "U")

        rows.append(
            {
                "gene_symbol": record["gene_symbol"],
                "transcript_id": transcript_id,
                "chromosome": chromosome,
                "strand": strand,
                "utr_length": len(sequence),
                "n_utr_segments": len(segments),
                "utr_3": sequence,
            }
        )

    return rows


def select_longest_per_gene(rows):
    best = {}

    for row in rows:
        gene = row["gene_symbol"]

        if gene not in best:
            best[gene] = row
            continue

        current = best[gene]

        if row["utr_length"] > current["utr_length"]:
            best[gene] = row
        elif (
            row["utr_length"] == current["utr_length"]
            and row["transcript_id"] < current["transcript_id"]
        ):
            best[gene] = row

    return [best[gene] for gene in sorted(best)]


def main() -> None:
    args = parse_args()

    genes_path = Path(args.genes)
    gtf_path = Path(args.gtf)
    genome_path = Path(args.genome)
    out_path = Path(args.out)
    missing_path = Path(args.out_missing)

    validate_file(genes_path, "gene-symbol list")
    validate_file(gtf_path, "GTF annotation")
    validate_file(genome_path, "reference genome")

    genes = read_gene_symbols(genes_path)
    requested = set(genes)

    print("Loading reference genome...")
    genome = read_fasta(genome_path)

    print("Reading annotated 3'UTR features...")
    transcripts = load_utr_features(
        gtf_path,
        requested,
        genome,
    )

    transcript_rows = assemble_transcript_utrs(
        transcripts,
        genome,
    )

    selected_rows = select_longest_per_gene(transcript_rows)

    recovered = {row["gene_symbol"] for row in selected_rows}
    missing_genes = [gene for gene in genes if gene not in recovered]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    missing_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selected_rows)

    with missing_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MISSING_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()

        for gene in missing_genes:
            writer.writerow(
                {
                    "gene_symbol": gene,
                    "status": "annotated_3UTR_not_recovered",
                }
            )

    print("===== MUSCULOSKELETAL 3'UTR EXTRACTION =====")
    print("Genes requested:", len(genes))
    print("Transcripts with annotated 3'UTR:", len(transcript_rows))
    print("Genes with recovered 3'UTR:", len(selected_rows))
    print("Genes without recovered 3'UTR:", len(missing_genes))

    if missing_genes:
        print("Missing:", ", ".join(missing_genes))

    print("Written:", out_path)
    print("Written:", missing_path)


if __name__ == "__main__":
    main()
