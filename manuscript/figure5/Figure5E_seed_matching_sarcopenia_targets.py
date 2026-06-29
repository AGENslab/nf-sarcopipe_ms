#!/usr/bin/env python3

"""
Figure 5E – Seed matching against sarcopenia targets

Description
-----------
This script performs seed-based miRNA target prediction against a
curated sarcopenia/exercise-related gene set using available 3'UTR
sequences. The same methodology is applied to BrumiR and miRDeep2
miRNAs.

Inputs
------
--genes
    Text file containing sarcopenia/exercise-related gene symbols.

--mirna
    TSV file containing miRNA seed summary information.

--utr
    TSV file containing gene_symbol and utr_3 columns.

--out
    Output TSV file with seed matches.

Outputs
-------
seed_matches_sarcopenia.tsv

Usage
-----
python3 Figure5E_seed_matching_sarcopenia_targets.py \\
  --genes input_sarcopenia_gene_set.txt \\
  --mirna all_miRNA_seed_summary.tsv \\
  --utr sarcopenia_gene_set_3UTR.tsv \\
  --out seed_matches_sarcopenia.tsv
"""

import argparse
import csv
from pathlib import Path


FIELDNAMES = [
    "gene_symbol",
    "renamed_miRNA",
    "source",
    "seed",
    "seed_rc",
    "n_sites",
    "site_positions",
    "log2FoldChange",
    "direction",
]


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


def revcomp_rna(seed: str) -> str:
    """Return the reverse complement of an RNA seed sequence."""
    complement = str.maketrans({"A": "U", "U": "A", "G": "C", "C": "G", "T": "A"})
    seed_rna = seed.upper().replace("T", "U")

    return seed_rna.translate(complement)[::-1]


def load_mirnas(mirna_file: Path) -> list:
    """Load miRNAs and seed information from the seed summary table."""
    mirnas = []

    with mirna_file.open(encoding="utf-8", errors="ignore") as handle:
        first = True

        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if not line:
                continue

            if first:
                first = False
                continue

            fields = line.split("\t")

            if len(fields) < 7:
                continue

            source, original_id, renamed_miRNA, seed, log2fc, padj, direction = fields[:7]

            seed = seed.strip().upper().replace("T", "U")

            if not seed:
                continue

            try:
                log2fc_val = float(log2fc)
            except ValueError:
                log2fc_val = 0.0

            mirnas.append({
                "source": source.strip(),
                "original_id": original_id.strip(),
                "renamed_miRNA": renamed_miRNA.strip(),
                "seed": seed,
                "seed_rc": revcomp_rna(seed),
                "log2FoldChange": log2fc_val,
                "padj": padj.strip(),
                "direction": direction.strip(),
            })

    return mirnas


def load_utrs(utr_file: Path) -> dict:
    """Load 3'UTR sequences by gene symbol."""
    utrs = {}

    with utr_file.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_cols = {"gene_symbol", "utr_3"}
        missing_cols = required_cols - set(reader.fieldnames or [])

        if missing_cols:
            raise SystemExit(
                "ERROR: UTR table missing required columns: "
                f"{', '.join(sorted(missing_cols))}"
            )

        for row in reader:
            gene = row["gene_symbol"].strip()
            sequence = row["utr_3"].strip().upper().replace("T", "U")

            if gene and sequence:
                utrs[gene] = sequence

    return utrs


def load_genes(genes_file: Path) -> set:
    """Load target gene symbols."""
    genes = set()

    with genes_file.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            gene = line.strip()

            if gene:
                genes.add(gene)

    return genes


def find_all_positions(sequence: str, motif: str) -> list:
    """Find all 1-based positions of a motif in a sequence."""
    positions = []
    start = 0

    while True:
        index = sequence.find(motif, start)

        if index == -1:
            break

        positions.append(index + 1)
        start = index + 1

    return positions


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run seed matching against sarcopenia/exercise-related 3'UTRs."
    )

    parser.add_argument(
        "--genes",
        required=True,
        help="Text file containing sarcopenia/exercise-related gene symbols.",
    )
    parser.add_argument(
        "--mirna",
        required=True,
        help="TSV file containing miRNA seed summary information.",
    )
    parser.add_argument(
        "--utr",
        required=True,
        help="TSV file containing gene_symbol and utr_3 columns.",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV file with seed matches.",
    )

    return parser.parse_args()


def main() -> None:
    """Run seed-based prediction against sarcopenia target 3'UTRs."""
    args = parse_args()

    genes_file = validate_file(args.genes, "Sarcopenia gene set")
    mirna_file = validate_file(args.mirna, "miRNA seed summary table")
    utr_file = validate_file(args.utr, "3'UTR table")
    out_file = prepare_output(args.out)

    mirnas = load_mirnas(mirna_file)
    utrs = load_utrs(utr_file)
    genes = load_genes(genes_file)

    output_rows = []

    for gene in sorted(genes):
        if gene not in utrs:
            continue

        sequence = utrs[gene]

        for mirna in mirnas:
            positions = find_all_positions(sequence, mirna["seed_rc"])

            if positions:
                output_rows.append({
                    "gene_symbol": gene,
                    "renamed_miRNA": mirna["renamed_miRNA"],
                    "source": mirna["source"],
                    "seed": mirna["seed"],
                    "seed_rc": mirna["seed_rc"],
                    "n_sites": len(positions),
                    "site_positions": ",".join(map(str, positions)),
                    "log2FoldChange": mirna["log2FoldChange"],
                    "direction": mirna["direction"],
                })

    with out_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print("miRNAs loaded:", len(mirnas))
    print("genes in sarcopenia set:", len(genes))
    print("genes with available 3'UTR:", len(utrs))
    print("matches found:", len(output_rows))
    print("output:", out_file)


if __name__ == "__main__":
    main()