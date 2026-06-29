#!/usr/bin/env python3

"""
Figure 5D – Prepare union of coherent target genes

Description
-----------
This script extracts the union of coherent target genes from the
seed-based miRNA–mRNA prediction table. It also writes separate
gene lists for BrumiR and miRDeep2 coherent target genes.

Inputs
------
--coherent
    TSV file containing coherent miRNA–mRNA seed-match pairs.

--out_union
    Output text file with the union of coherent target genes.

--out_brumir
    Output text file with BrumiR coherent target genes.

--out_mirdeep2
    Output text file with miRDeep2 coherent target genes.

Outputs
-------
coherent_target_genes_union.txt
coherent_target_genes_BrumiR.txt
coherent_target_genes_miRDeep2.txt

Usage
-----
python3 Figure5D_prepare_union_targets.py \\
  --coherent seed_matches_coherent.tsv \\
  --out_union coherent_target_genes_union.txt \\
  --out_brumir coherent_target_genes_BrumiR.txt \\
  --out_mirdeep2 coherent_target_genes_miRDeep2.txt
"""

import argparse
import csv
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


def write_gene_list(path: Path, genes: set[str]) -> None:
    """Write sorted gene symbols, one per line."""
    with path.open("w") as handle:
        for gene in sorted(genes):
            handle.write(f"{gene}\n")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract BrumiR, miRDeep2, and union coherent target gene lists."
    )

    parser.add_argument(
        "--coherent",
        required=True,
        help="TSV file containing coherent miRNA–mRNA seed-match pairs.",
    )
    parser.add_argument(
        "--out_union",
        required=True,
        help="Output text file with union coherent target genes.",
    )
    parser.add_argument(
        "--out_brumir",
        required=True,
        help="Output text file with BrumiR coherent target genes.",
    )
    parser.add_argument(
        "--out_mirdeep2",
        required=True,
        help="Output text file with miRDeep2 coherent target genes.",
    )

    return parser.parse_args()


def main() -> None:
    """Extract coherent target genes by source and union."""
    args = parse_args()

    coherent_file = validate_file(args.coherent, "Coherent seed-match table")
    out_union = prepare_output(args.out_union)
    out_brumir = prepare_output(args.out_brumir)
    out_mirdeep2 = prepare_output(args.out_mirdeep2)

    brumir_genes = set()
    mirdeep2_genes = set()
    union_genes = set()

    with coherent_file.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {"source", "gene_symbol"}
        missing_columns = required_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise SystemExit(
                "ERROR: Coherent table missing required columns: "
                f"{', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            source = row["source"].strip()
            gene = row["gene_symbol"].strip()

            union_genes.add(gene)

            if source == "BrumiR":
                brumir_genes.add(gene)
            elif source == "miRDeep2":
                mirdeep2_genes.add(gene)

    write_gene_list(out_brumir, brumir_genes)
    write_gene_list(out_mirdeep2, mirdeep2_genes)
    write_gene_list(out_union, union_genes)

    print("BrumiR coherent genes:", len(brumir_genes))
    print("miRDeep2 coherent genes:", len(mirdeep2_genes))
    print("Union coherent genes:", len(union_genes))
    print("Written:", out_union)


if __name__ == "__main__":
    main()
