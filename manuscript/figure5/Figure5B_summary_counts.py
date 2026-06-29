#!/usr/bin/env python3

"""
Figure 5B – Seed prediction summary counts

Description
-----------
This script summarizes seed-based miRNA–mRNA target prediction
results for BrumiR and miRDeep2. It counts evaluated miRNAs, genes
with seed matches, raw miRNA–mRNA pairs, and coherent inverse
miRNA–mRNA pairs.

Coherent pairs are defined as:
- miRNA Up_in_Active with gene Up_in_Sedentary
- miRNA Up_in_Sedentary with gene Up_in_Active

Inputs
------
--mirna_table
    TSV file containing miRNA seed summary information.

--gene_table
    CSV file containing DEG gene symbols and direction.

--brumir_match
    TSV file with BrumiR seed matches.

--mirdeep_match
    TSV file with miRDeep2 seed matches.

--out_summary
    Output TSV summary table.

--out_coherent
    Output TSV with coherent miRNA–mRNA pairs.

Outputs
-------
fig5a_summary.tsv
seed_matches_coherent.tsv

Usage
-----
python3 Figure5B_summary_counts.py \\
  --mirna_table all_miRNA_seed_summary.tsv \\
  --gene_table DEG_36_for_networks.csv \\
  --brumir_match seed_matches_36genes_BrumiR.tsv \\
  --mirdeep_match seed_matches_36genes_miRDeep2.tsv \\
  --out_summary fig5a_summary.tsv \\
  --out_coherent seed_matches_coherent.tsv
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


def read_mirnas(mirna_table: Path) -> dict:
    """Read evaluated miRNAs grouped by source."""
    mirnas = {"BrumiR": set(), "miRDeep2": set()}

    with mirna_table.open(encoding="utf-8", errors="ignore") as handle:
        next(handle)

        for line in handle:
            fields = line.rstrip("\n\r").split("\t")

            if len(fields) < 7:
                continue

            source = fields[0].strip()
            renamed = fields[2].strip()

            if source in mirnas and renamed:
                mirnas[source].add(renamed)

    return mirnas


def read_gene_direction(gene_table: Path) -> dict:
    """Read gene regulation direction by HGNC gene symbol."""
    gene_direction = {}

    with gene_table.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle)

        required_columns = {"gene_symbol", "direction"}
        missing_columns = required_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise SystemExit(
                "ERROR: Gene table missing required columns: "
                f"{', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            gene_direction[row["gene_symbol"].strip()] = row["direction"].strip()

    return gene_direction


def read_matches(path: Path) -> list:
    """Read seed match table as a list of dictionaries."""
    rows = []

    with path.open(encoding="utf-8", errors="ignore") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            rows.append({
                key.strip(): value.strip()
                for key, value in row.items()
            })

    return rows


def is_coherent(mirna_dir: str, gene_dir: str) -> bool:
    """Return True when miRNA and gene directions are inverse."""
    return (
        (mirna_dir == "Up_in_Active" and gene_dir == "Up_in_Sedentary") or
        (mirna_dir == "Up_in_Sedentary" and gene_dir == "Up_in_Active")
    )


def summarize(source: str, rows: list, gene_dir: dict, n_mirnas: int) -> tuple:
    """Summarize seed matches and extract coherent rows."""
    genes = set()
    coherent_rows = []

    for row in rows:
        gene = row["gene_symbol"]
        genes.add(gene)

        gene_direction = gene_dir.get(gene, "")

        if is_coherent(row["direction"], gene_direction):
            coherent_row = dict(row)
            coherent_row["gene_direction"] = gene_direction
            coherent_rows.append(coherent_row)

    summary = {
        "source": source,
        "n_miRNAs": n_mirnas,
        "n_genes": len(genes),
        "n_pairs": len(rows),
        "n_coherent": len(coherent_rows),
    }

    return summary, coherent_rows


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarize seed-based target prediction results for Figure 5B."
    )

    parser.add_argument(
        "--mirna_table",
        required=True,
        help="TSV file containing miRNA seed summary information.",
    )
    parser.add_argument(
        "--gene_table",
        required=True,
        help="CSV file containing DEG gene symbols and direction.",
    )
    parser.add_argument(
        "--brumir_match",
        required=True,
        help="TSV file with BrumiR seed matches.",
    )
    parser.add_argument(
        "--mirdeep_match",
        required=True,
        help="TSV file with miRDeep2 seed matches.",
    )
    parser.add_argument(
        "--out_summary",
        required=True,
        help="Output TSV summary table.",
    )
    parser.add_argument(
        "--out_coherent",
        required=True,
        help="Output TSV with coherent miRNA–mRNA pairs.",
    )

    return parser.parse_args()


def main() -> None:
    """Build Figure 5B summary tables."""
    args = parse_args()

    mirna_table = validate_file(args.mirna_table, "miRNA seed summary table")
    gene_table = validate_file(args.gene_table, "DEG gene table")
    brumir_match = validate_file(args.brumir_match, "BrumiR seed match table")
    mirdeep_match = validate_file(args.mirdeep_match, "miRDeep2 seed match table")

    out_summary = prepare_output(args.out_summary)
    out_coherent = prepare_output(args.out_coherent)

    mirnas = read_mirnas(mirna_table)
    gene_direction = read_gene_direction(gene_table)

    brumir_rows = read_matches(brumir_match)
    mirdeep_rows = read_matches(mirdeep_match)

    brumir_summary, brumir_coherent = summarize(
        "BrumiR",
        brumir_rows,
        gene_direction,
        len(mirnas["BrumiR"]),
    )

    mirdeep_summary, mirdeep_coherent = summarize(
        "miRDeep2",
        mirdeep_rows,
        gene_direction,
        len(mirnas["miRDeep2"]),
    )

    with out_summary.open("w", encoding="utf-8") as handle:
        handle.write("source\tn_miRNAs\tn_genes\tn_pairs\tn_coherent\n")

        for row in [brumir_summary, mirdeep_summary]:
            handle.write(
                f"{row['source']}\t"
                f"{row['n_miRNAs']}\t"
                f"{row['n_genes']}\t"
                f"{row['n_pairs']}\t"
                f"{row['n_coherent']}\n"
            )

    all_coherent = brumir_coherent + mirdeep_coherent

    if all_coherent:
        keys = list(all_coherent[0].keys())

        with out_coherent.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
            writer.writeheader()

            for row in all_coherent:
                writer.writerow(row)

    print("=== FIG 5B SUMMARY ===")
    print(brumir_summary)
    print(mirdeep_summary)
    print("Summary written:", out_summary)
    print("Coherent pairs written:", out_coherent)


if __name__ == "__main__":
    main()
