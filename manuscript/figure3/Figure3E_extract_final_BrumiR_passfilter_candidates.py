#!/usr/bin/env python3

"""
Figure3E_extract_final_BrumiR_passfilter_candidates.py

Extract the final differentially expressed BrumiR-RF candidates
from the BrumiR2Reference passfilter catalog.

Selection criteria
------------------
1. Candidate belongs to an athlete, sedentary, or shared core set.
2. Candidate is classified as either:
   - de novo exact, or
   - known exact.
3. Candidate is differentially expressed at padj < threshold.
4. One row is retained per candidate.
5. If multiple precursor loci exist, the locus with the most
   negative BrumiR2Reference MFE is selected.

The output contains the complete selected passfilter record plus
DESeq2 statistics.

Usage
-----
python3 Figure3E_extract_final_BrumiR_passfilter_candidates.py \
  --catalog denovo_passfilter_catalog.tsv \
  --deseq2 brumir_deseq2_results.csv \
  --out Figure3E_final_BrumiR_passfilter_candidates.tsv
"""

import argparse
import csv
from pathlib import Path


VALID_GROUPS = {"athlete", "sedentary", "shared"}


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--catalog",
        required=True,
        help="denovo_passfilter_catalog.tsv",
    )

    parser.add_argument(
        "--deseq2",
        required=True,
        help="BrumiR DESeq2 results CSV",
    )

    parser.add_argument(
        "--out",
        required=True,
        help="Output TSV",
    )

    parser.add_argument(
        "--padj",
        type=float,
        default=0.05,
        help="Adjusted p-value threshold (default: 0.05)",
    )

    return parser.parse_args()


def read_significant_deseq2(path, threshold):
    significant = {}

    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle)

        required = {
            "feature",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise SystemExit(
                "ERROR: DESeq2 table is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            try:
                padj = float(row["padj"])
                log2fc = float(row["log2FoldChange"])
            except (TypeError, ValueError):
                continue

            if padj < threshold:
                row["regulation"] = (
                    "Up in Active"
                    if log2fc > 0
                    else "Up in Sedentary"
                )

                significant[row["feature"].strip()] = row

    return significant


def locus_priority(row):
    """
    Prefer the locus with the most negative MFE.

    Ties are resolved deterministically using chromosome,
    start, stop, and precursor length.
    """
    return (
        float(row["brumir2ref_mfe"]),
        row["chr"],
        int(row["start"]),
        int(row["stop"]),
        int(row["precursor_len"]),
    )


def main():
    args = parse_args()

    catalog_path = Path(args.catalog)
    deseq2_path = Path(args.deseq2)
    output_path = Path(args.out)

    for path in (catalog_path, deseq2_path):
        if not path.is_file():
            raise SystemExit(f"ERROR: input file not found: {path}")

    significant = read_significant_deseq2(
        deseq2_path,
        args.padj,
    )

    best_by_candidate = {}
    catalog_columns = None

    with catalog_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        catalog_columns = list(reader.fieldnames or [])

        required = {
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
        }

        missing = required - set(catalog_columns)

        if missing:
            raise SystemExit(
                "ERROR: catalog is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            candidate = row["candidate"].strip()
            group = row["group"].strip().lower()

            if group not in VALID_GROUPS:
                continue

            if candidate not in significant:
                continue

            known = row["is_known_exact"].strip() == "1"
            denovo = row["is_denovo_exact"].strip() == "1"

            if not known and not denovo:
                continue

            if (
                candidate not in best_by_candidate
                or locus_priority(row)
                < locus_priority(best_by_candidate[candidate])
            ):
                best_by_candidate[candidate] = row

    rows = []

    for candidate, row in best_by_candidate.items():
        de = significant[candidate]

        output_row = dict(row)

        output_row["candidate_type"] = (
            "known_exact"
            if row["is_known_exact"] == "1"
            else "de_novo"
        )

        output_row["baseMean"] = de["baseMean"]
        output_row["log2FoldChange"] = de["log2FoldChange"]
        output_row["lfcSE"] = de["lfcSE"]
        output_row["stat"] = de["stat"]
        output_row["pvalue"] = de["pvalue"]
        output_row["padj"] = de["padj"]
        output_row["regulation"] = de["regulation"]

        rows.append(output_row)

    group_order = {
        "athlete": 0,
        "shared": 1,
        "sedentary": 2,
    }

    rows.sort(
        key=lambda row: (
            group_order.get(row["group"], 99),
            row["candidate_type"],
            int(row["candidate"].split("_")[-1]),
        )
    )

    output_columns = catalog_columns + [
        "candidate_type",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj",
        "regulation",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=output_columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)

    n_denovo = sum(
        row["candidate_type"] == "de_novo"
        for row in rows
    )

    n_known = sum(
        row["candidate_type"] == "known_exact"
        for row in rows
    )

    print("Figure 3E BrumiR-RF final passfilter catalog")
    print(f"Total candidates: {len(rows)}")
    print(f"De novo: {n_denovo}")
    print(f"Known exact: {n_known}")
    print(f"Written: {output_path}")

    if n_denovo != 7 or n_known != 1:
        raise SystemExit(
            "ERROR: expected 7 de novo and 1 known exact, "
            f"but obtained {n_denovo} de novo and {n_known} known."
        )


if __name__ == "__main__":
    main()
