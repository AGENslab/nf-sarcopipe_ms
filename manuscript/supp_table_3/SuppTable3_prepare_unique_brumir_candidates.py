#!/usr/bin/env python3

"""
SuppTable3_prepare_unique_brumir_candidates.py

Prepare Supplementary Table 3 from the complete
BrumiR2Reference candidate catalog.

Filtering:
1. Remove group == NA
2. Keep athlete, sedentary and shared candidates
3. Keep one row per candidate
4. If multiple loci exist, keep the locus with the
   lowest (most negative) BrumiR2Reference MFE.

Usage
-----
python3 SuppTable3_prepare_unique_brumir_candidates.py \
    --catalog denovo_passfilter_catalog.tsv \
    --out Supplementary_Table_3.tsv
"""

import argparse
import csv
import re
from pathlib import Path


parser = argparse.ArgumentParser()

parser.add_argument(
    "--catalog",
    required=True,
    help="denovo_passfilter_catalog.tsv"
)

parser.add_argument(
    "--out",
    required=True,
    help="Output TSV"
)

args = parser.parse_args()

catalog = Path(args.catalog)
outfile = Path(args.out)

best = {}

with catalog.open() as handle:

    reader = csv.DictReader(handle, delimiter="\t")

    for row in reader:

        if row["group"] == "NA":
            continue

        candidate = row["candidate"]
        mfe = float(row["brumir2ref_mfe"])

        if candidate not in best:
            best[candidate] = row
            continue

        old_mfe = float(best[candidate]["brumir2ref_mfe"])

        if mfe < old_mfe:

            best[candidate] = row

        elif mfe == old_mfe:

            old = best[candidate]

            tie_new = (
                row["chr"],
                int(row["start"]),
                int(row["stop"])
            )

            tie_old = (
                old["chr"],
                int(old["start"]),
                int(old["stop"])
            )

            if tie_new < tie_old:
                best[candidate] = row


def cluster_number(x):

    m = re.search(r"(\d+)$", x)

    if m:
        return int(m.group(1))

    return 999999


group_order = {
    "athlete": 0,
    "sedentary": 1,
    "shared": 2
}

rows = sorted(
    best.values(),
    key=lambda r: (
        group_order[r["group"]],
        cluster_number(r["candidate"])
    )
)

outfile.parent.mkdir(
    parents=True,
    exist_ok=True
)

columns = [
    "Candidate",
    "Group",
    "Chromosome",
    "Start",
    "End",
    "Mature sequence (5'->3')",
    "Mature length (nt)",
    "Precursor length (nt)",
    "BrumiR2Reference MFE (kcal/mol)",
    "Known exact match",
    "De novo exact",
]

with outfile.open(
    "w",
    newline=""
) as out:

    writer = csv.writer(
        out,
        delimiter="\t"
    )

    writer.writerow(columns)

    for r in rows:

        writer.writerow([
            r["candidate"],
            r["group"],
            r["chr"],
            r["start"],
            r["stop"],
            r["mature_seq"],
            r["mature_len"],
            r["precursor_len"],
            r["brumir2ref_mfe"],
            r["is_known_exact"],
            r["is_denovo_exact"],
        ])

known = sum(
    int(r["is_known_exact"])
    for r in rows
)

denovo = sum(
    int(r["is_denovo_exact"])
    for r in rows
)

print("------------------------------------------------")
print("Supplementary Table 3")
print("------------------------------------------------")
print(f"Unique candidates : {len(rows)}")
print(f"Known exact       : {known}")
print(f"De novo exact     : {denovo}")
print(f"Written           : {outfile}")
print("------------------------------------------------")
