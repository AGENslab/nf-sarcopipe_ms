#!/usr/bin/env python3
# ============================================================
# build_shared_miRGeneDB_table.py
# Description:
# Builds a supplementary table of shared miRNAs between
# BrumiR and miRDeep2 based on miRGeneDB annotation.
#
# Inputs:
#   - brumir_core220_vs_miRGeneDB.tsv
#   - md_core_vs_miRGeneDB.tsv
#   - shared_miRGeneDB_hits.txt
#
# Output:
#   - Supplementary_Table_S4_shared_miRNAs.tsv
# ============================================================

from pathlib import Path
import pandas as pd

base = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/overlap_miRGeneDB")

brumir_file = base / "brumir_core220_vs_miRGeneDB.tsv"
mirdeep_file = base / "md_core_vs_miRGeneDB.tsv"
shared_file = base / "shared_miRGeneDB_hits.txt"

out_file = base / "Supplementary_Table_S4_shared_miRNAs.tsv"

# Load shared miRGeneDB IDs
shared_ids = set(shared_file.read_text().splitlines())

# Load BLAST outputs
cols = ["query", "subject", "identity", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

br = pd.read_csv(brumir_file, sep="\t", names=cols)
md = pd.read_csv(mirdeep_file, sep="\t", names=cols)

# Filter only shared miRGeneDB hits
br_shared = br[br["subject"].isin(shared_ids)].copy()
md_shared = md[md["subject"].isin(shared_ids)].copy()

# Keep best hit per subject
br_shared = br_shared.sort_values("bitscore", ascending=False).drop_duplicates("subject")
md_shared = md_shared.sort_values("bitscore", ascending=False).drop_duplicates("subject")

# Merge both datasets
merged = pd.merge(
    br_shared,
    md_shared,
    on="subject",
    suffixes=("_BrumiR", "_miRDeep2")
)

# Select relevant columns
final = merged[[
    "subject",
    "query_BrumiR",
    "identity_BrumiR",
    "length_BrumiR",
    "query_miRDeep2",
    "identity_miRDeep2",
    "length_miRDeep2"
]]

final.columns = [
    "miRGeneDB_ID",
    "BrumiR_cluster",
    "BrumiR_identity",
    "BrumiR_alignment_length",
    "miRDeep2_miRNA",
    "miRDeep2_identity",
    "miRDeep2_alignment_length"
]

final.to_csv(out_file, sep="\t", index=False)

print("Written:", out_file)
print("Total shared miRNAs:", len(final))
