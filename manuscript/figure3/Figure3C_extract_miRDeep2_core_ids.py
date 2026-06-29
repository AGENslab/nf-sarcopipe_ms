#!/usr/bin/env python3
# ============================================================
# extract_mirdeep2_core941_ids.py
# Description:
# Extracts a unique list of miRDeep2 known miRNA IDs from the
# p0.8 core sets table, removing duplicated shared entries.
#
# Input:
#   - /mnt/beegfs/home/npoblete/sarcopipe/results/miRDeep2_known/miRDeep2.known.p08.core_sets.tsv
#
# Output:
#   - overlap_miRGeneDB/md_p08_unique_ids.txt
# ============================================================

from pathlib import Path

infile = Path("/mnt/beegfs/home/npoblete/sarcopipe/results/miRDeep2_known/miRDeep2.known.p08.core_sets.tsv")
outfile = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/overlap_miRGeneDB/md_p08_unique_ids.txt")

ids = set()
with infile.open() as f:
    next(f)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts and parts[0].strip():
            ids.add(parts[0].strip())

with outfile.open("w") as out:
    for x in sorted(ids):
        out.write(x + "\n")

print("Unique IDs:", len(ids))
print("Written:", outfile)
