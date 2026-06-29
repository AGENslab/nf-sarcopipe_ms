#!/usr/bin/env python3
# ============================================================
# extract_brumir_core220_fasta.py
# Description:
# Extracts the representative centroid sequences corresponding
# to the BrumiR-RF 0.95 core union cluster list and writes them
# to a FASTA file for downstream annotation against miRGeneDB.
#
# Inputs:
#   - /mnt/beegfs/home/npoblete/sarcopipe/results/clustering/0.95/brumir.0.95.core_union.txt
#   - /mnt/beegfs/home/npoblete/sarcopipe/results/clustering/0.95/all.candidates_clustered.cluster_named.fasta
#
# Output:
#   - overlap_miRGeneDB/brumir_core220_sequences.fa
# ============================================================

from pathlib import Path

core_ids = Path("/mnt/beegfs/home/npoblete/sarcopipe/results/clustering/0.95/brumir.0.95.core_union.txt")
cluster_fasta = Path("/mnt/beegfs/home/npoblete/sarcopipe/results/clustering/0.95/all.candidates_clustered.cluster_named.fasta")
out_fasta = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/overlap_miRGeneDB/brumir_core220_sequences.fa")

wanted = set(x.strip() for x in core_ids.read_text().splitlines() if x.strip())

header = None
seq = []

def flush(h, s, out):
    if h and h.split()[0] in wanted:
        out.write(f">{h}\n{s}\n")

with cluster_fasta.open() as f, out_fasta.open("w") as out:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                flush(header, "".join(seq), out)
            header = line[1:]
            seq = []
        else:
            seq.append(line)
    if header is not None:
        flush(header, "".join(seq), out)

print("Written:", out_fasta)
