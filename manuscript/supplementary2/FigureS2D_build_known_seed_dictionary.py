#!/usr/bin/env python3
# ============================================================
# build_known_seed_dictionary.py
# Description:
# Extracts canonical 7-mer seeds (positions 2-8) from known
# human miRNAs and builds a seed reference dictionary.
#
# Input:
#   - /mnt/beegfs/home/npoblete/databases/miRBase/hsa/hsa_mature.mirdeep2.OK.fa
#
# Output:
#   - kmer_seed_analysis/known_seed_dictionary.tsv
# ============================================================

from pathlib import Path

fasta = Path("/mnt/beegfs/home/npoblete/databases/miRBase/hsa/hsa_mature.mirdeep2.OK.fa")
out = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/kmer_seed_analysis/known_seed_dictionary.tsv")

header = None
seq = []

def family_from_id(mirna_id):
    # simple family approximation for summary
    x = mirna_id.replace("hsa-", "")
    return x

rows = []

with fasta.open() as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                s = "".join(seq).upper().replace("T", "U")
                if len(s) >= 8:
                    seed = s[1:8]   # canonical 2–8
                    rows.append((header, family_from_id(header), s, seed))
            header = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)
    if header is not None:
        s = "".join(seq).upper().replace("T", "U")
        if len(s) >= 8:
            seed = s[1:8]
            rows.append((header, family_from_id(header), s, seed))

with out.open("w") as o:
    o.write("miRNA_id\tfamily\tmature_seq\tcanonical_seed\n")
    for r in rows:
        o.write("\t".join(r) + "\n")

print("Written:", out)
print("Total known miRNAs:", len(rows))
