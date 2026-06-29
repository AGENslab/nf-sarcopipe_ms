#!/usr/bin/env python3
# ============================================================
# extract_mirdeep2_core941_fasta_from_miRGeneDB.py
# Description:
# Extracts mature sequences from the human miRGeneDB FASTA for
# the unique miRDeep2 p0.8 core miRNA IDs.
#
# Inputs:
#   - overlap_miRGeneDB/md_p08_unique_ids.txt
#   - /mnt/beegfs/home/npoblete/databases/miRBase/hsa/mirgenedb_hsa_mature.fa
#
# Output:
#   - overlap_miRGeneDB/md_p08_unique_sequences.fa
# ============================================================

from pathlib import Path

ids_file = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/overlap_miRGeneDB/md_p08_unique_ids.txt")
ref_fasta = Path("/mnt/beegfs/home/npoblete/databases/miRBase/hsa/mirgenedb_hsa_mature.fa")
out_fasta = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/overlap_miRGeneDB/md_p08_unique_sequences.fa")

wanted = set(x.strip() for x in ids_file.read_text().splitlines() if x.strip())

header = None
seq = []
keep = False
written = 0

with ref_fasta.open() as f, out_fasta.open("w") as out:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None and keep:
                out.write(f">{header}\n{''.join(seq)}\n")
                written += 1
            header = line[1:].split()[0]
            keep = header in wanted
            seq = []
        else:
            seq.append(line)
    if header is not None and keep:
        out.write(f">{header}\n{''.join(seq)}\n")
        written += 1

print("Requested IDs:", len(wanted))
print("Written sequences:", written)
print("Output:", out_fasta)
