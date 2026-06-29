#!/usr/bin/env python3
# ============================================================
# summarize_miRGeneDB_overlap.py
# Description:
# Builds a compact summary table for the exact-match overlap
# against human miRGeneDB for:
#   - BrumiR core220 sequences
#   - miRDeep2 core941 sequences
#
# Inputs:
#   - overlap_miRGeneDB/brumir_core220_sequences.fa
#   - overlap_miRGeneDB/md_p08_unique_sequences.fa
#   - overlap_miRGeneDB/brumir_core220_vs_miRGeneDB.tsv
#   - overlap_miRGeneDB/md_core_vs_miRGeneDB.tsv
#
# Output:
#   - overlap_miRGeneDB/overlap_miRGeneDB_summary.tsv
# ============================================================

from pathlib import Path

base = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/overlap_miRGeneDB")

def count_fasta_headers(path: Path) -> int:
    return sum(1 for line in path.open() if line.startswith(">"))

def count_blast_queries(path: Path) -> int:
    hits = set()
    if path.exists():
        with path.open() as f:
            for line in f:
                line = line.strip()
                if line:
                    hits.add(line.split("\t")[0])
    return len(hits)

br_fa = base / "brumir_core220_sequences.fa"
md_fa = base / "md_p08_unique_sequences.fa"

br_blast = base / "brumir_core220_vs_miRGeneDB.tsv"
md_blast = base / "md_core_vs_miRGeneDB.tsv"

out = base / "overlap_miRGeneDB_summary.tsv"

br_total = count_fasta_headers(br_fa)
md_total = count_fasta_headers(md_fa)

br_match = count_blast_queries(br_blast)
md_match = count_blast_queries(md_blast)

with out.open("w") as f:
    f.write("algorithm\ttotal\tmatched_100pct\tunmatched_putative_novel\n")
    f.write(f"BrumiR_core220\t{br_total}\t{br_match}\t{br_total - br_match}\n")
    f.write(f"miRDeep2_core941\t{md_total}\t{md_match}\t{md_total - md_match}\n")

print("Written:", out)
print(f"BrumiR_core220: total={br_total}, matched_100pct={br_match}, unmatched={br_total - br_match}")
print(f"miRDeep2_core941: total={md_total}, matched_100pct={md_match}, unmatched={md_total - md_match}")
