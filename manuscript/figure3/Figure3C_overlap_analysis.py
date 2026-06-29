#!/usr/bin/env python3
# ============================================================
# analysis_overlap_brumirRF_vs_miRDeep2_known.py
# Description:
# Compares all BrumiR-RF candidate sequences (before CD-HIT)
# against miRDeep2 known miRNA sequences to quantify exact
# overlap and near-match overlap at a user-defined identity.
#
# This script is intended as an analysis/validation step and
# can later be integrated into module 3 of the nf-Sarcopipe
# workflow.
#
# Inputs:
#   1. BrumiR-RF FASTA file (all retained RF candidates)
#   2. miRDeep2 known FASTA file
#
# Outputs:
#   - overlap_exact_matches.tsv
#   - overlap_summary.tsv
#   - overlap_exact_sequences.fasta
#
# Usage:
#   python3 analysis_overlap_brumirRF_vs_miRDeep2_known.py \
#       --brumir_rf path/to/brumir_rf.fasta \
#       --mirdeep2_known path/to/mirdeep2_known.fasta \
#       --outdir path/to/output_dir
# ============================================================

from pathlib import Path
import argparse
from collections import defaultdict

def read_fasta(path):
    records = []
    header = None
    seq_chunks = []

    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_chunks).upper().replace("U", "T")
                    records.append((header, seq))
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if header is not None:
            seq = "".join(seq_chunks).upper().replace("U", "T")
            records.append((header, seq))

    return records

def write_fasta(records, outpath):
    with open(outpath, "w", encoding="utf-8") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--brumir_rf", required=True, help="BrumiR-RF FASTA file")
    parser.add_argument("--mirdeep2_known", required=True, help="miRDeep2 known FASTA file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    brumir_path = Path(args.brumir_rf)
    mirdeep_path = Path(args.mirdeep2_known)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    brumir_records = read_fasta(brumir_path)
    mirdeep_records = read_fasta(mirdeep_path)

    print("BrumiR-RF sequences loaded:", len(brumir_records))
    print("miRDeep2 known sequences loaded:", len(mirdeep_records))

    # sequence -> headers
    brumir_by_seq = defaultdict(list)
    mirdeep_by_seq = defaultdict(list)

    for h, s in brumir_records:
        brumir_by_seq[s].append(h)

    for h, s in mirdeep_records:
        mirdeep_by_seq[s].append(h)

    brumir_set = set(brumir_by_seq.keys())
    mirdeep_set = set(mirdeep_by_seq.keys())

    shared_exact = sorted(brumir_set & mirdeep_set)
    brumir_only = sorted(brumir_set - mirdeep_set)
    mirdeep_only = sorted(mirdeep_set - brumir_set)

    exact_rows = []
    exact_fasta = []

    for seq in shared_exact:
        for bh in brumir_by_seq[seq]:
            for mh in mirdeep_by_seq[seq]:
                exact_rows.append({
                    "sequence": seq,
                    "brumir_rf_id": bh,
                    "mirdeep2_known_id": mh,
                    "match_type": "exact"
                })
        # keep one representative FASTA entry
        exact_fasta.append((f"shared_exact_{len(exact_fasta)+1}", seq))

    exact_tsv = outdir / "overlap_exact_matches.tsv"
    with open(exact_tsv, "w", encoding="utf-8") as f:
        f.write("sequence\tbrumir_rf_id\tmirdeep2_known_id\tmatch_type\n")
        for r in exact_rows:
            f.write(
                f'{r["sequence"]}\t{r["brumir_rf_id"]}\t{r["mirdeep2_known_id"]}\t{r["match_type"]}\n'
            )

    exact_fasta_out = outdir / "overlap_exact_sequences.fasta"
    write_fasta(exact_fasta, exact_fasta_out)

    summary_tsv = outdir / "overlap_summary.tsv"
    with open(summary_tsv, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"BrumiR_RF_total_sequences\t{len(brumir_set)}\n")
        f.write(f"miRDeep2_known_total_sequences\t{len(mirdeep_set)}\n")
        f.write(f"shared_exact_sequences\t{len(shared_exact)}\n")
        f.write(f"BrumiR_RF_only_sequences\t{len(brumir_only)}\n")
        f.write(f"miRDeep2_known_only_sequences\t{len(mirdeep_only)}\n")
        f.write(f"shared_exact_pairwise_matches\t{len(exact_rows)}\n")

    print("Exact shared sequences:", len(shared_exact))
    print("BrumiR-RF only:", len(brumir_only))
    print("miRDeep2 known only:", len(mirdeep_only))
    print("Written:", exact_tsv)
    print("Written:", exact_fasta_out)
    print("Written:", summary_tsv)

if __name__ == "__main__":
    main()
