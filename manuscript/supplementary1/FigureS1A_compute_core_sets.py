#!/usr/bin/env python3
"""
cdhit_cluster_presence.py

Generate per-cluster sample presence tables from a CD-HIT .clstr file
and a FASTA containing the corresponding sequence identifiers.

Purpose
-------
This script summarizes how many distinct samples are represented in each
CD-HIT cluster and exports:

1) a detailed per-sequence table
2) a summary table counting how many clusters appear in exactly N samples

Expected FASTA header format
----------------------------
>SAMPLE|SEQID optional_description

Only the first token of the header is used as the sequence ID.

Expected CD-HIT member format in .clstr
---------------------------------------
>SRR11363967|95718... *

Outputs
-------
- Detailed CSV:
    status, identity, cluster_id, sequence_id, representative, sample, count, sequence
- Summary TSV:
    status, identity, count_in_n_samples, n_clusters
"""

import argparse
import csv
import re
from collections import defaultdict


def read_fasta_as_dict(fasta_path: str):
    """
    Read a FASTA file into a dictionary: sequence_id -> sequence.

    Assumes header format:
      >SAMPLE|SEQID optional_description

    Only the first token after '>' is stored as the record ID.
    """
    seqs = {}
    cur_id = None
    cur_seq = []

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())

        if cur_id is not None:
            seqs[cur_id] = "".join(cur_seq)

    return seqs


def parse_clstr(clstr_path: str):
    """
    Parse a CD-HIT .clstr file.

    Returns
    -------
    list of tuples:
        (cluster_id:int, seq_id:str, is_rep:bool)

    Example line:
      0  25nt, >SRR11363967|95718... *
    """
    entries = []
    cluster_id = None

    pat = re.compile(r">(.+?)\.\.\.")

    with open(clstr_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">Cluster"):
                cluster_id = int(line.split()[-1])
                continue

            m = pat.search(line)
            if m and cluster_id is not None:
                sid = m.group(1).strip()
                is_rep = line.endswith("*")
                entries.append((cluster_id, sid, is_rep))

    return entries


def split_sample(seq_id: str):
    """
    Extract the sample prefix from a sequence ID.

    Expected format:
      SAMPLE|SEQID

    If the separator is absent, returns 'NA'.
    """
    if "|" in seq_id:
        sample = seq_id.split("|", 1)[0]
    else:
        sample = "NA"
    return sample


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Generate per-cluster sample presence from a CD-HIT .clstr file "
            "and a FASTA with matching sequence identifiers."
        )
    )
    ap.add_argument("--status", required=True, help="Run label, e.g. athlete / sedentary / all")
    ap.add_argument("--identity", required=True, help="Identity label, e.g. 0.85")
    ap.add_argument("--clstr", required=True, help="Path to CD-HIT .clstr")
    ap.add_argument(
        "--fasta",
        required=True,
        help="FASTA containing the same sequence IDs referenced in the .clstr"
    )
    ap.add_argument("--out_csv", required=True, help="Output detailed CSV path")
    ap.add_argument("--out_summary", required=True, help="Output summary TSV path")
    args = ap.parse_args()

    fasta_seqs = read_fasta_as_dict(args.fasta)
    entries = parse_clstr(args.clstr)

    # cluster -> set of distinct samples represented in that cluster
    cluster_samples = defaultdict(set)
    for cid, sid, _ in entries:
        cluster_samples[cid].add(split_sample(sid))

    # cluster -> number of unique samples
    cluster_count = {cid: len(sample_set) for cid, sample_set in cluster_samples.items()}

    # Detailed per-sequence output
    with open(args.out_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow([
            "status",
            "identity",
            "cluster_id",
            "sequence_id",
            "representative",
            "sample",
            "count",
            "sequence",
        ])

        for cid, sid, is_rep in entries:
            sample = split_sample(sid)
            seq = fasta_seqs.get(sid, "")
            writer.writerow([
                args.status,
                args.identity,
                cid,
                sid,
                int(is_rep),
                sample,
                cluster_count.get(cid, 0),
                seq,
            ])

    # Summary: how many clusters appear in exactly N samples
    freq = defaultdict(int)
    for _, cnt in cluster_count.items():
        freq[cnt] += 1

    with open(args.out_summary, "w", newline="") as out:
        out.write("status\tidentity\tcount_in_n_samples\tn_clusters\n")
        for cnt in sorted(freq.keys()):
            out.write(f"{args.status}\t{args.identity}\t{cnt}\t{freq[cnt]}\n")


if __name__ == "__main__":
    main()
    