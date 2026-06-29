#!/usr/bin/env python3
"""
brumir_cdhit_counts_matrix.py

Build a cluster-level BrumiR count matrix from:
1) CD-HIT cluster memberships (.clstr), and
2) per-sample BrumiR-RF high-confidence FASTA headers containing KC counts.

Purpose
-------
This script reconstructs cluster abundance after CD-HIT clustering by summing
KC values from BrumiR-RF FASTA headers for all members belonging to each cluster.

It also generates:
- a counts matrix (cluster x sample; KC sums)
- a binary presence/absence matrix
- group-specific core sets based on a minimum within-group presence fraction

Expected sample metadata CSV
----------------------------
Must contain at least:
- ID
- status

Typical statuses:
- active or athlete
- sedentary

Expected BrumiR-RF FASTA header example
---------------------------------------
>7712 RANK=1 KM=1371228 KC=12341052 LN=22

Expected CD-HIT member example in .clstr
----------------------------------------
>SRR11363967|7712... *

Outputs
-------
- counts TSV
- binary TSV
- core sets TSV
"""

import argparse
import csv
import math
import re
from collections import defaultdict, OrderedDict


def read_samples_csv(csv_path: str):
    """
    Read the sample metadata CSV.

    Expected header:
      ID, fastq_raw, status

    Returns
    -------
    statuses : dict
        sample_id -> status
    groups : dict
        status -> list(sample_id), preserving CSV order
    """
    statuses = {}
    groups = defaultdict(list)

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sid = row["ID"]
            status = row["status"]
            statuses[sid] = status
            groups[status].append(sid)

    return statuses, groups


def parse_rf_fasta_kc(fasta_path: str, sample_id: str):
    """
    Parse one BrumiR-RF high-confidence FASTA and extract KC counts.

    Example header:
      >7712 RANK=1 KM=1371228 KC=12341052 LN=22

    Parameters
    ----------
    fasta_path : str
        Path to one per-sample RF FASTA
    sample_id : str
        Sample ID inferred from filename

    Returns
    -------
    dict
        "SAMPLE|SEQID" -> KC (int)
    """
    kc_map = {}
    kc_re = re.compile(r"\bKC=(\d+)\b")
    seqid_re = re.compile(r"^>(\S+)")

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                m_id = seqid_re.match(line)
                if not m_id:
                    continue

                seqid = m_id.group(1)
                m_kc = kc_re.search(line)

                if not m_kc:
                    kc = 0
                else:
                    kc = int(m_kc.group(1))

                key = f"{sample_id}|{seqid}"
                kc_map[key] = kc

    return kc_map


def parse_rf_fastas_kc(rf_fastas):
    """
    Parse multiple BrumiR-RF FASTAs and collect KC values.

    FASTA filenames are expected to look like:
      <SAMPLE>.brumir_rf.high_confidence.fasta

    Sample ID is inferred from the basename prefix before the first dot.

    Returns
    -------
    kc_map : dict
        "SAMPLE|SEQID" -> KC
    sample_ids_found : list
        Sample IDs in discovery order
    """
    import os

    kc_map = {}
    sample_ids_found = []

    for fp in rf_fastas:
        base = os.path.basename(fp)
        sample_id = base.split(".", 1)[0]
        sample_ids_found.append(sample_id)
        kc_map.update(parse_rf_fasta_kc(fp, sample_id))

    return kc_map, sample_ids_found


def parse_clstr(clstr_path: str):
    """
    Parse a CD-HIT .clstr file.

    Returns
    -------
    list of tuples:
        (cluster_id, member_id, is_rep)

    Example extracted member_id:
      SRR11363967|7712
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
                member_id = m.group(1).strip()
                is_rep = line.endswith("*")
                entries.append((cluster_id, member_id, is_rep))

    return entries


def write_matrix_tsv(out_path: str, header_samples, feature_to_counts):
    """
    Write a matrix-style TSV with features as rows and samples as columns.
    """
    with open(out_path, "w", newline="") as out:
        out.write("miRNA\t" + "\t".join(header_samples) + "\n")
        for feat, counts in feature_to_counts.items():
            row = [feat] + [str(counts.get(s, 0)) for s in header_samples]
            out.write("\t".join(row) + "\n")


def write_core_sets(out_path: str, athlete_core, sedentary_core):
    """
    Write group-specific core feature sets as a two-column TSV.
    """
    with open(out_path, "w", newline="") as out:
        out.write("set\tmiRNA\n")
        for feat in athlete_core:
            out.write(f"athlete_core\t{feat}\n")
        for feat in sedentary_core:
            out.write(f"sedentary_core\t{feat}\n")


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Build a BrumiR cluster-level count matrix using KC values from "
            "BrumiR-RF FASTA headers and CD-HIT cluster assignments."
        )
    )
    ap.add_argument("--csv", required=True, help="Sample metadata CSV with columns ID, fastq_raw, status")
    ap.add_argument("--identity", required=True, help="Identity label for provenance (e.g. 0.85)")
    ap.add_argument("--clstr", required=True, help="CD-HIT .clstr file")
    ap.add_argument("--clustered_fasta", required=True, help="CD-HIT clustered FASTA (kept for provenance)")
    ap.add_argument("--rf_fastas", required=True, nargs="+", help="Per-sample BrumiR-RF high-confidence FASTAs")
    ap.add_argument("--min_presence_frac", type=float, default=0.8, help="Minimum within-group presence fraction for core sets")
    ap.add_argument("--out_counts", required=True, help="Output counts TSV")
    ap.add_argument("--out_binary", required=True, help="Output binary TSV")
    ap.add_argument("--out_core_sets", required=True, help="Output core sets TSV")
    args = ap.parse_args()

    statuses, groups = read_samples_csv(args.csv)

    kc_map, rf_samples = parse_rf_fastas_kc(args.rf_fastas)

    # Use CSV order for columns, keeping only samples present in RF FASTAs
    rf_sample_set = set(rf_samples)
    all_samples_ordered = [sid for sid in statuses.keys() if sid in rf_sample_set]

    athlete_samples = [sid for sid in groups.get("athlete", []) if sid in rf_sample_set]
    sedentary_samples = [sid for sid in groups.get("sedentary", []) if sid in rf_sample_set]

    if not all_samples_ordered:
        raise SystemExit("ERROR: No sample IDs from CSV were found among --rf_fastas filenames.")

    entries = parse_clstr(args.clstr)

    # cluster -> sample -> KC sum
    cluster_to_counts = defaultdict(lambda: defaultdict(int))
    cluster_to_binary = defaultdict(lambda: defaultdict(int))

    for cid, member_id, _is_rep in entries:
        if "|" not in member_id:
            continue

        sample_id = member_id.split("|", 1)[0]
        if sample_id not in rf_sample_set:
            continue

        kc = kc_map.get(member_id, 0)
        feat = f"cluster_{cid}"

        if kc > 0:
            cluster_to_counts[feat][sample_id] += kc
            cluster_to_binary[feat][sample_id] = 1

    features_sorted = sorted(cluster_to_counts.keys(), key=lambda x: int(x.split("_")[1]))
    counts_ordered = OrderedDict((f, cluster_to_counts[f]) for f in features_sorted)
    binary_ordered = OrderedDict((f, cluster_to_binary[f]) for f in features_sorted)

    write_matrix_tsv(args.out_counts, all_samples_ordered, counts_ordered)
    write_matrix_tsv(args.out_binary, all_samples_ordered, binary_ordered)

    def core_features(sample_list):
        """
        Return the set of features present in at least ceil(min_presence_frac * n)
        samples within a group.
        """
        if not sample_list:
            return set()

        n = len(sample_list)
        min_n = int(math.ceil(args.min_presence_frac * n))
        core = set()

        for feat, pres in binary_ordered.items():
            present_n = sum(1 for s in sample_list if pres.get(s, 0) == 1)
            if present_n >= min_n:
                core.add(feat)

        return core

    athlete_core = sorted(core_features(athlete_samples), key=lambda x: int(x.split("_")[1]))
    sedentary_core = sorted(core_features(sedentary_samples), key=lambda x: int(x.split("_")[1]))

    write_core_sets(args.out_core_sets, athlete_core, sedentary_core)


if __name__ == "__main__":
    main()
