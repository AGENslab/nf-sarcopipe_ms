#!/usr/bin/env python3
"""
denovo_passfilter_catalog.py

Description
-----------
This script builds a catalog of BrumiR core candidates that are supported by
BrumiR2Reference passfilter results, and classifies them as known or de novo
based on exact mature-sequence matching against miRBase and MirGeneDB.

It produces:
1. A candidate-level catalog of supported mature sequences found inside
   BrumiR2Reference precursor sequences
2. A summary table of unique supported candidates
3. A group-level count table for downstream plotting

Inputs
------
--brumir_core_fasta   FASTA with BrumiR core mature sequences
--passfilter_tsv      BrumiR2Reference passfilter table
--core_sets_tsv       Core sets TSV from Module II
--mirbase_mature      miRBase mature FASTA
--mirgenedb_mature    MirGeneDB mature FASTA
--out_catalog         Output TSV with supported candidate catalog
--out_summary         Output TSV with summary metrics
--out_counts          Output TSV with group-level counts

Outputs
-------
1. denovo_passfilter_catalog.tsv
2. denovo_passfilter_summary.tsv
3. status_counts_by_group.tsv
"""

import argparse
import re
import sys
from pathlib import Path


def read_fasta(path):
    """
    Read FASTA into a dict: header_id -> sequence (RNA alphabet, uppercase).
    """
    seqs = {}
    name, buf = None, []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper().replace("T", "U")
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(re.sub(r"\s+", "", line))

        if name is not None:
            seqs[name] = "".join(buf).upper().replace("T", "U")

    return seqs


def fasta_values(path):
    """
    Return all sequence values from a FASTA file, or an empty list if missing.
    """
    if not path or not Path(path).exists():
        return []
    return list(read_fasta(path).values())


def load_known(mirbase_path, mirgenedb_path):
    """
    Load known mature sequences from miRBase and MirGeneDB into one set.
    """
    known = set()

    for s in fasta_values(mirbase_path):
        if s:
            known.add(s.upper().replace("T", "U"))

    for s in fasta_values(mirgenedb_path):
        if s:
            known.add(s.upper().replace("T", "U"))

    return known


def parse_core_sets(path):
    """
    Read core sets TSV into:
      candidate_id -> set(core_set_labels)

    Expected format:
      set<TAB>miRNA
    """
    mapping = {}

    with open(path) as f:
        _header = f.readline()

        for line in f:
            if not line.strip():
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue

            core_set, cid = parts[0], parts[1]
            mapping.setdefault(cid, set()).add(core_set)

    return mapping


def core_class(sets):
    """
    Convert core-set labels into a simplified group label.
    """
    if not sets:
        return "NA"

    s = set(sets)
    a = "athlete_core" in s
    d = "sedentary_core" in s

    if a and d:
        return "shared"
    if a:
        return "athlete"
    if d:
        return "sedentary"
    if "core_union" in s:
        return "core_union"
    return "other"


def main():
    ap = argparse.ArgumentParser(
        description="Build a catalog of BrumiR passfilter-supported core candidates and classify them as known or de novo."
    )
    ap.add_argument("--brumir_core_fasta", required=True, help="BrumiR core mature FASTA")
    ap.add_argument("--passfilter_tsv", required=True, help="BrumiR2Reference passfilter TSV")
    ap.add_argument("--core_sets_tsv", required=True, help="BrumiR core sets TSV")
    ap.add_argument("--mirbase_mature", required=True, help="miRBase mature FASTA")
    ap.add_argument("--mirgenedb_mature", required=True, help="MirGeneDB mature FASTA")
    ap.add_argument("--out_catalog", default="denovo_passfilter_catalog.tsv", help="Output catalog TSV")
    ap.add_argument("--out_summary", default="denovo_passfilter_summary.tsv", help="Output summary TSV")
    ap.add_argument("--out_counts", default="status_counts_by_group.tsv", help="Output count TSV")
    args = ap.parse_args()

    core = read_fasta(args.brumir_core_fasta)       # candidate -> mature sequence
    known = load_known(args.mirbase_mature, args.mirgenedb_mature)
    core_sets = parse_core_sets(args.core_sets_tsv) # candidate -> set(core_set)

    # Read passfilter table
    with open(args.passfilter_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        cols = {c: i for i, c in enumerate(header)}

        required = ["miRNA", "chr", "start", "stop", "MFE", "Precursor_Seq"]
        for k in required:
            if k not in cols:
                raise SystemExit(f"[ERROR] passfilter missing column '{k}'. Found: {header}")

        pass_rows = []
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            pass_rows.append({
                "passfilter_miRNA": parts[cols["miRNA"]],
                "chr": parts[cols["chr"]],
                "start": parts[cols["start"]],
                "stop": parts[cols["stop"]],
                "brumir2ref_mfe": parts[cols["MFE"]],
                "precursor_seq": parts[cols["Precursor_Seq"]].upper().replace("T", "U"),
            })

    # Index BrumiR core mature sequences by length for substring matching
    by_len = {}
    for cid, mature in core.items():
        by_len.setdefault(len(mature), []).append((cid, mature))

    catalog = []
    supported_candidates = set()  # (candidate, mature_seq)

    for row in pass_rows:
        precursor = row["precursor_seq"]
        hit = None

        # Search mature sequence inside precursor sequence
        for L in sorted(by_len.keys()):
            for cid, mature in by_len[L]:
                pos = precursor.find(mature)
                if pos != -1:
                    hit = (cid, mature, pos)
                    break
            if hit:
                break

        if not hit:
            continue

        cid, mature, pos0 = hit
        supported_candidates.add((cid, mature))

        sets = core_sets.get(cid, set())
        cclass = core_class(sets)
        is_known = "1" if mature in known else "0"
        is_denovo = "1" if mature not in known else "0"

        catalog.append({
            "candidate": cid,
            "group": cclass,
            "core_sets": ",".join(sorted(sets)) if sets else "NA",
            "passfilter_miRNA": row["passfilter_miRNA"],
            "chr": row["chr"],
            "start": row["start"],
            "stop": row["stop"],
            "brumir2ref_mfe": row["brumir2ref_mfe"],
            "precursor_len": str(len(precursor)),
            "mature_len": str(len(mature)),
            "mature_start": str(pos0 + 1),
            "mature_end": str(pos0 + len(mature)),
            "is_known_exact": is_known,
            "is_denovo_exact": is_denovo,
            "mature_seq": mature,
            "precursor_seq": precursor,
        })

    def all_core_in_group(group):
        return {
            (cid, mature)
            for cid, mature in core.items()
            if core_class(core_sets.get(cid, set())) == group
        }

    groups = ["athlete", "sedentary", "shared"]
    out_counts = []
    supported_set = set(supported_candidates)

    for g in groups:
        core_g = all_core_in_group(g)
        total = len(core_g)
        supported = len(core_g & supported_set)
        not_supported = total - supported
        denovo_core = {(cid, m) for (cid, m) in core_g if m not in known}
        denovo_supported = len(denovo_core & supported_set)

        out_counts.append((g, "core_total", total))
        out_counts.append((g, "supported", supported))
        out_counts.append((g, "not_supported", not_supported))
        out_counts.append((g, "denovo_supported_exact", denovo_supported))

    total_supported = len(supported_set)
    denovo_supported = len({k for k in supported_set if k[1] not in known})
    known_supported = total_supported - denovo_supported

    # Write catalog
    cols = [
        "candidate", "group", "core_sets", "passfilter_miRNA", "chr", "start", "stop",
        "brumir2ref_mfe", "precursor_len", "mature_len", "mature_start", "mature_end",
        "is_known_exact", "is_denovo_exact", "mature_seq", "precursor_seq"
    ]

    with open(args.out_catalog, "w") as w:
        w.write("\t".join(cols) + "\n")
        for row in catalog:
            w.write("\t".join(row[c] for c in cols) + "\n")

    with open(args.out_summary, "w") as w:
        w.write("metric\tvalue\n")
        w.write(f"supported_unique_candidates\t{total_supported}\n")
        w.write(f"known_exact_supported\t{known_supported}\n")
        w.write(f"denovo_exact_supported\t{denovo_supported}\n")

    with open(args.out_counts, "w") as w:
        w.write("group\tmetric\tcount\n")
        for g, m, c in out_counts:
            w.write(f"{g}\t{m}\t{c}\n")

    print("=== SUMMARY (core ∩ passfilter) ===", file=sys.stderr)
    print(f"Supported unique candidates: {total_supported}", file=sys.stderr)
    print(f"Known exact supported (miRBase ∪ MirGeneDB): {known_supported}", file=sys.stderr)
    print(f"De novo exact supported: {denovo_supported}", file=sys.stderr)
    print(f"[OK] Wrote: {args.out_catalog}, {args.out_summary}, {args.out_counts}", file=sys.stderr)


if __name__ == "__main__":
    main()
