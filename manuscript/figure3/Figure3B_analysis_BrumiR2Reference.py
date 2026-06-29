#!/usr/bin/env python3

"""
analysis_normalize_brumir2ref.py

Description
-----------
This script integrates CD-HIT cluster membership, BrumiR2Reference mapping support,
and sample metadata to generate cluster-level support summaries for downstream analysis.

It computes:
- number of passfilter, nonpassfilter, and unmapped members per cluster
- sample and group presence per cluster
- optional core membership flags using p0.8 and/or p1.0 core sets

Inputs
------
--csv            Sample metadata CSV with columns: ID,status
--clstr          CD-HIT .clstr file
--passfilter     BrumiR2Reference passfilter table
--nonpassfilter  BrumiR2Reference nonpassfilter table
--core_p08       Optional core sets TSV for p=0.8
--core_p10       Optional core sets TSV for p=1.0
--out_prefix     Output prefix

Outputs
-------
<prefix>.cluster_support.tsv
<prefix>.core_summary.tsv
<prefix>.core_membership.tsv (optional)
"""

"""
analysis_normalize_brumir2ref.py

Purpose
-------
Normalize and reconcile identifiers across:
  1) CD-HIT cluster definitions (cluster_X; defined by the .clstr file)
  2) BrumiR-RF member identifiers (SAMPLE|SEQID; embedded in CD-HIT .clstr)
  3) BrumiR2Reference mapping outputs (passfilter/nonpassfilter tables, where
     the miRNA ID column contains SAMPLE|SEQID)

This script produces cluster-centric reports that are directly usable for
Module III (Analysis), enabling:
  - cluster-level genome support statistics (pass/non-pass/unmapped)
  - sample/group presence summaries (athlete vs sedentary)
  - optional annotation of "core" clusters under two thresholds (e.g., 0.8 and 1.0)

Inputs
------
- CD-HIT .clstr file:
    all.candidates_clustered.fasta.clstr
  Defines cluster -> list of member IDs (SAMPLE|SEQID) and representative (*).

- BrumiR2Reference tables:
    brumir2ref.<identity>.passfilter.txt
    brumir2ref.<identity>.nonpassfilter.txt
  These contain an ID column (miRNA; column 2) matching SAMPLE|SEQID.

- Sample sheet CSV:
    files_to_process.csv
  Must contain columns: ID, status (status expected: athlete/sedentary)

- Optional core set tables (from Module II):
    brumir.<identity>.core_sets.tsv        (e.g., p=0.8)
    brumir.<identity>.p10.core_sets.tsv    (e.g., p=1.0)
  Format: set<TAB>miRNA where miRNA == cluster_X.

Outputs
-------
Given an output prefix (e.g., results/analysis/brumir095), the script writes:
  - <prefix>.cluster_support.tsv
      One row per cluster with genome support + sample/group presence + core flags.

  - <prefix>.core_summary.tsv
      High-level summaries per core set (if provided).

  - <prefix>.core_membership.tsv
      Optional: detailed listing of clusters that belong to each core set.

Notes
-----
- The .clstr file is the source of truth for cluster membership.
- passfilter/nonpassfilter IDs are matched against cluster member IDs (SAMPLE|SEQID).
- Some cluster members may be absent from both passfilter and nonpassfilter
  (e.g., not reported by BrumiR2Reference). These are counted as "unmapped".
"""

import argparse
import csv
import os
import re
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional


# Capture SAMPLE|SEQID from lines like: ">SRR11363960|30072... *"
CLSTR_MEMBER_RE = re.compile(r">(.+?)\.\.\.")


def read_samples_csv(csv_path: str) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Read sample metadata CSV.

    Required columns:
      - ID
      - status (expected values: athlete/sedentary)

    Returns
    -------
    statuses : dict
        sample_id -> status (lowercase)
    groups : dict
        status -> [sample_ids] preserving CSV order
    """
    statuses: Dict[str, str] = {}
    groups: Dict[str, List[str]] = defaultdict(list)

    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise SystemExit("ERROR: Sample CSV is empty or missing a header row.")

        required = {"ID", "status"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit(f"ERROR: Sample CSV missing required columns: {', '.join(sorted(missing))}")

        for row in reader:
            sid = (row.get("ID") or "").strip()
            st = (row.get("status") or "").strip().lower()
            if not sid:
                continue
            statuses[sid] = st
            groups[st].append(sid)

    return statuses, groups


def parse_clstr(clstr_path: str) -> Tuple[Dict[int, List[str]], Dict[int, str]]:
    """
    Parse a CD-HIT .clstr file.

    Returns
    -------
    cluster_members : dict
        cluster_id -> [member_id (SAMPLE|SEQID), ...]
    cluster_rep : dict
        cluster_id -> representative member_id (marked by '*'; fallback to first member)
    """
    cluster_members: Dict[int, List[str]] = defaultdict(list)
    cluster_rep: Dict[int, str] = {}

    current: Optional[int] = None

    with open(clstr_path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue

            if line.startswith(">Cluster"):
                current = int(line.split()[-1])
                continue

            if current is None:
                continue

            m = CLSTR_MEMBER_RE.search(line)
            if not m:
                continue

            member_id = m.group(1).strip()  # SAMPLE|SEQID
            cluster_members[current].append(member_id)

            if line.endswith("*"):
                cluster_rep[current] = member_id

    # Fallback representative: first member if none marked with '*'
    for cid, mems in cluster_members.items():
        if cid not in cluster_rep and mems:
            cluster_rep[cid] = mems[0]

    return dict(cluster_members), cluster_rep


def sample_from_member(member_id: str) -> str:
    """
    Extract sample ID from a cluster member ID.

    member_id format: SAMPLE|SEQID
    """
    if "|" in member_id:
        return member_id.split("|", 1)[0]
    return "NA"


def parse_brumir2ref_ids(tab_path: str) -> Set[str]:
    """
    Parse BrumiR2Reference passfilter/nonpassfilter tables.

    Expected format:
      - First line is a header beginning with '#'
      - Column 2 (index 1) contains miRNA ID == SAMPLE|SEQID

    Returns
    -------
    ids : set
        Set of SAMPLE|SEQID values observed in the table.
    """
    ids: Set[str] = set()

    with open(tab_path, "r", encoding="utf-8", errors="replace") as fh:
        for i, raw in enumerate(fh):
            line = raw.rstrip("\n")
            if not line:
                continue
            if i == 0 and line.startswith("#"):
                continue
            if line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            mid = parts[1].strip()
            if mid:
                ids.add(mid)

    return ids


def read_core_sets(core_tsv: str) -> Tuple[Set[str], Set[str]]:
    """
    Read core set file produced by Module II.

    Format:
      set<TAB>miRNA
    where:
      set is athlete_core or sedentary_core
      miRNA is cluster_X

    Returns
    -------
    athlete_core : set[str]
    sedentary_core : set[str]
    """
    athlete: Set[str] = set()
    sedentary: Set[str] = set()

    with open(core_tsv, "r", encoding="utf-8", errors="replace") as fh:
        _header = fh.readline()
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            sname, item = line.split("\t", 1)
            sname = sname.strip()
            item = item.strip()

            if sname == "athlete_core":
                athlete.add(item)
            elif sname == "sedentary_core":
                sedentary.add(item)

    return athlete, sedentary


def safe_div(a: int, b: int) -> float:
    return (a / b) if b else 0.0


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Normalize CD-HIT cluster membership against BrumiR2Reference IDs and generate analysis-ready reports."
    )
    ap.add_argument("--csv", required=True, help="Sample sheet CSV (must contain columns: ID,status).")
    ap.add_argument("--clstr", required=True, help="CD-HIT .clstr file defining cluster membership.")
    ap.add_argument("--passfilter", required=True, help="BrumiR2Reference passfilter table (tab-delimited).")
    ap.add_argument("--nonpassfilter", required=True, help="BrumiR2Reference nonpassfilter table (tab-delimited).")
    ap.add_argument("--core_p08", default=None, help="Optional core sets TSV (e.g., p=0.8).")
    ap.add_argument("--core_p10", default=None, help="Optional core sets TSV (e.g., p=1.0).")
    ap.add_argument("--out_prefix", required=True, help="Output prefix (e.g., results/analysis/brumir095).")
    args = ap.parse_args()

    statuses, _groups = read_samples_csv(args.csv)
    cluster_members, cluster_rep = parse_clstr(args.clstr)

    pass_ids = parse_brumir2ref_ids(args.passfilter)
    nonpass_ids = parse_brumir2ref_ids(args.nonpassfilter)

    # Optional core sets
    core08_a: Set[str] = set()
    core08_s: Set[str] = set()
    core10_a: Set[str] = set()
    core10_s: Set[str] = set()

    if args.core_p08:
        core08_a, core08_s = read_core_sets(args.core_p08)
    if args.core_p10:
        core10_a, core10_s = read_core_sets(args.core_p10)

    out_support = f"{args.out_prefix}.cluster_support.tsv"
    out_summary = f"{args.out_prefix}.core_summary.tsv"
    out_membership = f"{args.out_prefix}.core_membership.tsv"

    os.makedirs(os.path.dirname(out_support) or ".", exist_ok=True)

    # Build cluster-centric metrics
    rows = []
    for cid in sorted(cluster_members.keys()):
        mems = cluster_members[cid]
        rep = cluster_rep.get(cid, mems[0] if mems else "")
        mem_set = set(mems)

        n_members = len(mems)
        n_pass = len(mem_set & pass_ids)
        n_nonpass = len(mem_set & nonpass_ids)

        # Members not reported by BrumiR2Reference at all
        n_unmapped = n_members - len(mem_set & (pass_ids | nonpass_ids))

        # Sample/group presence based on member IDs
        sample_set = {sample_from_member(m) for m in mems}
        sample_set.discard("NA")

        a_samples = {s for s in sample_set if statuses.get(s, "") == "athlete"}
        s_samples = {s for s in sample_set if statuses.get(s, "") == "sedentary"}

        feat = f"cluster_{cid}"

        rows.append({
            "cluster": feat,
            "cluster_id": cid,
            "rep_member": rep,
            "n_members": n_members,
            "n_pass": n_pass,
            "n_nonpass": n_nonpass,
            "n_unmapped": n_unmapped,
            "pass_frac": f"{safe_div(n_pass, n_members):.4f}",
            "samples_total": len(sample_set),
            "samples_athlete": len(a_samples),
            "samples_sedentary": len(s_samples),
            "is_core08_athlete": int(feat in core08_a),
            "is_core08_sedentary": int(feat in core08_s),
            "is_core10_athlete": int(feat in core10_a),
            "is_core10_sedentary": int(feat in core10_s),
        })

    # Write cluster_support.tsv
    header = [
        "cluster", "cluster_id", "rep_member",
        "n_members", "n_pass", "n_nonpass", "n_unmapped", "pass_frac",
        "samples_total", "samples_athlete", "samples_sedentary",
        "is_core08_athlete", "is_core08_sedentary",
        "is_core10_athlete", "is_core10_sedentary",
    ]
    with open(out_support, "w", newline="") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(str(r[h]) for h in header) + "\n")

    # Core summaries (only if core sets are provided)
    def summarize_core(name: str, core_set: Set[str]) -> Dict[str, str]:
        if not core_set:
            return {
                "set": name,
                "n_clusters": "0",
                "clusters_with_any_pass": "0",
                "clusters_all_pass": "0",
                "mean_pass_frac": "NA",
            }

        sub = [r for r in rows if r["cluster"] in core_set]
        n = len(sub)

        any_pass = sum(1 for r in sub if int(r["n_pass"]) > 0)

        # All members mapped and all mapped members are passfilter
        all_pass = sum(
            1 for r in sub
            if int(r["n_unmapped"]) == 0
            and int(r["n_nonpass"]) == 0
            and int(r["n_pass"]) == int(r["n_members"])
        )

        mean_pass = sum(float(r["pass_frac"]) for r in sub) / n if n else 0.0

        return {
            "set": name,
            "n_clusters": str(n),
            "clusters_with_any_pass": str(any_pass),
            "clusters_all_pass": str(all_pass),
            "mean_pass_frac": f"{mean_pass:.4f}",
        }

    summaries = []
    if args.core_p08:
        summaries.append(summarize_core("core08_athlete", core08_a))
        summaries.append(summarize_core("core08_sedentary", core08_s))
        summaries.append(summarize_core("core08_union", core08_a | core08_s))
    if args.core_p10:
        summaries.append(summarize_core("core10_athlete", core10_a))
        summaries.append(summarize_core("core10_sedentary", core10_s))
        summaries.append(summarize_core("core10_union", core10_a | core10_s))

    with open(out_summary, "w", newline="") as out:
        out.write("set\tn_clusters\tclusters_with_any_pass\tclusters_all_pass\tmean_pass_frac\n")
        for s in summaries:
            out.write(
                f"{s['set']}\t{s['n_clusters']}\t{s['clusters_with_any_pass']}\t{s['clusters_all_pass']}\t{s['mean_pass_frac']}\n"
            )

    # Optional: detailed membership per core set (traceability/debug)
    if args.core_p08 or args.core_p10:
        with open(out_membership, "w", newline="") as out:
            out.write("core_set\tcluster\tcluster_id\trep_member\tn_members\tn_pass\tn_nonpass\tn_unmapped\tpass_frac\n")

            def emit_set(label: str, core_set: Set[str]) -> None:
                for r in rows:
                    if r["cluster"] in core_set:
                        out.write(
                            f"{label}\t{r['cluster']}\t{r['cluster_id']}\t{r['rep_member']}\t"
                            f"{r['n_members']}\t{r['n_pass']}\t{r['n_nonpass']}\t{r['n_unmapped']}\t{r['pass_frac']}\n"
                        )

            if args.core_p08:
                emit_set("core08_athlete", core08_a)
                emit_set("core08_sedentary", core08_s)
            if args.core_p10:
                emit_set("core10_athlete", core10_a)
                emit_set("core10_sedentary", core10_s)

    print("OK")
    print(f"  wrote: {out_support}")
    print(f"  wrote: {out_summary}")
    if args.core_p08 or args.core_p10:
        print(f"  wrote: {out_membership}")


if __name__ == "__main__":
    main()