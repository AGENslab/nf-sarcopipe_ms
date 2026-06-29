#!/usr/bin/env python3

"""
Figure 3B – BrumiR2Reference cluster support analysis

Description
-----------
This script integrates CD-HIT cluster membership, BrumiR2Reference
mapping support, and sample metadata to generate cluster-level support
summaries for downstream manuscript analysis.

It computes:
- number of passfilter, nonpassfilter, and unmapped members per cluster
- sample and group presence per cluster
- optional core membership flags using p=0.8 and/or p=1.0 core sets

Inputs
------
--csv
    Sample metadata CSV with columns: ID,status

--clstr
    CD-HIT .clstr file defining cluster membership.

--passfilter
    BrumiR2Reference passfilter table.

--nonpassfilter
    BrumiR2Reference nonpassfilter table.

--core_p08
    Optional core sets TSV for p=0.8.

--core_p10
    Optional core sets TSV for p=1.0.

--out_prefix
    Output prefix.

Outputs
-------
<out_prefix>.cluster_support.tsv
<out_prefix>.core_summary.tsv
<out_prefix>.core_membership.tsv, if core sets are provided

Usage
-----
python3 Figure3B_analysis_BrumiR2Reference.py \\
  --csv files_to_process.csv \\
  --clstr all.candidates_clustered.fasta.clstr \\
  --passfilter brumir2ref.passfilter.txt \\
  --nonpassfilter brumir2ref.nonpassfilter.txt \\
  --core_p08 brumir.core_sets.tsv \\
  --core_p10 brumir.p10.core_sets.tsv \\
  --out_prefix results/analysis/brumir095

Notes
-----
- The .clstr file is the source of truth for cluster membership.
- passfilter/nonpassfilter IDs are matched against cluster member IDs.
- Unreported members are counted as unmapped.
"""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


CLSTR_MEMBER_RE = re.compile(r">(.+?)\.\.\.")


def validate_file(path: str, label: str) -> Path:
    """Validate that an input file exists."""
    file_path = Path(path)

    if not file_path.is_file():
        raise SystemExit(f"ERROR: {label} not found: {file_path}")

    return file_path


def read_samples_csv(csv_path: Path) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
    """
    Read sample metadata CSV.

    Required columns:
    - ID
    - status
    """
    statuses: Dict[str, str] = {}
    groups: Dict[str, List[str]] = defaultdict(list)

    with csv_path.open(newline="") as fh:
        reader = csv.DictReader(fh)

        if not reader.fieldnames:
            raise SystemExit("ERROR: Sample CSV is empty or missing a header row.")

        required_columns = {"ID", "status"}
        missing_columns = required_columns - set(reader.fieldnames)

        if missing_columns:
            raise SystemExit(
                "ERROR: Sample CSV missing required columns: "
                f"{', '.join(sorted(missing_columns))}"
            )

        for row in reader:
            sample_id = (row.get("ID") or "").strip()
            status = (row.get("status") or "").strip().lower()

            if not sample_id:
                continue

            statuses[sample_id] = status
            groups[status].append(sample_id)

    return statuses, groups


def parse_clstr(clstr_path: Path) -> Tuple[Dict[int, List[str]], Dict[int, str]]:
    """
    Parse a CD-HIT .clstr file.

    Returns:
    - cluster_members: cluster_id -> member IDs
    - cluster_rep: cluster_id -> representative member ID
    """
    cluster_members: Dict[int, List[str]] = defaultdict(list)
    cluster_rep: Dict[int, str] = {}

    current_cluster: Optional[int] = None

    with clstr_path.open("r") as fh:
        for raw_line in fh:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[-1])
                continue

            if current_cluster is None:
                continue

            match = CLSTR_MEMBER_RE.search(line)
            if not match:
                continue

            member_id = match.group(1).strip()
            cluster_members[current_cluster].append(member_id)

            if line.endswith("*"):
                cluster_rep[current_cluster] = member_id

    for cluster_id, members in cluster_members.items():
        if cluster_id not in cluster_rep and members:
            cluster_rep[cluster_id] = members[0]

    return dict(cluster_members), cluster_rep


def sample_from_member(member_id: str) -> str:
    """Extract sample ID from a member ID formatted as SAMPLE|SEQID."""
    if "|" in member_id:
        return member_id.split("|", 1)[0]

    return "NA"


def parse_brumir2ref_ids(table_path: Path) -> Set[str]:
    """
    Parse BrumiR2Reference passfilter/nonpassfilter tables.

    Column 2 is expected to contain miRNA IDs matching SAMPLE|SEQID.
    """
    ids: Set[str] = set()

    with table_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line_number, raw_line in enumerate(fh):
            line = raw_line.rstrip("\n")

            if not line:
                continue

            if line_number == 0 and line.startswith("#"):
                continue

            if line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 2:
                continue

            mirna_id = fields[1].strip()

            if mirna_id:
                ids.add(mirna_id)

    return ids


def read_core_sets(core_tsv: Path) -> Tuple[Set[str], Set[str]]:
    """
    Read core set file produced by Module II.

    Expected format:
    set<TAB>miRNA
    """
    athlete_core: Set[str] = set()
    sedentary_core: Set[str] = set()

    with core_tsv.open("r", encoding="utf-8", errors="replace") as fh:
        _header = fh.readline()

        for raw_line in fh:
            line = raw_line.strip()

            if not line:
                continue

            set_name, item = line.split("\t", 1)
            set_name = set_name.strip()
            item = item.strip()

            if set_name == "athlete_core":
                athlete_core.add(item)
            elif set_name == "sedentary_core":
                sedentary_core.add(item)

    return athlete_core, sedentary_core


def safe_div(numerator: int, denominator: int) -> float:
    """Divide safely, returning 0.0 when denominator is zero."""
    return (numerator / denominator) if denominator else 0.0


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Normalize CD-HIT cluster membership against BrumiR2Reference IDs "
            "and generate analysis-ready reports."
        )
    )

    parser.add_argument(
        "--csv",
        required=True,
        help="Sample sheet CSV. Must contain columns: ID,status.",
    )
    parser.add_argument(
        "--clstr",
        required=True,
        help="CD-HIT .clstr file defining cluster membership.",
    )
    parser.add_argument(
        "--passfilter",
        required=True,
        help="BrumiR2Reference passfilter table.",
    )
    parser.add_argument(
        "--nonpassfilter",
        required=True,
        help="BrumiR2Reference nonpassfilter table.",
    )
    parser.add_argument(
        "--core_p08",
        default=None,
        help="Optional core sets TSV, for example p=0.8.",
    )
    parser.add_argument(
        "--core_p10",
        default=None,
        help="Optional core sets TSV, for example p=1.0.",
    )
    parser.add_argument(
        "--out_prefix",
        required=True,
        help="Output prefix, for example results/analysis/brumir095.",
    )

    return parser.parse_args()


def main() -> None:
    """Run BrumiR2Reference cluster support analysis."""
    args = parse_args()

    csv_path = validate_file(args.csv, "Sample CSV")
    clstr_path = validate_file(args.clstr, "CD-HIT .clstr file")
    passfilter_path = validate_file(args.passfilter, "Passfilter table")
    nonpassfilter_path = validate_file(args.nonpassfilter, "Nonpassfilter table")

    core_p08_path = validate_file(args.core_p08, "Core p=0.8 TSV") if args.core_p08 else None
    core_p10_path = validate_file(args.core_p10, "Core p=1.0 TSV") if args.core_p10 else None

    out_prefix = Path(args.out_prefix)
    out_support = Path(f"{out_prefix}.cluster_support.tsv")
    out_summary = Path(f"{out_prefix}.core_summary.tsv")
    out_membership = Path(f"{out_prefix}.core_membership.tsv")

    out_support.parent.mkdir(parents=True, exist_ok=True)

    statuses, _groups = read_samples_csv(csv_path)
    cluster_members, cluster_rep = parse_clstr(clstr_path)

    pass_ids = parse_brumir2ref_ids(passfilter_path)
    nonpass_ids = parse_brumir2ref_ids(nonpassfilter_path)

    core08_athlete: Set[str] = set()
    core08_sedentary: Set[str] = set()
    core10_athlete: Set[str] = set()
    core10_sedentary: Set[str] = set()

    if core_p08_path:
        core08_athlete, core08_sedentary = read_core_sets(core_p08_path)

    if core_p10_path:
        core10_athlete, core10_sedentary = read_core_sets(core_p10_path)

    rows = []

    for cluster_id in sorted(cluster_members.keys()):
        members = cluster_members[cluster_id]
        representative = cluster_rep.get(cluster_id, members[0] if members else "")
        member_set = set(members)

        n_members = len(members)
        n_pass = len(member_set & pass_ids)
        n_nonpass = len(member_set & nonpass_ids)
        n_unmapped = n_members - len(member_set & (pass_ids | nonpass_ids))

        sample_set = {sample_from_member(member) for member in members}
        sample_set.discard("NA")

        athlete_samples = {
            sample for sample in sample_set
            if statuses.get(sample, "") == "athlete"
        }
        sedentary_samples = {
            sample for sample in sample_set
            if statuses.get(sample, "") == "sedentary"
        }

        cluster_name = f"cluster_{cluster_id}"

        rows.append({
            "cluster": cluster_name,
            "cluster_id": cluster_id,
            "rep_member": representative,
            "n_members": n_members,
            "n_pass": n_pass,
            "n_nonpass": n_nonpass,
            "n_unmapped": n_unmapped,
            "pass_frac": f"{safe_div(n_pass, n_members):.4f}",
            "samples_total": len(sample_set),
            "samples_athlete": len(athlete_samples),
            "samples_sedentary": len(sedentary_samples),
            "is_core08_athlete": int(cluster_name in core08_athlete),
            "is_core08_sedentary": int(cluster_name in core08_sedentary),
            "is_core10_athlete": int(cluster_name in core10_athlete),
            "is_core10_sedentary": int(cluster_name in core10_sedentary),
        })

    header = [
        "cluster",
        "cluster_id",
        "rep_member",
        "n_members",
        "n_pass",
        "n_nonpass",
        "n_unmapped",
        "pass_frac",
        "samples_total",
        "samples_athlete",
        "samples_sedentary",
        "is_core08_athlete",
        "is_core08_sedentary",
        "is_core10_athlete",
        "is_core10_sedentary",
    ]

    with out_support.open("w", newline="") as out:
        out.write("\t".join(header) + "\n")

        for row in rows:
            out.write("\t".join(str(row[column]) for column in header) + "\n")

    def summarize_core(name: str, core_set: Set[str]) -> Dict[str, str]:
        """Summarize BrumiR2Reference support within a core set."""
        if not core_set:
            return {
                "set": name,
                "n_clusters": "0",
                "clusters_with_any_pass": "0",
                "clusters_all_pass": "0",
                "mean_pass_frac": "NA",
            }

        subset = [row for row in rows if row["cluster"] in core_set]
        n_clusters = len(subset)

        clusters_with_any_pass = sum(
            1 for row in subset
            if int(row["n_pass"]) > 0
        )

        clusters_all_pass = sum(
            1 for row in subset
            if int(row["n_unmapped"]) == 0
            and int(row["n_nonpass"]) == 0
            and int(row["n_pass"]) == int(row["n_members"])
        )

        mean_pass_frac = (
            sum(float(row["pass_frac"]) for row in subset) / n_clusters
            if n_clusters else 0.0
        )

        return {
            "set": name,
            "n_clusters": str(n_clusters),
            "clusters_with_any_pass": str(clusters_with_any_pass),
            "clusters_all_pass": str(clusters_all_pass),
            "mean_pass_frac": f"{mean_pass_frac:.4f}",
        }

    summaries = []

    if core_p08_path:
        summaries.append(summarize_core("core08_athlete", core08_athlete))
        summaries.append(summarize_core("core08_sedentary", core08_sedentary))
        summaries.append(summarize_core("core08_union", core08_athlete | core08_sedentary))

    if core_p10_path:
        summaries.append(summarize_core("core10_athlete", core10_athlete))
        summaries.append(summarize_core("core10_sedentary", core10_sedentary))
        summaries.append(summarize_core("core10_union", core10_athlete | core10_sedentary))

    with out_summary.open("w", newline="") as out:
        out.write(
            "set\tn_clusters\tclusters_with_any_pass\t"
            "clusters_all_pass\tmean_pass_frac\n"
        )

        for summary in summaries:
            out.write(
                f"{summary['set']}\t"
                f"{summary['n_clusters']}\t"
                f"{summary['clusters_with_any_pass']}\t"
                f"{summary['clusters_all_pass']}\t"
                f"{summary['mean_pass_frac']}\n"
            )

    if core_p08_path or core_p10_path:
        with out_membership.open("w", newline="") as out:
            out.write(
                "core_set\tcluster\tcluster_id\trep_member\tn_members\t"
                "n_pass\tn_nonpass\tn_unmapped\tpass_frac\n"
            )

            def emit_set(label: str, core_set: Set[str]) -> None:
                for row in rows:
                    if row["cluster"] in core_set:
                        out.write(
                            f"{label}\t{row['cluster']}\t{row['cluster_id']}\t"
                            f"{row['rep_member']}\t{row['n_members']}\t"
                            f"{row['n_pass']}\t{row['n_nonpass']}\t"
                            f"{row['n_unmapped']}\t{row['pass_frac']}\n"
                        )

            if core_p08_path:
                emit_set("core08_athlete", core08_athlete)
                emit_set("core08_sedentary", core08_sedentary)

            if core_p10_path:
                emit_set("core10_athlete", core10_athlete)
                emit_set("core10_sedentary", core10_sedentary)

    print("OK")
    print(f"  wrote: {out_support}")
    print(f"  wrote: {out_summary}")

    if core_p08_path or core_p10_path:
        print(f"  wrote: {out_membership}")


if __name__ == "__main__":
    main()