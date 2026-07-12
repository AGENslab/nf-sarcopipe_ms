#!/usr/bin/env python3
"""
summarize_mirdeep2_raw_predictions.py

Purpose
-------
Summarize the raw per-sample miRDeep2 predictions produced after processing
the 18–28 nt small-RNA reads.

For each sample, the script counts the number of entries in:

  known_miRNAs_<sample>.tsv
  novel_miRNAs_<sample>.tsv

These are raw miRDeep2 per-library results. The counts are calculated before:

  - presence filtering (for example, p >= 0.8 or p >= 0.3)
  - core-set selection
  - differential-expression analysis
  - overlap or redundancy filtering across samples
  - downstream target integration

Therefore, the output is suitable for preprocessing/QC tables reporting the
number of known and novel miRNAs detected by miRDeep2 in each dataset.

Inputs
------
--known_dir
    Directory containing known_miRNAs_<sample>.tsv files.

--novel_dir
    Directory containing novel_miRNAs_<sample>.tsv files.

--out
    Output TSV file.

Outputs
-------
A TSV table with these columns:

  sample_id
  mirdeep2_known_raw
  mirdeep2_novel_raw
  known_file
  novel_file

A final TOTAL row is added with the accumulated known and novel counts.

Counting rule
-------------
The first non-empty line of each file is treated as the header and excluded.
Every remaining non-empty line is counted as one miRDeep2 prediction.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Optional


KNOWN_RE = re.compile(r"^known_miRNAs_(.+)\.tsv$")
NOVEL_RE = re.compile(r"^novel_miRNAs_(.+)\.tsv$")


def collect_files(directory: Path, pattern: re.Pattern[str]) -> Dict[str, Path]:
    """Return sample_id -> file for all matching non-empty TSV files."""
    files: Dict[str, Path] = {}

    if not directory.is_dir():
        raise SystemExit(f"ERROR: directory does not exist: {directory}")

    for path in sorted(directory.glob("*.tsv")):
        match = pattern.match(path.name)
        if not match:
            continue

        if path.stat().st_size == 0:
            print(f"WARNING: empty file ignored: {path}", file=sys.stderr)
            continue

        sample_id = match.group(1)
        files[sample_id] = path

    return files


def count_predictions(path: Optional[Path]) -> int:
    """Count non-empty data rows after the first non-empty line (header)."""
    if path is None:
        return 0

    nonempty_lines = 0

    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.strip():
                nonempty_lines += 1

    return max(nonempty_lines - 1, 0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Count raw known and novel miRDeep2 predictions per sample."
    )

    parser.add_argument(
        "--known_dir",
        required=True,
        help="Directory with known_miRNAs_<sample>.tsv files.",
    )
    parser.add_argument(
        "--novel_dir",
        required=True,
        help="Directory with novel_miRNAs_<sample>.tsv files.",
    )
    parser.add_argument(
        "--out",
        default="mirdeep2_raw_predictions_summary.tsv",
        help="Output TSV path.",
    )

    args = parser.parse_args()

    known_dir = Path(args.known_dir)
    novel_dir = Path(args.novel_dir)
    out_path = Path(args.out)

    known_files = collect_files(known_dir, KNOWN_RE)
    novel_files = collect_files(novel_dir, NOVEL_RE)

    samples = sorted(set(known_files) | set(novel_files))

    if not samples:
        raise SystemExit(
            "ERROR: no known_miRNAs_*.tsv or novel_miRNAs_*.tsv files found."
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_known = 0
    total_novel = 0

    with out_path.open("w", encoding="utf-8") as output:
        output.write(
            "sample_id\tmirdeep2_known_raw\tmirdeep2_novel_raw\t"
            "known_file\tnovel_file\n"
        )

        for sample_id in samples:
            known_path = known_files.get(sample_id)
            novel_path = novel_files.get(sample_id)

            known_count = count_predictions(known_path)
            novel_count = count_predictions(novel_path)

            total_known += known_count
            total_novel += novel_count

            output.write(
                f"{sample_id}\t{known_count}\t{novel_count}\t"
                f"{known_path if known_path else 'MISSING'}\t"
                f"{novel_path if novel_path else 'MISSING'}\n"
            )

        output.write(f"TOTAL\t{total_known}\t{total_novel}\tNA\tNA\n")

    print(f"Samples summarized: {len(samples)}")
    print(f"Total raw known miRNAs: {total_known}")
    print(f"Total raw novel miRNAs: {total_novel}")
    print(f"[OK] Wrote: {out_path}")


if __name__ == "__main__":
    main()
