#!/usr/bin/env python3
"""
venn_core_2set.py

Create a 2-set Venn diagram comparing athlete_core vs sedentary_core
from a core-sets TSV file.

Purpose
-------
This script reads a table of core miRNA assignments and generates a PNG
showing the overlap between:

- athlete_core
- sedentary_core

Expected input format
---------------------
A TSV file with two columns:

    set    miRNA

Example:

    athlete_core    hsa-miR-1-3p
    athlete_core    hsa-miR-133a-3p
    sedentary_core  hsa-miR-1-3p

Expected set names
------------------
The script expects the following set labels:
- athlete_core
- sedentary_core

Output
------
A PNG image containing a 2-set Venn diagram.

Notes
-----
- The script first tries to use `matplotlib_venn`.
- If that package is not available, it falls back to a simple circle-based plot.
"""

import argparse


def read_core_sets(path: str):
    """
    Read a core-sets TSV file with format:

        set<TAB>miRNA

    Returns
    -------
    dict
        set_name -> set(miRNA)
    """
    sets = {}

    with open(path, "r") as fh:
        header = fh.readline()
        for line in fh:
            line = line.strip()
            if not line:
                continue
            sname, item = line.split("\t", 1)
            sets.setdefault(sname, set()).add(item.strip())

    return sets


def main():
    ap = argparse.ArgumentParser(
        description="Create a 2-set Venn plot for core sets (athlete_core vs sedentary_core)."
    )
    ap.add_argument("--core_sets", required=True, help="Input core_sets.tsv with columns: set, miRNA")
    ap.add_argument("--out_png", required=True, help="Output PNG path")
    args = ap.parse_args()

    sets = read_core_sets(args.core_sets)
    A = sets.get("athlete_core", set())
    B = sets.get("sedentary_core", set())

    onlyA = len(A - B)
    onlyB = len(B - A)
    inter = len(A & B)

    import matplotlib.pyplot as plt

    # Preferred option: use matplotlib_venn if available
    try:
        from matplotlib_venn import venn2  # type: ignore

        fig = plt.figure(figsize=(5, 5))
        venn2(subsets=(onlyA, onlyB, inter), set_labels=("athlete_core", "sedentary_core"))
        plt.title("Core miRNAs (presence-based)")
        plt.tight_layout()
        plt.savefig(args.out_png, dpi=200)
        plt.close(fig)
        return

    except Exception:
        pass

    # Fallback: draw two simple circles if matplotlib_venn is unavailable
    from matplotlib.patches import Circle

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.add_patch(Circle((0.45, 0.5), 0.30, fill=False, linewidth=2))
    ax.add_patch(Circle((0.65, 0.5), 0.30, fill=False, linewidth=2))

    ax.text(0.32, 0.50, str(onlyA), ha="center", va="center", fontsize=14)
    ax.text(0.78, 0.50, str(onlyB), ha="center", va="center", fontsize=14)
    ax.text(0.55, 0.50, str(inter), ha="center", va="center", fontsize=14)

    ax.text(0.35, 0.18, "athlete_core", ha="center", va="center", fontsize=10)
    ax.text(0.75, 0.18, "sedentary_core", ha="center", va="center", fontsize=10)

    ax.set_title("Core miRNAs (presence-based)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(args.out_png, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
