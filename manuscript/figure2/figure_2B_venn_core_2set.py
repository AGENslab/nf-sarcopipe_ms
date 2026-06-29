#!/usr/bin/env python3
"""
Figure 2B – Core miRNA overlap.

Create a 2-set Venn diagram comparing active and sedentary miRNA core sets
from a TSV file.

Expected input
--------------
A TSV file with two columns:

    set    miRNA

Expected set labels
-------------------
- athlete_core
- sedentary_core

Output
------
PNG image containing a 2-set Venn diagram.
"""

from pathlib import Path
import argparse


def read_core_sets(path: Path) -> dict[str, set[str]]:
    """Read a two-column core-set TSV file: set<TAB>miRNA."""
    sets: dict[str, set[str]] = {}

    with path.open("r") as handle:
        header = handle.readline()

        for line_number, line in enumerate(handle, start=2):
            line = line.strip()
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Invalid line {line_number} in {path}: expected at least 2 tab-separated columns."
                )

            set_name, item = fields[0], fields[1]
            sets.setdefault(set_name, set()).add(item.strip())

    return sets


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a 2-set Venn plot for active and sedentary miRNA core sets."
    )
    parser.add_argument(
        "--core_sets",
        required=True,
        type=Path,
        help="Input TSV with columns: set and miRNA.",
    )
    parser.add_argument(
        "--out_png",
        required=True,
        type=Path,
        help="Output PNG path.",
    )
    args = parser.parse_args()

    sets = read_core_sets(args.core_sets)

    required_sets = {"athlete_core", "sedentary_core"}
    missing = required_sets - set(sets)
    if missing:
        raise ValueError(
            f"Missing expected set labels in {args.core_sets}: {', '.join(sorted(missing))}"
        )

    active_core = sets["athlete_core"]
    sedentary_core = sets["sedentary_core"]

    only_active = len(active_core - sedentary_core)
    only_sedentary = len(sedentary_core - active_core)
    overlap = len(active_core & sedentary_core)

    args.out_png.parent.mkdir(parents=True, exist_ok=True)

    import matplotlib.pyplot as plt

    try:
        from matplotlib_venn import venn2  # type: ignore

        fig = plt.figure(figsize=(5, 5))
        venn2(
            subsets=(only_active, only_sedentary, overlap),
            set_labels=("Active core", "Sedentary core"),
        )
        plt.title("Core miRNAs")
        plt.tight_layout()
        plt.savefig(args.out_png, dpi=300)
        plt.close(fig)
        return

    except ImportError:
        pass

    from matplotlib.patches import Circle

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.add_patch(Circle((0.45, 0.5), 0.30, fill=False, linewidth=2))
    ax.add_patch(Circle((0.65, 0.5), 0.30, fill=False, linewidth=2))

    ax.text(0.32, 0.50, str(only_active), ha="center", va="center", fontsize=14)
    ax.text(0.78, 0.50, str(only_sedentary), ha="center", va="center", fontsize=14)
    ax.text(0.55, 0.50, str(overlap), ha="center", va="center", fontsize=14)

    ax.text(0.35, 0.18, "Active core", ha="center", va="center", fontsize=10)
    ax.text(0.75, 0.18, "Sedentary core", ha="center", va="center", fontsize=10)

    ax.set_title("Core miRNAs")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(args.out_png, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
