#!/usr/bin/env python3

"""
============================================================
summarize_fastqc_test_brumir.py
============================================================

Purpose
-------
Summarize the FastQC quality-control reports generated for the
input small RNA sequencing reads used in the BrumiR workflow.

The script reads all `*_fastqc.zip` files within a directory,
extracts the `fastqc_data.txt` report from each archive, and
parses the main FastQC modules to generate a single tab-separated
summary table containing one row per sample.

Reported metrics
----------------
For each sample, the script reports:

- Basic sequencing statistics (total reads, sequence length, GC content)
- Mean per-sequence Phred quality score
- Percentage of sequences with quality scores ≥ Q20 and ≥ Q30
- Minimum and final median per-base quality scores
- Position with the lowest median base quality
- Maximum adapter contamination, including adapter identity and position
- Number and abundance of overrepresented sequences
- Estimated sequence duplication percentage
- PASS/WARN/FAIL status for key FastQC quality modules.

Input
-----
Directory containing FastQC compressed reports:

    *_fastqc.zip

Output
------
A tab-separated summary table containing one row per sample,
providing a compact overview of sequencing quality prior to
BrumiR-based de novo miRNA discovery.
"""

import csv
import re
import sys
import zipfile
from pathlib import Path


def read_fastqc_zip(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        target = next(
            name for name in z.namelist()
            if name.endswith("fastqc_data.txt")
        )
        return z.read(target).decode("utf-8", errors="replace")


def parse_modules(text):
    statuses = {}
    contents = {}
    current = None

    for line in text.splitlines():
        if line.startswith(">>END_MODULE"):
            current = None
            continue

        if line.startswith(">>"):
            parts = line[2:].split("\t")
            if len(parts) >= 2:
                current = parts[0]
                statuses[current] = parts[1].upper()
                contents[current] = []
            continue

        if current:
            contents[current].append(line)

    return statuses, contents


def parse_basic(lines):
    data = {}

    for line in lines:
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t", 1)

        if len(parts) == 2:
            data[parts[0]] = parts[1]

    return data


def parse_per_sequence_quality(lines):
    total = 0.0
    weighted = 0.0
    q20 = 0.0
    q30 = 0.0

    for line in lines:
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")

        if len(parts) < 2:
            continue

        quality = float(parts[0])
        count = float(parts[1])

        total += count
        weighted += quality * count

        if quality >= 20:
            q20 += count

        if quality >= 30:
            q30 += count

    if total == 0:
        return "NA", "NA", "NA"

    return (
        f"{weighted / total:.2f}",
        f"{100 * q20 / total:.2f}",
        f"{100 * q30 / total:.2f}",
    )


def parse_per_base_quality(lines):
    values = []

    for line in lines:
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")

        if len(parts) >= 3:
            position = parts[0]
            median = float(parts[2])
            values.append((position, median))

    if not values:
        return "NA", "NA", "NA"

    worst_position, minimum = min(values, key=lambda x: x[1])
    final_median = values[-1][1]

    return (
        f"{minimum:.2f}",
        f"{final_median:.2f}",
        worst_position,
    )


def parse_adapter_content(lines):
    header = []
    maximum = None
    maximum_name = "NA"
    maximum_position = "NA"

    for line in lines:
        if line.startswith("#"):
            header = line.lstrip("#").split("\t")
            continue

        if not line:
            continue

        parts = line.split("\t")

        for index, value in enumerate(parts[1:], start=1):
            percentage = float(value)

            if maximum is None or percentage > maximum:
                maximum = percentage
                maximum_position = parts[0]

                if index < len(header):
                    maximum_name = header[index]

    return (
        f"{maximum:.2f}" if maximum is not None else "NA",
        maximum_name,
        maximum_position,
    )


def parse_overrepresented(lines):
    count = 0
    maximum = None
    maximum_sequence = "NA"
    maximum_source = "NA"

    for line in lines:
        if not line or line.startswith("#"):
            continue

        parts = line.split("\t")

        if len(parts) < 3:
            continue

        count += 1
        percentage = float(parts[2])

        if maximum is None or percentage > maximum:
            maximum = percentage
            maximum_sequence = parts[0]
            maximum_source = parts[3] if len(parts) >= 4 else "NA"

    return (
        count,
        f"{maximum:.2f}" if maximum is not None else "NA",
        maximum_sequence,
        maximum_source,
    )


def parse_duplication(lines):
    for line in lines:
        if line.startswith("#Total Deduplicated Percentage"):
            parts = line.split("\t")

            if len(parts) >= 2:
                deduplicated = float(parts[1])
                return f"{100 - deduplicated:.2f}"

    return "NA"


def sample_from_filename(filename):
    name = Path(filename).name

    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            name = name[:-len(suffix)]

    return name


def parse_report(zip_path):
    text = read_fastqc_zip(zip_path)
    statuses, modules = parse_modules(text)

    basic = parse_basic(modules.get("Basic Statistics", []))

    mean_q, q20, q30 = parse_per_sequence_quality(
        modules.get("Per sequence quality scores", [])
    )

    min_median, final_median, worst_position = parse_per_base_quality(
        modules.get("Per base sequence quality", [])
    )

    max_adapter, adapter_name, adapter_position = parse_adapter_content(
        modules.get("Adapter Content", [])
    )

    (
        overrep_count,
        overrep_max,
        overrep_sequence,
        overrep_source,
    ) = parse_overrepresented(
        modules.get("Overrepresented sequences", [])
    )

    filename = basic.get("Filename", zip_path.stem)

    return {
        "sample": sample_from_filename(filename),
        "total_sequences": basic.get("Total Sequences", "NA"),
        "poor_quality_sequences": basic.get(
            "Sequences flagged as poor quality",
            "NA",
        ),
        "sequence_length": basic.get("Sequence length", "NA"),
        "gc_percent": basic.get("%GC", "NA"),
        "mean_sequence_quality": mean_q,
        "sequences_q20_percent": q20,
        "sequences_q30_percent": q30,
        "minimum_median_base_quality": min_median,
        "final_median_base_quality": final_median,
        "worst_quality_position": worst_position,
        "maximum_adapter_percent": max_adapter,
        "maximum_adapter_name": adapter_name,
        "maximum_adapter_position": adapter_position,
        "overrepresented_sequence_count": overrep_count,
        "maximum_overrepresented_percent": overrep_max,
        "maximum_overrepresented_sequence": overrep_sequence,
        "maximum_overrepresented_source": overrep_source,
        "estimated_duplication_percent": parse_duplication(
            modules.get("Sequence Duplication Levels", [])
        ),
        "status_per_base_sequence_quality": statuses.get(
            "Per base sequence quality",
            "NA",
        ),
        "status_per_sequence_quality_scores": statuses.get(
            "Per sequence quality scores",
            "NA",
        ),
        "status_sequence_length_distribution": statuses.get(
            "Sequence Length Distribution",
            "NA",
        ),
        "status_adapter_content": statuses.get(
            "Adapter Content",
            "NA",
        ),
        "status_overrepresented_sequences": statuses.get(
            "Overrepresented sequences",
            "NA",
        ),
        "status_sequence_duplication_levels": statuses.get(
            "Sequence Duplication Levels",
            "NA",
        ),
    }


def main():
    if len(sys.argv) != 3:
        raise SystemExit(
            "Uso: python3 summarize_fastqc_test.py "
            "<directorio_fastqc> <salida.tsv>"
        )

    input_dir = Path(sys.argv[1])
    output_file = Path(sys.argv[2])

    zip_files = sorted(input_dir.glob("*_fastqc.zip"))

    if not zip_files:
        raise SystemExit(
            f"ERROR: no se encontraron *_fastqc.zip en {input_dir}"
        )

    rows = [parse_report(path) for path in zip_files]

    with output_file.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=rows[0].keys(),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"FastQC procesados: {len(rows)}")
    print(f"Tabla creada: {output_file}")


if __name__ == "__main__":
    main()
