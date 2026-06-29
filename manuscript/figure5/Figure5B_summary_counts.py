#!/usr/bin/env python3
from pathlib import Path
import csv

BASE = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2")
INPUT = BASE / "input" / "input_mRNA_36"
OUTPUT = BASE / "output"

MIRNA_TABLE = INPUT / "all_miRNA_seed_summary.tsv"
GENE_TABLE = INPUT / "DEG_36_for_networks.csv"
BRUMIR_MATCH = OUTPUT / "seed_matches_36genes_BrumiR.tsv"
MIRDEEP_MATCH = OUTPUT / "seed_matches_36genes_miRDeep2.tsv"

OUT_SUMMARY = OUTPUT / "fig5a_summary.tsv"
OUT_COHERENT = OUTPUT / "seed_matches_coherent.tsv"

def read_mirnas():
    d = {"BrumiR": set(), "miRDeep2": set()}
    with open(MIRNA_TABLE, encoding="utf-8", errors="ignore") as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 7:
                continue
            source = parts[0].strip()
            renamed = parts[2].strip()
            if source in d and renamed:
                d[source].add(renamed)
    return d

def read_gene_direction():
    d = {}
    with open(GENE_TABLE, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        for r in reader:
            d[r["gene_symbol"].strip()] = r["direction"].strip()
    return d

def read_matches(path):
    rows = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            rows.append({k.strip(): v.strip() for k, v in r.items()})
    return rows

def is_coherent(mirna_dir, gene_dir):
    return (
        (mirna_dir == "Up_in_Active" and gene_dir == "Up_in_Sedentary") or
        (mirna_dir == "Up_in_Sedentary" and gene_dir == "Up_in_Active")
    )

def summarize(source, rows, gene_dir, n_mirnas):
    genes = set()
    coherent_rows = []

    for r in rows:
        gene = r["gene_symbol"]
        genes.add(gene)
        gdir = gene_dir.get(gene, "")
        if is_coherent(r["direction"], gdir):
            rr = dict(r)
            rr["gene_direction"] = gdir
            coherent_rows.append(rr)

    summary = {
        "source": source,
        "n_miRNAs": n_mirnas,
        "n_genes": len(genes),
        "n_pairs": len(rows),
        "n_coherent": len(coherent_rows),
    }
    return summary, coherent_rows

def main():
    mirnas = read_mirnas()
    gene_dir = read_gene_direction()

    br_rows = read_matches(BRUMIR_MATCH)
    md_rows = read_matches(MIRDEEP_MATCH)

    br_sum, br_coh = summarize("BrumiR", br_rows, gene_dir, len(mirnas["BrumiR"]))
    md_sum, md_coh = summarize("miRDeep2", md_rows, gene_dir, len(mirnas["miRDeep2"]))

    with open(OUT_SUMMARY, "w", encoding="utf-8") as f:
        f.write("source\tn_miRNAs\tn_genes\tn_pairs\tn_coherent\n")
        for row in [br_sum, md_sum]:
            f.write(
                f'{row["source"]}\t{row["n_miRNAs"]}\t{row["n_genes"]}\t{row["n_pairs"]}\t{row["n_coherent"]}\n'
            )

    all_coh = br_coh + md_coh
    if all_coh:
        keys = list(all_coh[0].keys())
        with open(OUT_COHERENT, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
            writer.writeheader()
            for r in all_coh:
                writer.writerow(r)

    print("=== FIG 5A SUMMARY ===")
    print(br_sum)
    print(md_sum)
    print("Summary written:", OUT_SUMMARY)
    print("Coherent pairs written:", OUT_COHERENT)

if __name__ == "__main__":
    main()
