#!/usr/bin/env python3
# ============================================================
# scan_candidate_7mer_space.py
# Description:
# Generates all possible 7-mers from selected candidate miRNA
# sequences and compares them against a dictionary of known
# canonical miRNA seeds.
#
# Inputs:
#   - input/brumir_DE_seeds_para_analisis_kmers.tsv
#   - kmer_seed_analysis/known_seed_dictionary.tsv
#
# Outputs:
#   - kmer_seed_analysis/candidate_7mer_long.tsv
#   - kmer_seed_analysis/candidate_7mer_summary.tsv
# ============================================================

from pathlib import Path
import csv
from collections import defaultdict

inp = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/input/brumir_DE_seeds_para_analisis_kmers.tsv")
seed_db = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/kmer_seed_analysis/known_seed_dictionary.tsv")

out_long = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/kmer_seed_analysis/candidate_7mer_long.tsv")
out_summary = Path("/mnt/beegfs/home/npoblete/sarcopipe/bin/module_3_analysis/v2/scripts/kmer_seed_analysis/candidate_7mer_summary.tsv")

# load known seeds
seed_to_ids = defaultdict(list)
seed_to_families = defaultdict(set)

with seed_db.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        seed = r["canonical_seed"].strip().upper().replace("T", "U")
        seed_to_ids[seed].append(r["miRNA_id"])
        seed_to_families[seed].add(r["family"])

# define controls
control_map = {
    "hsa_miR-660-5p": "known_miRNA_control",
    "known_seed_A": "known_seed_control",
    "hsa-miR-novel_A": "novel_candidate",
    "hsa-miR-novel_B": "novel_candidate",
    "hsa-miR-novel_C": "novel_candidate",
    "hsa-miR-novel_D": "novel_candidate",
}

long_rows = []
summary_rows = []

with inp.open() as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        name, mature_seq, canonical_seed = line.split("\t")
        mature_seq = mature_seq.upper().replace("T", "U")
        canonical_seed = canonical_seed.upper().replace("T", "U")

        total_7mers = 0
        matched_positions = []
        matched_seeds = set()
        matched_families = set()

        for i in range(len(mature_seq) - 6):
            kmer = mature_seq[i:i+7]
            pos1 = i + 1  # 1-based position
            is_canonical = (pos1 == 2)
            match_status = "known_seed_match" if kmer in seed_to_ids else "novel_seed"
            matched_ids = ",".join(sorted(seed_to_ids.get(kmer, [])))
            families = ",".join(sorted(seed_to_families.get(kmer, [])))

            if kmer in seed_to_ids:
                matched_positions.append(str(pos1))
                matched_seeds.add(kmer)
                matched_families.update(seed_to_families[kmer])

            long_rows.append({
                "candidate": name,
                "control_type": control_map.get(name, "candidate"),
                "mature_seq": mature_seq,
                "kmer_start_pos_1based": pos1,
                "kmer_7mer": kmer,
                "is_canonical_seed_2_8": "yes" if is_canonical else "no",
                "kmer_status": match_status,
                "matched_known_miRNAs": matched_ids if matched_ids else "-",
                "matched_families": families if families else "-"
            })
            total_7mers += 1

        canonical_status = "known_seed" if canonical_seed in seed_to_ids else "novel_seed"
        canonical_families = ",".join(sorted(seed_to_families.get(canonical_seed, []))) if canonical_seed in seed_to_ids else "-"

        summary_rows.append({
            "candidate": name,
            "control_type": control_map.get(name, "candidate"),
            "mature_seq": mature_seq,
            "canonical_seed_2_8": canonical_seed,
            "canonical_seed_status": canonical_status,
            "canonical_seed_matched_families": canonical_families,
            "n_total_7mers": total_7mers,
            "n_matching_known_7mers": len(matched_seeds),
            "matching_positions": ",".join(matched_positions) if matched_positions else "-",
            "matched_known_7mers": ",".join(sorted(matched_seeds)) if matched_seeds else "-",
            "matched_families_all_positions": ",".join(sorted(matched_families)) if matched_families else "-"
        })

# write long
with out_long.open("w") as o:
    writer = csv.DictWriter(o, fieldnames=list(long_rows[0].keys()), delimiter="\t")
    writer.writeheader()
    writer.writerows(long_rows)

# write summary
with out_summary.open("w") as o:
    writer = csv.DictWriter(o, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
    writer.writeheader()
    writer.writerows(summary_rows)

print("Written:", out_long)
print("Written:", out_summary)
