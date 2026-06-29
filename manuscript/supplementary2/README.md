# Supplementary Figure 2 – Structural and sequence-based characterization of de novo miRNA candidates identified by nf-Sarcopipe.

This directory contains the scripts used to generate Supplementary Figure 2 of the manuscript.

---

## Supplementary Figure S2A – Sequence length distribution

### Script

- FigureS2A_sequence_length_distribution.R

### Output

- Comparison of mature miRNA length distributions between BrumiR-RF core candidates and miRDeep2 known miRNAs.

---

## Supplementary Figure S2B – BrumiR2Reference top candidate summary

### Scripts

- FigureS2B_filter_top_BrumiR2Reference_candidates.sh

### Output

- Selection of the top BrumiR2Reference-supported candidates used to assemble the summary table.

**Note**

The final publication table was manually assembled from the filtered BrumiR2Reference output.

---

## Supplementary Figure S2C – Secondary structure visualization of the top BrumiR candidates

### Scripts

- FigureS2C_run_RNAfold_denovo.py
- FigureS2C_parse_RNAfold_denovo.py

### Output

- RNAfold secondary structures.
- Minimum free energy (MFE).
- Dot-bracket secondary structure annotation.

**Note**

Final precursor structure images were rendered with VARNA using the selected candidates and manually assembled for the manuscript.

---

## Supplementary Figure S2D – Seed-based classification of final BrumiR candidates

### Scripts

- FigureS2D_build_known_seed_dictionary.py
- FigureS2D_scan_candidate_7mer_space.py

### Output

- Classification of candidate miRNAs according to canonical and novel 7-mer seed sequences.

---

## Supplementary Figure S2E – Classification of final de novo candidates

### Output

- Summary table describing the final BrumiR de novo candidates and their classification.

**Note**

This table was manually assembled from the seed classification results generated in Supplementary Figure S2D.

---

## Supplementary Figure S2F – 7-mer seed-space heatmap

### Script

- FigureS2F_candidate_seed_heatmap.R

### Output

- Heatmap showing the distribution of canonical and non-canonical 7-mer seed matches across BrumiR candidate miRNAs.

---

## Supplementary Figure S2G – Complete 7-mer seed analysis

### Output

- Complete seed-space profiling table for all evaluated candidate miRNAs.

**Note**

This table was manually assembled from the comprehensive outputs generated during the seed-space analysis.

---

## Notes

Supplementary Figures S2D–S2G correspond to a single seed-space characterization workflow.

The workflow consists of:

1. Construction of a reference dictionary of known miRNA seeds.
2. Exhaustive scanning of all possible candidate 7-mer seeds.
3. Classification of canonical and novel seed sequences.
4. Visualization of seed-space profiles.
5. Manual assembly of publication-ready summary tables from the generated outputs.

Likewise, Supplementary Figures S2B and S2C correspond to the structural validation stage of BrumiR candidates. RNA secondary structures and minimum free energy values were computed automatically, whereas the final publication figures were assembled manually after candidate selection.
