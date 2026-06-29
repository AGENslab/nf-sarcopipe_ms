# Supplementary Figure 4 – Extended analysis of predicted miRNA–mRNA interactions and pathway enrichment.

This directory contains the scripts used to generate Supplementary Figure 4 of the manuscript.

---

## Supplementary Figure S4A – KEGG enrichment of coherent target genes

### Scripts

- FigureS4A_union_KEGG_enrichment.R

### Output

- KEGG pathway enrichment analysis performed using the union of coherent predicted target genes.

**Note**

This panel presents the complete KEGG enrichment results obtained from the union of coherent targets. Because only a limited number of pathways reached statistical significance, this analysis is provided as supplementary material.

---

## Supplementary Figure S4B – Prioritized predicted target genes

### Scripts

- FigureS4B_prepare_prioritized_targets.py
- FigureS4B_plot_prioritized_targets.R

### Output

- Bubble plot showing prioritized predicted target genes ranked according to interaction strength (number of seed-matched sites) and prediction support.

Bubble size represents the total number of predicted seed-matched interaction sites, while colors indicate whether targets were supported by both miRNA discovery approaches or uniquely identified by miRDeep2.

---

## Notes

Supplementary Figure 4 summarizes complementary downstream analyses of the coherent miRNA–mRNA interaction network. These analyses provide additional biological context and prioritize candidate target genes for future experimental validation.
