# Figure 5 –  Integrated miRNA–mRNA regulatory landscape combining dataset-derived and biologically prioritized predicted targets.

This directory contains the scripts used to generate the analyses presented in Figure 5 of the manuscript.

---

## Figure 5A – Coherent miRNA–mRNA regulation

This panel is a schematic representation of coherent miRNA–mRNA regulation.

The diagram was prepared manually and illustrates the expected inverse relationship between miRNA and mRNA expression after seed-based target prediction.

---

## Figure 5B – Seed-based target prediction summary

### Scripts

- Figure5B_summary_counts.py
- Figure5B_seed_prediction_summary.R

### Output

- Summary of:
  - miRNAs evaluated.
  - Genes with seed matches.
  - Raw predicted interactions.
  - Final coherent miRNA–mRNA pairs.

---

## Figure 5C – Coherent miRNA–mRNA pairs

### Scripts

- Figure5C_prepare_coherent_pairs.py
- Figure5C_plot_coherent_pairs_lollipop.R

### Output

- Lollipop plot summarizing coherent miRNA–mRNA interactions selected after inverse expression filtering.

---

## Figure 5D – Functional enrichment of coherent target genes

### Scripts

- Figure5D_prepare_union_targets.py
- Figure5D_GO_KEGG_enrichment.R
- Figure5D_plot_filtered_GO.R

### Output

- Gene Ontology (GO) and KEGG enrichment analyses of coherent predicted target genes.
- Filtered GO enrichment plot presented in the manuscript.

---

## Figure 5E – Biological prioritization of BrumiR candidate targets in sarcopenia 

### Scripts

- Figure5E_get_sarcopenia_3UTR.R
- Figure5E_seed_matching_sarcopenia_targets.py
- Figure5E_annotate_sarcopenia_targets.py
- Figure5E_prepare_brumir_biological_targets.py
- Figure5E_bubbleplot_biological_targets.R

### Output

- Prioritized biological target genes associated with sarcopenia-related processes.
- Bubble plot summarizing predicted targets grouped by biological category.

---

## Figure 5F – miRNA–mRNA interaction network

### Scripts

- Figure5F_build_network_table.py
- Figure5F_build_node_table.py

### Output

- Cytoscape edge table.
- Cytoscape node table.

The final interaction network was assembled in Cytoscape using the generated node and edge tables. Node layout, colors and graphical adjustments were performed manually for publication.

---

## Notes

This directory contains the downstream scripts used for seed-based target prediction, coherent miRNA–mRNA integration, functional enrichment, biological prioritization and network generation.

The complete workflow is implemented in the nf-Sarcopipe pipeline. This repository contains only the scripts and supporting files required to reproduce the manuscript figures.
