# Figure 4 – Correction of tissue-driven confounding reveals a robust transcriptional signature associated with physical activity.

This directory contains the scripts used to generate the analyses presented in Figure 4 of the manuscript.

---

## Figure 4A – Public RNA-seq dataset

This panel is a schematic representation of the public RNA sequencing datasets used in this study.

The diagram was prepared manually and summarizes:

- Active group (whole blood samples from skaters).
- Sedentary group (vastus lateralis muscle samples from healthy young women).
- RNA-seq processing using the nf-core/rnaseq pipeline.
- Transcriptomic analysis performed to provide functional context for the miRNA analyses.

---

## Figure 4B – Linear model corrected by PC1 and PC2

### Scripts

- Figure4B_linear_model_PC1_PC2.R
- Figure4B_volcano_PC1_PC2.R

### Output

- Linear model correcting for the first two principal components.
- Volcano plot showing differentially expressed genes after confounding correction.

---

## Figure 4C – Progressive reduction of differentially expressed genes

### Script

- Figure4C_gene_reduction_barplot.R

### Output

- Bar plot summarizing the number of genes after each filtering step:
  - Initial differential expression.
  - Significant genes.
  - Final gene set after PC1 + PC2 correction.

---

## Figure 4D – Heatmap of corrected genes

### Script

- Figure4D_heatmap_36_genes.R

### Output

- Heatmap showing the expression profile of the 36 corrected genes.

---

## Figure 4E – Direction and magnitude of differential expression

### Script

- Figure4E_barplot_36_DE_genes.R

### Output

- Bar plot showing log2 fold change and regulation direction of the 36 corrected genes.

---

## Figure 4F – Functional enrichment analysis

### Scripts

- Figure4F_GO_KEGG_enrichment.R
- Figure4F_GO_KEGG_manual_curation.R
- Figure4F_GO_KEGG_curated_final.R

### Output

- Gene Ontology (GO) enrichment analysis.
- KEGG pathway enrichment analysis.
- Curated functional enrichment figure presented in the manuscript.

---

## Notes

RNA-seq preprocessing, alignment and quantification were performed using the nf-core/rnaseq pipeline.

The scripts included in this directory correspond to the downstream statistical analyses, visualization and functional interpretation used to generate Figure 4 of the manuscript.
