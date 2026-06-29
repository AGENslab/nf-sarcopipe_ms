# Supplementary Figure 3 – Evidence of confounding effects and their correction in RNA-seq analysis.

This directory contains the scripts used to generate Supplementary Figure 3 of the manuscript.

---

## Supplementary Figure S3A – PCA of the 5000 most variable genes

### Script

- FigureS3AB_raw_mRNAseq_DESeq2_PCA_volcano.R

### Output

- Principal component analysis (PCA) of the 5000 most variable genes before confounding correction.

---

## Supplementary Figure S3B – Differential expression without confounding correction

### Script

- FigureS3AB_raw_mRNAseq_DESeq2_PCA_volcano.R

### Output

- Volcano plot of differentially expressed genes obtained from the initial DESeq2 analysis prior to correction for tissue-driven confounding.

---

## Supplementary Figure S3C – PCA of the corrected 36-gene signature

### Script

- FigureS3C_PCA_36_corrected_genes.R

### Output

- PCA performed using the final set of 36 genes retained after linear model correction with PC1 and PC2.

---

## Supplementary Figure S3D – Sensitivity analysis of principal component correction

### Scripts

- FigureS3D_model_PC1_PC2.R
- FigureS3D_model_PC1_PC2_PC3.R
- FigureS3D_model_PC1_PC2_PC3_PC4.R
- FigureS3D_sensitivity_barplot.R

### Output

- Comparison of the number of significantly differentially expressed genes after progressively including additional principal components in the linear model.

---

## Notes

Supplementary Figures S3A and S3B correspond to the initial exploratory RNA-seq analysis performed directly after the nf-core/rnaseq pipeline and DESeq2 differential expression analysis, before correction for tissue-driven confounding.

Supplementary Figures S3C and S3D document the sensitivity analyses performed after identifying tissue as the major source of variation. These analyses support the selection of the PC1 + PC2 correction model used throughout the manuscript by demonstrating that adding additional principal components substantially reduces the number of significant genes without providing additional biological interpretability.
