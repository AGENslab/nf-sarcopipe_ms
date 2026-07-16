# Supplementary Figure 4 – Extended pathway enrichment and validation of musculoskeletal-associated miRNA targets

This directory contains the scripts used to generate Supplementary Figure 4 of the manuscript.

---

## Supplementary Figure S4A – KEGG pathway enrichment

### Scripts

- FigureS4AB_prepare_source_specific_enrichment.R
- FigureS4AB_plot_source_specific_KEGG_Reactome_enrichment.R

### Output

KEGG pathway enrichment analysis performed separately for validated target genes regulated by de novo and annotated miRNAs.

This panel provides the complete KEGG enrichment results complementary to the Gene Ontology analysis presented in Figure 5D.

---

## Supplementary Figure S4B – Reactome pathway enrichment

### Scripts

- FigureS4AB_prepare_source_specific_enrichment.R
- FigureS4AB_plot_source_specific_KEGG_Reactome_enrichment.R

### Output

Reactome pathway enrichment analysis performed separately for validated target genes regulated by de novo and annotated miRNAs.

Although no pathways remained significant after multiple-testing correction, the highest-ranking pathways provide additional biological context for the coherent regulatory networks reconstructed by nf-Sarcopipe.

---

## Supplementary Figure S4C – Validated interactions with sarcopenia-associated genes

### Scripts

- FigureS4C_prepare_validated_sarcopenia_pairs.py
- FigureS4C_plot_validated_sarcopenia_pairs.R

### Output

Summary of validated miRNA–mRNA interactions involving a curated set of genes associated with sarcopenia and skeletal muscle biology.

Interactions are grouped according to their biological function (muscle growth/atrophy, senescence/damage and extracellular matrix remodeling) and display:

- de novo and annotated miRNAs;
- absolute miRNA log2 fold change;
- number of predicted seed-matched sites;
- validation support by miRanda or multiMiR.

---

## Notes

Supplementary Figure 4 provides complementary analyses supporting the integrated miRNA–mRNA regulatory framework presented in Figure 5.

The KEGG and Reactome enrichment analyses were generated using the same source-specific enrichment workflow used for Figure 5D, whereas Supplementary Figure 4C summarizes validated interactions obtained from the sarcopenia-focused analysis presented in Figure 5F.
