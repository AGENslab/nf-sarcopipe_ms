#!/usr/bin/env Rscript

############################################################
# Figure 5D – Source-specific functional enrichment
#
# Description
# -----------
# Performs GO Biological Process, KEGG, and Reactome
# enrichment separately for:
#
# 1. BrumiR de novo miRNA target genes validated by miRanda
# 2. miRDeep2 annotated miRNA target genes supported by multiMiR
# 3. The integrated union of all validated target genes
#
# The script reads the final validated Cytoscape edge table.
# No genes, pathways, candidate IDs, or counts are hardcoded.
#
# Bubbleplots
# -----------
# Each enrichment database is shown in a single faceted plot:
#
# - BrumiR: purple
# - miRDeep2: fuchsia
# - Integrated: dark grey
#
# Point size represents the number of genes contributing to
# the term. The x-axis represents GeneRatio.
#
# Significant terms (BH-adjusted p < 0.05) are shown at full
# opacity. When a source has no significant terms, its best
# nominal terms are retained at reduced opacity and explicitly
# labeled as non-significant in the plotting table.
#
# Usage
# -----
# Rscript Figure5D_source_specific_enrichment_and_bubbleplots.R \
#   <cytoscape_edges.tsv> \
#   <top_n_per_source> \
#   <output_directory>
############################################################

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    paste0(
      "Usage:\n",
      "Rscript Figure5D_source_specific_enrichment_and_bubbleplots.R ",
      "<cytoscape_edges.tsv> <top_n_per_source> <output_directory>\n"
    ),
    call. = FALSE
  )
}

input_edges <- args[1]
top_n <- suppressWarnings(as.integer(args[2]))
output_dir <- args[3]

if (!file.exists(input_edges)) {
  stop(
    "Input Cytoscape edge table not found: ",
    input_edges,
    call. = FALSE
  )
}

if (is.na(top_n) || top_n < 1) {
  stop(
    "top_n_per_source must be a positive integer.",
    call. = FALSE
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Helper functions
# ============================================================

require_columns <- function(data, required, label) {
  missing_columns <- setdiff(required, colnames(data))

  if (length(missing_columns) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}


safe_ratio_to_numeric <- function(ratio) {
  ratio <- as.character(ratio)

  vapply(
    ratio,
    function(value) {
      parts <- strsplit(value, "/", fixed = TRUE)[[1]]

      if (length(parts) != 2) {
        return(NA_real_)
      }

      numerator <- suppressWarnings(as.numeric(parts[1]))
      denominator <- suppressWarnings(as.numeric(parts[2]))

      if (
        is.na(numerator) ||
        is.na(denominator) ||
        denominator == 0
      ) {
        return(NA_real_)
      }

      numerator / denominator
    },
    numeric(1)
  )
}


map_symbols <- function(symbols) {
  symbols <- unique(
    trimws(
      as.character(symbols)
    )
  )

  symbols <- symbols[
    !is.na(symbols) &
      symbols != ""
  ]

  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c(
      "SYMBOL",
      "ENTREZID"
    )
  ) %>%
    as_tibble() %>%
    filter(
      !is.na(SYMBOL),
      !is.na(ENTREZID),
      SYMBOL != "",
      ENTREZID != ""
    ) %>%
    distinct(
      SYMBOL,
      ENTREZID
    ) %>%
    arrange(
      SYMBOL,
      ENTREZID
    )

  list(
    symbols = symbols,
    mapping = mapping,
    entrez = unique(
      as.character(mapping$ENTREZID)
    ),
    unmapped = setdiff(
      symbols,
      unique(mapping$SYMBOL)
    )
  )
}


empty_enrichment_table <- function() {
  tibble(
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    Count = numeric()
  )
}


result_to_table <- function(result) {
  if (is.null(result)) {
    return(empty_enrichment_table())
  }

  result_df <- as.data.frame(result)

  if (nrow(result_df) == 0) {
    return(empty_enrichment_table())
  }

  as_tibble(result_df)
}


run_go <- function(entrez_ids) {
  if (length(entrez_ids) == 0) {
    return(empty_enrichment_table())
  }

  result <- tryCatch(
    enrichGO(
      gene = entrez_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(error) {
      warning(
        "GO enrichment failed: ",
        conditionMessage(error)
      )
      NULL
    }
  )

  result_to_table(result)
}


run_kegg <- function(entrez_ids) {
  if (length(entrez_ids) == 0) {
    return(empty_enrichment_table())
  }

  result <- tryCatch(
    {
      enriched <- enrichKEGG(
        gene = entrez_ids,
        organism = "hsa",
        keyType = "ncbi-geneid",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )

      if (
        !is.null(enriched) &&
        nrow(as.data.frame(enriched)) > 0
      ) {
        enriched <- setReadable(
          enriched,
          OrgDb = org.Hs.eg.db,
          keyType = "ENTREZID"
        )
      }

      enriched
    },
    error = function(error) {
      warning(
        "KEGG enrichment failed: ",
        conditionMessage(error)
      )
      NULL
    }
  )

  result_to_table(result)
}


run_reactome <- function(entrez_ids) {
  if (length(entrez_ids) == 0) {
    return(empty_enrichment_table())
  }

  result <- tryCatch(
    enrichPathway(
      gene = entrez_ids,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      readable = TRUE
    ),
    error = function(error) {
      warning(
        "Reactome enrichment failed: ",
        conditionMessage(error)
      )
      NULL
    }
  )

  result_to_table(result)
}


add_metadata <- function(
  result_table,
  source_name,
  database_name,
  n_input_genes
) {
  if (nrow(result_table) == 0) {
    return(
      result_table %>%
        mutate(
          Source = character(),
          Database = character(),
          n_input_genes = integer(),
          GeneRatio_numeric = numeric(),
          significance = character(),
          minus_log10_padj = numeric()
        )
    )
  }

  result_table %>%
    mutate(
      Source = source_name,
      Database = database_name,
      n_input_genes = n_input_genes,
      GeneRatio_numeric = safe_ratio_to_numeric(GeneRatio),

      significance = case_when(
        !is.na(p.adjust) &
          p.adjust < 0.05 ~
          "BH-adjusted p < 0.05",

        TRUE ~
          "Not significant"
      ),

      minus_log10_padj = -log10(
        pmax(
          p.adjust,
          .Machine$double.xmin,
          na.rm = TRUE
        )
      )
    )
}


select_terms_for_plot <- function(data, n_terms) {
  if (nrow(data) == 0) {
    return(data)
  }

  significant <- data %>%
    filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    arrange(
      p.adjust,
      desc(Count),
      Description
    ) %>%
    slice_head(
      n = n_terms
    )

  if (nrow(significant) > 0) {
    return(significant)
  }

  data %>%
    arrange(
      p.adjust,
      pvalue,
      desc(Count),
      Description
    ) %>%
    slice_head(
      n = n_terms
    )
}


prepare_plot_table <- function(data, n_terms) {
  sources <- c(
    "BrumiR",
    "miRDeep2",
    "Integrated"
  )

  selected <- lapply(
    sources,
    function(source_name) {
      data %>%
        filter(
          Source == source_name
        ) %>%
        select_terms_for_plot(
          n_terms
        )
    }
  )

  bind_rows(selected) %>%
    mutate(
      Source = factor(
        Source,
        levels = sources
      )
    )
}


make_bubbleplot <- function(
  plot_data,
  title,
  subtitle,
  output_prefix
) {
  if (nrow(plot_data) == 0) {
    warning(
      "No enrichment terms available for plot: ",
      title
    )
    return(invisible(NULL))
  }

  plot_data <- plot_data %>%
    mutate(
      Description_wrapped = str_wrap(
        Description,
        width = 48
      ),

      Description_key = paste(
        Source,
        Description_wrapped,
        sep = "__"
      ),

      plot_alpha = ifelse(
        significance == "BH-adjusted p < 0.05",
        1,
        0.42
      )
    ) %>%
    arrange(
      Source,
      GeneRatio_numeric,
      p.adjust
    )

  description_levels <- unique(
    plot_data$Description_key
  )

  plot_data <- plot_data %>%
    mutate(
      Description_key = factor(
        Description_key,
        levels = description_levels
      )
    )

  description_labels <- setNames(
    plot_data$Description_wrapped,
    plot_data$Description_key
  )

  source_colors <- c(
    "BrumiR" = "#6A3D9A",
    "miRDeep2" = "#E7298A",
    "Integrated" = "#4D4D4D"
  )

  p <- ggplot(
    plot_data,
    aes(
      x = GeneRatio_numeric,
      y = Description_key,
      size = Count,
      color = Source,
      alpha = plot_alpha
    )
  ) +
    geom_point(
      stroke = 0.3
    ) +
    facet_grid(
      Source ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    scale_y_discrete(
      labels = description_labels
    ) +
    scale_x_continuous(
      labels = percent_format(
        accuracy = 1
      ),
      expand = expansion(
        mult = c(0.02, 0.12)
      )
    ) +
    scale_color_manual(
      values = source_colors,
      guide = "none"
    ) +
    scale_alpha_identity() +
    scale_size_continuous(
      range = c(3.5, 10),
      breaks = pretty_breaks(
        n = 4
      ),
      name = "Target genes"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10.5
      ),
      strip.text = element_text(
        face = "bold",
        size = 11
      ),
      strip.background = element_rect(
        fill = "grey96",
        color = "grey60"
      ),
      axis.text.y = element_text(
        size = 9
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  plot_height <- max(
    7,
    2.5 + 0.33 * nrow(plot_data)
  )

  ggsave(
    filename = paste0(
      output_prefix,
      ".png"
    ),
    plot = p,
    width = 11.5,
    height = plot_height,
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = paste0(
      output_prefix,
      ".pdf"
    ),
    plot = p,
    width = 11.5,
    height = plot_height,
    bg = "white"
  )
}


# ============================================================
# Read validated regulatory network
# ============================================================

edges <- read_tsv(
  input_edges,
  show_col_types = FALSE,
  progress = FALSE
)

require_columns(
  edges,
  c(
    "target_node",
    "miRNA_source",
    "validated_by"
  ),
  "Cytoscape edge table"
)

validated_edges <- edges %>%
  filter(
    validated_by %in% c(
      "miRanda",
      "multiMiR"
    )
  )

if (nrow(validated_edges) == 0) {
  stop(
    "No miRanda- or multiMiR-validated interactions found.",
    call. = FALSE
  )
}


# ============================================================
# Construct source-specific gene sets
# ============================================================

brumir_genes <- validated_edges %>%
  filter(
    miRNA_source == "BrumiR",
    validated_by == "miRanda"
  ) %>%
  pull(target_node) %>%
  unique()

mirdeep2_genes <- validated_edges %>%
  filter(
    miRNA_source == "miRDeep2",
    validated_by == "multiMiR"
  ) %>%
  pull(target_node) %>%
  unique()

integrated_genes <- validated_edges %>%
  pull(target_node) %>%
  unique()

gene_sets <- list(
  BrumiR = brumir_genes,
  miRDeep2 = mirdeep2_genes,
  Integrated = integrated_genes
)


# ============================================================
# Map genes and export source-specific lists
# ============================================================

mapping_results <- list()
mapping_summary <- list()

for (source_name in names(gene_sets)) {
  mapped <- map_symbols(
    gene_sets[[source_name]]
  )

  mapping_results[[source_name]] <- mapped

  write_tsv(
    tibble(
      SYMBOL = sort(
        mapped$symbols
      )
    ),
    file.path(
      output_dir,
      paste0(
        "Figure5D_",
        source_name,
        "_target_genes.tsv"
      )
    )
  )

  write_tsv(
    mapped$mapping,
    file.path(
      output_dir,
      paste0(
        "Figure5D_",
        source_name,
        "_gene_mapping.tsv"
      )
    )
  )

  mapping_summary[[source_name]] <- tibble(
    Source = source_name,
    input_genes = length(
      mapped$symbols
    ),
    mapped_genes = length(
      unique(mapped$mapping$SYMBOL)
    ),
    unmapped_genes = length(
      mapped$unmapped
    ),
    unmapped_symbols = ifelse(
      length(mapped$unmapped) > 0,
      paste(
        mapped$unmapped,
        collapse = ";"
      ),
      ""
    )
  )
}

mapping_summary_table <- bind_rows(
  mapping_summary
)

write_tsv(
  mapping_summary_table,
  file.path(
    output_dir,
    "Figure5D_source_specific_gene_summary.tsv"
  )
)


# ============================================================
# Run enrichment independently for each source
# ============================================================

go_tables <- list()
kegg_tables <- list()
reactome_tables <- list()

for (source_name in names(mapping_results)) {
  mapped <- mapping_results[[source_name]]

  go_tables[[source_name]] <- add_metadata(
    run_go(
      mapped$entrez
    ),
    source_name,
    "GO Biological Process",
    length(mapped$symbols)
  )

  kegg_tables[[source_name]] <- add_metadata(
    run_kegg(
      mapped$entrez
    ),
    source_name,
    "KEGG",
    length(mapped$symbols)
  )

  reactome_tables[[source_name]] <- add_metadata(
    run_reactome(
      mapped$entrez
    ),
    source_name,
    "Reactome",
    length(mapped$symbols)
  )
}

go_combined <- bind_rows(
  go_tables
)

kegg_combined <- bind_rows(
  kegg_tables
)

reactome_combined <- bind_rows(
  reactome_tables
)


# ============================================================
# Export complete combined enrichment tables
# ============================================================

go_output <- file.path(
  output_dir,
  "Figure5D_GO_source_specific.tsv"
)

kegg_output <- file.path(
  output_dir,
  "FigureS4A_KEGG_source_specific.tsv"
)

reactome_output <- file.path(
  output_dir,
  "FigureS4B_Reactome_source_specific.tsv"
)

write_tsv(
  go_combined,
  go_output
)

write_tsv(
  kegg_combined,
  kegg_output
)

write_tsv(
  reactome_combined,
  reactome_output
)


# ============================================================
# Prepare top terms and plotting tables
# ============================================================

go_plot_table <- prepare_plot_table(
  go_combined,
  top_n
)

kegg_plot_table <- prepare_plot_table(
  kegg_combined,
  top_n
)

reactome_plot_table <- prepare_plot_table(
  reactome_combined,
  top_n
)

write_tsv(
  go_plot_table,
  file.path(
    output_dir,
    "Figure5D_GO_source_specific_plotting_table.tsv"
  )
)

write_tsv(
  kegg_plot_table,
  file.path(
    output_dir,
    "FigureS4A_KEGG_source_specific_plotting_table.tsv"
  )
)

write_tsv(
  reactome_plot_table,
  file.path(
    output_dir,
    "FigureS4B_Reactome_source_specific_plotting_table.tsv"
  )
)


# ============================================================
# Generate bubbleplots
# ============================================================

make_bubbleplot(
  go_plot_table,
  title = "Source-specific GO Biological Process enrichment",
  subtitle = paste0(
    "Validated targets of de novo, annotated, and integrated ",
    "miRNA regulatory networks"
  ),
  output_prefix = file.path(
    output_dir,
    "Figure5D_GO_source_specific_bubbleplot"
  )
)

make_bubbleplot(
  kegg_plot_table,
  title = "Source-specific KEGG pathway enrichment",
  subtitle = paste0(
    "Validated targets of de novo, annotated, and integrated ",
    "miRNA regulatory networks"
  ),
  output_prefix = file.path(
    output_dir,
    "FigureS4A_KEGG_source_specific_bubbleplot"
  )
)

make_bubbleplot(
  reactome_plot_table,
  title = "Source-specific Reactome pathway enrichment",
  subtitle = paste0(
    "Validated targets of de novo, annotated, and integrated ",
    "miRNA regulatory networks"
  ),
  output_prefix = file.path(
    output_dir,
    "FigureS4B_Reactome_source_specific_bubbleplot"
  )
)


# ============================================================
# Final summary
# ============================================================

count_summary <- bind_rows(
  go_combined %>%
    mutate(
      Analysis = "GO BP"
    ),

  kegg_combined %>%
    mutate(
      Analysis = "KEGG"
    ),

  reactome_combined %>%
    mutate(
      Analysis = "Reactome"
    )
) %>%
  group_by(
    Analysis,
    Source
  ) %>%
  summarise(
    total_terms = n(),
    significant_terms = sum(
      !is.na(p.adjust) &
        p.adjust < 0.05
    ),
    .groups = "drop"
  )

write_tsv(
  count_summary,
  file.path(
    output_dir,
    "Figure5D_source_specific_enrichment_summary.tsv"
  )
)

cat("===== SOURCE-SPECIFIC ENRICHMENT =====\n")
print(mapping_summary_table)

cat("\n===== ENRICHMENT TERM SUMMARY =====\n")
print(count_summary)

cat("\nOutput directory:", output_dir, "\n")
