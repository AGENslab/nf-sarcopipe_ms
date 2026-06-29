#!/usr/bin/env Rscript

# ------------------------------------------------------------
# This script summarizes differentially expressed miRNAs from
# DESeq2 results (BrumiR and miRDeep2) and generates a stacked
# bar plot showing the number of miRNAs upregulated in Active
# and Sedentary groups.
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(readr)

# input args
args <- commandArgs(trailingOnly = TRUE)

brumir_file <- args[1]
mirdeep_file <- args[2]
out_plot <- args[3]

# load data
brumir <- read_csv(brumir_file)
mirdeep <- read_csv(mirdeep_file)

# summarize DE results
summarize_de <- function(df, label){
  df %>%
    filter(!is.na(padj), padj < 0.05) %>%
    mutate(direction = case_when(
      log2FoldChange > 0 ~ "Up in Active",
      log2FoldChange < 0 ~ "Up in Sedentary"
    )) %>%
    count(direction) %>%
    mutate(method = label)
}

# process both datasets
brumir_sum <- summarize_de(brumir, "BrumiR")
mirdeep_sum <- summarize_de(mirdeep, "miRDeep2")

df_plot <- bind_rows(brumir_sum, mirdeep_sum)

# plot
p <- ggplot(df_plot, aes(x=method, y=n, fill=direction)) +
  geom_bar(stat="identity") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white")
  ) +
  scale_fill_manual(values=c(
    "Up in Active"="blue",
    "Up in Sedentary"="red"
  )) +
  labs(
    title="Differentially expressed miRNAs",
    x="Method",
    y="Count",
    fill="miRNA regulation"
  )

# save plot with white background
ggsave(out_plot, p, width=6, height=5, bg = "white")
