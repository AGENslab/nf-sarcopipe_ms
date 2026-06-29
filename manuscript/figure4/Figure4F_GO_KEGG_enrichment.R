# ================================
# LIBRERÍAS
# ================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)

# ================================
# CARGAR CSV (CLAVE)
# ================================

df <- read.csv("DEG_36_for_networks.csv")

# ================================
# Separar grupos
# ================================
active <- df %>% filter(direction == "Up_in_Active")
sedentary <- df %>% filter(direction == "Up_in_Sedentary")

# Quitar NA
active_genes <- unique(active$ENTREZID[!is.na(active$ENTREZID)])
sedentary_genes <- unique(sedentary$ENTREZID[!is.na(sedentary$ENTREZID)])

# ================================
# GO ACTIVE
# ================================
ego_active <- enrichGO(
  gene = active_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.1,
  readable = TRUE
)

# ================================
# GO SEDENTARY
# ================================
ego_sedentary <- enrichGO(
  gene = sedentary_genes,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.1,
  readable = TRUE
)

# ================================
# PLOTS
# ================================
dotplot(ego_active, showCategory=6, title="Active (immune signature)")

dotplot(ego_sedentary, showCategory=6, title="Sedentary (metabolic signature)")