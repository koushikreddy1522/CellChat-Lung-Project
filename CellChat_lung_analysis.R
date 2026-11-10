library(CellChat)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
getwd()
cellchat <- readRDS("CellChat_lung_results.rds")
length(levels(cellchat@idents))
length(cellchat@netP$pathways)
sum(cellchat@net$count)

# ------------------------------------------------------------
# Summary and interpretation of results
# ------------------------------------------------------------

# The analysis included 9 different cell types.
# A balanced subset of 1,800 nuclei was used,
# with 200 nuclei selected from each cell type.

# CellChat identified 50 significant signalling pathways
# and 2,013 predicted cell-cell interactions.

# These results suggest that several cell populations
# may communicate through different signalling pathways.

# Important:
# These are computational predictions based on gene
# expression data. They do not prove that communication
# actually occurs in the tissue.

# The analysis used a subset of the available data.
# Results may differ if all 39,778 nuclei are analysed.