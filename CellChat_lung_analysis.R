
# ============================================================
# CellChat analysis of human lung single-cell RNA-seq data
# Dataset: Mendeley Data, DOI: 10.17632/7r2cwbw44m.2
# Extension: Human lung dataset, with 200 cells sampled per cell type
# ============================================================

# 1. Load required packages
library(CellChat)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)

# 2. Create folders for results
dir.create("results", showWarnings = FALSE)
dir.create("project", showWarnings = FALSE)

# 3. Load cell metadata
metadata <- read.csv(
  "project/Metadata_lung_cells_update.csv",
  check.names = TRUE
)

# Check the available cell types
print(table(metadata$Cell.type))

# 4. Select a balanced sample of up to 200 cells per cell type
set.seed(123)

sampled_metadata <- do.call(
  rbind,
  lapply(
    split(metadata, metadata$Cell.type),
    function(x) x[sample(nrow(x), min(200, nrow(x))), ]
  )
)

# 5. Read the count-file header without loading the full dataset
counts_header <- fread(
  "project/Counts_lung_cells.csv",
  nrows = 0,
  check.names = FALSE
)

count_names <- names(counts_header)[-1]

# Match count-file column names to metadata cell IDs
count_ids <- sub(
  "\\.1$", "-1",
  sub("^X", "", count_names)
)

selected_columns <- c(
  "V1",
  count_names[match(sampled_metadata$X, count_ids)]
)

if (anyNA(selected_columns)) {
  stop("Some metadata cell IDs could not be matched to the count file.")
}

# 6. Load only the selected cells from the large count file
counts_subset <- fread(
  "project/Counts_lung_cells.csv",
  select = selected_columns,
  data.table = FALSE
)

# Confirm that the cell IDs match
stopifnot(
  all(
    sub("\\.1$", "-1", sub("^X", "", names(counts_subset)[-1])) ==
      sampled_metadata$X
  )
)

# 7. Create the gene expression matrix
counts_mat <- as.matrix(counts_subset[, -1])

rownames(counts_mat) <- make.unique(
  as.character(counts_subset[[1]])
)

colnames(counts_mat) <- sampled_metadata$X

# Check the matrix
print(dim(counts_mat))
print(sum(is.na(counts_mat)))

# 8. Prepare cell metadata for CellChat
cellchat_meta <- data.frame(
  labels = sampled_metadata$Cell.type,
  samples = sampled_metadata$ID,
  row.names = sampled_metadata$X
)

print(table(cellchat_meta$labels))

# 9. Load the human CellChat database
data("CellChatDB.human")
CellChatDB <- CellChatDB.human

# 10. Normalise expression data and create CellChat object
data.input <- normalizeData(counts_mat)

cellchat <- createCellChat(
  object = data.input,
  meta = cellchat_meta,
  group.by = "labels"
)

cellchat@DB <- CellChatDB.human

# 11. Identify overexpressed genes and interactions
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(
  cellchat,
  do.fast = FALSE
)

# 12. Calculate communication probabilities
cellchat <- computeCommunProb(cellchat)

# 13. Filter low-cell-number communications
cellchat <- filterCommunication(
  cellchat,
  min.cells = 10
)

# 14. Calculate pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)

# 15. Aggregate the communication network
cellchat <- aggregateNet(cellchat)

# 16. Calculate network centrality
cellchat <- netAnalysis_computeCentrality(
  cellchat,
  slot.name = "netP"
)

# 17. Print key results
cat("Number of cell types:",
    length(levels(cellchat@idents)), "\n")

cat("Number of signalling pathways:",
    length(cellchat@netP$pathways), "\n")

cat("Total predicted interactions:",
    sum(cellchat@net$count), "\n")

# 18. Overall cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))

png(
  "results/CellChat_overall_network.png",
  width = 2400,
  height = 2400,
  res = 300
)

netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)

dev.off()

# 19. Signalling role scatter plot
p_role <- netAnalysis_signalingRole_scatter(cellchat)

ggsave(
  "results/CellChat_signaling_role_scatter.png",
  plot = p_role,
  width = 10,
  height = 7,
  dpi = 300
)

# 20. Outgoing signalling heatmap
p_out <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing"
)

png(
  "results/CellChat_outgoing_heatmap_clear.png",
  width = 3000,
  height = 2400,
  res = 300
)

ComplexHeatmap::draw(p_out)
dev.off()

# 21. Incoming signalling heatmap
p_in <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "incoming"
)

png(
  "results/CellChat_incoming_heatmap_clear.png",
  width = 3000,
  height = 2400,
  res = 300
)

ComplexHeatmap::draw(p_in)
dev.off()

# 22. Investigate COLLAGEN signalling
if ("COLLAGEN" %in% cellchat@netP$pathways) {
  
  # COLLAGEN contribution
  p <- netAnalysis_contribution(
    cellchat,
    signaling = "COLLAGEN"
  )
  
  ggsave(
    "results/COLLAGEN_contribution.png",
    plot = p,
    width = 12,
    height = 8,
    units = "in",
    dpi = 300
  )
  
  # COLLAGEN communication network
  png(
    "results/COLLAGEN_signalling.png",
    width = 2400,
    height = 2400,
    res = 300
  )
  
  netVisual_aggregate(
    cellchat,
    signaling = "COLLAGEN",
    layout = "circle"
  )
  
  dev.off()
  
} else {
  message("COLLAGEN pathway was not found in this analysis.")
}

# 23. Save the final CellChat object
saveRDS(
  cellchat,
  file = "CellChat_lung_results.rds"
)

# 24. Save key results to a text file
sink("CellChat_key_results.txt")

cat("CellChat human lung analysis\n")
cat("Number of cell types:",
    length(levels(cellchat@idents)), "\n")
cat("Number of signalling pathways:",
    length(cellchat@netP$pathways), "\n")
cat("Total predicted interactions:",
    sum(cellchat@net$count), "\n")

sink()

cat("CellChat analysis script completed.\n")
