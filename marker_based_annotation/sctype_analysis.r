## ============================================================
##  scType — marker-based annotation using ScTypeDB
##  + integration with our Approach A consensus
## ============================================================
##
##  Paper: Ianevski et al., Nature Communications 2022
##  GitHub: github.com/IanevskiAleksandr/sc-type
##
##  What scType does differently from our Approach A:
##    - Uses BOTH positive AND negative markers
##    - Scores the full expression matrix (not just DGE genes)
##    - Score = sum(positive marker expr) - sum(negative marker expr)
##    - Has its own curated DB: ScTypeDB (tissue-organised Excel file)
##
##  Install:
##    install.packages(c("Seurat","HGNChelper","openxlsx","dplyr","ggplot2"))
## ============================================================

setwd("C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\sctype")

library(Seurat)
library(HGNChelper)    # for gene name standardisation
library(openxlsx)      # to read ScTypeDB Excel file
library(dplyr)
library(ggplot2)

#paths
cellmarker_db_path = "C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\marker_data\\cellmarker/Cell_marker_Human.csv"
panglao_db_path= "C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\marker_data\\PanglaoDB_markers_27_Mar_2020.tsv"


## ============================================================
## BLOCK 1 — Load scType functions directly from GitHub
## ============================================================
##
##  scType is not on CRAN or Bioconductor — you source the
##  functions directly from the GitHub repository

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# This loads two functions:
#   gene_sets_prepare()  — loads and processes ScTypeDB for your tissue
#   sctype_score()       — scores your Seurat object against the DB


## ============================================================
## BLOCK 2 — Load ScTypeDB
## ============================================================
##
##  ScTypeDB is an Excel file hosted on GitHub
##  Columns: tissueType | CellName | geneSymbolmore1 (positive) |
##           geneSymbolmore2 (negative)
##
##  Available tissue types in ScTypeDB:
##    Immune system · Pancreas · Liver · Eye · Kidney
##    Brain · Lung · Adrenal · Heart · Intestine
##    Muscle · Placenta · Spleen · Stomach · Thymus

db_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"

# Download locally (recommended — avoids repeated GitHub calls)
if (!file.exists("ScTypeDB_full.xlsx")) {
  download.file(db_url, "ScTypeDB_full.xlsx", mode = "wb")
}

# Check available tissues in the DB
db_raw <- read.xlsx("ScTypeDB_full.xlsx", sheet = 1)
cat("Available tissues in ScTypeDB:\n")
print(unique(db_raw$tissueType))


## ============================================================
## BLOCK 3 — Prepare gene sets for your tissue
## ============================================================

# Set your tissue type — must match exactly what's in ScTypeDB
# Options: "Immune system", "Pancreas", "Liver", "Eye",
#          "Kidney", "Brain", "Lung", "Adrenal",
#          "Heart", "Intestine", "Muscle", "Placenta",
#          "Spleen", "Stomach", "Thymus"

tissue_type <- "Immune system"   # ← change for your tissue

gs_list <- gene_sets_prepare(
  path_to_db_file = "ScTypeDB_full.xlsx",
  cell_type       = tissue_type
)

# gs_list has two components:
cat("Cell types in positive marker list:\n")
print(names(gs_list$gs_positive))

cat("\nSample positive markers for Plasma B cells:\n")
print(gs_list$gs_positive[["Plasma B cells"]])

cat("\nSample negative markers for Plasma B cells:\n")
print(gs_list$gs_negative[["Plasma B cells"]])


## ============================================================
## BLOCK 4 — Score clusters with scType
## ============================================================
##
##  sctype_score() needs:
##    scRNAseqData : scaled expression matrix (genes x cells)
##    scaled       : TRUE if already scaled
##    gs           : positive marker gene sets
##    gs2          : negative marker gene sets

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/PriyabrataPanigrahi/Downloads/ai/singlecell/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Add %mitochondrial gene metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Filtering
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
pbmc <- NormalizeData(pbmc)

# Find top variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scaling genes
pbmc <- ScaleData(pbmc)

# PCA analysis
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc)

# UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

seurat_obj = pbmc
#rm(pbmc)
rm(pbmc.data)

# Extract scaled expression matrix from Seurat
scaled_matrix <- GetAssayData(
  seurat_obj,
  assay = "RNA",
  layer  = "scale.data"    # must be scaled data
)

dim(scaled_matrix)

# Run scType scoring
es_max <- sctype_score(
  scRNAseqData = scaled_matrix,
  scaled       = TRUE,
  gs           = gs_list$gs_positive,
  gs2          = gs_list$gs_negative
)

# es_max is a matrix: cell types x cells
# Each value is the scType score for that cell type in that cell
dim(es_max)      # rows = cell types, cols = cells
head(es_max[,1:5])


## ============================================================
## BLOCK 5 — Assign cell types per cluster
## ============================================================
##
##  For each Seurat cluster:
##    - Sum scType scores across all cells in the cluster
##    - Cell type with highest sum = cluster annotation
##    - Also compute a confidence: top score / (top score + 2nd score)

clusters <- sort(unique(seurat_obj$seurat_clusters))
print(clusters)

sctype_results <- lapply(clusters, function(cl) {
  
  # Get cell barcodes for this cluster
  cell_ids <- WhichCells(seurat_obj,
                         idents = cl)
  
  # Sum scores across cells in cluster
  if (length(cell_ids) == 1) {
    cluster_scores <- es_max[, cell_ids, drop = FALSE]
    scores_sum     <- rowSums(cluster_scores)
  } else {
    cluster_scores <- es_max[, cell_ids]
    scores_sum     <- rowSums(cluster_scores)
  }
  
  # Sort and take top 3
  top3 <- sort(scores_sum, decreasing = TRUE)[1:3]
  
  # Confidence = gap between 1st and 2nd normalised by total
  confidence <- ifelse(
    length(top3) >= 2 && top3[2] > 0,
    round((top3[1] - top3[2]) / abs(top3[1]) * 100, 1),
    100
  )
  
  data.frame(
    cluster      = cl,
    rank1_type   = names(top3)[1],
    rank1_score  = round(top3[1], 2),
    rank2_type   = names(top3)[2],
    rank2_score  = round(top3[2], 2),
    rank3_type   = names(top3)[3],
    rank3_score  = round(top3[3], 2),
    confidence_gap = confidence,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

cat("\n========== scType results per cluster ==========\n")
View(sctype_results)


## ============================================================
## BLOCK 6 — Add scType labels to Seurat object
## ============================================================

# Build label map from top-1 result
sctype_label_map <- setNames(
  sctype_results$rank1_type,
  as.character(sctype_results$cluster)
)

sctype_label_map

# Add labes in seurat metadata
seurat_obj@meta.data$sctype_label <- sctype_label_map[
  as.character(seurat_obj$seurat_clusters)
]

seurat_obj@meta.data$sctype_confidence <- sctype_results$confidence_gap[
  match(seurat_obj$seurat_clusters, sctype_results$cluster)
]


## ============================================================
## BLOCK 7 — Visualize scType results
## ============================================================

# UMAP with scType labels
DimPlot(
  seurat_obj,
  group.by   = "sctype_label",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle(paste("scType annotation —", tissue_type)) +
  theme(legend.position = "right")


# Score heatmap — shows all cell type scores per cluster
library(pheatmap)

# Build cluster-level score matrix
cluster_score_mat <- sapply(clusters, function(cl) {
  cell_ids <- WhichCells(seurat_obj, idents = cl)
  if (length(cell_ids) == 1) {
    rowSums(es_max[, cell_ids, drop = FALSE])
  } else {
    rowSums(es_max[, cell_ids])
  }
})

dim(cluster_score_mat)
colnames(cluster_score_mat) <- paste0("C", clusters)
cluster_score_mat[1:3,1:3]

# Keep only top scoring cell types for readability
top_types <- names(sort(rowMeans(cluster_score_mat),
                        decreasing = TRUE))[1:15]

pheatmap(
  cluster_score_mat[top_types, ],
  main          = "scType scores per cluster",
  color         = colorRampPalette(c("#f7f7f7","#2166ac"))(50),
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  display_numbers = FALSE,
  fontsize_row  = 9,
  fontsize_col  = 10
)

# Bubble plot — score strength per cluster
score_long <- as.data.frame(cluster_score_mat[top_types,]) %>%
  tibble::rownames_to_column("cell_type") %>%
  tidyr::pivot_longer(-cell_type,
                      names_to  = "cluster",
                      values_to = "score") %>%
  filter(score > 0)

ggplot(score_long,
       aes(x = cluster, y = cell_type,
           size = score, color = score)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(1, 8)) +
  scale_color_gradient(low = "#d9f0d3", high = "#1a6634") +
  labs(title = "scType score per cluster",
       x = "Cluster", y = "Cell type",
       size = "Score", color = "Score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ============================================================
## BLOCK 8 — Use ScTypeDB with a CUSTOM marker list
## ============================================================
##
##  You are NOT limited to the built-in ScTypeDB.
##  You can pass any marker list in the same format.
##  This lets you use PanglaoDB / CellMarker / your own list
##  with the scType SCORING ENGINE (pos - neg scoring).
##
##  Format needed:
##    gs_positive: named list of positive marker vectors
##    gs_negative: named list of negative marker vectors (can be empty)

# ── Option A: Use PanglaoDB markers with scType scoring ──────

# Load PanglaoDB (already downloaded)
panglaodb <- read.delim(
  panglao_db_path,
  stringsAsFactors = FALSE
) %>%
  filter(grepl("Hs", species),
         organ %in% c("Blood", "Immune system")) %>%
  dplyr::select(gene = official.gene.symbol,
         cell_type = cell.type,
         spec = specificity_human)

# Build positive marker list from PanglaoDB
panglao_positive <- panglaodb %>%
  group_by(cell_type) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

gs_panglao_positive <- setNames(
  panglao_positive$genes,
  panglao_positive$cell_type
)

head(gs_panglao_positive)

# No negative markers from PanglaoDB (it doesn't have them)
# Create empty negative list of same length
gs_panglao_negative <- lapply(gs_panglao_positive, function(x) c())

# Score with PanglaoDB markers using scType engine
es_panglao <- sctype_score(
  scRNAseqData = scaled_matrix,
  scaled       = TRUE,
  gs           = gs_panglao_positive,
  gs2          = gs_panglao_negative
)

cat("scType scoring with PanglaoDB markers — done\n")

seurat_obj = pbmc

sctype_results_pangalodb <- lapply(clusters, function(cl) {
  
  # Get cell barcodes for this cluster
  cell_ids <- WhichCells(seurat_obj,
                         idents = cl)
  
  # Sum scores across cells in cluster
  if (length(cell_ids) == 1) {
    cluster_scores <- es_panglao[, cell_ids, drop = FALSE]
    scores_sum     <- rowSums(cluster_scores)
  } else {
    cluster_scores <- es_panglao[, cell_ids]
    scores_sum     <- rowSums(cluster_scores)
  }
  
  # Sort and take top 3
  top3 <- sort(scores_sum, decreasing = TRUE)[1:3]
  
  # Confidence = gap between 1st and 2nd normalised by total
  confidence <- ifelse(
    length(top3) >= 2 && top3[2] > 0,
    round((top3[1] - top3[2]) / abs(top3[1]) * 100, 1),
    100
  )
  
  data.frame(
    cluster      = cl,
    rank1_type   = names(top3)[1],
    rank1_score  = round(top3[1], 2),
    rank2_type   = names(top3)[2],
    rank2_score  = round(top3[2], 2),
    rank3_type   = names(top3)[3],
    rank3_score  = round(top3[3], 2),
    confidence_gap = confidence,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

# Build label map from top-1 result
sctype_label_map <- setNames(
  sctype_results_pangalodb$rank1_type,
  as.character(sctype_results_pangalodb$cluster)
)

sctype_label_map

# Add labes in seurat metadata
seurat_obj@meta.data$sctype_label <- sctype_label_map[
  as.character(seurat_obj$seurat_clusters)
]

seurat_obj@meta.data$sctype_confidence <- sctype_results_pangalodb$confidence_gap[
  match(seurat_obj$seurat_clusters, sctype_results_pangalodb$cluster)
]


## ============================================================
## BLOCK 7 — Visualize scType results
## ============================================================

# UMAP with scType labels
DimPlot(
  seurat_obj,
  group.by   = "sctype_label",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle(paste("scType annotation PangaloDB")) +
  theme(legend.position = "right")


# ── Option B: Use CellMarker markers with scType scoring ─────

cellmarker <- read.csv(
  cellmarker_db_path,
  sep = ",", stringsAsFactors = FALSE
) %>%
  dplyr::select(tissue = tissue_type,
         cell_type = cell_name,
         gene = Symbol) %>%
  dplyr::filter(grepl("blood", tissue, ignore.case = TRUE)) %>%
  dplyr::filter(!is.na(gene), gene != "", gene != "NA") %>%
  dplyr::mutate(gene = strsplit(gene, "[,;]+")) %>%
  tidyr::unnest(gene) %>%
  dplyr::mutate(gene = trimws(gene)) %>%
  dplyr::filter(gene != "")

dim(cellmarker)

gs_cm_positive <- cellmarker %>%
  group_by(cell_type) %>%
  summarise(genes = list(unique(gene)), .groups = "drop") %>%
  { setNames(.$genes, .$cell_type) }

gs_cm_negative <- lapply(gs_cm_positive, function(x) c())

print(length(gs_cm_positive))
print(length(gs_cm_negative))

es_cellmarker <- sctype_score(
  scRNAseqData = scaled_matrix,
  scaled       = TRUE,
  gs           = gs_cm_positive,
  gs2          = gs_cm_negative
)

cat("scType scoring with CellMarker markers — done\n")

seurat_obj = pbmc
sctype_results_cellmarker <- lapply(clusters, function(cl) {
  
  # Get cell barcodes for this cluster
  cell_ids <- WhichCells(seurat_obj,
                         idents = cl)
  
  # Sum scores across cells in cluster
  if (length(cell_ids) == 1) {
    cluster_scores <- es_cellmarker[, cell_ids, drop = FALSE]
    scores_sum     <- rowSums(cluster_scores)
  } else {
    cluster_scores <- es_cellmarker[, cell_ids]
    scores_sum     <- rowSums(cluster_scores)
  }
  
  # Sort and take top 3
  top3 <- sort(scores_sum, decreasing = TRUE)[1:3]
  
  # Confidence = gap between 1st and 2nd normalised by total
  confidence <- ifelse(
    length(top3) >= 2 && top3[2] > 0,
    round((top3[1] - top3[2]) / abs(top3[1]) * 100, 1),
    100
  )
  
  data.frame(
    cluster      = cl,
    rank1_type   = names(top3)[1],
    rank1_score  = round(top3[1], 2),
    rank2_type   = names(top3)[2],
    rank2_score  = round(top3[2], 2),
    rank3_type   = names(top3)[3],
    rank3_score  = round(top3[3], 2),
    confidence_gap = confidence,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

# Build label map from top-1 result
sctype_label_map <- setNames(
  sctype_results_cellmarker$rank1_type,
  as.character(sctype_results_cellmarker$cluster)
)

sctype_label_map

# Add labes in seurat metadata
seurat_obj@meta.data$sctype_label <- sctype_label_map[
  as.character(seurat_obj$seurat_clusters)
]

seurat_obj@meta.data$sctype_confidence <- sctype_results_cellmarker$confidence_gap[
  match(seurat_obj$seurat_clusters, sctype_results_cellmarker$cluster)
]


## ============================================================
## BLOCK 7 — Visualize scType results
## ============================================================

# UMAP with scType labels
DimPlot(
  seurat_obj,
  group.by   = "sctype_label",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle(paste("scType annotation Cellmarker")) +
  theme(legend.position = "right")

## ============================================================
## BLOCK 9 — Get cluster-level top hits for each source
## ============================================================

seurat_obj = pbmc

get_top1_per_cluster <- function(es_matrix, clusters, seurat_obj,
                                 source_name) {
  lapply(clusters, function(cl) {
    cell_ids <- WhichCells(seurat_obj, idents = cl)
    if (length(cell_ids) == 1) {
      scores <- rowSums(es_matrix[, cell_ids, drop = FALSE])
    } else {
      scores <- rowSums(es_matrix[, cell_ids])
    }
    top1 <- sort(scores, decreasing = TRUE)[1]
    data.frame(
      cluster    = cl,
      source     = source_name,
      cell_type  = names(top1),
      score      = round(top1, 2),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()
}

top1_sctype    <- get_top1_per_cluster(es_max,        clusters,
                                       seurat_obj, "ScTypeDB")
top1_panglao   <- get_top1_per_cluster(es_panglao,    clusters,
                                       seurat_obj, "PanglaoDB_sctype")
top1_cellmarker <- get_top1_per_cluster(es_cellmarker, clusters,
                                        seurat_obj, "CellMarker_sctype")

# Combine
all_sctype <- bind_rows(top1_sctype, top1_panglao, top1_cellmarker)

# Wide view
sctype_wide <- all_sctype %>%
  dplyr::select(cluster, source, cell_type) %>%
  tidyr::pivot_wider(names_from  = source,
                     values_from = cell_type) %>%
  arrange(as.numeric(as.character(cluster)))

cat("\n========== scType — all sources side by side ==========\n")
print(sctype_wide, n = Inf)
View(sctype_wide)

