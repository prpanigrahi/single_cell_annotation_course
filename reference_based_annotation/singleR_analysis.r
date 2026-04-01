#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("scrapper")
# install.packages("viridis")

library(celldex)
library(SingleR)
library(dplyr)
library(Seurat)

# Metadata survey for all available reference datasets in the celldex package.
surveyReferences()

# List the available reference datasets and the associated versions in celldex.
listReferences()

# A SummarizedExperiment object 
ref <- HumanPrimaryCellAtlasData()

# gene*sample. 19363 713 
ref

ref@colData
# DataFrame with 713 rows and 3 columns
# label.main             label.fine   label.ont

ref@colData$label.main %>% table() %>% sort() %>% View()

# 
ref@assays@data$logcounts[1:3,1:3]

## ============================================================
## Run Seurat Pipeline to get top markers of each cluster
## ============================================================

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
rm(pbmc.data)

# Convert Seurat object to SingleCellExperiment for SingleR
sce <- as.SingleCellExperiment(pbmc)

# Run SingleR (Labeling by Cluster)
# clusters = pbmc$seurat_clusters tells SingleR to give one label per cluster
predictions <- SingleR(test = sce, 
                       ref = ref, 
                       labels = ref$label.main, 
                       clusters = pbmc$seurat_clusters)


# View the labels assigned to each cluster
print(predictions$labels)

table(predictions$labels)

plotScoreHeatmap(predictions)
plotScoreHeatmap(predictions, normalize = F)

plotDeltaDistribution(predictions, ncol = 3, show = "delta.next")
plotDeltaDistribution(predictions, ncol = 3, show = "delta.med")

# plotScoreDistribution(results = predictions)

# Create a named vector for mapping
new.cluster.ids <- predictions$labels
names(new.cluster.ids) <- levels(pbmc)

# Rename the clusters in the Seurat object
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Add to metadata for permanent storage
pbmc$singleR_labels <- Idents(pbmc)

# Visualize
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
