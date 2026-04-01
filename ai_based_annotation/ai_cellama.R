# AI based
setwd("C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\cellama")


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

DimPlot(seurat_obj, reduction = "umap", label =  T)

# Clusters
levels(seurat_obj)

### Find DGE markers of all clusters
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos    = TRUE,     # positive markers only for annotation
  min.pct     = 0.25,     # expressed in at least 25% of cluster cells
  logfc.threshold = 0.25, # minimum log2 fold change
  test.use    = "wilcox"  # Wilcoxon rank-sum (default, robust)
)

# Quick look at output structure
head(all_markers)

# Total cluster wise marker genes detected
all_markers %>% group_by(cluster) %>% summarise(total_marker_genes = n()) %>% View()


#============================================
# https://github.com/CelVoxes/ceLLama
# 2024 featured on Nature
#============================================

# Installation
#install.packages("devtools")
#devtools::install_github("CelVoxes/ceLLama")

library(ceLLama)
library(tidyverse)
library(httr)

DimPlot(seurat_obj, label = T, label.size = 3) + theme_void() + theme(aspect.ratio = 1)

# split into a lists per cluster
top_markers <- all_markers %>%
  filter(
    p_val_adj  < 0.05,
    avg_log2FC > 0.5,
    pct.1      > 0.25,
    pct.2      < 0.50       # specificity filter — key
  )

pbmc.markers.list <- split(top_markers, top_markers$cluster)

# Increase temperature to diversify outputs. 
# Set different base_prompt to customize annotations.
res <- ceLLama(pbmc.markers.list, 
        n_genes = 20, 
        seed = 101, 
        model = "llama3.2", 
        base_prompt =
        "Act like an expert immunologist and give me the cell type annotation for this cluster. Please, reply the final cell type association and nothing else! If you're not sure just label it as 'unsure'.",
        get_reason = F,
        url = "http://localhost:11434/api/generate",
        temperature = 0.1
        )

res_reason <- ceLLama(pbmc.markers.list, 
               n_genes = 20, 
               seed = 101, 
               model = "llama3.2", 
               base_prompt =
                 "Act like an expert immunologist and give me the cell type annotation for this cluster. Please, reply the final cell type association and nothing else! If you're not sure just label it as 'unsure'.",
               get_reason = T,
               url = "http://localhost:11434/api/generate",
               temperature = 0.1
)

#openai_key="abc"
res1 = ceLLama(pbmc.markers.list, 
               n_genes = 20, 
               seed = 101, 
               use_openai = T, 
               openai_api_key = openai_key,
               base_prompt =
                 "Act like an expert immunologist and give me the cell type annotation for this cluster. Please, reply the final cell type association and nothing else! If you're not sure just label it as 'unsure'.",
               get_reason = F,
               #url = "http://localhost:11434/api/generate",
               temperature = 0.1,
               model = "gpt-4o-mini"
)

  

# transfer the labels
annotations <- lapply(res, FUN = function(x){x[[1]]}) %>% unlist

names(annotations) <-  levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, annotations)
DimPlot(seurat_obj, label = T, label.size = 3) + theme_void() + theme(aspect.ratio = 1)

# These creates html report in the current directory
generate_report_md(res)
create_html_report()

