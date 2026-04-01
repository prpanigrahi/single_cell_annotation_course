# Install from CRAN (recommended)
#install.packages("mLLMCelltype")

# Or install development version from GitHub
#devtools::install_github("cafferychen777/mLLMCelltype", subdir = "R")

setwd("C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\mllmcelltype")

# openrouter key
# okey="sk-or-v1-"

# Load required packages
library(mLLMCelltype)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot) # Added for plot_grid

DimPlot(pbmc, label = T, label.size = 3)


# Set up cache directory to speed up processing
cache_dir <- "./mllmcelltype_cache"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

top_markers %>% colnames()

c1_marker = top_markers %>% filter(cluster == 0) %>% head(20)
 
# Run LLMCelltype annotation with multiple LLM models
consensus_results <- interactive_consensus_annotation(
  input = c1_marker,
  tissue_name = "human PBMC",  # provide tissue context
  models = c(
 #   "deepseek/deepseek-r1",
#    "openai/gpt-oss-20b:free"
    "qwen/qwen3.6-plus-preview:free",
#    "minimax/minimax-m2.5:free",
"qwen/qwen3-next-80b-a3b-instruct:free"
      ),
  api_keys = list(
    openrouter = okey,
    openai = "your-openai-key",
    gemini = "your-google-key",
    qwen = "your-qwen-key"
  ),
  top_gene_count = 10,
  controversy_threshold = 1.0,
  entropy_threshold = 1.0,
  cache_dir = cache_dir
)

print("Available fields in consensus_results:")
print(names(consensus_results))
consensus_results$initial_results
consensus_results$final_annotations
consensus_results$discussion_logs
consensus_results$discussion_results

cluster_to_celltype_map <- consensus_results$final_annotations
