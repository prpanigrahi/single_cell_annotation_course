## ============================================================
##  EasyCellType — Marker-based Cell Type Annotation
##
##  Paper : Li et al., Bioinformatics Advances 2023
##          doi: 10.1093/bioadv/vbad029
##  Package: Bioconductor
##  Shiny  : biostatistics.mdanderson.org/shinyapps/EasyCellType
##
##  What it does:
##    Takes your FindAllMarkers() output, converts gene symbols
##    to Entrez IDs, then runs either GSEA or Fisher's exact test
##    against CellMarker / PanglaoDB / clustermole databases.
##
##    Two output modes:
##      Hard classification — one best label per cluster
##      Soft classification — top N candidate labels per cluster
##
##    Built-in plots:
##      Dot plot  — top 5 candidates across all clusters
##      Bar plots — ranked candidates within each cluster
##
##  Databases supported:
##      "cellmarker"  — recommended, best accuracy for most tissues
##      "panglaodb"   — good for tissues not in CellMarker
##      "clustermole" — broader coverage fallback
##
##  Install:
##    if (!require("BiocManager")) install.packages("BiocManager")
##    BiocManager::install("EasyCellType")
##    BiocManager::install("org.Hs.eg.db")   # human gene ID mapping
##    BiocManager::install("org.Mm.eg.db")   # mouse (if needed)
## ============================================================

setwd("C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\easycelltype")
library(EasyCellType)
library(org.Hs.eg.db)      # for gene symbol → Entrez ID conversion
library(AnnotationDbi)
library(Seurat)
library(dplyr)
library(ggplot2)


## ============================================================
## BLOCK 1 — Prepare input from FindAllMarkers output
## ============================================================
##
##  EasyCellType requires a very specific input format:
##    Column 1: gene     — Entrez ID (NOT gene symbol)
##    Column 2: cluster  — cluster number
##    Column 3: score    — expression score (avg_log2FC)
##
##  Within each cluster, genes MUST be sorted by score descending.
##  The function uses rank order — sorting is not optional.
##
##  Important: EasyCellType accepts gene symbols OR Entrez IDs
##  depending on the gene_type parameter. We show both approaches.

# ── Option A: Use gene symbols directly (simplest) ───────────
# gene_type = "symbol" — no conversion needed

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "C:/Users/PriyabrataPanigrahi/Downloads/ai/singlecell/pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc)

# PCA analysis
# Default npcs = 50
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc@meta.data$mt_class = pbmc@meta.data$percent.mt > 2
### 
pbmc@graphs
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# Computing nearest neighbor graph
# Computing SNN
pbmc@graphs
pbmc <- FindClusters(pbmc)
head(pbmc@meta.data)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

summary(pbmc@meta.data$seurat_clusters)

# UMAP
pbmc@reductions
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

pbmc@reductions

seurat_obj = pbmc
rm(pbmc.data)

### Find DGE markers of all clusters
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos    = TRUE,     # positive markers only for annotation
  min.pct     = 0.25,     # expressed in at least 25% of cluster cells
  logfc.threshold = 0.25, # minimum log2 fold change
  test.use    = "wilcox"  # Wilcoxon rank-sum (default, robust)
)

# Start from your filtered markers
top_markers_ect <- all_markers %>%
  filter(
    p_val_adj  < 0.05,
    avg_log2FC > 0.25,    # EasyCellType uses a slightly broader
    pct.1      > 0.10     # threshold than our Approach A
  ) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()

# Build input data frame — symbol version
input_symbol <- data.frame(
  gene    = top_markers_ect$gene,
  cluster = top_markers_ect$cluster,
  score   = top_markers_ect$avg_log2FC,
  stringsAsFactors = FALSE
)

dim(input_symbol)

# Confirm structure
head(input_symbol)
cat("Clusters:", unique(input_symbol$cluster), "\n")
cat("Total genes:", nrow(input_symbol), "\n")


# ── Option B: Convert to Entrez IDs (more accurate matching) ─
# gene_type = "Entrezid" — required for GSEA method

# Map gene symbols to Entrez IDs using org.Hs.eg.db
entrez_map <- mapIds(
  org.Hs.eg.db,
  keys      = top_markers_ect$gene,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"         # take first match if multiple
)

# Add Entrez IDs to marker table
top_markers_ect$entrezid <- entrez_map[top_markers_ect$gene]

# Check conversion rate
n_mapped   <- sum(!is.na(top_markers_ect$entrezid))
n_total    <- nrow(top_markers_ect)
cat(sprintf("Gene symbol → Entrez mapping: %d/%d (%.1f%%)\n",
            n_mapped, n_total, n_mapped/n_total*100))

# Genes that didn't map — usually non-coding or pseudogenes
unmapped <- top_markers_ect %>%
  filter(is.na(entrezid)) %>%
  pull(gene) %>%
  unique()
cat("Unmapped genes (sample):", head(unmapped, 10), "\n")

# Build input data frame — Entrez version
# Remove unmapped genes, sort by score within cluster
input_entrez <- top_markers_ect %>%
  filter(!is.na(entrezid)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup() 

input_entrez = input_entrez %>% dplyr::select(
    gene    = entrezid,   # Entrez ID in column 1
    cluster = cluster,    # cluster in column 2
    score   = avg_log2FC  # expression score in column 3
  )

input_entrez <- as.data.frame(input_entrez)

# Confirm — must be data.frame not tibble
class(input_entrez)
head(input_entrez)


## ============================================================
## BLOCK 2 — Check available tissues in each database
## ============================================================
##
##  Before running, check which tissue names are valid.
##  You must pass tissue names EXACTLY as they appear in the DB.

# CellMarker tissue options
data(cellmarker_tissue)
cat("\n=== CELLMARKER TISSUES ===\n")
print(cellmarker_tissue)

# PanglaoDB tissue options
data(panglao_tissue)
cat("\n=== PANGLAODB TISSUES ===\n")
print(panglao_tissue)

# clustermole tissue options
data(clustermole_tissue)
cat("\n=== CLUSTERMOLE TISSUES ===\n")
print(clustermole_tissue)

# For PBMC — relevant CellMarker tissue options:
# "Blood", "Peripheral blood", "Blood vessel",
# "Umbilical cord blood", "Venous blood"


## ============================================================
## BLOCK 3 — Run annotation — Method 1: GSEA
## ============================================================
##
##  GSEA approach:
##    Ranks your marker genes by avg_log2FC within each cluster,
##    then runs GSEA against database gene sets.
##    Tests whether high-scoring genes are enriched in a
##    cell type's signature gene set.
##
##  Parameters:
##    db        : "cellmarker", "panglaodb", or "clustermole"
##    species   : "Human" or "Mouse"
##    tissue    : vector of tissue names from the DB
##    p_cut     : p-value cutoff for reporting results (0.3 default)
##    test      : "GSEA" or "fisher"

cat("Running EasyCellType GSEA — CellMarker...\n")

annot_GSEA <- easyct(
  data  = input_entrez,          # Entrez ID input
  db       = "cellmarker",          # database
  species  = "Human",
  tissue   = c("Blood",
               "Peripheral blood",
               "Blood vessel",
               "Umbilical cord blood",
               "Venous blood"),
  p_cut    = 0.3,                   # relaxed threshold for exploration
  test     = "GSEA"
)

# Output is a list with results per cluster
class(annot_GSEA)
length(annot_GSEA)

# Peek at results for cluster 0
View(annot_GSEA[[3]]@result)

# Columns in GSEA results:
#   ID          : cell type name
#   enrichmentScore : GSEA enrichment score
#   NES         : normalised enrichment score
#   pvalue      : raw p-value
#   p.adjust    : BH-adjusted p-value
#   core_enrichment : marker genes driving the enrichment


## ============================================================
## BLOCK 4 — Run annotation — Method 2: Fisher's exact test
## ============================================================
##
##  Fisher approach:
##    For each cluster × cell type pair: builds a 2×2
##    contingency table of your markers vs DB markers.
##    Tests whether the overlap is greater than chance.
##    More intuitive — similar logic to our Approach A
##    but adds statistical significance.
##
##  Note: Fisher test works with gene symbols (gene_type="symbol")
##        which avoids the Entrez conversion step.

cat("Running EasyCellType Fisher — CellMarker...\n")

annot_Fisher <- easyct(
  data   = input_symbol,         # symbol input works for Fisher
  db        = "cellmarker",
  species   = "Human",
  tissue    = c("Blood",
                "Peripheral blood",
                "Blood vessel",
                "Umbilical cord blood",
                "Venous blood"),
  p_cut     = 0.3,
  test      = "fisher",
  genetype = "symbol"              # specify symbol input
)

# Peek at Fisher results for cluster 0
head(annot_Fisher[[1]])
View(annot_Fisher[[3]])

# Columns in Fisher results:
#   ID          : cell type name
#   p.value     : raw Fisher p-value
#   p.adjust    : BH-adjusted p-value
#   overlap     : number of overlapping genes
#   overlap_gene: actual overlapping gene symbols


## ============================================================
## BLOCK 5 — Run against all 3 databases for comparison
## ============================================================

# ── PanglaoDB ────────────────────────────────────────────────
# First check which PanglaoDB tissue to use
# data(panglao_tissue) — look for "Blood", "Immune system" etc.

cat("Running EasyCellType Fisher — PanglaoDB...\n")

annot_panglao <- easyct(
  data   = input_symbol,
  db        = "panglao",
  species   = "Human",
  tissue    = c("Blood", "Immune system"),  # PanglaoDB tissue names
  p_cut     = 0.3,
  test      = "fisher",
  genetype = "symbol"
)

# ── clustermole ──────────────────────────────────────────────
cat("Running EasyCellType Fisher — clustermole...\n")

annot_clustermole <- easyct(
  data   = input_symbol,
  db        = "clustermole",
  species   = "Human",
  tissue    = NULL,       # clustermole has no tissue filter in EasyCellType
  p_cut     = 0.3,
  test      = "fisher",
  genetype = "symbol"
)


## ============================================================
## BLOCK 6 — Extract results: hard and soft classification
## ============================================================
##
##  Hard classification — best label per cluster (top-1)
##  Soft classification — top N candidates per cluster

# ── Hard classification (top-1 per cluster) ──────────────────
get_hard_labels <- function(annot_result, source_name) {
  lapply(seq_along(annot_result), function(i) {
    res <- annot_result[[i]]
    if (is.null(res) || nrow(res) == 0) {
      return(data.frame(cluster = i - 1, cell_type = "Unknown",
                        score = NA, source = source_name))
    }
    # Top hit — lowest p.adjust or highest NES
    top1 <- if ("NES" %in% names(res)) {
      res %>% arrange(desc(NES)) %>% slice(1)
    } else {
      res %>% arrange(p_adjust) %>% dplyr::slice(1)
    }
    data.frame(
      cluster   = i - 1,
      cell_type = top1$cellName,
      p_adj     = round(top1$p_adjust, 4),
      source    = source_name,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
}

hard_cellmarker  <- get_hard_labels(annot_Fisher,     "CellMarker")
hard_panglao     <- get_hard_labels(annot_panglao,    "PanglaoDB")
hard_clustermole <- get_hard_labels(annot_clustermole,"clustermole")

# Side by side comparison
hard_comparison <- hard_cellmarker %>%
  dplyr::rename(CellMarker = cell_type, CM_padj = p_adj) %>%
  dplyr::left_join(hard_panglao %>%
              dplyr::rename(PanglaoDB = cell_type, PG_padj = p_adj),
            by = "cluster") %>%
  dplyr::left_join(hard_clustermole %>%
                     dplyr::rename(clustermole = cell_type, CL_padj = p_adj),
            by = "cluster") %>%
  dplyr::select(cluster, CellMarker, PanglaoDB, clustermole,
         CM_padj, PG_padj, CL_padj)

cat("\n=== HARD CLASSIFICATION — ALL 3 DATABASES ===\n")
View(hard_comparison)

# ── Soft classification (top 5 per cluster) ──────────────────
get_soft_labels <- function(annot_result, source_name, top_n = 5) {
  lapply(seq_along(annot_result), function(i) {
    res <- annot_result[[i]]
    if (is.null(res) || nrow(res) == 0) return(NULL)
    top <- if ("NES" %in% names(res)) {
      res %>% dplyr::arrange(desc(NES)) %>% dplyr::slice_head(n = top_n)
    } else {
      res %>% dplyr::arrange(p_adjust) %>% dplyr::slice_head(n = top_n)
    }
    top$cluster <- i - 1
    top$source  <- source_name
    top
  }) %>% bind_rows()
}

soft_cellmarker <- get_soft_labels(annot_Fisher, "CellMarker", top_n = 5)

cat("\n=== SOFT CLASSIFICATION — TOP 5 PER CLUSTER (CellMarker) ===\n")
soft_cellmarker %>%
        dplyr::select(cluster, cellName, p_adjust, source) %>% View()


## ============================================================
## BLOCK 7 — Built-in EasyCellType plots
## ============================================================
##
##  EasyCellType generates two built-in plot types:
##    plot_dot()  — dot plot: top 5 candidates across all clusters
##    plot_bar()  — bar plots: ranked candidates within each cluster
##
##  These use the raw annotation result objects directly.

# ── Plot 1: Dot plot — CellMarker GSEA results ───────────────
#
#  X = clusters, Y = cell type candidates
#  Dot size = number of marker genes in gene set
#  Dot colour = normalised enrichment score (NES) or p-value

p_dot_GSEA <- plot_dot(annot_GSEA, test = "GSEA")
print(p_dot_GSEA)
ggsave("easycelltype_dot_GSEA_cellmarker.pdf",
       p_dot_GSEA, width = 12, height = 8)

# Dot plot — Fisher results
p_dot_Fisher <- plot_dot(annot_Fisher, test = "fisher")
print(p_dot_Fisher)
ggsave("easycelltype_dot_Fisher_cellmarker.pdf",
       p_dot_Fisher, width = 12, height = 8)


# ── Plot 2: Bar plots — one per cluster ──────────────────────
#
#  Each bar plot shows ranked cell type candidates for one cluster
#  Bar length = significance score
#  Shows top 5 candidates

p_bar <- plot_bar(annot_Fisher, test = "fisher")
print(p_bar)
ggsave("easycelltype_bar_Fisher_cellmarker.pdf",
       p_bar, width = 14, height = 10)


## ============================================================
## BLOCK 8 — Custom downstream plots
## ============================================================

# ── Plot 3: Heatmap of top candidate p-values ────────────────
library(pheatmap)

# Build matrix: clusters × top cell types (p.adjust values)
top_types <- soft_cellmarker %>%
  count(cellName, sort = TRUE) %>%
  slice_head(n = 15) %>%
  pull(cellName)

top_types

pval_mat <- soft_cellmarker %>%
  dplyr::filter(cellName %in% top_types) %>%
  dplyr::select(cluster, cellName, p_adjust) %>%
  dplyr::mutate(neg_log_p = -log10(p_adjust + 1e-10)) %>%
  dplyr::select(cluster, cellName, neg_log_p) %>%
  tidyr::pivot_wider(names_from  = cluster,
                     values_from = neg_log_p,
                     values_fill = 0) %>%
  tibble::column_to_rownames("cellName") %>%
  as.matrix()

View(pval_mat)

pheatmap(
  pval_mat,
  main          = "EasyCellType — -log10(p.adjust) per cluster",
  color         = colorRampPalette(c("white","#d73027","#a50026"))(50),
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  display_numbers = TRUE,
  number_format = "%.1f",
  fontsize_row  = 9,
  fontsize_col  = 9
)


# ── Plot 4: Confidence plot — top hit significance per cluster ─
hard_sig <- hard_cellmarker %>%
  dplyr::mutate(
    neg_log_padj = -log10(p_adj + 1e-10),
    confidence   = case_when(
      p_adj < 0.01  ~ "High",
      p_adj < 0.05  ~ "Medium",
      p_adj < 0.3   ~ "Low",
      TRUE          ~ "Not significant"
    )
  )

ggplot(hard_sig,
       aes(x    = as.factor(cluster),
           y    = neg_log_padj,
           fill = confidence)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = cell_type),
            angle = 45, hjust = 0,
            vjust = 0.5, size = 3,
            color = "black") +
  scale_fill_manual(values = c(
    "High"            = "#1a9641",
    "Medium"          = "#a6d96a",
    "Low"             = "#fdae61",
    "Not significant" = "#d7191c"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(
    title = "EasyCellType — annotation confidence per cluster",
    x     = "Cluster",
    y     = "-log10(p.adjust)",
    fill  = "Confidence"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ── Plot 5: Comparison dot plot across 3 databases ────────────
# Shows how consistent predictions are across CellMarker,
# PanglaoDB and clustermole for each cluster

compare_long <- bind_rows(
  hard_cellmarker  %>% dplyr::mutate(source = "CellMarker"),
  hard_panglao     %>% dplyr::mutate(source = "PanglaoDB"),
  hard_clustermole %>% dplyr::mutate(source = "clustermole")
) %>%
  dplyr::mutate(
    neg_log_padj = -log10(p_adj + 1e-10),
    cluster      = as.factor(cluster)
  )

View(compare_long)

ggplot(compare_long,
       aes(x     = cluster,
           y     = source,
           size  = neg_log_padj,
           color = cell_type)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = cell_type),
            size = 2.5, vjust = -1.2,
            color = "black", check_overlap = TRUE) +
  scale_size_continuous(range = c(3, 10),
                        name  = "-log10(p.adj)") +
  labs(
    title  = "EasyCellType — top prediction per cluster per database",
    x      = "Cluster",
    y      = "Database",
    color  = "Cell type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.text     = element_text(size = 8)
  )


## ============================================================
## BLOCK 9 — Apply labels to Seurat object
## ============================================================

# Use CellMarker Fisher results as primary labels
# (EasyCellType recommends CellMarker for best accuracy)

label_map <- setNames(
  hard_cellmarker$cell_type,
  as.character(hard_cellmarker$cluster)
)

Idents(seurat_obj) <- "seurat_clusters"
seurat_obj <- RenameIdents(seurat_obj, label_map)
seurat_obj[["cell_type_easycelltype"]] <- Idents(seurat_obj)

# Store p.adjust as confidence metric
seurat_obj@meta.data$ect_confidence <-
  hard_cellmarker$p_adj[
    match(seurat_obj@meta.data$seurat_clusters,
          hard_cellmarker$cluster)
  ]

# Final UMAP
DimPlot(
  seurat_obj,
  group.by   = "cell_type_easycelltype",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle("EasyCellType annotation — CellMarker Fisher test")


## ============================================================
## BLOCK 10 — Save outputs
## ============================================================

write.csv(hard_comparison,
          "easycelltype_hard_classification.csv",  row.names = FALSE)
write.csv(soft_cellmarker,
          "easycelltype_soft_classification.csv",  row.names = FALSE)
saveRDS(seurat_obj, "seurat_easycelltype_annotated.rds")

cat("\nOutputs saved:\n")
cat("  easycelltype_hard_classification.csv\n")
cat("  easycelltype_soft_classification.csv\n")
cat("  easycelltype_dot_GSEA_cellmarker.pdf\n")
cat("  easycelltype_dot_Fisher_cellmarker.pdf\n")
cat("  easycelltype_bar_Fisher_cellmarker.pdf\n")
cat("  seurat_easycelltype_annotated.rds\n")

## ============================================================
## SHINY APP — no-code alternative
## ============================================================
##
##  If you prefer a point-and-click interface:
##  biostatistics.mdanderson.org/shinyapps/EasyCellType/
##
##  Steps:
##    1. Export your markers as CSV:
##       write.csv(input_symbol, "markers_for_easycelltype.csv",
##                 row.names = FALSE)
##    2. Upload to the Shiny app
##    3. Select database, species, tissue, test method
##    4. Download dot plot and bar plots directly
##
##  Or launch locally:
##  shiny::runApp(system.file("app", package = "EasyCellType"))

## ============================================================
## END
## ============================================================