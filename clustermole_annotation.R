## ============================================================
##  clustermole — Marker-based Cell Type Annotation
##  Queries 7 databases simultaneously in one function call
##
##  Package: clustermole (CRAN)
##  Author:  Igor Dolgalev
##  Ref:     https://igordot.github.io/clustermole
##
##  Two annotation methods covered here:
##    Method 1 — clustermole_overlaps()
##               Input : your DGE marker gene list
##               Test  : hypergeometric test (overrepresentation)
##               Output: ranked cell types with p-values
##
##    Method 2 — clustermole_enrichment()
##               Input : average expression matrix per cluster
##               Test  : ssGSEA / GSVA / singscore
##               Output: enrichment scores per cell type per cluster
##
##  Databases queried automatically:
##    ARCHS4     — RNA-seq co-expression derived markers
##    CellMarker — manually curated literature markers
##    MSigDB     — C8 cell type signature gene sets
##    PanglaoDB  — crowdsourced + scRNA-seq derived markers
##    SaVanT     — tissue-specific gene expression signatures
##    TISSUES    — tissue expression database
##    xCell      — cell type enrichment signatures
##
##  Install:
##    install.packages("clustermole")
##    install.packages(c("Seurat","dplyr","ggplot2","tidyr","pheatmap"))
## ============================================================

setwd("C:\\Users\\PriyabrataPanigrahi\\Downloads\\ai\\singlecell\\final\\clustermole")


library(clustermole)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)


## ============================================================
## Explore the clustermole database
## ============================================================
##
##  Before running annotation it is worth understanding
##  what is inside the clustermole database.

# Load all markers (422,000+ rows — all databases, human + mouse)
all_markers_db <- clustermole_markers(species = "hs")

cat("Total entries (human):", nrow(all_markers_db), "\n")

# Structure of the database
# Columns:
#   celltype_full : full name including database source
#   db            : source database (ARCHS4, CellMarker, etc.)
#   species       : hs or mm
#   organ         : tissue / organ
#   celltype      : short cell type name
#   n_genes       : number of marker genes for this cell type
#   gene          : individual gene symbol (one per row)

glimpse(all_markers_db)
dim(all_markers_db)

# Count cell types per database
all_markers_db %>%
  distinct(db, celltype_full) %>%
  count(db, name = "n_cell_types") %>%
  arrange(desc(n_cell_types)) %>%
  print()

# Available databases:
# ARCHS4     — derived from RNA-seq co-expression patterns
# CellMarker — literature curated (same as our Approach A source)
# MSigDB     — C8 collection of scRNA-seq cell type signatures
# PanglaoDB  — same as our Approach A source
# SaVanT     — tissue/cell type gene expression signatures
# TISSUES    — tissue-based expression signatures
# xCell      — 64 immune and stromal cell type signatures

# Browse cell types in a specific database
all_markers_db %>%
  filter(db == "PanglaoDB") %>%
  distinct(celltype, organ) %>%
  arrange(organ, celltype) %>%
  print(n = 30)

# Filter to immune / blood cell types across all databases
immune_types <- all_markers_db %>%
  filter(
    grepl("blood|immune|lymph|T cell|B cell|NK|monocyte|macrophage",
          organ, ignore.case = TRUE) |
    grepl("blood|immune|lymph|T cell|B cell|NK|monocyte|macrophage",
          celltype, ignore.case = TRUE)
  ) %>%
  distinct(db, celltype, organ)

cat("\nImmune-related cell types across databases:",
    nrow(immune_types), "\n")

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

#=====================
# Top ranked marker genes
#=====================
# Filter and rank DGE output w.r.t rank_score
# The filtering logic in Step 1 is the most important thing to get right. 
# Most people filter only on p_val_adj and avg_log2FC, 
# but pct.2 is the real specificity filter. 
# A gene with pct.1 = 0.8 and pct.2 = 0.75 is expressed nearly everywhere — 
# useless for annotation even with a good fold change. 
# The rank_score = avg_log2FC * (pct.1 - pct.2) formula 
# combines both signal strength and specificity into one sortable number.


# Filter and rank 
top_markers <- all_markers %>%
  filter(
    p_val_adj  < 0.05,
    avg_log2FC > 0.5,
    pct.1      > 0.25,
    pct.2      < 0.50          # specificity filter
  ) %>%
  mutate(
    specificity_score = pct.1 - pct.2,
    rank_score        = avg_log2FC * specificity_score
  ) %>%
  group_by(cluster) %>%
  arrange(desc(rank_score), .by_group = TRUE) %>%
  ungroup()

clusters <- sort(unique(top_markers$cluster))
cat("Clusters to annotate:", paste(clusters, collapse = ", "), "\n")

## ============================================================
## METHOD 1 — clustermole_overlaps()
## ============================================================
##
##  What it does:
##    Takes your top N marker genes for one cluster,
##    runs a hypergeometric test against EVERY cell type
##    in ALL 7 databases simultaneously.
##
##  Hypergeometric test asks:
##    "Given the size of your gene list, the size of the
##     reference cell type gene set, and the size of the
##     total gene universe — is the overlap between them
##     larger than expected by chance?"
##
##  Returns: p-value + FDR for every cell type in every DB
##  Lower p-value = better match
##
##  This is DIFFERENT from our Approach A Jaccard:
##    Approach A — Jaccard similarity (no statistics)
##    clustermole — hypergeometric test (p-value based)

# ── Run for all clusters ─────────────────────────────────────

overlaps_all <- lapply(clusters, function(cl) {

  # Get top 20 marker genes for this cluster (ranked)
  genes <- top_markers %>%
    filter(cluster == cl) %>%
    slice_max(rank_score, n = 20) %>%
    pull(gene)

  cat("Cluster", cl, "— querying", length(genes),
      "genes across all databases...\n")

  # Run clustermole_overlaps — queries ALL 7 databases at once
  result <- tryCatch(
    clustermole_overlaps(genes   = genes,
                         species = "hs"),   # hs = human
    error = function(e) {
      warning("Cluster ", cl, " failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(result) || nrow(result) == 0) return(NULL)

  # Add cluster column
  result$cluster <- cl
  result

}) %>%
  bind_rows()

# Quick look at the output structure
View(overlaps_all)

# Columns returned by clustermole_overlaps():
#   celltype_full : full label including DB name
#   db            : source database
#   organ         : tissue / organ
#   celltype      : short cell type name
#   n_genes       : size of reference gene set
#   n_overlap     : number of overlapping genes
#   p_value       : hypergeometric test p-value
#   fdr           : FDR-adjusted p-value (Benjamini-Hochberg)
#   cluster       : your cluster number (added above)

cat("\nTotal overlap results:", nrow(overlaps_all), "\n")


## ============================================================
## BLOCK 3 — Filter and summarise overlaps results
## ============================================================

# Filter to significant results only (FDR < 0.05)
overlaps_sig <- overlaps_all %>%
  dplyr::filter(fdr < 0.05) %>%
  dplyr::arrange(cluster, fdr)

cat("Significant results (FDR < 0.05):", nrow(overlaps_sig), "\n")

# Top 5 hits per cluster across ALL databases
top5_all_db <- overlaps_all %>%
  dplyr::filter(!is.na(fdr)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_min(fdr, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(cluster, db, celltype, organ, overlap, p_value, fdr)

cat("\n========== TOP 5 HITS PER CLUSTER (ALL DBS) ==========\n")
print(top5_all_db, n = Inf)
View(top5_all_db)

# Top 3 hits per cluster PER database
top3_per_db <- overlaps_all %>%
  dplyr::filter(!is.na(fdr)) %>%
  dplyr::group_by(cluster, db) %>%
  dplyr::slice_min(fdr, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::arrange(cluster, db, fdr)

cat("\n========== TOP 3 HITS PER CLUSTER PER DATABASE ==========\n")
print(top3_per_db %>%
        dplyr::select(cluster, db, celltype, overlap, fdr),
  n = Inf)

View(top3_per_db)

# Best hit (top-1) per cluster per database — for consensus
top1_per_db <- overlaps_all %>%
  dplyr::filter(!is.na(fdr), fdr < 0.05) %>%
  dplyr::group_by(cluster, db) %>%
  dplyr::slice_min(fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(cluster, db, celltype, fdr) %>%
  dplyr::arrange(cluster, fdr)

# Wide view — one column per database
top1_wide <- top1_per_db %>%
  dplyr::select(cluster, db, celltype) %>%
  pivot_wider(names_from  = db,
              values_from = celltype,
              values_fill = "—") %>%
  dplyr::arrange(as.numeric(as.character(cluster)))

cat("\n========== TOP-1 PER CLUSTER PER DATABASE (WIDE) ==========\n")
print(top1_wide, n = Inf)
View(top1_wide)

## ============================================================
## Consensus from clustermole overlaps
## ============================================================
##
##  Strategy: across all 7 databases, what is the most
##  common top-hit cell type for each cluster?
##  More databases agreeing = higher confidence.

# Normalize cell type names — same function as Approach A
normalize_cell_type <- function(x) {
  
  # ── Layer 1: Exact match ─────────────────────────────────
  # Maps verbose clustermole labels directly to broad categories.
  # Built from real PBMC clustermole output (47 unique labels).
  
  exact_map <- c(
    
    # T cell — generic
    "Cd3plus t Cell"                            = "T cell",
    "Tlymphocyte"                               = "T cell",
    "t Memory Cells"                            = "T cell",
    "Activated t Cell"                          = "T cell",
    "Hay Bone Marrow Naive t Cell"              = "T cell",
    "Fan Embryonic Ctx Brain Naive Like t Cell" = "T cell",
    "Imgn t 8eff Sp Ot1 D10 Lisova"            = "T cell",
    "Gamma Delta t Cells"                       = "T cell",
    
    # T cell — CD4
    "Naive Cd4plus t Cell"                      = "CD4 T cell",
    "Cd4plus t Cells Hpca 2"                    = "CD4 T cell",
    "Cd4plus Cytotoxic t Cell"                  = "CD4 T cell",
    
    # T cell — CD8
    "Cd8plus Tem Blueprint 3"                   = "CD8 T cell",
    
    # Treg
    "Treg"                                      = "Treg",
    
    # B cell
    "b Lymphocyte Cell Line"                    = "B cell",
    "b Cell (Renal Cell Carcinoma)"             = "B cell",
    "b Cells Naive"                             = "B cell",
    "Hpca b Cells"                              = "B cell",
    "Cd19plus b Cells"                          = "B cell",
    "b Cells Novershtern 2"                     = "B cell",
    "Durante Adult Olfactory Neuroepithelium b Cells" = "B cell",
    
    # NK cell
    "NK cell"                                   = "NK cell",
    "Aizarani Liver C1 Nk Nkt Cells 1"         = "NK cell",
    "Hba Cd56plus Nkcells"                      = "NK cell",
    "Hay Bone Marrow Nk Cells"                  = "NK cell",
    "Large Granular Lymphocyte"                 = "NK cell",
    
    # Monocyte — classical
    "Monocyte"                                  = "Monocyte",
    "Monocytes Blueprint 3"                     = "Monocyte",
    "Monocytes Fantom 1"                        = "Monocyte",
    "Monocytes Novershtern 1"                   = "Monocyte",
    "Hay Bone Marrow Monocyte"                  = "Monocyte",
    
    # Monocyte — non-classical
    "Hpca Monocytes Cd16pluscd14minus"          = "Non-Classical Monocyte",
    
    # Macrophage
    "Macrophage"                                = "Macrophage",
    
    # DC — broad
    "DC"                                        = "DC",
    "Dc Iris 3"                                 = "DC",
    "Durante Adult Olfactory Neuroepithelium Dendritic Cells" = "DC",
    "Hay Bone Marrow Dendritic Cell"            = "DC",
    "Langerhans Cells"                          = "DC",
    
    # DC — myeloid specific
    "Cd1c Cd141 Dendritic Cell"                 = "DC myeloid",
    "Hpca Dendritic Cells Bdca3plus"            = "DC myeloid",
    
    # Platelet / Megakaryocyte
    "Platelet"                                  = "Platelet",
    "Hpca Platelets"                            = "Platelet",
    "Platelets Hpca 1"                          = "Platelet",
    "Zheng Cord Blood C1 Putative Megakaryocyte Progenitor" = "Platelet",
    
    # Too broad to vote with — excluded from consensus
    "Blood"                                     = "Unresolved",
    "Blood Plasma"                              = "Unresolved",
    "Peripheral Blood"                          = "Unresolved",
    "Immune System"                             = "Unresolved"
  )
  
  # ── Layer 2: Keyword / regex fallback ────────────────────
  # Applied only when label is NOT in exact_map above.
  # Handles new labels from future DB updates or other tissues.
  # Order matters — specific patterns checked before broad ones.
  
  keyword_fallback <- function(label) {
    l <- tolower(label)
    
    if (grepl("treg|regulatory t|foxp3",              l)) return("Treg")
    if (grepl("cd8|cytotoxic t|\\bctl\\b|\\btem\\b|
               effector.*t cell",                     l)) return("CD8 T cell")
    if (grepl("cd4|helper t|\\bth1\\b|\\bth2\\b|
               \\bth17\\b|naive.*cd4",                l)) return("CD4 T cell")
    if (grepl("t.cell|t.lymph|tlymph|t.memory|
               t.naive|t.activ|gamma.delta|\\bcd3\\b", l)) return("T cell")
    if (grepl("b.cell|b.lymph|blymph|\\bcd19\\b|
               plasma.cell|plasmablast",              l)) return("B cell")
    if (grepl("\\bnk\\b|natural.killer|\\bcd56\\b|
               \\bgnly\\b|large.granular",            l)) return("NK cell")
    if (grepl("cd16.*mono|non.*classical|cd14minus",  l)) return("Non-Classical Monocyte")
    if (grepl("monocyte|monocytes|\\bcd14\\b",        l)) return("Monocyte")
    if (grepl("macrophage|kupffer|alveolar.macro",    l)) return("Macrophage")
    if (grepl("plasmacytoid|\\bpdc\\b",               l)) return("Plasmacytoid DC")
    if (grepl("\\bdc\\b|dendritic|langerhans",        l)) return("DC")
    if (grepl("platelet|megakaryocyte|\\bppbp\\b",    l)) return("Platelet")
    
    return("Unresolved")
  }
  
  # Apply both layers — exact first, regex fallback if not found
  sapply(x, function(label) {
    if (label %in% names(exact_map)) {
      exact_map[[label]]
    } else {
      keyword_fallback(label)
    }
  }, USE.NAMES = FALSE)
}

normalize_cell_type_old_dontuse <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[_\\-]", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- gsub("\\+", "plus", x)

  mapping <- c(
    # T cell family
    "t cell"="T cell", "t cells"="T cell",
    "cd4 t cell"="CD4 T cell", "cd4 t cells"="CD4 T cell",
    "cd4plus t cell"="CD4 T cell", "cd4plus alpha beta t cell"="CD4 T cell",
    "cd8 t cell"="CD8 T cell", "cd8 t cells"="CD8 T cell",
    "cd8plus t cell"="CD8 T cell", "cd8plus alpha beta t cell"="CD8 T cell",
    "cytotoxic t cell"="CD8 T cell", "cytotoxic t lymphocyte"="CD8 T cell",
    "naive t cell"="Naive T cell", "naive t cells"="Naive T cell",
    "memory t cell"="Memory T cell", "memory t cells"="Memory T cell",
    "regulatory t cell"="Treg", "regulatory t cells"="Treg",
    "treg"="Treg", "t regulatory cell"="Treg",
    "effector t cell"="Effector T cell",
    "exhausted t cell"="Exhausted T cell",

    # B cell family
    "b cell"="B cell", "b cells"="B cell",
    "naive b cell"="B cell naive", "memory b cell"="B cell memory",
    "plasma cell"="Plasma cell", "plasma cells"="Plasma cell",
    "plasmablast"="Plasma cell",

    # NK
    "nk cell"="NK cell", "nk cells"="NK cell",
    "natural killer cell"="NK cell",
    "natural killer cells"="NK cell",

    # Monocyte / macrophage
    "monocyte"="Monocyte", "monocytes"="Monocyte",
    "classical monocyte"="Classical Monocyte",
    "classical monocytes"="Classical Monocyte",
    "cd14plus monocyte"="Classical Monocyte",
    "nonclassical monocyte"="Non-Classical Monocyte",
    "non classical monocyte"="Non-Classical Monocyte",
    "cd16plus monocyte"="Non-Classical Monocyte",
    "macrophage"="Macrophage", "macrophages"="Macrophage",

    # DC
    "dendritic cell"="DC", "dendritic cells"="DC",
    "myeloid dendritic cell"="DC myeloid",
    "plasmacytoid dendritic cell"="Plasmacytoid DC",
    "pdc"="Plasmacytoid DC",

    # Other immune
    "mast cell"="Mast cell", "mast cells"="Mast cell",
    "neutrophil"="Neutrophil", "neutrophils"="Neutrophil",
    "platelet"="Platelet", "platelets"="Platelet",
    "megakaryocyte"="Platelet",
    "basophil"="Basophil", "eosinophil"="Eosinophil",

    # Structural
    "epithelial cell"="Epithelial", "epithelial cells"="Epithelial",
    "fibroblast"="Fibroblast", "fibroblasts"="Fibroblast",
    "endothelial cell"="Endothelial", "endothelial cells"="Endothelial"
  )

  ifelse(x %in% names(mapping),
         mapping[x],
         tools::toTitleCase(x))
}

# Apply normalization
top1_per_db_norm <- top1_per_db %>%
  mutate(celltype_norm = normalize_cell_type(celltype))

# Majority vote across all databases
consensus_clustermole <- top1_per_db_norm %>%
  group_by(cluster) %>%
  summarise(
    # Count votes per normalized cell type
    vote_result  = {
      tbl <- sort(table(celltype_norm), decreasing = TRUE)
      names(tbl)[1]
    },
    n_agree      = {
      tbl <- sort(table(celltype_norm), decreasing = TRUE)
      as.integer(tbl[1])
    },
    n_db_total   = n_distinct(db),  # how many DBs had a sig hit
    dbs_agreeing = {
      tbl <- sort(table(celltype_norm), decreasing = TRUE)
      top_type <- names(tbl)[1]
      paste(db[celltype_norm == top_type], collapse = ", ")
    },
    mean_fdr     = mean(fdr[celltype_norm == vote_result],
                        na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_agree  = round(n_agree / n_db_total * 100, 0),
    confidence = case_when(
      n_agree >= 5                    ~ "Very High",
      n_agree >= 4                    ~ "High",
      n_agree >= 3 & mean_fdr < 0.01 ~ "High",
      n_agree >= 3                    ~ "Medium",
      n_agree == 2                    ~ "Low",
      n_agree == 1                    ~ "Conflict",
      TRUE                            ~ "Unknown"
    )
  ) %>%
  arrange(as.numeric(as.character(cluster)))

cat("\n========== CLUSTERMOLE CONSENSUS ==========\n")
print(consensus_clustermole %>%
        dplyr::select(cluster, vote_result, n_agree, n_db_total,
         pct_agree, mean_fdr, confidence, dbs_agreeing),
  n = Inf)

View(consensus_clustermole)


## ============================================================
## Apply labels to Seurat object
## ============================================================
seurat_obj = pbmc

# Use overlap-based consensus as primary label
label_map <- setNames(
  consensus_clustermole$vote_result,
  as.character(consensus_clustermole$cluster)
)

View(label_map)

Idents(seurat_obj) <- "seurat_clusters"
seurat_obj <- RenameIdents(seurat_obj, label_map)
seurat_obj[["cell_type_clustermole"]] <- Idents(seurat_obj)

seurat_obj@meta.data$clustermole_confidence <-
  consensus_clustermole$confidence[
    match(seurat_obj@meta.data$seurat_clusters,
          consensus_clustermole$cluster)
  ]

seurat_obj@meta.data$clustermole_n_agree <-
  consensus_clustermole$n_agree[
    match(seurat_obj@meta.data$seurat_clusters,
          consensus_clustermole$cluster)
  ]

# UMAP
DimPlot(
  seurat_obj,
  group.by   = "cell_type_clustermole",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3.5,
  pt.size    = 0.4
) +
  ggtitle("clustermole annotation — 7 databases consensus")


