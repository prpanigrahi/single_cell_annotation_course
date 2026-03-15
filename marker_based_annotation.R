## ============================================================
##  APPROACH A — Overlap scoring across 3 marker databases
##  Sources: Canonical (own list) | PanglaoDB | CellMarker 2.0
##  Output:  Per-cluster consensus annotation with confidence
## ============================================================
##
##  Files needed before running:
##    1. PanglaoDB_markers_27_Mar_2020.tsv  → panglaodb.se/markers.html
##    2. Human_cell_markers.txt             → cellmarker.bioinf.sdu.edu.cn
##
##  Packages:
##    install.packages(c("dplyr","tibble","tidyr","ggplot2","pheatmap"))
## ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


## ============================================================
## BLOCK 1 — Your DGE markers (already filtered and ranked)
## If not below snippet shows some steps of seurat pipeline
## ============================================================
## If you haven't run FindAllMarkers yet:
##
## all_markers <- FindAllMarkers(
##   seurat_obj, only.pos = TRUE,
##   min.pct = 0.25, logfc.threshold = 0.25
## )

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
rm(pbmc)
rm(pbmc.data)

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

# Filter and rank DGE output
top_markers <- all_markers %>%
  filter(
    p_val_adj  < 0.05,
    avg_log2FC > 0.5,
    pct.1      > 0.25,
    pct.2      < 0.50       # specificity filter — key
  ) %>%
  mutate(
    specificity_score = pct.1 - pct.2,
    rank_score        = avg_log2FC * specificity_score
  ) %>%
  group_by(cluster) %>%
  arrange(desc(rank_score), .by_group = TRUE) %>%
  ungroup()

clusters <- sort(unique(top_markers$cluster))
cat("Clusters found:", paste(clusters, collapse = ", "), "\n")

dim(top_markers)


## ============================================================
## BLOCK 2 — Source 1: Canonical marker list (your own)
## ============================================================

canonical_markers <- list(

  # T cells
  "T_cell"                = c("CD3D","CD3E","CD3G","TRAC"),
  "CD4_T"                 = c("CD4","IL7R","MAL","LEF1"),
  "CD8_T"                 = c("CD8A","CD8B","GZMK","GZMA"),
  "T_naive"               = c("CCR7","SELL","TCF7","LEF1"),
  "T_effector_memory"     = c("GZMB","PRF1","IFNG","FASLG"),
  "T_exhausted"           = c("PDCD1","HAVCR2","LAG3","TIGIT","TOX"),
  "Treg"                  = c("FOXP3","IL2RA","CTLA4","TIGIT"),

  # B cells
  "B_cell"                = c("CD19","MS4A1","CD79A","CD79B","IGHM"),
  "Plasma_cell"           = c("MZB1","IGHG1","SDC1","PRDM1","XBP1"),

  # NK
  "NK_cell"               = c("GNLY","NKG7","KLRD1","NCAM1","FCGR3A"),

  # Myeloid
  "Monocyte_classical"    = c("CD14","LYZ","S100A8","S100A9","FCN1"),
  "Monocyte_nonclassical" = c("FCGR3A","MS4A7","CDKN1C","LST1"),
  "Macrophage"            = c("CD68","MRC1","CSF1R","C1QA","C1QB"),
  "DC_myeloid"            = c("FCER1A","CLEC10A","CD1C","CLEC9A"),
  "DC_plasmacytoid"       = c("LILRA4","IL3RA","CLEC4C","GZMB"),
  "Mast_cell"             = c("TPSAB1","CPA3","KIT","MS4A2"),
  "Neutrophil"            = c("FCGR3B","CSF3R","CXCR2","S100A8"),
  "Platelet"              = c("PPBP","PF4","GP1BA","ITGA2B"),

  # Structural / stromal
  "Epithelial"            = c("EPCAM","KRT8","KRT18","KRT19","CDH1"),
  "Fibroblast"            = c("COL1A1","COL1A2","DCN","LUM","FAP"),
  "Endothelial"           = c("PECAM1","CDH5","VWF","ENG","CLDN5"),
  "Pericyte"              = c("ACTA2","RGS5","PDGFRB","MYH11"),

  # Neural
  "Neuron"                = c("RBFOX3","MAP2","SYP","SNAP25","TUBB3"),
  "Astrocyte"             = c("GFAP","AQP4","S100B","ALDH1L1"),
  "Microglia"             = c("CX3CR1","P2RY12","TMEM119","HEXB"),
  "Oligodendrocyte"       = c("MBP","MOG","PLP1","MAG"),

  # Hepatic / other
  "Hepatocyte"            = c("ALB","APOB","TTR","CYP3A4","HNF4A"),
  "Adipocyte"             = c("ADIPOQ","FABP4","LEP","PLIN1","PPARG")
)


## ============================================================
## BLOCK 3 — Source 2: PanglaoDB
## ============================================================

panglaodb_raw <- read.delim(
  "PanglaoDB_markers_27_Mar_2020.tsv",
  stringsAsFactors = FALSE
)

panglaodb_human <- panglaodb_raw %>%
  filter(grepl("Hs", species)) %>%
  select(
    gene      = official.gene.symbol,
    cell_type = cell.type,
    organ     = organ,
    sens      = sensitivity_human,
    spec      = specificity_human
  ) %>%
  filter(!is.na(gene), gene != "")

# Convert to named list — filter by organ relevant to your tissue
# Common organ values: "Blood", "Immune system", "Liver", "Lung", "Brain"
# Set organ_filter = NULL to include all organs

build_panglaodb_list <- function(db,
                                  organ_filter = c("Blood","Immune system"),
                                  min_spec     = 0.0) {
  if (!is.null(organ_filter)) {
    db <- db %>% filter(organ %in% organ_filter)
  }
  if (min_spec > 0) {
    db <- db %>% filter(!is.na(spec), spec >= min_spec)
  }
  marker_list <- db %>%
    group_by(cell_type) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  setNames(marker_list$genes, marker_list$cell_type)
}

panglao_list <- build_panglaodb_list(
  panglaodb_human,
  organ_filter = c("Blood", "Immune system"),  # ← adjust for your tissue
  min_spec     = 0.0
)

cat("PanglaoDB cell types loaded:", length(panglao_list), "\n")


## ============================================================
## BLOCK 4 — Source 3: CellMarker 2.0
## ============================================================

cellmarker_raw <- read.csv(
  "cellmarker/Cell_marker_Human.csv",
  sep = ",", stringsAsFactors = FALSE
)

cellmarker_clean <- cellmarker_raw %>%
  select(tissue = tissue_type, cell_type = cell_name, gene = Symbol) %>%
  filter(!is.na(gene), gene != "", gene != "NA") %>%
  mutate(gene = strsplit(gene, "[,;]+")) %>%
  unnest(gene) %>%
  mutate(gene = trimws(gene)) %>%
  filter(gene != "")

# Convert to named list — filter by tissue
# Common tissue values: "Blood", "Peripheral blood", "Liver", "Lung", "Brain"
# Set tissue_filter = NULL to include all tissues

build_cellmarker_list <- function(db, tissue_filter = "blood") {
  if (!is.null(tissue_filter)) {
    db <- db %>%
      filter(grepl(tissue_filter, tissue, ignore.case = TRUE))
  }
  marker_list <- db %>%
    group_by(cell_type) %>%
    summarise(genes = list(unique(gene)), .groups = "drop")
  setNames(marker_list$genes, marker_list$cell_type)
}

cellmarker_list <- build_cellmarker_list(
  cellmarker_clean,
  tissue_filter = "blood"   # ← adjust for your tissue
)

cat("CellMarker cell types loaded:", length(cellmarker_list), "\n")


## ============================================================
## BLOCK 5 — Core overlap function (used for all 3 sources)
## ============================================================
##
##  For each cell type in the database:
##    - counts how many cluster marker genes match
##    - computes pct_ref_covered  : % of DB markers found in cluster
##    - computes pct_query_covered: % of cluster genes that are DB markers
##    - computes jaccard_index    : overlap / union (best single metric)
##
##  jaccard_index = |A ∩ B| / |A ∪ B|
##  Ranges 0–1. Higher = more similar. Use this for ranking.

overlap_score <- function(cluster_genes,
                           marker_db,
                           top_n = 20) {

  # Use only top N ranked genes from the cluster
  cluster_genes <- head(cluster_genes, top_n)
  n_query       <- length(cluster_genes)

  results <- lapply(names(marker_db), function(cell_type) {
    ref_genes <- unique(marker_db[[cell_type]])
    n_ref     <- length(ref_genes)
    overlap   <- intersect(cluster_genes, ref_genes)
    n_overlap <- length(overlap)

    if (n_overlap == 0) return(NULL)

    union_size   <- length(union(cluster_genes, ref_genes))
    jaccard      <- round(n_overlap / union_size, 4)
    pct_ref      <- round(n_overlap / n_ref    * 100, 1)
    pct_query    <- round(n_overlap / n_query  * 100, 1)

    data.frame(
      cell_type         = cell_type,
      n_overlap         = n_overlap,
      n_ref_markers     = n_ref,
      pct_ref_covered   = pct_ref,    # how much of DB marker set you captured
      pct_query_covered = pct_query,  # how much of your cluster genes match
      jaccard_index     = jaccard,    # best single ranking metric
      matched_genes     = paste(overlap, collapse = ", "),
      stringsAsFactors  = FALSE
    )
  }) %>%
    bind_rows() %>%
    arrange(desc(jaccard_index))

  return(results)
}


## ============================================================
## BLOCK 6 — Run overlap scoring against all 3 sources
## ============================================================

run_overlap_all_clusters <- function(top_markers,
                                      clusters,
                                      marker_db,
                                      source_label,
                                      top_n = 20) {
  cat("Scoring against", source_label, "...\n")

  lapply(clusters, function(cl) {
    genes <- top_markers %>%
      filter(cluster == cl) %>%
      pull(gene)

    result <- overlap_score(genes, marker_db, top_n = top_n)

    if (is.null(result) || nrow(result) == 0) {
      return(data.frame(
        cluster           = cl,
        source            = source_label,
        cell_type         = "Unknown",
        n_overlap         = 0,
        n_ref_markers     = NA,
        pct_ref_covered   = 0,
        pct_query_covered = 0,
        jaccard_index     = 0,
        matched_genes     = "",
        stringsAsFactors  = FALSE
      ))
    }

    result$cluster <- cl
    result$source  <- source_label
    result
  }) %>%
    bind_rows()
}

# Run all three
results_canonical  <- run_overlap_all_clusters(
  top_markers, clusters,
  canonical_markers, "Canonical", top_n = 20
)

results_panglaodb  <- run_overlap_all_clusters(
  top_markers, clusters,
  panglao_list, "PanglaoDB", top_n = 20
)

results_cellmarker <- run_overlap_all_clusters(
  top_markers, clusters,
  cellmarker_list, "CellMarker", top_n = 20
)

# Combine all results into one long table
all_results <- bind_rows(
  results_canonical,
  results_panglaodb,
  results_cellmarker
)

# Save full results
write.csv(all_results, "overlap_all_sources_full.csv", row.names = FALSE)
cat("Full results saved → overlap_all_sources_full.csv\n")


## ============================================================
## BLOCK 7 — Top 3 hits per cluster per source (for inspection)
## ============================================================

top3_per_source <- all_results %>%
  filter(cell_type != "Unknown") %>%
  group_by(cluster, source) %>%
  slice_max(jaccard_index, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, source, desc(jaccard_index))

cat("\n========== TOP 3 HITS PER CLUSTER PER SOURCE ==========\n")
print(top3_per_source %>%
  select(cluster, source, cell_type, n_overlap,
         pct_ref_covered, jaccard_index, matched_genes),
  n = Inf)


## ============================================================
## BLOCK 8 — Consensus: top-1 hit from each source per cluster
## ============================================================
##
##  Strategy:
##    1. Take the best hit (highest jaccard) from each source
##    2. Normalize cell type names for comparison
##       (PanglaoDB and CellMarker use different naming)
##    3. Check agreement across sources → confidence score
##    4. Majority vote → consensus label

# Step 8a: pull top-1 per cluster per source
top1_per_source <- all_results %>%
  filter(cell_type != "Unknown", jaccard_index > 0) %>%
  group_by(cluster, source) %>%
  slice_max(jaccard_index, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, source, cell_type, jaccard_index, matched_genes)


# Step 8b: normalize cell type names
##
##  PanglaoDB uses: "T cells", "B cells", "NK cells"
##  CellMarker uses: "T Cell", "B Cell", "Natural killer cell"
##  Canonical uses:  "T_cell", "B_cell", "NK_cell"
##
##  Map everything to a shared vocabulary for comparison.
##  Extend this lookup table for your cell types.

normalize_cell_type <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)

  # Broad grouping map — add entries as needed
  mapping <- c(
    # T cell family
    "t cell"                 = "T cell",
    "t cells"                = "T cell",
    "cd4 t"                  = "CD4 T cell",
    "cd4 t cell"             = "CD4 T cell",
    "cd4+ t cell"            = "CD4 T cell",
    "cd4 t naive"            = "CD4 T naive",
    "t naive"                = "CD4 T naive",
    "naive t cell"           = "CD4 T naive",
    "cd8 t"                  = "CD8 T cell",
    "cd8 t cell"             = "CD8 T cell",
    "cd8+ t cell"            = "CD8 T cell",
    "cytotoxic t cell"       = "CD8 T cell",
    "regulatory t cell"      = "Treg",
    "treg"                   = "Treg",
    "t effector memory"      = "T effector",
    "exhausted t cell"       = "T exhausted",
    "t exhausted"            = "T exhausted",

    # B cell family
    "b cell"                 = "B cell",
    "b cells"                = "B cell",
    "plasma cell"            = "Plasma cell",
    "plasmablast"            = "Plasma cell",

    # NK
    "nk cell"                = "NK cell",
    "nk cells"               = "NK cell",
    "natural killer cell"    = "NK cell",
    "natural killer cells"   = "NK cell",

    # Monocyte / macrophage
    "monocyte"               = "Monocyte",
    "monocytes"              = "Monocyte",
    "monocyte classical"     = "Monocyte classical",
    "classical monocyte"     = "Monocyte classical",
    "monocyte nonclassical"  = "Monocyte nonclassical",
    "non-classical monocyte" = "Monocyte nonclassical",
    "macrophage"             = "Macrophage",
    "macrophages"            = "Macrophage",

    # DC
    "dendritic cell"         = "DC",
    "dendritic cells"        = "DC",
    "dc myeloid"             = "DC myeloid",
    "myeloid dendritic cell" = "DC myeloid",
    "plasmacytoid dendritic cell" = "DC plasmacytoid",
    "dc plasmacytoid"        = "DC plasmacytoid",
    "pdc"                    = "DC plasmacytoid",

    # Other immune
    "mast cell"              = "Mast cell",
    "mast cells"             = "Mast cell",
    "neutrophil"             = "Neutrophil",
    "neutrophils"            = "Neutrophil",
    "platelet"               = "Platelet",
    "platelets"              = "Platelet",
    "megakaryocyte"          = "Platelet",
    "megakaryocytes"          = "Platelet",

    # Stromal
    "epithelial cell"        = "Epithelial",
    "epithelial"             = "Epithelial",
    "fibroblast"             = "Fibroblast",
    "fibroblasts"            = "Fibroblast",
    "endothelial cell"       = "Endothelial",
    "endothelial"            = "Endothelial",
    "pericyte"               = "Pericyte",

    # Neural
    "neuron"                 = "Neuron",
    "neurons"                = "Neuron",
    "astrocyte"              = "Astrocyte",
    "microglia"              = "Microglia",
    "oligodendrocyte"        = "Oligodendrocyte"
  )

  ifelse(x %in% names(mapping), mapping[x], tools::toTitleCase(x))
}

top1_normalized <- top1_per_source %>%
  mutate(cell_type_norm = normalize_cell_type(cell_type))


# Step 8c: pivot wide — one row per cluster, one column per source
consensus_wide <- top1_normalized %>%
  select(cluster, source, cell_type_norm, jaccard_index) %>%
  pivot_wider(
    names_from  = source,
    values_from = c(cell_type_norm, jaccard_index),
    names_glue  = "{source}_{.value}"
  ) %>%
  arrange(as.numeric(as.character(cluster)))

# Rename for readability
consensus_wide <- consensus_wide %>%
  rename_with(~ gsub("_cell_type_norm", "", .x)) %>%
  rename_with(~ gsub("_jaccard_index", "_jaccard", .x))


# Step 8d: majority vote + confidence
source_cols <- c("Canonical", "PanglaoDB", "CellMarker")

consensus_final <- consensus_wide %>%
  rowwise() %>%
  mutate(
    # collect all three calls into a vector
    all_calls = list(c_across(all_of(source_cols))),

    # majority vote (most common label; ties broken by first)
    vote_table   = list(sort(table(unlist(all_calls)),
                             decreasing = TRUE)),
    consensus    = names(vote_table)[1],
    n_agree      = as.integer(vote_table[1]),

    # mean jaccard across sources that returned a result
    mean_jaccard = mean(c_across(ends_with("_jaccard")),
                        na.rm = TRUE),

    # confidence label
    confidence = case_when(
      n_agree == 3                          ~ "High",
      n_agree == 2 & mean_jaccard >= 0.10  ~ "Medium",
      n_agree == 2 & mean_jaccard <  0.10  ~ "Low",
      n_agree == 1                          ~ "Conflict",
      TRUE                                  ~ "Unknown"
    )
  ) %>%
  select(-vote_table, -all_calls) %>%
  ungroup()


## ============================================================
## BLOCK 9 — Print the final consensus table
## ============================================================

cat("\n========== CONSENSUS ANNOTATION TABLE ==========\n")
consensus_final %>%
  select(
    cluster,
    Canonical, PanglaoDB, CellMarker,   # individual calls
    consensus, n_agree,                  # vote result
    mean_jaccard, confidence             # quality metrics
  ) %>%
  print(n = Inf)

# Flag conflicts for manual review
conflicts <- consensus_final %>%
  filter(confidence %in% c("Conflict", "Low", "Unknown"))

if (nrow(conflicts) > 0) {
  cat("\n⚠ Clusters needing manual review:\n")
  print(conflicts %>% select(cluster, Canonical, PanglaoDB,
                              CellMarker, confidence))
} else {
  cat("\n✓ All clusters reached consensus.\n")
}


## ============================================================
## BLOCK 10 — Visualization
## ============================================================

# --- Plot 1: Jaccard heatmap (top hit per cluster per source) ---
library(pheatmap)

jaccard_mat <- top1_normalized %>%
  select(cluster, source, jaccard_index) %>%
  pivot_wider(names_from = source, values_from = jaccard_index,
              values_fill = 0) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

pheatmap(
  jaccard_mat,
  main          = "Jaccard index — top hit per source",
  color         = colorRampPalette(c("white","#2166ac"))(50),
  display_numbers = TRUE,
  number_format = "%.2f",
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  fontsize       = 10
)


# --- Plot 2: Dot plot — confidence per cluster ---
ggplot(
  consensus_final,
  aes(x = as.factor(cluster),
      y = mean_jaccard,
      color = confidence,
      size  = n_agree)
) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c(
    "High"     = "#1a9641",
    "Medium"   = "#a6d96a",
    "Low"      = "#fdae61",
    "Conflict" = "#d7191c",
    "Unknown"  = "#999999"
  )) +
  scale_size_continuous(range = c(4, 10), breaks = c(1, 2, 3)) +
  labs(
    title  = "Annotation confidence per cluster",
    x      = "Cluster",
    y      = "Mean Jaccard index",
    color  = "Confidence",
    size   = "Sources agree"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# --- Plot 3: Stacked bar — what each source called each cluster ---
source_calls_long <- top1_normalized %>%
  select(cluster, source, cell_type_norm)

ggplot(source_calls_long,
       aes(x = as.factor(cluster), fill = cell_type_norm)) +
  geom_bar(position = "dodge") +
  facet_wrap(~ source, ncol = 1) +
  labs(
    title = "Top cell type call per source per cluster",
    x     = "Cluster", y = "Count", fill = "Cell type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )


## ============================================================
## BLOCK 11 — Apply consensus labels to Seurat object
## ============================================================

# Build the label map from consensus table
label_map <- setNames(
  consensus_final$consensus,
  as.character(consensus_final$cluster)
)

# Apply to Seurat object
Idents(seurat_obj) <- "seurat_clusters"
seurat_obj <- RenameIdents(seurat_obj, label_map)
seurat_obj[["cell_type_consensus"]] <- Idents(seurat_obj)

# Store confidence in metadata too
seurat_obj@meta.data$annotation_confidence <-
  consensus_final$confidence[
    match(seurat_obj@meta.data$seurat_clusters,
          consensus_final$cluster)
  ]

# Final UMAP
DimPlot(
  seurat_obj,
  group.by   = "cell_type_consensus",
  label      = TRUE,
  label.size = 3.5,
  repel      = TRUE,
  pt.size    = 0.4
) +
  ggtitle("Consensus annotation — Canonical + PanglaoDB + CellMarker") +
  theme(legend.position = "right")


## ============================================================
## BLOCK 12 — Save outputs
## ============================================================

write.csv(consensus_final,
          "consensus_annotation_final.csv",
          row.names = FALSE)

write.csv(top3_per_source,
          "top3_hits_per_source.csv",
          row.names = FALSE)

saveRDS(seurat_obj, "seurat_annotated.rds")

cat("\nOutputs saved:\n")
cat("  consensus_annotation_final.csv\n")
cat("  top3_hits_per_source.csv\n")
cat("  overlap_all_sources_full.csv\n")
cat("  seurat_annotated.rds\n")

## ============================================================
## END
## ============================================================
