# Single Cell RNA-seq — Marker-Based Cell Type Annotation

A reproducible pipeline for annotating single-cell clusters using marker gene overlap scoring across three curated databases, with consensus voting and confidence scoring.

---

## Table of Contents

- [Overview](#overview)
- [Why Marker-Based Annotation](#why-marker-based-annotation)
- [Databases Used](#databases-used)
- [Pipeline Overview](#pipeline-overview)
- [Jaccard Index — The Ranking Metric](#jaccard-index--the-ranking-metric)
- [Consensus Strategy](#consensus-strategy)
- [Prerequisites](#prerequisites)
- [Input Requirements](#input-requirements)
- [Code Walkthrough](#code-walkthrough)
  - [Block 1 — Filter and Rank DGE Markers](#block-1--filter-and-rank-dge-markers)
  - [Block 2 — Source 1: Canonical Marker List](#block-2--source-1-canonical-marker-list)
  - [Block 3 — Source 2: PanglaoDB](#block-3--source-2-panglaodb)
  - [Block 4 — Source 3: CellMarker 2.0](#block-4--source-3-cellmarker-20)
  - [Block 5 — Core Overlap Function](#block-5--core-overlap-function)
  - [Block 6 — Run Across All Clusters and Sources](#block-6--run-across-all-clusters-and-sources)
  - [Block 7 — Top 3 Hits Per Cluster Per Source](#block-7--top-3-hits-per-cluster-per-source)
  - [Block 8 — Normalize Names and Majority Vote](#block-8--normalize-names-and-majority-vote)
  - [Block 9 — Print Consensus Table](#block-9--print-consensus-table)
  - [Block 10 — Visualization](#block-10--visualization)
  - [Block 11 — Apply Labels to Seurat Object](#block-11--apply-labels-to-seurat-object)
  - [Block 12 — Save Outputs](#block-12--save-outputs)
- [Output Files](#output-files)
- [Interpreting Results](#interpreting-results)
- [Confidence Levels](#confidence-levels)
- [Adapting for Non-PBMC Tissues](#adapting-for-non-pbmc-tissues)
- [References](#references)

---

## Overview

After clustering single-cell RNA-seq data with Seurat and running differential gene expression (DGE) with `FindAllMarkers()`, each cluster has a number — not a biological identity. The goal of annotation is to assign a meaningful cell type label to each cluster.

This pipeline takes your `FindAllMarkers` output and annotates every cluster by comparing its top marker genes against three independent marker databases:

1. **Canonical markers** — your own curated list of well-established marker genes
2. **PanglaoDB** — a crowdsourced database of single-cell marker genes with sensitivity and specificity scores
3. **CellMarker 2.0** — a manually curated literature-based marker database with 36,000+ entries

For each cluster, every database produces a ranked list of candidate cell types. The pipeline then normalises naming differences across databases and runs a **majority vote** to produce a single consensus label per cluster, along with a confidence score.

The key design principle is that **no single database is trusted alone**. When all three sources agree, annotation confidence is high. When they disagree, the cluster is flagged for manual review.

---

## Why Marker-Based Annotation

| Approach | What it needs | Best when |
|---|---|---|
| Marker-based (this pipeline) | DGE output + marker databases | You want interpretable, auditable results |
| Reference-based (SingleR, Azimuth) | A labelled reference atlas | A well-matched atlas exists for your tissue |
| AI-assisted (GPTCelltype, CASSIA) | Marker list + LLM API | You want a quick second opinion |

Marker-based annotation is the most interpretable approach — you can see exactly which genes drove each decision through the `matched_genes` column. It does not require a reference dataset, making it applicable to any tissue.

---

## Databases Used

### Canonical Marker List

A manually curated list maintained in the script covering ~30 immune, stromal, neural, and hepatic cell types. Each cell type has 4–7 well-established marker genes. This list is fully customisable — add or remove cell types to match your tissue.

### PanglaoDB

> Franzén O, Gan LM, Björkegren JLM. PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. *Database* 2019. doi: 10.1093/database/baz046

- Sourced from crowdsourced contributions and automated extraction from published scRNA-seq datasets
- Includes per-marker sensitivity and specificity scores for human
- Download: [panglaodb.se/markers.html](https://panglaodb.se/markers.html)

### CellMarker 2.0

> Zhang X et al. CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools. *Nucleic Acids Research* 2023. doi: 10.1093/nar/gkac947

- Manually curated from published literature
- 36,000+ entries covering human and mouse
- Organised by tissue type and cancer type
- Download: [cellmarker.bioinf.sdu.edu.cn](http://cellmarker.bioinf.sdu.edu.cn/CellMarker2/)

---

## Pipeline Overview

```
FindAllMarkers() output
        │
        ▼
  Filter + Rank DGE
  (p_val_adj, log2FC, pct.1, pct.2)
        │
        ├─────────────────────────────────────┐
        │                                     │
   Top 20 genes                          Top 20 genes
   per cluster                           per cluster
        │                                     │
        ▼                                     ▼
  overlap_score()              overlap_score()        overlap_score()
  vs Canonical                 vs PanglaoDB            vs CellMarker
        │                           │                       │
        └───────────────────────────┴───────────────────────┘
                                    │
                              all_results
                         (long table, all scores)
                                    │
                         ┌──────────┴──────────┐
                         ▼                     ▼
                   Top 3 per source      Top 1 per source
                   (inspection)          (for voting)
                                               │
                                    normalize_cell_type()
                                    (shared vocabulary)
                                               │
                                        Majority vote
                                               │
                                    consensus + confidence
                                               │
                              ┌────────────────┴──────────────────┐
                              ▼                ▼                   ▼
                     Seurat object      CSV outputs           Visualizations
                     (labelled)
```

---

## Jaccard Index — The Ranking Metric

The Jaccard index measures the similarity between your cluster's marker genes and a database's reference gene set for a given cell type:

```
         | cluster genes ∩ reference genes |
J =  ─────────────────────────────────────────
         | cluster genes ∪ reference genes |
```

**Why Jaccard over simple overlap count?**

A raw overlap count rewards large gene sets unfairly. Jaccard normalises by the union — penalising both cases where your query is too broad and where the reference set is very large.

**Practical interpretation for scRNA-seq:**

| Jaccard | Interpretation | Action |
|---|---|---|
| > 0.25 | Strong match | Accept confidently |
| 0.15 – 0.25 | Good match — typical for curated DBs | Accept |
| 0.10 – 0.15 | Moderate signal | Cross-check with FeaturePlot |
| < 0.10 | Weak — possible noise | Do not accept without visual check |
| 0.00 | No overlap | Definite non-match |

> **Note:** Jaccard scores in scRNA-seq annotation are typically in the 0.10–0.35 range — this is normal. Databases include protein markers and bulk RNA markers that may not appear due to dropout. Scores of 0.8+ are not expected.

---

## Consensus Strategy

For each cluster, each database independently nominates a top cell type. The pipeline then:

1. Normalises all cell type names to a shared vocabulary (e.g. "T cells", "T Cell", "T_cell" → "T cell")
2. Performs a majority vote across the three sources
3. Computes `n_agree` (how many sources agree) and `mean_jaccard` (average signal strength)
4. Assigns a confidence level based on those two values

```
Cluster 3 example:
  Canonical  → "CD8 T cell"   ✓
  PanglaoDB  → "CD8 T cell"   ✓
  CellMarker → "CD8 T cell"   ✓
  consensus  = "CD8 T cell"   n_agree = 3   confidence = High

Cluster 0 example:
  Canonical  → "CD4 T cell"   ✓
  PanglaoDB  → "Naive T cell" (subtype disagreement)
  CellMarker → "CD4 T cell"   ✓
  consensus  = "CD4 T cell"   n_agree = 2   confidence = Medium
  → check naive markers (CCR7, SELL, TCF7) manually
```

---

## Prerequisites

```r
# CRAN packages
install.packages(c(
  "Seurat",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "pheatmap"
))
```

**Data files required:**

| File | Source |
|---|---|
| `PanglaoDB_markers_27_Mar_2020.tsv` | [panglaodb.se/markers.html](https://panglaodb.se/markers.html) |
| `Human_cell_markers.txt` | [cellmarker.bioinf.sdu.edu.cn](http://cellmarker.bioinf.sdu.edu.cn/CellMarker2/) → Download → Human cell markers |

---

## Input Requirements

The pipeline starts from `FindAllMarkers()` output. If you have not run it yet:

```r
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)
```

Expected columns in `all_markers`:

| Column | Description |
|---|---|
| `p_val_adj` | Adjusted p-value |
| `avg_log2FC` | Average log2 fold change |
| `pct.1` | Fraction of cells in cluster expressing the gene |
| `pct.2` | Fraction of cells outside cluster expressing the gene |
| `cluster` | Cluster identity |
| `gene` | Gene symbol |

---

## Code Walkthrough

### Source 1: Canonical Marker List

A hand-curated named list in R covering ~30 cell types. Each entry maps a cell type name to a vector of well-established marker genes. This list is tissue-agnostic — remove entries irrelevant to your tissue (e.g. remove `Hepatocyte` when analysing PBMC).

```r
canonical_markers <- list(
  "T_cell"                = c("CD3D","CD3E","CD3G","TRAC"),
  "CD4_T"                 = c("CD4","IL7R","MAL","LEF1"),
  "CD8_T"                 = c("CD8A","CD8B","GZMK","GZMA"),
  "T_naive"               = c("CCR7","SELL","TCF7","LEF1"),
  "T_effector_memory"     = c("GZMB","PRF1","IFNG","FASLG"),
  "T_exhausted"           = c("PDCD1","HAVCR2","LAG3","TIGIT","TOX"),
  "Treg"                  = c("FOXP3","IL2RA","CTLA4","TIGIT"),
  "B_cell"                = c("CD19","MS4A1","CD79A","CD79B","IGHM"),
  "Plasma_cell"           = c("MZB1","IGHG1","SDC1","PRDM1","XBP1"),
  "NK_cell"               = c("GNLY","NKG7","KLRD1","NCAM1","FCGR3A"),
  "Monocyte_classical"    = c("CD14","LYZ","S100A8","S100A9","FCN1"),
  "Monocyte_nonclassical" = c("FCGR3A","MS4A7","CDKN1C","LST1"),
  "Macrophage"            = c("CD68","MRC1","CSF1R","C1QA","C1QB"),
  "DC_myeloid"            = c("FCER1A","CLEC10A","CD1C","CLEC9A"),
  "DC_plasmacytoid"       = c("LILRA4","IL3RA","CLEC4C","GZMB"),
  "Mast_cell"             = c("TPSAB1","CPA3","KIT","MS4A2"),
  "Neutrophil"            = c("FCGR3B","CSF3R","CXCR2","S100A8"),
  "Platelet"              = c("PPBP","PF4","GP1BA","ITGA2B"),
  "Epithelial"            = c("EPCAM","KRT8","KRT18","KRT19","CDH1"),
  "Fibroblast"            = c("COL1A1","COL1A2","DCN","LUM","FAP"),
  "Endothelial"           = c("PECAM1","CDH5","VWF","ENG","CLDN5"),
  "Pericyte"              = c("ACTA2","RGS5","PDGFRB","MYH11"),
  "Neuron"                = c("RBFOX3","MAP2","SYP","SNAP25","TUBB3"),
  "Astrocyte"             = c("GFAP","AQP4","S100B","ALDH1L1"),
  "Microglia"             = c("CX3CR1","P2RY12","TMEM119","HEXB"),
  "Oligodendrocyte"       = c("MBP","MOG","PLP1","MAG"),
  "Hepatocyte"            = c("ALB","APOB","TTR","CYP3A4","HNF4A"),
  "Adipocyte"             = c("ADIPOQ","FABP4","LEP","PLIN1","PPARG")
)
```

---

### Source 2: PanglaoDB

PanglaoDB is loaded from the downloaded TSV file, filtered to human markers, and converted to the same named list format as the canonical markers. The `organ_filter` argument restricts the database to cell types relevant to your tissue.

```r
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

# Adjust organ_filter for your tissue:
# Blood / Immune system | Liver | Lung | Brain | Kidney
panglao_list <- build_panglaodb_list(
  panglaodb_human,
  organ_filter = c("Blood", "Immune system"),
  min_spec     = 0.0
)
```

---

### Source 3: CellMarker 2.0

CellMarker is loaded from the downloaded text file. Gene entries that are comma-separated within a single field are split into individual rows. The `tissue_filter` argument restricts to your tissue type.

```r
cellmarker_raw <- read.csv(
  "Human_cell_markers.txt",
  sep = "\t", stringsAsFactors = FALSE
)

cellmarker_clean <- cellmarker_raw %>%
  select(tissue = tissueType, cell_type = cellName, gene = geneSymbol) %>%
  filter(!is.na(gene), gene != "", gene != "NA") %>%
  mutate(gene = strsplit(gene, "[,;]+")) %>%
  unnest(gene) %>%
  mutate(gene = trimws(gene)) %>%
  filter(gene != "")

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

# Adjust tissue_filter: "blood" | "liver" | "lung" | "brain"
cellmarker_list <- build_cellmarker_list(
  cellmarker_clean,
  tissue_filter = "blood"
)
```

---

### Filter and Rank DGE Markers

Before querying any database, we filter the DGE output to keep only meaningful, specific markers.

**Filtering criteria:**

| Filter | Threshold | Reason |
|---|---|---|
| `p_val_adj` | < 0.05 | Statistical significance |
| `avg_log2FC` | > 0.5 | Minimum expression difference |
| `pct.1` | > 0.25 | Expressed in at least 25% of cluster cells |
| `pct.2` | < 0.50 | **Specificity filter** — not expressed in majority of other cells |

The `pct.2` filter is the most important and most commonly overlooked. A gene with high fold change but `pct.2 = 0.70` is expressed almost everywhere — it is useless for annotation.

**Ranking:** Markers are ranked by `rank_score = avg_log2FC × (pct.1 - pct.2)`, which combines expression strength with specificity into a single sortable value.

```r
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

top_markers <- all_markers %>%
  filter(
    p_val_adj  < 0.05,
    avg_log2FC > 0.5,
    pct.1      > 0.25,
    pct.2      < 0.50
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
```

---

### Core Overlap Function

`overlap_score()` is the heart of the pipeline. It takes one cluster's top marker genes and one database, then for every cell type in that database computes four metrics:

| Metric | Formula | Meaning |
|---|---|---|
| `n_overlap` | \|A ∩ B\| | Raw count of shared genes |
| `pct_ref_covered` | \|A ∩ B\| / \|B\| × 100 | % of DB markers found in cluster |
| `pct_query_covered` | \|A ∩ B\| / \|A\| × 100 | % of cluster genes matching DB |
| `jaccard_index` | \|A ∩ B\| / \|A ∪ B\| | Primary ranking metric |

Results are returned sorted by `jaccard_index` descending — the best-matching cell type is always at row 1.

```r
overlap_score <- function(cluster_genes,
                           marker_db,
                           top_n = 20) {

  cluster_genes <- head(cluster_genes, top_n)
  n_query       <- length(cluster_genes)

  results <- lapply(names(marker_db), function(cell_type) {
    ref_genes <- unique(marker_db[[cell_type]])
    n_ref     <- length(ref_genes)
    overlap   <- intersect(cluster_genes, ref_genes)
    n_overlap <- length(overlap)

    if (n_overlap == 0) return(NULL)

    union_size <- length(union(cluster_genes, ref_genes))
    jaccard    <- round(n_overlap / union_size, 4)
    pct_ref    <- round(n_overlap / n_ref   * 100, 1)
    pct_query  <- round(n_overlap / n_query * 100, 1)

    data.frame(
      cell_type         = cell_type,
      n_overlap         = n_overlap,
      n_ref_markers     = n_ref,
      pct_ref_covered   = pct_ref,
      pct_query_covered = pct_query,
      jaccard_index     = jaccard,
      matched_genes     = paste(overlap, collapse = ", "),
      stringsAsFactors  = FALSE
    )
  }) %>%
    bind_rows() %>%
    arrange(desc(jaccard_index))

  return(results)
}
```

---

### Run Across All Clusters and Sources

`run_overlap_all_clusters()` is a wrapper that loops `overlap_score()` over every cluster in your dataset for a given database. The three calls produce `results_canonical`, `results_panglaodb`, and `results_cellmarker`, which are then combined with `bind_rows()` into a single long table: `all_results`.

```r
run_overlap_all_clusters <- function(top_markers, clusters,
                                      marker_db, source_label,
                                      top_n = 20) {
  lapply(clusters, function(cl) {
    genes  <- top_markers %>% filter(cluster == cl) %>% pull(gene)
    result <- overlap_score(genes, marker_db, top_n = top_n)
    if (is.null(result) || nrow(result) == 0) {
      return(data.frame(cluster=cl, source=source_label,
                        cell_type="Unknown", n_overlap=0,
                        jaccard_index=0, matched_genes=""))
    }
    result$cluster <- cl
    result$source  <- source_label
    result
  }) %>% bind_rows()
}

results_canonical  <- run_overlap_all_clusters(
  top_markers, clusters, canonical_markers, "Canonical")
results_panglaodb  <- run_overlap_all_clusters(
  top_markers, clusters, panglao_list, "PanglaoDB")
results_cellmarker <- run_overlap_all_clusters(
  top_markers, clusters, cellmarker_list, "CellMarker")

# Combine into one long table
all_results <- bind_rows(results_canonical,
                         results_panglaodb,
                         results_cellmarker)

write.csv(all_results, "overlap_all_sources_full.csv", row.names = FALSE)
```

**Structure of `all_results`:**

| Column | Description |
|---|---|
| `cluster` | Seurat cluster number |
| `source` | Database name (Canonical / PanglaoDB / CellMarker) |
| `cell_type` | Candidate cell type from database |
| `n_overlap` | Number of shared genes |
| `jaccard_index` | Similarity score (primary ranking) |
| `matched_genes` | Comma-separated list of matching genes |

---

### Top 3 Hits Per Cluster Per Source

`all_results` can contain hundreds of rows per cluster. Block 7 filters it down to the top 3 candidates per cluster per source, ranked by Jaccard. This is used for inspection — the 2nd and 3rd hits tell you about annotation ambiguity.

```r
top3_per_source <- all_results %>%
  filter(cell_type != "Unknown") %>%
  group_by(cluster, source) %>%
  slice_max(jaccard_index, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, source, desc(jaccard_index))
```

---

### Normalize Names and Majority Vote

**Step 8a — Pull top-1 per source:**

```r
top1_per_source <- all_results %>%
  filter(cell_type != "Unknown", jaccard_index > 0) %>%
  group_by(cluster, source) %>%
  slice_max(jaccard_index, n = 1, with_ties = FALSE) %>%
  ungroup()
```

**Step 8b — Normalize cell type names:**

Each database uses different naming conventions for the same cell type. `normalize_cell_type()` maps all variants to a shared vocabulary. **This function must be extended for any cell types specific to your tissue.**

```r
normalize_cell_type <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("_", " ", x)
  x <- gsub("\\s+", " ", x)

  mapping <- c(
    "t cells"="T cell", "t cell"="T cell",
    "cd4+ t cells"="CD4 T cell", "cd4 t"="CD4 T cell",
    "cd8+ t cells"="CD8 T cell", "cd8 t"="CD8 T cell",
    "naive t cells"="Naive T cell", "naive t cell"="Naive T cell",
    "regulatory t cells"="Treg", "regulatory t cell"="Treg",
    "nk cells"="NK cell", "natural killer cell"="NK cell",
    "b cells"="B cell", "b cell"="B cell",
    "plasma cells"="Plasma cell", "plasma cell"="Plasma cell",
    "classical monocytes"="Classical Monocyte",
    "classical monocyte"="Classical Monocyte",
    "nonclassical monocytes"="Non-Classical Monocyte",
    "non-classical monocyte"="Non-Classical Monocyte",
    "dendritic cells"="DC myeloid", "myeloid dc"="DC myeloid",
    "plasmacytoid dc"="Plasmacytoid DC",
    "megakaryocytes"="Platelet", "platelet"="Platelet",
    "macrophages"="Macrophage", "macrophage"="Macrophage"
    # Add more mappings as needed for your tissue
  )
  ifelse(x %in% names(mapping), mapping[x], tools::toTitleCase(x))
}

top1_normalized <- top1_per_source %>%
  mutate(cell_type_norm = normalize_cell_type(cell_type))
```

**Steps 8c–8d — Pivot wide and vote:**

```r
source_cols <- c("Canonical", "PanglaoDB", "CellMarker")

consensus_final <- consensus_wide %>%
  rowwise() %>%
  mutate(
    all_calls    = list(c_across(all_of(source_cols))),
    vote_table   = list(sort(table(unlist(all_calls)), decreasing=TRUE)),
    consensus    = names(vote_table)[1],
    n_agree      = as.integer(vote_table[1]),
    mean_jaccard = mean(c_across(ends_with("_jaccard")), na.rm=TRUE),
    confidence   = case_when(
      n_agree == 3                         ~ "High",
      n_agree == 2 & mean_jaccard >= 0.10 ~ "Medium",
      n_agree == 2 & mean_jaccard <  0.10 ~ "Low",
      n_agree == 1                         ~ "Conflict",
      TRUE                                 ~ "Unknown"
    )
  ) %>%
  select(-vote_table, -all_calls) %>%
  ungroup()
```

---

### Print Consensus Table

```r
consensus_final %>%
  select(cluster, Canonical, PanglaoDB, CellMarker,
         consensus, n_agree, mean_jaccard, confidence) %>%
  print(n = Inf)

# Flag any clusters needing manual review
conflicts <- consensus_final %>%
  filter(confidence %in% c("Conflict", "Low", "Unknown"))

if (nrow(conflicts) > 0) {
  cat("⚠ Clusters needing manual review:\n")
  print(conflicts %>% select(cluster, Canonical, PanglaoDB,
                              CellMarker, confidence))
}
```

---

### Visualization

Three plots are generated to visualise annotation quality:

**Plot 1 — Jaccard heatmap:** Rows = clusters, columns = sources, values = Jaccard score. Shows signal strength across all clusters at a glance.

**Plot 2 — Confidence dot plot:** X = cluster, Y = mean Jaccard, dot size = n_agree, colour = confidence level. Most useful single visualisation for a presentation or paper.

**Plot 3 — Source comparison bar chart:** Faceted by source — shows what each database called each cluster independently. Useful for spotting disagreements visually.

```r
library(pheatmap)

# Plot 1: Jaccard heatmap
jaccard_mat <- top1_normalized %>%
  select(cluster, source, jaccard_index) %>%
  pivot_wider(names_from=source, values_from=jaccard_index,
              values_fill=0) %>%
  column_to_rownames("cluster") %>%
  as.matrix()

pheatmap(jaccard_mat,
         main="Jaccard index per source",
         color=colorRampPalette(c("white","#2166ac"))(50),
         display_numbers=TRUE, number_format="%.2f",
         cluster_rows=FALSE, cluster_cols=FALSE)

# Plot 2: Confidence dot plot
ggplot(consensus_final,
  aes(x=as.factor(cluster), y=mean_jaccard,
      color=confidence, size=n_agree)) +
  geom_point(alpha=0.85) +
  scale_color_manual(values=c(
    "High"="#1a9641", "Medium"="#a6d96a",
    "Low"="#fdae61",  "Conflict"="#d7191c")) +
  labs(title="Annotation confidence per cluster",
       x="Cluster", y="Mean Jaccard") +
  theme_minimal()
```

---

### Apply Labels to Seurat Object

```r
label_map <- setNames(
  consensus_final$consensus,
  as.character(consensus_final$cluster)
)

Idents(seurat_obj) <- "seurat_clusters"
seurat_obj <- RenameIdents(seurat_obj, label_map)
seurat_obj[["cell_type_consensus"]] <- Idents(seurat_obj)

seurat_obj@meta.data$annotation_confidence <-
  consensus_final$confidence[
    match(seurat_obj@meta.data$seurat_clusters,
          consensus_final$cluster)
  ]

DimPlot(seurat_obj, group.by="cell_type_consensus",
        label=TRUE, repel=TRUE)
```

---

### Save Outputs

```r
write.csv(consensus_final,
          "consensus_annotation_final.csv", row.names=FALSE)
write.csv(top3_per_source,
          "top3_hits_per_source.csv",        row.names=FALSE)
saveRDS(seurat_obj, "seurat_annotated.rds")
```

---

## Output Files

| File | Contents | Use |
|---|---|---|
| `consensus_annotation_final.csv` | One row per cluster — all source calls, consensus label, confidence | Share in paper supplementary |
| `top3_hits_per_source.csv` | Top 3 candidates per cluster per source | Inspect ambiguous clusters |
| `overlap_all_sources_full.csv` | All overlap scores, all clusters, all cell types | Deep debugging and exploration |
| `seurat_annotated.rds` | Seurat object with `cell_type_consensus` in metadata | All downstream analyses |

---

## Interpreting Results

**High confidence — all 3 sources agree:**
Accept the label. Document the three matching sources in your methods.

**Medium confidence — 2/3 sources agree, mean Jaccard ≥ 0.10:**
Accept with a note. Run `FeaturePlot()` for the top matched genes to visually confirm.

**Low confidence — 2/3 agree but mean Jaccard < 0.10:**
Weak signal overall. Check `matched_genes` column — is the overlap biologically meaningful even if numerically small? Run visual checks before labelling.

**Conflict — only 1 source has a hit:**
Investigate manually. Common causes:
- Rare or tissue-specific cell type not well represented in all three databases
- Cluster is a transitional state between two cell types
- Cluster contains doublets — check with `scDblFinder`
- Batch effect creating a spurious cluster

**Always label genuine unknowns as `"Unknown"` — never force a label on a low-confidence cluster.**

---

## Confidence Levels

```
n_agree = 3  →  High      ████████  Accept directly
n_agree = 2, J ≥ 0.10  →  Medium   ██████    Accept, visual check recommended
n_agree = 2, J < 0.10  →  Low      ████      Visual check required
n_agree = 1  →  Conflict  ██        Manual review required
```

---

## Adapting for Non-PBMC Tissues

The pipeline works for any tissue — only two things need changing:

**1. Filter parameters in Blocks 3 and 4:**

```r
# Lung dataset
panglao_list    <- build_panglaodb_list(panglaodb_human,
                     organ_filter = c("Lung"))
cellmarker_list <- build_cellmarker_list(cellmarker_clean,
                     tissue_filter = "lung")

# Liver dataset
panglao_list    <- build_panglaodb_list(panglaodb_human,
                     organ_filter = c("Liver"))
cellmarker_list <- build_cellmarker_list(cellmarker_clean,
                     tissue_filter = "liver")

# All tissues (no filter)
panglao_list    <- build_panglaodb_list(panglaodb_human,
                     organ_filter = NULL)
cellmarker_list <- build_cellmarker_list(cellmarker_clean,
                     tissue_filter = NULL)
```

**2. Extend the `normalize_cell_type()` mapping table in Block 8** to include cell types specific to your tissue (hepatocytes, pneumocytes, neurons, etc.).

---

## References

1. Franzén O, Gan LM, Björkegren JLM. **PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data.** *Database* 2019. [doi:10.1093/database/baz046](https://doi.org/10.1093/database/baz046)

2. Zhang X et al. **CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools.** *Nucleic Acids Research* 2023. [doi:10.1093/nar/gkac947](https://doi.org/10.1093/nar/gkac947)

3. Hao Y et al. **Integrated analysis of multimodal single-cell data.** *Cell* 2021. [doi:10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048) — Seurat v4

4. Jaccard P. **Étude comparative de la distribution florale dans une portion des Alpes et des Jura.** *Bulletin de la Société Vaudoise des Sciences Naturelles* 1901. — Original Jaccard index

---

*Part of a YouTube series on single-cell RNA-seq annotation. See the full series for reference-based annotation (SingleR, Azimuth, CellTypist) and AI-assisted annotation (CASSIA, GPTCelltype).  https://www.youtube.com/playlist?list=PLAbpA2ThJEBSG-E1L0K9zwvlJB1vg4gbL*