# Single Cell RNA-seq — clustermole-Based Cell Type Annotation

A pipeline for annotating single-cell clusters using **clustermole**, an R package that queries 7 marker databases simultaneously and returns ranked cell type candidates with statistical significance scores.

> **Part of the single-cell annotation course.**
> This builds on the [marker-based annotation pipeline](../marker_based_annotation/marker_based_annotation.md) (Approach A) where we compared 3 databases using Jaccard overlap scoring. clustermole queries 7 databases in a single function call and uses a statistical test (hypergeometric) rather than Jaccard similarity.

---

## Table of Contents

- [Overview](#overview)
- [How clustermole Differs from Approach A](#how-clustermole-differs-from-approach-a)
- [Databases Queried](#databases-queried)
- [Two Annotation Methods](#two-annotation-methods)
- [The Verbose Label Problem and the Fix](#the-verbose-label-problem-and-the-fix)
- [Prerequisites](#prerequisites)
- [Input Requirements](#input-requirements)
- [Code Walkthrough](#code-walkthrough)
  - [Block 1 — Explore the Database](#block-1--explore-the-database)
  - [Block 2 — Prepare DGE Input](#block-2--prepare-dge-input)
  - [Method 1 — clustermole_overlaps](#method-1--clustermole_overlaps)
  - [Block 3 — Filter and Summarise Results](#block-3--filter-and-summarise-results)
  - [Block 4 — Normalize and Consensus Vote](#block-4--normalize-and-consensus-vote)
  - [Method 2 — clustermole_enrichment](#method-2--clustermole_enrichment)
  - [Block 5 — Compare Method 1 vs Method 2](#block-5--compare-method-1-vs-method-2)
  - [Block 6 — Visualization](#block-6--visualization)
  - [Block 7 — Apply Labels to Seurat Object](#block-7--apply-labels-to-seurat-object)
  - [Block 8 — Meta-Consensus with Approach A](#block-8--meta-consensus-with-approach-a)
  - [Block 9 — Save Outputs](#block-9--save-outputs)
- [Output Files](#output-files)
- [Interpreting Results](#interpreting-results)
- [Confidence Levels](#confidence-levels)
- [Adapting for Non-PBMC Tissues](#adapting-for-non-pbmc-tissues)
- [References](#references)

---

## Overview

clustermole is an R package that takes your cluster marker gene list and queries 7 curated databases in a single function call. For each cluster it returns a ranked list of candidate cell types from every database with a statistical significance score (FDR-adjusted p-value from a hypergeometric test).

The package supports two fundamentally different annotation modes:

1. **`clustermole_overlaps()`** — takes your top marker genes per cluster, runs a hypergeometric overrepresentation test against every cell type in all 7 databases. Returns a p-value and FDR per candidate.
2. **`clustermole_enrichment()`** — takes the full average expression matrix per cluster, runs ssGSEA enrichment scoring. Does not require you to define marker genes explicitly.

---

## How clustermole Differs from Approach A

| Aspect | Approach A (Jaccard) | clustermole |
|---|---|---|
| Databases | 3 (Canonical, PanglaoDB, CellMarker) | 7 simultaneously |
| Scoring metric | Jaccard index (0–1 similarity) | Hypergeometric p-value or ssGSEA score |
| Statistical test | None — similarity only | Yes — FDR-adjusted p-value |
| Needs DGE output | Yes | Method 1: Yes / Method 2: No |
| Needs Seurat object | No (works from CSV) | Method 2: Yes |
| External files needed | PanglaoDB TSV + CellMarker TXT | None — all DBs bundled |

**Key difference — statistical testing vs similarity scoring:**
Approach A computes a Jaccard index and asks *"how similar are these two gene lists?"* clustermole runs a hypergeometric test and asks *"is this overlap larger than expected by chance?"* The hypergeometric test accounts for the size of both gene lists when evaluating significance — a Jaccard of 0.15 may be meaningful or noise depending on list sizes, but a p-value of 0.001 always carries the same meaning.

---

## Databases Queried

All 7 databases are bundled inside the clustermole package — no separate downloads required.

| Database | Source | Cell Types |
|---|---|---|
| **ARCHS4** | Derived from RNA-seq co-expression patterns | 108 |
| **CellMarker** | Manually curated from literature | 692 |
| **MSigDB** | C8 scRNA-seq cell type signatures | 295 |
| **PanglaoDB** | Crowdsourced + scRNA-seq derived | 322 |
| **SaVanT** | Tissue/cell type gene expression signatures | 619 |
| **TISSUES** | Tissue-based expression signatures | 537 |
| **xCell** | 64 well-validated immune and stromal signatures | 466 |

> CellMarker and PanglaoDB are the same databases used in Approach A. Seeing them agree across both pipelines gives strong cross-method validation.

---

## Two Annotation Methods

```
METHOD 1 — clustermole_overlaps()
──────────────────────────────────
Input  : Top N marker genes per cluster (from FindAllMarkers)
Test   : Hypergeometric overrepresentation test
Output : FDR p-value per cell type per database
Use    : When you have clean DGE markers

  Cluster genes ──► overlap with every DB gene set
                              │
                    hypergeometric test
                              │
                   p-value + FDR per cell type


METHOD 2 — clustermole_enrichment()
─────────────────────────────────────
Input  : Average expression matrix (genes × clusters)
Test   : ssGSEA / GSVA / singscore
Output : Enrichment score per cell type per cluster
Use    : Independent of marker definition, uses all genes

  Average expression ──► ssGSEA against all DB gene sets
  per cluster                        │
                          enrichment score per cell type
```

---

## The Verbose Label Problem and the Fix

**This is the most important concept to understand before running clustermole.**

clustermole returns dataset-specific, verbose labels. For example, all of the following mean T cell:

```
"Fan Embryonic Ctx Brain Naive Like t Cell"
"Imgn t 8eff Sp Ot1 D10 Lisova"
"Hay Bone Marrow Naive t Cell"
"Cd3plus t Cell"
"Tlymphocyte"
"t Memory Cells"
```

If you vote directly on raw labels, each database returns a different string — all counts are 1 — and the "winner" is picked alphabetically. This is meaningless.

**The fix is `normalize_cell_type()`** — a two-layer function that maps all verbose labels to broad PBMC categories before voting:

```
BEFORE — every count = 1, vote picks alphabetically ✗

  ARCHS4     → "Fan Embryonic Ctx Brain Naive Like t Cell"
  CellMarker → "Naive Cd4plus t Cell"
  MSigDB     → "Cd3plus t Cell"
  PanglaoDB  → "t Memory Cells"
  SaVanT     → "Tlymphocyte"
  TISSUES    → "Blood"
  xCell      → "Hay Bone Marrow Naive t Cell"


AFTER normalize_cell_type() — vote works correctly ✓

  ARCHS4     → "T cell"           ✓
  CellMarker → "CD4 T cell"       ✓
  MSigDB     → "T cell"           ✓
  PanglaoDB  → "T cell"           ✓
  SaVanT     → "T cell"           ✓
  TISSUES    → "Unresolved"       ← excluded from vote
  xCell      → "T cell"           ✓

  Winner: T cell   n_agree = 5/6   confidence = Very High
```

---

## Prerequisites

```r
install.packages("clustermole")
install.packages(c("Seurat", "dplyr", "tidyr", "ggplot2", "pheatmap"))
```

No external database files needed.

---

## Input Requirements

The pipeline starts from `FindAllMarkers()` output. If not yet run:

```r
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  test.use        = "wilcox"
)
```

For Method 2 you additionally need the Seurat object to compute `AverageExpression()`.

---

## Code Walkthrough

### Block 1 — Explore the Database

Before running annotation, inspect what is available inside clustermole.

```r
library(clustermole)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# Load all human markers
all_markers_db <- clustermole_markers(species = "hs")

cat("Total entries (human):", nrow(all_markers_db), "\n")
glimpse(all_markers_db)
dim(all_markers_db)
```

**Columns in the clustermole marker database:**

| Column | Description |
|---|---|
| `celltype_full` | Full name including source dataset |
| `db` | Source database (ARCHS4, CellMarker, etc.) |
| `species` | `hs` (human) or `mm` (mouse) |
| `organ` | Tissue / organ |
| `celltype` | Short cell type name |
| `n_genes` | Number of marker genes for this cell type |
| `gene` | Individual gene symbol — one per row |

```r
# Count cell types per database
all_markers_db %>%
  distinct(db, celltype_full) %>%
  count(db, name = "n_cell_types") %>%
  arrange(desc(n_cell_types))

# Browse cell types in PanglaoDB
all_markers_db %>%
  filter(db == "PanglaoDB") %>%
  distinct(celltype, organ) %>%
  arrange(organ, celltype) %>%
  print(n = 30)

# Filter to immune / blood-related entries
immune_types <- all_markers_db %>%
  filter(
    grepl("blood|immune|lymph|T cell|B cell|NK|monocyte|macrophage",
          organ, ignore.case = TRUE) |
    grepl("blood|immune|lymph|T cell|B cell|NK|monocyte|macrophage",
          celltype, ignore.case = TRUE)
  ) %>%
  distinct(db, celltype, organ)

cat("Immune-related cell types:", nrow(immune_types), "\n")
dim(immune_types)
```

---

### Block 2 — Prepare DGE Input

Same filtering and ranking as Approach A. The key filter is `pct.2 < 0.50` — the specificity filter that removes genes expressed broadly across all clusters.

```r
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
```

---

### Method 1 — clustermole_overlaps()

`clustermole_overlaps()` takes your top marker genes for one cluster and runs a hypergeometric test against every cell type in all 7 databases simultaneously.

**How the hypergeometric test works:**

```
Given:
  n_query   = size of your gene list             (e.g. 20)
  n_ref     = size of reference gene set         (e.g. 50)
  overlap   = genes in common                    (e.g. 5)
  N         = background gene universe           (~20,000)

The test asks: what is the probability of seeing this many
overlapping genes purely by chance, given these sizes?

Low p-value → overlap is not random → strong candidate match
```

```r
overlaps_all <- lapply(clusters, function(cl) {

  # Top 20 marker genes ranked by rank_score
  genes <- top_markers %>%
    filter(cluster == cl) %>%
    slice_max(rank_score, n = 20) %>%
    pull(gene)

  cat("Cluster", cl, "— querying", length(genes),
      "genes across all databases...\n")

  # Single call queries ALL 7 databases
  result <- tryCatch(
    clustermole_overlaps(genes   = genes,
                         species = "hs"),
    error = function(e) {
      warning("Cluster ", cl, " failed: ", e$message)
      return(NULL)
    }
  )

  if (is.null(result) || nrow(result) == 0) return(NULL)
  result$cluster <- cl
  result

}) %>%
  bind_rows()

glimpse(overlaps_all)
cat("Total overlap results:", nrow(overlaps_all), "\n")
```

**Columns returned by `clustermole_overlaps()`:**

| Column | Description |
|---|---|
| `celltype_full` | Full label including source dataset name |
| `db` | Database (ARCHS4, CellMarker, MSigDB, etc.) |
| `organ` | Tissue context |
| `celltype` | Short cell type name |
| `n_genes` | Size of reference gene set |
| `overlap` | Number of overlapping genes |
| `p_value` | Raw hypergeometric p-value |
| `fdr` | BH-adjusted FDR — **primary ranking metric** |
| `cluster` | Cluster number (added in the loop above) |

> **Note:** The overlap count column is named `overlap` (not `n_overlap`). All filtering and selection steps in this script use this exact column name.

---

### Block 3 — Filter and Summarise Results

```r
# Significant hits only (FDR < 0.05)
overlaps_sig <- overlaps_all %>%
  filter(fdr < 0.05) %>%
  arrange(cluster, fdr)

cat("Significant results (FDR < 0.05):", nrow(overlaps_sig), "\n")

# Top 5 per cluster across ALL databases combined
top5_all_db <- overlaps_all %>%
  filter(!is.na(fdr)) %>%
  group_by(cluster) %>%
  slice_min(fdr, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, db, celltype, organ, overlap, p_value, fdr)

cat("\n========== TOP 5 HITS PER CLUSTER (ALL DBS) ==========\n")
print(top5_all_db, n = Inf)

# Top 3 per cluster per database
top3_per_db <- overlaps_all %>%
  filter(!is.na(fdr)) %>%
  group_by(cluster, db) %>%
  slice_min(fdr, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, db, fdr)

cat("\n========== TOP 3 HITS PER CLUSTER PER DATABASE ==========\n")
print(top3_per_db %>%
  select(cluster, db, celltype, overlap, fdr),
  n = Inf)

# Top-1 per database — used for consensus voting
top1_per_db <- overlaps_all %>%
  filter(!is.na(fdr), fdr < 0.05) %>%
  group_by(cluster, db) %>%
  slice_min(fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, db, celltype, fdr) %>%
  arrange(cluster, fdr)

# Wide view — one column per database
top1_wide <- top1_per_db %>%
  select(cluster, db, celltype) %>%
  pivot_wider(names_from  = db,
              values_from = celltype,
              values_fill = "—") %>%
  arrange(as.numeric(as.character(cluster)))

cat("\n========== TOP-1 PER CLUSTER PER DATABASE (WIDE) ==========\n")
print(top1_wide, n = Inf)
```

---

### Block 4 — Normalize and Consensus Vote

#### The `normalize_cell_type()` function

This is the critical function that solves the verbose label problem. It uses two layers:

- **Layer 1 — Exact match lookup:** 47 specific labels observed in real PBMC clustermole output, mapped directly to broad categories
- **Layer 2 — Keyword / regex fallback:** For any label not in the exact table — handles new database entries or non-PBMC tissue labels

```r
normalize_cell_type <- function(x) {

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
    # Platelet
    "Platelet"                                  = "Platelet",
    "Hpca Platelets"                            = "Platelet",
    "Platelets Hpca 1"                          = "Platelet",
    "Zheng Cord Blood C1 Putative Megakaryocyte Progenitor" = "Platelet",
    # Too broad to vote with
    "Blood"                                     = "Unresolved",
    "Blood Plasma"                              = "Unresolved",
    "Peripheral Blood"                          = "Unresolved",
    "Immune System"                             = "Unresolved"
  )

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

  sapply(x, function(label) {
    if (label %in% names(exact_map)) exact_map[[label]]
    else keyword_fallback(label)
  }, USE.NAMES = FALSE)
}
```

> **Note:** The script also retains `normalize_cell_type_old_dontuse()` — the previous version — for reference only. Do not use it. It cannot handle verbose clustermole labels and defaults to alphabetical sort when counts tie.

#### Applying normalization and voting

```r
# Apply normalization to top-1 labels
top1_per_db_norm <- top1_per_db %>%
  mutate(celltype_norm = normalize_cell_type(celltype))

# Majority vote across databases
consensus_clustermole <- top1_per_db_norm %>%
  group_by(cluster) %>%
  summarise(
    vote_result  = {
      tbl <- sort(table(celltype_norm), decreasing = TRUE)
      names(tbl)[1]
    },
    n_agree      = {
      tbl <- sort(table(celltype_norm), decreasing = TRUE)
      as.integer(tbl[1])
    },
    n_db_total   = n_distinct(db),
    dbs_agreeing = {
      tbl  <- sort(table(celltype_norm), decreasing = TRUE)
      top  <- names(tbl)[1]
      paste(db[celltype_norm == top], collapse = ", ")
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
  select(cluster, vote_result, n_agree, n_db_total,
         pct_agree, mean_fdr, confidence, dbs_agreeing),
  n = Inf)
```

**Key columns in `consensus_clustermole`:**

| Column | Description |
|---|---|
| `vote_result` | Winning broad category after normalization |
| `n_agree` | Databases voting for the winner |
| `n_db_total` | Total databases with a significant hit |
| `pct_agree` | `n_agree / n_db_total × 100` |
| `dbs_agreeing` | Names of agreeing databases |
| `mean_fdr` | Mean FDR across agreeing databases |
| `confidence` | Confidence level — see table below |

---

### Method 2 — clustermole_enrichment()

Method 2 uses the full average expression profile per cluster instead of DGE markers. ssGSEA enrichment is computed against every cell type gene set in all databases.

**When to prefer Method 2:**
- DGE markers are noisy or very few per cluster
- You want annotation independent of clustering resolution
- You want an independent cross-validation of Method 1 results

```r
# Average expression per cluster
avg_exp     <- AverageExpression(seurat_obj,
                 assays = "RNA", group.by = "seurat_clusters",
                 slot = "data")   # log-normalised
avg_exp_mat <- avg_exp$RNA

cat("Matrix:", nrow(avg_exp_mat), "genes ×",
    ncol(avg_exp_mat), "clusters\n")

# Note: if slot="counts" was used, log-transform before passing:
# avg_exp_mat <- log1p(avg_exp_mat)

cat("Running enrichment — this may take a few minutes...\n")

enrich_results <- clustermole_enrichment(
  expr_mat = avg_exp_mat,
  species  = "hs",
  method   = "ssgsea"     # options: ssgsea, gsva, singscore, all
)

glimpse(enrich_results)

# Top 3 enrichment hits per cluster per database
top_enrich <- enrich_results %>%
  group_by(cluster, db) %>%
  slice_max(score, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, db, desc(score))

cat("\n========== TOP 3 ENRICHMENT HITS PER CLUSTER PER DB ==========\n")
print(top_enrich %>% select(cluster, db, celltype, score), n = Inf)

# Best enrichment hit per cluster
best_enrich <- enrich_results %>%
  group_by(cluster) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(celltype_norm = normalize_cell_type(celltype)) %>%
  select(cluster, celltype_norm, db, score) %>%
  arrange(as.numeric(as.character(cluster)))

cat("\n========== BEST ENRICHMENT HIT PER CLUSTER ==========\n")
print(best_enrich, n = Inf)
```

> `normalize_cell_type()` is applied to Method 2 results here too — so both methods produce labels in the same broad category vocabulary, making the Block 5 comparison straightforward.

---

### Block 5 — Compare Method 1 vs Method 2

```r
comparison <- consensus_clustermole %>%
  select(cluster,
         overlap_label      = vote_result,
         overlap_confidence = confidence,
         overlap_n_agree    = n_agree) %>%
  left_join(
    best_enrich %>%
      select(cluster,
             enrichment_label = celltype_norm,
             enrichment_db    = db,
             enrichment_score = score),
    by = "cluster"
  ) %>%
  mutate(
    methods_agree = (overlap_label == enrichment_label),
    final_label   = overlap_label,   # overlap-based label is primary
    note = case_when(
      methods_agree  ~ "Both methods agree",
      !methods_agree ~ paste("Review — Overlap:", overlap_label,
                             "| Enrichment:", enrichment_label),
      TRUE           ~ ""
    )
  )

cat("\n========== METHOD 1 vs METHOD 2 COMPARISON ==========\n")
print(comparison %>%
  select(cluster, overlap_label, enrichment_label,
         methods_agree, overlap_confidence, note),
  n = Inf)
```

---

### Block 6 — Visualization

**Plot 1 — FDR bubble plot** (significance per database per cluster)

X = cluster, Y = database, dot size and colour = `-log10(FDR)`. Shows which database has the most significant hit for each cluster.

```r
plot_data <- top1_per_db %>%
  mutate(neg_log_fdr = -log10(fdr + 1e-10),
         cluster     = as.factor(cluster))

ggplot(plot_data,
       aes(x = cluster, y = db,
           size = neg_log_fdr, color = neg_log_fdr)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2, 10), name = "-log10(FDR)") +
  scale_color_gradient(low = "#d9f0d3", high = "#1a6634",
                       name = "-log10(FDR)") +
  labs(title = "Annotation significance per cluster per DB",
       x = "Cluster", y = "Database") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

**Plot 2 — Label matrix** (what each database called each cluster)

Printed as a table — shows normalized label from each database side by side.

```r
label_mat <- top1_per_db %>%
  mutate(celltype_short = normalize_cell_type(celltype)) %>%
  select(cluster, db, celltype_short) %>%
  pivot_wider(names_from  = db,
              values_from = celltype_short,
              values_fill = "—") %>%
  column_to_rownames("cluster")

cat("\n========== LABEL MATRIX (cluster × database) ==========\n")
print(label_mat)
```

**Plot 3 — Confidence bar chart** (agreement per cluster)

Bar height = `pct_agree`, fill = confidence level, label = `n_agree / n_db_total`.

```r
conf_levels <- c("Very High","High","Medium","Low","Conflict","Unknown")
conf_colors <- c("#1a9641","#a6d96a","#ffffbf","#fdae61","#d7191c","#999999")

consensus_clustermole %>%
  mutate(confidence = factor(confidence, levels = conf_levels)) %>%
  ggplot(aes(x = as.factor(cluster), y = pct_agree,
             fill = confidence)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(n_agree, "/", n_db_total)),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = setNames(conf_colors, conf_levels),
                    drop = FALSE) +
  labs(title = "Database agreement per cluster",
       x = "Cluster", y = "% databases agreeing") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

**Plot 4 — Enrichment score heatmap** (Method 2 only)

Top 20 cell types by maximum enrichment score, row-scaled for contrast.

```r
top_types_enrich <- enrich_results %>%
  group_by(celltype) %>%
  summarise(max_score = max(score), .groups = "drop") %>%
  slice_max(max_score, n = 20) %>%
  pull(celltype)

enrich_mat <- enrich_results %>%
  filter(celltype %in% top_types_enrich) %>%
  select(celltype, cluster, score) %>%
  pivot_wider(names_from = cluster, values_from = score,
              values_fill = 0) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

pheatmap(enrich_mat,
         main         = "clustermole enrichment scores (ssGSEA)",
         color        = colorRampPalette(c("#f7f7f7","#2166ac"))(50),
         cluster_rows = TRUE, cluster_cols = FALSE,
         fontsize_row = 8, fontsize_col = 9,
         scale        = "row")
```

---

### Block 7 — Apply Labels to Seurat Object

Labels from the overlap-based consensus are applied to the Seurat object. Confidence and n_agree are also stored in metadata for downstream quality filtering.

```r
label_map <- setNames(
  consensus_clustermole$vote_result,
  as.character(consensus_clustermole$cluster)
)

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

DimPlot(seurat_obj,
        group.by = "cell_type_clustermole",
        label = TRUE, repel = TRUE,
        label.size = 3.5, pt.size = 0.4) +
  ggtitle("clustermole annotation — 7 databases consensus")
```

---

### Block 8 — Meta-Consensus with Approach A

Merges results from the Approach A Jaccard pipeline and clustermole. Where both pipelines agree you have very high cross-method confidence. Where they disagree the cluster needs manual investigation.

```r
approachA <- read.csv("consensus_annotation_final.csv")

meta_consensus <- approachA %>%
  select(cluster,
         approachA_label      = consensus,
         approachA_confidence = confidence) %>%
  left_join(
    consensus_clustermole %>%
      select(cluster,
             clustermole_label      = vote_result,
             clustermole_confidence = confidence,
             clustermole_n_agree    = n_agree,
             clustermole_n_db       = n_db_total),
    by = "cluster"
  ) %>%
  mutate(
    approachA_norm     = normalize_cell_type(approachA_label),
    clustermole_norm   = normalize_cell_type(clustermole_label),
    pipelines_agree    = (approachA_norm == clustermole_norm),
    final_label        = ifelse(pipelines_agree,
                                approachA_norm,
                                paste0(approachA_norm, "* (review)")),
    overall_confidence = case_when(
      pipelines_agree & approachA_confidence == "High"   ~ "Very High",
      pipelines_agree & approachA_confidence == "Medium" ~ "High",
      pipelines_agree                                    ~ "Medium",
      !pipelines_agree                                   ~ "Review"
    )
  )

cat("\n========== META-CONSENSUS: APPROACH A + CLUSTERMOLE ==========\n")
meta_consensus %>%
  select(cluster, approachA_norm, clustermole_norm,
         pipelines_agree, overall_confidence, final_label) %>%
  print(n = Inf)

# Flag disagreements
disagreements <- meta_consensus %>% filter(!pipelines_agree)

if (nrow(disagreements) > 0) {
  cat("\n⚠ Clusters where pipelines disagree — manual review needed:\n")
  print(disagreements %>%
    select(cluster, approachA_norm, clustermole_norm,
           approachA_confidence, clustermole_confidence))
} else {
  cat("\n✓ Both pipelines agree on all clusters.\n")
}
```

---

### Block 9 — Save Outputs

```r
write.csv(overlaps_all,
          "clustermole_overlaps_full.csv",            row.names = FALSE)
write.csv(top3_per_db,
          "clustermole_top3_per_db.csv",              row.names = FALSE)
write.csv(consensus_clustermole,
          "clustermole_consensus.csv",                row.names = FALSE)
write.csv(meta_consensus,
          "meta_consensus_approachA_clustermole.csv", row.names = FALSE)
saveRDS(seurat_obj, "seurat_clustermole_annotated.rds")
```

---

## Output Files

| File | Contents | Use |
|---|---|---|
| `clustermole_overlaps_full.csv` | All overlap results, all clusters, all databases | Deep debugging |
| `clustermole_top3_per_db.csv` | Top 3 hits per cluster per database | Inspect candidates |
| `clustermole_consensus.csv` | Final consensus with confidence per cluster | Main result |
| `meta_consensus_approachA_clustermole.csv` | Merged Approach A + clustermole results | Cross-pipeline validation |
| `seurat_clustermole_annotated.rds` | Seurat object with `cell_type_clustermole` in metadata | Downstream analyses |

---

## Interpreting Results

**Very High / High — 4–7 databases agree:**
Accept directly. Consistent signal across diverse database types.

**Medium — 3 databases agree:**
Accept with a note. Check `dbs_agreeing` — agreement across diverse database types (e.g. CellMarker + ARCHS4 + xCell) is stronger evidence than agreement within similar types (e.g. CellMarker + PanglaoDB only, which share similar curation origins).

**Low — 2 databases agree:**
Run `FeaturePlot()` for the top matched genes before labelling. Weak signal overall — may be a rare cell type, a transitional state, or a cluster with very few marker genes.

**Conflict — only 1 database has a significant hit:**
Investigate manually. Common causes:
- Cell type is rare or poorly represented across databases
- Cluster contains doublets — check with `scDblFinder`
- Cluster reflects a cell state (e.g. cycling cells) rather than a stable cell type identity
- Too few marker genes — consider adjusting clustering resolution

---

## Confidence Levels

```
n_agree ≥ 5              →  Very High  ██████  Accept directly
n_agree = 4              →  High       █████   Accept directly
n_agree = 3, FDR < 0.01  →  High       █████   Accept
n_agree = 3              →  Medium     ████    Visual check recommended
n_agree = 2              →  Low        ███     Visual check required
n_agree = 1              →  Conflict   █       Manual review required
```

---

## Adapting for Non-PBMC Tissues

Only `normalize_cell_type()` needs updating.

**Step 1:** Run the pipeline once. Look at the raw `celltype` column in `clustermole_top3_per_db.csv`. Find any labels still showing as `"Unresolved"` that should map to something specific.

**Step 2:** Add those labels to `exact_map` inside `normalize_cell_type()`:

```r
# Example — liver-specific labels
"Aizarani Liver C3 Hepatocytes 1"            = "Hepatocyte",
"Aizarani Liver C5 Hepatocytes 2"            = "Hepatocyte",
"Aizarani Liver C9 Kupffer Cells"            = "Macrophage",
"Aizarani Liver C25 Portal Endothelial Cells" = "Endothelial",
"Aizarani Liver C31 Cholangiocytes 1"        = "Cholangiocyte",
```

**Step 3:** For truly novel labels the regex cannot catch, add a line to `keyword_fallback()`:

```r
if (grepl("hepatocyte|\\balb\\b|\\bhnf4",  l)) return("Hepatocyte")
if (grepl("podocyte|\\bnphs1",             l)) return("Podocyte")
if (grepl("club cell|\\bscgb1a1",          l)) return("Club cell")
```

---

## References

1. Dolgalev I. **clustermole: Cell Type Identification from Single-Cell Transcriptomics Data.** CRAN. [igordot.github.io/clustermole](https://igordot.github.io/clustermole)

2. Hao Y et al. **Integrated analysis of multimodal single-cell data.** *Cell* 2021. [doi:10.1016/j.cell.2021.04.048](https://doi.org/10.1016/j.cell.2021.04.048) — Seurat v4

3. Zhang X et al. **CellMarker 2.0: an updated database of manually curated cell markers.** *Nucleic Acids Research* 2023. [doi:10.1093/nar/gkac947](https://doi.org/10.1093/nar/gkac947)

4. Franzén O et al. **PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data.** *Database* 2019. [doi:10.1093/database/baz046](https://doi.org/10.1093/database/baz046)

5. Liberzon A et al. **The Molecular Signatures Database Hallmark Gene Set Collection.** *Cell Systems* 2015. — MSigDB C8 collection

6. Aran D et al. **xCell: digitally portraying the tissue cellular heterogeneity landscape.** *Genome Biology* 2017. [doi:10.1186/s13059-017-1349-1](https://doi.org/10.1186/s13059-017-1349-1)

---

*Part of a YouTube series on single-cell RNA-seq annotation. This is Module 1 — marker-based annotation. See the repository for reference-based annotation (SingleR, Azimuth, CellTypist) and AI-assisted annotation (CASSIA, GPTCelltype) modules.*