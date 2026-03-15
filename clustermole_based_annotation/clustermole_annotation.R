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

library(clustermole)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)


## ============================================================
## BLOCK 1 — Explore the clustermole database
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

dim(immune_types)

## ============================================================
## BLOCK 2 — Prepare your DGE input
## ============================================================
##
##  clustermole_overlaps() needs a simple character vector
##  of gene symbols for each cluster — same as our Approach A.
##
##  Start from FindAllMarkers() output filtered and ranked.

# If you haven't run FindAllMarkers yet:
# all_markers <- FindAllMarkers(
#   seurat_obj, only.pos = TRUE,
#   min.pct = 0.25, logfc.threshold = 0.25
# )

# Filter and rank — same logic as Approach A
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
glimpse(overlaps_all)

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
  filter(fdr < 0.05) %>%
  arrange(cluster, fdr)

cat("Significant results (FDR < 0.05):", nrow(overlaps_sig), "\n")

# Top 5 hits per cluster across ALL databases
top5_all_db <- overlaps_all %>%
  filter(!is.na(fdr)) %>%
  group_by(cluster) %>%
  slice_min(fdr, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, db, celltype, organ, overlap, p_value, fdr)

cat("\n========== TOP 5 HITS PER CLUSTER (ALL DBS) ==========\n")
print(top5_all_db, n = Inf)

# Top 3 hits per cluster PER database
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

# Best hit (top-1) per cluster per database — for consensus
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


## ============================================================
## BLOCK 4 — Consensus from clustermole overlaps
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
  select(cluster, vote_result, n_agree, n_db_total,
         pct_agree, mean_fdr, confidence, dbs_agreeing),
  n = Inf)


## ============================================================
## METHOD 2 — clustermole_enrichment()
## ============================================================
##
##  This method is fundamentally different from Method 1.
##  Instead of using your DGE marker genes, it uses the
##  FULL average expression profile per cluster.
##
##  What it does:
##    - Takes average expression per cluster (genes × clusters)
##    - Runs ssGSEA / GSVA / singscore for every cell type
##      signature in the database
##    - Returns an enrichment SCORE not a p-value
##    - Higher score = stronger cell type signature activity
##
##  Advantage over Method 1:
##    - Does not depend on how you defined markers
##    - Uses all genes, not just top N
##    - Less sensitive to clustering resolution
##
##  Disadvantage:
##    - Slower (running enrichment across thousands of gene sets)
##    - Harder to interpret which genes drove the result

# ── Step 1: Compute average expression per cluster ───────────

# AverageExpression returns a list; we want the RNA assay matrix
avg_exp <- AverageExpression(
  seurat_obj,
  assays   = "RNA",
  group.by = "seurat_clusters",
  slot     = "data"          # log-normalised values
)

# Extract the matrix (genes × clusters)
avg_exp_mat <- avg_exp$RNA

cat("Expression matrix dimensions:",
    nrow(avg_exp_mat), "genes ×",
    ncol(avg_exp_mat), "clusters\n")

# ── Step 2: Log-transform (required by clustermole) ──────────
# AverageExpression with slot="data" already returns log1p values
# If you used slot="counts", apply log1p here:
# avg_exp_mat <- log1p(avg_exp_mat)

# ── Step 3: Run clustermole_enrichment ───────────────────────
##
##  method options:
##    "ssgsea"   — single-sample GSEA (default, recommended)
##    "gsva"     — GSVA (Hanzelmann et al 2013)
##    "singscore" — rank-based single-sample scoring
##    "all"      — run all three and return combined results

cat("Running enrichment — this may take a few minutes...\n")

enrich_results <- clustermole_enrichment(
  expr_mat = avg_exp_mat,
  species  = "hs",
  method   = "ssgsea"        # ssgsea is fastest and robust
)

# Structure of enrichment results:
# celltype_full : full cell type label
# db            : source database
# organ         : tissue context
# celltype      : short name
# cluster       : your cluster (column name from matrix)
# score         : enrichment score (higher = more enriched)

glimpse(enrich_results)

# ── Step 4: Top hits from enrichment per cluster ─────────────

top_enrich <- enrich_results %>%
  group_by(cluster, db) %>%
  slice_max(score, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, db, desc(score))

cat("\n========== TOP 3 ENRICHMENT HITS PER CLUSTER PER DB ==========\n")
print(top_enrich %>%
  select(cluster, db, celltype, score),
  n = Inf)

# ── Step 5: Best enrichment hit per cluster ──────────────────
best_enrich <- enrich_results %>%
  group_by(cluster) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(celltype_norm = normalize_cell_type(celltype)) %>%
  select(cluster, celltype_norm, db, score) %>%
  arrange(as.numeric(as.character(cluster)))

cat("\n========== BEST ENRICHMENT HIT PER CLUSTER ==========\n")
print(best_enrich, n = Inf)


## ============================================================
## BLOCK 5 — Compare Method 1 vs Method 2
## ============================================================

# Merge overlap-based and enrichment-based results
comparison <- consensus_clustermole %>%
  select(cluster,
         overlap_label = vote_result,
         overlap_confidence = confidence,
         overlap_n_agree = n_agree) %>%
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
    final_label   = ifelse(methods_agree,
                           overlap_label,
                           overlap_label),  # overlap wins if disagree
    note = case_when(
      methods_agree                ~ "Both methods agree",
      !methods_agree               ~ paste("Review —",
                                           "Overlap:", overlap_label,
                                           "| Enrichment:", enrichment_label),
      TRUE                         ~ ""
    )
  )

cat("\n========== METHOD 1 vs METHOD 2 COMPARISON ==========\n")
print(comparison %>%
  select(cluster, overlap_label, enrichment_label,
         methods_agree, overlap_confidence, note),
  n = Inf)


## ============================================================
## BLOCK 6 — Visualization
## ============================================================

# ── Plot 1: FDR bubble plot across databases ─────────────────
# Shows significance of top hit per cluster per database
# Smaller FDR = more significant = darker/larger bubble

plot_data <- top1_per_db %>%
  mutate(
    neg_log_fdr = -log10(fdr + 1e-10),
    cluster     = as.factor(cluster)
  )

ggplot(plot_data,
       aes(x = cluster, y = db,
           size = neg_log_fdr,
           color = neg_log_fdr)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(2, 10),
                        name  = "-log10(FDR)") +
  scale_color_gradient(low  = "#d9f0d3",
                       high = "#1a6634",
                       name = "-log10(FDR)") +
  labs(
    title = "clustermole — annotation significance per cluster per DB",
    x     = "Cluster",
    y     = "Database"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ── Plot 2: Heatmap of top cell type labels ──────────────────
# Shows which cell type each database predicts for each cluster

label_mat <- top1_per_db %>%
  mutate(celltype_short = normalize_cell_type(celltype)) %>%
  select(cluster, db, celltype_short) %>%
  pivot_wider(names_from  = db,
              values_from = celltype_short,
              values_fill = "—") %>%
  column_to_rownames("cluster")

# Convert labels to numeric for heatmap coloring
# (just show as annotation table instead)
cat("\n========== LABEL MATRIX (cluster × database) ==========\n")
print(label_mat)


# ── Plot 3: Confidence summary bar ───────────────────────────
conf_levels <- c("Very High","High","Medium","Low","Conflict","Unknown")
conf_colors <- c("#1a9641","#a6d96a","#ffffbf","#fdae61","#d7191c","#999999")

consensus_clustermole %>%
  mutate(confidence = factor(confidence, levels = conf_levels)) %>%
  ggplot(aes(x = as.factor(cluster),
             y = pct_agree,
             fill = confidence)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(n_agree, "/", n_db_total)),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = setNames(conf_colors, conf_levels),
                    drop = FALSE) +
  labs(
    title = "clustermole — database agreement per cluster",
    x     = "Cluster",
    y     = "% databases agreeing",
    fill  = "Confidence"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ── Plot 4: Enrichment score heatmap (Method 2) ──────────────
# Top 20 cell types by enrichment score across all clusters

top_types_enrich <- enrich_results %>%
  group_by(celltype) %>%
  summarise(max_score = max(score), .groups = "drop") %>%
  slice_max(max_score, n = 20) %>%
  pull(celltype)

enrich_mat <- enrich_results %>%
  filter(celltype %in% top_types_enrich) %>%
  select(celltype, cluster, score) %>%
  pivot_wider(names_from  = cluster,
              values_from = score,
              values_fill = 0) %>%
  column_to_rownames("celltype") %>%
  as.matrix()

pheatmap(
  enrich_mat,
  main         = "clustermole enrichment scores (ssGSEA)",
  color        = colorRampPalette(c("#f7f7f7","#2166ac"))(50),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 9,
  scale        = "row"   # scale per cell type for better contrast
)


## ============================================================
## BLOCK 7 — Apply labels to Seurat object
## ============================================================

# Use overlap-based consensus as primary label
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


## ============================================================
## BLOCK 8 — Merge with Approach A for final meta-consensus
## ============================================================
##
##  Now you have two independent pipelines:
##    Approach A : Jaccard overlap — Canonical + PanglaoDB + CellMarker
##    clustermole: Hypergeometric  — 7 databases simultaneously
##
##  Where they agree = very high confidence
##  Where they disagree = investigate those clusters

# Read Approach A results (from our previous pipeline)
approachA <- read.csv("consensus_annotation_final.csv")

# Merge
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
    # Normalize both labels to same vocabulary
    approachA_norm     = normalize_cell_type(approachA_label),
    clustermole_norm   = normalize_cell_type(clustermole_label),
    pipelines_agree    = (approachA_norm == clustermole_norm),

    # Final label — if both agree use it, else flag
    final_label        = ifelse(pipelines_agree,
                                approachA_norm,
                                paste0(approachA_norm,
                                       "* (review)")),

    # Overall confidence
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
disagreements <- meta_consensus %>%
  filter(!pipelines_agree)

if (nrow(disagreements) > 0) {
  cat("\n⚠ Clusters where pipelines disagree — manual review needed:\n")
  print(disagreements %>%
    select(cluster, approachA_norm, clustermole_norm,
           approachA_confidence, clustermole_confidence))
} else {
  cat("\n✓ Both pipelines agree on all clusters.\n")
}


## ============================================================
## BLOCK 9 — Save outputs
## ============================================================

write.csv(overlaps_all,
          "clustermole_overlaps_full.csv",       row.names = FALSE)
write.csv(top3_per_db,
          "clustermole_top3_per_db.csv",         row.names = FALSE)
write.csv(consensus_clustermole,
          "clustermole_consensus.csv",            row.names = FALSE)
write.csv(meta_consensus,
          "meta_consensus_approachA_clustermole.csv",
                                                  row.names = FALSE)
saveRDS(seurat_obj, "seurat_clustermole_annotated.rds")

cat("\nOutputs saved:\n")
cat("  clustermole_overlaps_full.csv\n")
cat("  clustermole_top3_per_db.csv\n")
cat("  clustermole_consensus.csv\n")
cat("  meta_consensus_approachA_clustermole.csv\n")
cat("  seurat_clustermole_annotated.rds\n")

## ============================================================
## END
## ============================================================