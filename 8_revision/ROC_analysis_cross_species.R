# =============================================================================
# ROC analysis: mesenchymal cell marker evaluation across species
#
# Evaluates single-gene and paired-gene markers for distinguishing mesenchymal
# (SMC) cells from other cell types across five datasets (human carotid,
# human coronary, pig aorta, pig coronary, mouse aorta).
#
# Dependencies: MarkerPairEval_v3.R, SingleMarkerEval_v3.R (sourced below)
#
# Input:  whole-plaque Seurat objects with all cell types (.rds)
# Output: paired-gene ROC results (.rds + .xlsx), single-gene results,
#         threshold comparison table
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(pROC)
  library(scCustomize)
  library(purrr)
  library(patchwork)
  library(ggplot2)
})

source("MarkerPairEval_v3.R")
source("SingleMarkerEval_v3.R")

# =============================================================================
# Helper: assign broad cell-type labels ("superduperclusters")
# =============================================================================

assign_broad_labels <- function(obj, cluster_col, label_map) {
  Idents(obj) <- cluster_col
  names(label_map) <- levels(Idents(obj))
  obj <- RenameIdents(obj, label_map)
  obj$superduperclusters <- Idents(obj)
  obj
}

# =============================================================================
# 1. Load whole-plaque Seurat objects (all cell types)
# =============================================================================

human <- readRDS("~/THOR/diana/integration_pigs_with_WT/humans/Carotids_Dani_human/coronaries_carotids_processed_2.rds")

human_carotid  <- human[, human$author != "wirka"]
human_coronary <- human[, human$author == "wirka"]
rm(human)

DefaultAssay(human_carotid)  <- "RNA"
human_carotid  <- NormalizeData(human_carotid)

DefaultAssay(human_coronary) <- "RNA"
human_coronary <- NormalizeData(human_coronary)

pig_aorta    <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SC.Analysis.V4.ThirdPass.RNA.Aorta.All.seuratSC.Final.rds")
pig_coronary <- readRDS("~/THOR/diana/integration_pigs_with_WT/Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SC.Analysis.V4.ThirdPass.RNA.RCA.All.seuratSC.Final.rds")
mouse        <- readRDS("~/THOR/diana/integration_pigs_with_WT/hdfc/SC.Analysis.V4.FirstPass.RNA.W4.MT15NbW034.W4.NoDoublets.Final.rds")

# =============================================================================
# 2. Assign broad cell-type labels per dataset
#    Each vector maps cluster levels (in order) to one of:
#    Lymphocytes, Mesenchymal, Macrophages, ECs
# =============================================================================

human_carotid <- assign_broad_labels(
  human_carotid, "superclusters",
  c("Lymphocytes", "Mesenchymal", "Macrophages", "Macrophages", "ECs",
    "Mesenchymal", "ECs", "Mesenchymal", "Lymphocytes", "Lymphocytes",
    "Lymphocytes", "Mesenchymal", "Lymphocytes", "Lymphocytes", "Macrophages", "ECs")
)

human_coronary <- assign_broad_labels(
  human_coronary, "superclusters",
  c("Lymphocytes", "Mesenchymal", "Macrophages", "Macrophages", "ECs",
    "Mesenchymal", "ECs", "Mesenchymal", "Lymphocytes", "Lymphocytes",
    "Lymphocytes", "Mesenchymal", "Lymphocytes", "Lymphocytes", "Macrophages", "ECs")
)

pig_aorta <- assign_broad_labels(
  pig_aorta, "res.Final",
  c("Macrophages", "Mesenchymal", "Mesenchymal", "Lymphocytes", "Macrophages",
    "Lymphocytes", "ECs", "Lymphocytes", "Macrophages", "ECs", "Mesenchymal")
)

pig_coronary <- assign_broad_labels(
  pig_coronary, "res.Final",
  c("Mesenchymal", "ECs", "Macrophages", "Mesenchymal", "Lymphocytes",
    "ECs", "Macrophages", "ECs", "ECs", "Lymphocytes")
)

mouse <- assign_broad_labels(
  mouse, "res.Final",
  c("Mesenchymal", "Mesenchymal", "Mesenchymal", "Macrophages",
    "Mesenchymal", "ECs", "Mesenchymal")
)

# =============================================================================
# 3. Paired-gene ROC analysis across all datasets
#    For each dataset: find mesenchymal DE markers, then evaluate all gene
#    pairs (including ACTA2 + LUM) for classification performance.
# =============================================================================

my_ds <- list(
  h_car = human_carotid,
  h_cor = human_coronary,
  p_ao  = pig_aorta,
  p_cor = pig_coronary,
  mouse = mouse
)

my_res <- list()

for (k in names(my_ds)) {
  message("Processing: ", k)
  seurat_obj <- my_ds[[k]]
  Idents(seurat_obj) <- "superduperclusters"

  # Find mesenchymal marker genes (DE vs all other cell types)
  tmp <- FindMarkers(
    seurat_obj,
    ident.1  = "Mesenchymal",
    only.pos = TRUE,
    assay    = "RNA",
    group.by = "superduperclusters"
  )
  tmp <- tmp[tmp$avg_log2FC > 0.6 & tmp$pct.1 > 0.3, ]
  tmp <- tmp[order(tmp$avg_log2FC, decreasing = TRUE), ]

  # Further filter by mean expression in mesenchymal cells
  mes_cells <- seurat_obj[, seurat_obj$superduperclusters == "Mesenchymal"]
  tmp$exp <- matrixStats::rowMeans2(
    as.matrix(expm1(mes_cells[["RNA"]]@data[rownames(tmp), ]))
  )
  tmp <- tmp[tmp$exp > 2, ]

  # Evaluate all gene pairs (always include ACTA2 and LUM)
  genes <- unique(c("ACTA2", "LUM", rownames(tmp)))
  best_genes <- MarkerPairEval(
    genes             = genes,
    labels            = seurat_obj$superduperclusters == "Mesenchymal",
    expression_matrix = seurat_obj@assays$RNA@data,
    combination       = "OR",
    hard_threshold    = 0
  )

  my_res[[k]] <- best_genes
  message("  Found ", nrow(best_genes), " gene pairs")
}

saveRDS(my_res, "lum.acta2.claude.thr0.OR.4.rds")
openxlsx::write.xlsx(x = my_res, file = "lum.acta2.claude.thr0.OR.4.xlsx")
message("Paired-gene results saved.")

# =============================================================================
# 4. Single-gene ROC analysis (mouse only)
# =============================================================================

Idents(mouse) <- "superduperclusters"

markers <- FindMarkers(
  mouse,
  ident.1  = "Mesenchymal",
  only.pos = TRUE,
  assay    = "RNA"
)
markers <- markers[markers$avg_log2FC > 0.6 & markers$pct.1 > 0.3, ]

single_gene_results <- SingleMarkerEval(
  genes             = rownames(markers),
  labels            = mouse$superduperclusters == "Mesenchymal",
  expression_matrix = mouse@assays$RNA@data,
  hard_threshold    = 0
)

saveRDS(single_gene_results, "mouse_single_gene_metrics.rds")
openxlsx::write.xlsx(x = single_gene_results, file = "mouse_single_gene_metrics.xlsx")
message("Single-gene results saved.")

# =============================================================================
# 5. Compare ROC-optimal vs hard thresholds for LUM / ACTA2
#    Pre-computed ROC-optimal thresholds per dataset
# =============================================================================

thresholds <- c(h_car = 1.74, h_cor = 1.92, p_ao = 2.3,
                p_cor = 2.99, mouse = 2.88)

genes_to_check <- c("LUM", "ACTA2")
results <- data.frame()

for (k in names(my_ds)) {
  expr_matrix   <- my_ds[[k]]@assays$RNA@data
  roc_threshold <- thresholds[k]

  for (gene in genes_to_check) {
    if (!gene %in% rownames(expr_matrix)) next

    expr    <- as.numeric(expr_matrix[gene, ])
    n_total <- length(expr)
    n_pos_roc  <- sum(expr > roc_threshold)
    n_pos_hard <- sum(expr > 0)

    results <- rbind(results, data.frame(
      Dataset        = k,
      Gene           = gene,
      n_total        = n_total,
      ROC.Threshold  = roc_threshold,
      n_pos_ROC      = n_pos_roc,
      pct_pos_ROC    = round(n_pos_roc  / n_total * 100, 2),
      Hard.Threshold = 0,
      n_pos_Hard     = n_pos_hard,
      pct_pos_Hard   = round(n_pos_hard / n_total * 100, 2)
    ))
  }
}

print(results)
