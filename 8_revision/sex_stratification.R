# =============================================================================
# Sex stratification analysis of SMC phenotypes across species
#
# Infers biological sex per donor via pseudobulk expression of X/Y-linked
# marker genes, then visualises SMC phenotype UMAPs split by sex and runs
# propeller cell-type proportion testing.
#
# Input:  per-species SMC Seurat objects (.rds)
# Output: sex-split UMAP PDFs, intermediate Seurat objects with sex labels
# =============================================================================

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC")
set.seed(1984)

# --- Libraries ---------------------------------------------------------------
library(Seurat)
library(scCustomize)
library(sctransform)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(Matrix)
library(MatrixGenerics)
library(openxlsx)
library(RColorBrewer)
library(MetBrewer)
library(MAST)
library(DESeq2)
library(psych)
library(data.table)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(mclust)
library(speckle)
library(limma)

# =============================================================================
# 1. Load per-species SMC Seurat objects
# =============================================================================

mouse <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/smc_plaque_Alencar_seurat.rds")
human <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/human/athero_comb_smc_seurat.rds")
pig   <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SMC/bengal/pig_plaque_seurat.rds")

human_cor <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/human_coronaries_seurat.rds")
pig_cor   <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/pig_coronaries_seurat.rds")

# =============================================================================
# 2. Annotate known sex metadata
# =============================================================================

# Mouse: all male (Alencar et al.)
mouse$sex <- "Male"

# Pig: sex stored in Gender column
pig$sex     <- pig$Gender
pig_cor$sex <- pig_cor$Gender

# Human coronary (Wirka et al.): per-sample sex from GEO metadata
human_cor_sex_map <- c(
  GSM3819856 = "Male",   GSM3819857 = "Male",
  GSM3819858 = "Male",   GSM3819859 = "Male",
  GSM3819860 = "Male",   GSM3819861 = "Male",
  GSM3819862 = "Male",   GSM3819863 = "Female"
)

# Human carotid (Pan et al.): partially known
human_car_sex_map <- c(
  GSM4705589 = "Male",   GSM4705590 = "Male",
  GSM4705591 = "Female",
  GSM4837523 = NA,        GSM4837525 = NA,        GSM4837527 = NA
)

# Map sex to cells via sample_id (proper per-cell lookup, not recycled ifelse)
human$sex     <- human_car_sex_map[human$sample_id]
human_cor$sex <- human_cor_sex_map[human_cor$sample_id]

# =============================================================================
# 3. Standardise sample_id and species columns
# =============================================================================

mouse$sample_id <- mouse$orig.ident
pig$sample_id   <- pig$Sample
pig_cor$sample_id <- pig_cor$Sample

mouse$species   <- "mouse"
pig$species     <- "pig"
pig_cor$species <- "pig"
human$species   <- "human"
human_cor$species <- "human"

# =============================================================================
# 4. Sex-linked marker genes (per species)
# =============================================================================

sex_markers <- list(
  human = list(
    X = c("XIST"),
    Y = c("RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y", "ZFY", "TXLNGY")
  ),
  mouse = list(
    X = c("Xist"),
    Y = c("Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Usp9y", "Zfy1", "Zfy2", "Uba1y", "Rps4y1")
  ),
  pig = list(
    X = c("XIST"),
    Y = c("RPS4Y1", "DDX3Y", "EIF1AY", "KDM5D", "UTY", "USP9Y", "ZFY",
           "LOC100624590", "LOC100624149")
  )
)

# --- Helper: case-insensitive gene matching ----------------------------------
.match_genes <- function(genes, feature_names) {
  feature_names[tolower(feature_names) %in% tolower(genes)]
}

# =============================================================================
# 5. Infer sex per donor via pseudobulk X/Y scores
# =============================================================================

infer_sex_by_pseudobulk <- function(
    seu,
    donor_col       = "sample_id",
    species_col     = "species",
    assay           = "RNA",
    use_cells       = NULL,
    min_cells_per_donor = 100,
    use_mclust      = FALSE
) {
  stopifnot(donor_col   %in% colnames(seu@meta.data),
            species_col %in% colnames(seu@meta.data))

  if (!is.null(use_cells)) {
    seu <- subset(seu, cells = intersect(colnames(seu), use_cells))
  }

  DefaultAssay(seu) <- assay
  counts <- GetAssayData(seu, assay = assay, slot = "counts")
  md     <- seu@meta.data

  # Filter donors with too few cells

  donor_tab   <- table(md[[donor_col]])
  donors_keep <- names(donor_tab)[donor_tab >= min_cells_per_donor]
  md     <- md[md[[donor_col]] %in% donors_keep, , drop = FALSE]
  counts <- counts[, rownames(md), drop = FALSE]

  # Pseudobulk: sum counts per gene per donor
  split_cells <- split(rownames(md), md[[donor_col]])
  pb_list <- lapply(split_cells, function(cells) {
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  pb <- do.call(cbind, pb_list)
  colnames(pb) <- names(split_cells)

  # CPM normalisation
  libsize <- Matrix::colSums(pb)
  cpm     <- t(t(pb) / pmax(libsize, 1)) * 1e6

  # Species per donor
  donor_species <- vapply(names(split_cells), function(d) {
    unique(md[md[[donor_col]] == d, species_col])[1]
  }, character(1))

  # Compute XIST and ChrY scores per donor
  donors <- colnames(cpm)

  get_marker_score <- function(gset, M) {
    g <- .match_genes(gset, rownames(M))
    if (length(g) == 0) return(rep(0, ncol(M)))
    colMeans(log1p(M[g, , drop = FALSE]))
  }

  xist_score <- numeric(length(donors))
  y_score    <- numeric(length(donors))

  for (i in seq_along(donors)) {
    d  <- donors[i]
    sp <- tolower(donor_species[[d]])
    sp <- ifelse(sp %in% names(sex_markers), sp, "human")
    xist_score[i] <- get_marker_score(sex_markers[[sp]]$X, cpm)[i]
    y_score[i]    <- get_marker_score(sex_markers[[sp]]$Y, cpm)[i]
  }

  # Donor-level summary table
  df <- data.frame(
    donor      = donors,
    species    = donor_species[donors],
    n_cells    = as.integer(donor_tab[donors]),
    xist_score = xist_score,
    y_score    = y_score,
    sex_score  = scale(xist_score)[, 1] - scale(y_score)[, 1],
    stringsAsFactors = FALSE
  )

  # --- Classification --------------------------------------------------------
  if (use_mclust && requireNamespace("mclust", quietly = TRUE)) {
    fit <- mclust::Mclust(cbind(df$xist_score, df$y_score), G = 2)
    cl  <- fit$classification
    cl_means    <- aggregate(df[, c("xist_score", "y_score")], list(cluster = cl), mean)
    fem_cluster <- with(cl_means, cluster[which.max(xist_score - y_score)])
    df$sex_inferred <- ifelse(cl == fem_cluster, "female", "male")
  } else {
    # Rule-based: quadrant logic on medians
    xm <- stats::median(df$xist_score, na.rm = TRUE)
    ym <- stats::median(df$y_score, na.rm = TRUE)
    df$sex_inferred <- ifelse(df$xist_score > xm & df$y_score < ym, "female",
                       ifelse(df$y_score > ym & df$xist_score < xm, "male",
                              "ambiguous"))
  }

  # Map inferred sex back to cells
  sex_map <- setNames(df$sex_inferred, df$donor)
  seu$sex_inferred <- unname(sex_map[seu@meta.data[[donor_col]]])

  list(donor_table = df, pseudobulk_cpm = cpm, seurat = seu)
}

# =============================================================================
# 6. Run sex inference on all datasets
# =============================================================================

res <- infer_sex_by_pseudobulk(human, use_mclust = TRUE)
human$sex_inferred <- res$seurat$sex_inferred[colnames(human)]

res <- infer_sex_by_pseudobulk(human_cor, use_mclust = TRUE)
human_cor$sex_inferred <- res$seurat$sex_inferred[colnames(human_cor)]

res <- infer_sex_by_pseudobulk(pig, use_mclust = TRUE)
pig$sex_inferred <- res$seurat$sex_inferred[colnames(pig)]

res <- infer_sex_by_pseudobulk(pig_cor, use_mclust = TRUE)
pig_cor$sex_inferred <- res$seurat$sex_inferred[colnames(pig_cor)]

res <- infer_sex_by_pseudobulk(mouse, use_mclust = TRUE)
mouse$sex_inferred <- res$seurat$sex_inferred[colnames(mouse)]

# Sanity checks
head(res$donor_table)
table(res$donor_table$sex_inferred)

table(human$sex_inferred,     human$orig.ident)
table(human_cor$sex_inferred, human_cor$sample_id)
table(pig$sex_inferred,       pig$id)
table(pig_cor$sex_inferred,   pig_cor$id)
table(mouse$sex_inferred,     mouse$orig.ident)

# =============================================================================
# 7. Propeller: test cell-type proportion differences by sex (pig example)
# =============================================================================

cell_info <- data.frame(
  clusters = pig$smc_phenotypes,
  sample   = pig$Sample,
  group    = pig$sex_inferred
)

propeller_res <- propeller(
  clusters  = cell_info$clusters,
  sample    = cell_info$sample,
  group     = cell_info$group,
  transform = "logit"
)

# Save intermediate objects with sex annotations
saveRDS(mouse, "mouse.clusters.rds")
saveRDS(pig,   "pig.clusters.rds")
saveRDS(human, "human.clusters.rds")

# =============================================================================
# 8. Sex-split UMAP plots
# =============================================================================

# --- Shared colour palette (elastic arteries: 3-species phenotypes) ----------
pal_elastic <- c(
  '#C1C8D9', '#8491B4df', '#4f576c',
  pal_npg("nrc")(10)[c(1, 2, 3, 10, 5, 7, 8, 9)]
)
names(pal_elastic) <- c(
  'Contractile', 'Transitional', 'Terminal',
  'Pericytes_Human_pig_specific', 'DLX5_SMC_Human_specific',
  'Terminal_Mouse_specific', 'Proliferating',
  'Contractile_Human_pig_specific', 'Too_small_cluster'
)

# --- Shared colour palette (coronaries: 2-species phenotypes) ----------------
pal_coronary <- c(
  '#C1C8D9', '#8491B4df', '#4f576c',
  pal_npg("nrc")(10)[c(1, 2, 3, 10, 5, 7, 8, 9)]
)
names(pal_coronary) <- c(
  'Contractile', 'Transitional', 'Terminal',
  'Pericytes_Human-specific', 'LBH_SMC_Human-specific', 'Fibroblasts',
  'Terminal_IGFBP2_Pig_specific', 'IFN_SMC', 'Transitional_2_Pig-specific'
)

# --- Shared ggplot theme for clean UMAPs -------------------------------------
theme_umap_clean <- theme(
  plot.title   = element_blank(),
  legend.position = 'none',
  axis.line    = element_blank(),
  axis.text    = element_blank(),
  axis.ticks   = element_blank(),
  axis.title   = element_blank()
)

# --- Helper: sex-split UMAP for one dataset ----------------------------------
plot_sex_umap <- function(obj, palette, filename) {

  # Subset palette to levels present in this object
  obj$shared_clusters <- factor(obj$shared_clusters)
  pal_use <- palette[levels(obj$shared_clusters)]
  names(pal_use) <- levels(obj$shared_clusters)

  obj$sex_inferred <- factor(tolower(obj$sex_inferred), levels = c("female", "male"))

  p <- DimPlot_scCustom(
    obj,
    reduction  = 'umap',
    group.by   = 'shared_clusters',
    split.by   = "sex_inferred",
    pt.size    = 0.5,
    shuffle    = TRUE,
    label      = TRUE,
    label.box  = TRUE,
    colors_use = pal_use,
    repel      = TRUE,
    label.size = 2.7,
    aspect_ratio = 1
  ) +
    scale_color_manual(values = pal_use) +
    theme_umap_clean

  print(p)
  ggsave(filename = filename, plot = p)
  invisible(p)
}

# --- Elastic artery UMAPs (mouse, pig, human) --------------------------------
plot_sex_umap(mouse, pal_elastic, "sex.split.mouse.pdf")
plot_sex_umap(pig,   pal_elastic, "sex.split.pig.pdf")
plot_sex_umap(human, pal_elastic, "sex.split.human.pdf")

# --- Coronary artery UMAPs (pig_cor, human_cor) ------------------------------
plot_sex_umap(pig_cor,   pal_coronary, "sex.split.pig_cor.pdf")
plot_sex_umap(human_cor, pal_coronary, "sex.split.human_cor.pdf")

# =============================================================================
# 9. ACTA2 / DLX5 co-expression blend plot (human example)
# =============================================================================

find_gene <- function(seu, gene) {
  hit <- rownames(seu)[tolower(rownames(seu)) == tolower(gene)]
  if (length(hit) == 0) stop(paste("Gene not found:", gene))
  hit[1]
}

obj <- human
g1  <- find_gene(obj, "ACTA2")
g2  <- find_gene(obj, "DLX5")

p_blend <- FeaturePlot(
  obj, features = c(g1, g2), reduction = "umap",
  blend = TRUE, blend.threshold = 0.1, pt.size = 0.4,
  order = TRUE, cols = c('magenta', 'green')
) + ggtitle(paste(g1, "&", g2, "co-expression"))

p_blend
