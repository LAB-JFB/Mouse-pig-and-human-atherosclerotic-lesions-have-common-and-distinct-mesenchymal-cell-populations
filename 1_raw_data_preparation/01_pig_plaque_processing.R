# ==============================================================================
# 01_pig_plaque_processing.R
#
# Pipeline: Load pig aorta scRNA-seq data -> Subset SMCs -> QC filter ->
#           Normalize & cluster -> Annotate phenotypes -> Add cross-species
#           metadata -> Find markers -> Generate plots
# ==============================================================================

# --- Libraries ----------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(scCustomize)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(sctransform)
  library(MatrixGenerics)
  library(RColorBrewer)
  library(MAST)
  library(openxlsx)
  library(cowplot)
  library(ggpubr)
  library(MetBrewer)
  library(ggrepel)
})

set.seed(1984)

# --- Paths --------------------------------------------------------------------

base_dir    <- "/faststorage/project/THOR/diana/integration_pigs_with_WT"
smc_dir     <- file.path(base_dir, "Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SMC")
output_dir  <- file.path(base_dir, "integrated_pig_mice_hum/closer/more_closer/pig")
palette_dir <- file.path(base_dir, "integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load all-cell-type pig aorta data -------------------------------------

pigs_all <- readRDS(file.path(smc_dir,
  "../SC.Analysis.V4.ThirdPass.RNA.Aorta.All.seuratSC.Final.rds"))

# --- 2. Subset SMCs -----------------------------------------------------------
# Keep ACTA2+/MYH11+ clusters (C1, C2, C10), apply QC thresholds,
# remove immune marker-positive cells and ribosomal protein genes

smc_data <- subset(pigs_all,
  res.Final %in% c("C1", "C2", "C10") &
    nFeature_RNA > 900 &
    nCount_RNA > 1200)
smc_data <- subset(smc_data, PTPRC == 0 & TYROBP == 0)
smc_data <- smc_data[!grepl("^RP[SL]", rownames(smc_data)), ]

rm(pigs_all)

# Remove genes expressed in fewer than 5 cells
counts <- GetAssayData(smc_data, slot = "counts", assay = "RNA")
keep_genes <- Matrix::rowSums(counts > 0) >= 5
smc_data <- CreateSeuratObject(counts[keep_genes, ], meta.data = smc_data@meta.data)
smc_data$Condition_Sample <- paste(smc_data$Condition, smc_data$Sample, sep = ".")

# --- 3. Normalize, reduce dimensions, cluster --------------------------------

all_genes <- rownames(smc_data)
smc_data <- NormalizeData(smc_data, normalization.method = "LogNormalize", scale.factor = 10000)
smc_data <- FindVariableFeatures(smc_data, selection.method = "vst", nfeatures = 2000)
smc_data <- ScaleData(smc_data, features = all_genes)

smc_data <- RunPCA(smc_data, npcs = 40, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:22, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:28, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5, 0.7, 1),
               verbose = FALSE)

# --- 4. Annotate SMC phenotypes at multiple resolutions -----------------------

smc_data$res.final.0.1 <- paste0("C", smc_data$RNA_snn_res.0.1)
smc_data$res.final.0.2 <- paste0("C", smc_data$RNA_snn_res.0.15)
smc_data$res.final.0.3 <- paste0("C", smc_data$RNA_snn_res.0.3)

# Resolution 0.3 annotations
Idents(smc_data) <- "res.final.0.3"
smc_labels_0.3 <- c("Transitional", "Contractile", "Terminal", "Pericytes",
                     "DES_SMC", "OGN_SMC", "IFN_SMC")
names(smc_labels_0.3) <- levels(Idents(smc_data))
smc_data <- RenameIdents(smc_data, smc_labels_0.3)
smc_data$smc_0.3 <- Idents(smc_data)

# Resolution 0.15 annotations
Idents(smc_data) <- "res.final.0.2"
smc_labels_0.15 <- c("Modulated", "Contractile", "Pericytes", "Medial",
                      "Proliferative", "OGN_SMC")
names(smc_labels_0.15) <- levels(Idents(smc_data))
smc_data <- RenameIdents(smc_data, smc_labels_0.15)
smc_data$smc_0.15 <- Idents(smc_data)

# Resolution 0.1 annotations
Idents(smc_data) <- "res.final.0.1"
smc_labels_0.1 <- c("Modulated", "Contractile", "Pericytes", "OGN_SMC")
names(smc_labels_0.1) <- levels(Idents(smc_data))
smc_data <- RenameIdents(smc_data, smc_labels_0.1)
smc_data$smc_0.1 <- Idents(smc_data)

saveRDS(smc_data, file.path(smc_dir, "SMC.aorta.processed.RNA.norm.rds"))

# --- 5. Prepare plaque subset (HFD only) with cross-species metadata ---------

pig_plaque <- smc_data[, smc_data$Condition == "HFD"]

# Metadata for cross-species integration
pig_plaque$fine_smc       <- pig_plaque$smc_0.3
pig_plaque$rough_smc      <- pig_plaque$smc_0.15
pig_plaque$dataset         <- pig_plaque$Sample
pig_plaque$species         <- "sscrofa"
pig_plaque$Origin          <- "Pig"
pig_plaque$Dataset_2       <- "pig_plaque"
pig_plaque$Artery_type     <- "Aorta_plaque"
pig_plaque$dataset_author  <- "Bentzon, pig"
pig_plaque$rough_org       <- paste("P", pig_plaque$rough_smc, sep = "_")
pig_plaque$fine_org        <- paste("P", pig_plaque$fine_smc, sep = "_")
pig_plaque$percent.mt      <- pig_plaque$subsets_MT_percent

# Rename phenotypes for cross-species labels
DefaultAssay(pig_plaque) <- "RNA"
Idents(pig_plaque) <- "smc_0.3"
phenotype_labels <- c("Transitional", "Contractile", "Terminal", "Pericytes",
                       "Medial", "Fibroblasts", "IFN_SMC")
names(phenotype_labels) <- levels(Idents(pig_plaque))
pig_plaque <- RenameIdents(pig_plaque, phenotype_labels)
pig_plaque$smc_phenotypes <- Idents(pig_plaque)

saveRDS(pig_plaque, file.path(output_dir, "pig_plaque_seurat.rds"))

# --- 6. Find markers and generate plots --------------------------------------

pal_celltypes_pig <- c("#E18727FF", "#0072B5FF", "#BC3C29FF", "#20854EFF",
                       "#7876B1FF", "#FFDC91FF", "#EE4C97FF")

Idents(pig_plaque) <- "smc_phenotypes"
markers <- FindAllMarkers(pig_plaque,
  assay = "RNA", logfc.threshold = 0.25, test.use = "MAST",
  min.pct = 0.3, only.pos = TRUE, max.cells.per.ident = 1000)
markers <- markers[markers$p_val_adj < 0.05, ]

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 6) %>%
  pull(gene) %>%
  unique()

# DotPlot of top markers
p_dot <- DotPlot(pig_plaque,
  features = top_markers, assay = "RNA", scale = TRUE,
  group.by = "smc_phenotypes") +
  ggtitle("Top markers for SMC phenotypes — pig aorta plaque") +
  scale_colour_gradientn(
    name = "log2 (count + 1)",
    colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme_bw() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("")

ggsave(file.path(output_dir, "dotplot_pig_smc_phenotypes.pdf"),
       plot = p_dot, device = "pdf", width = 8, height = 4, units = "cm", scale = 3)

# UMAP DimPlot
p_umap <- DimPlot_scCustom(pig_plaque,
  reduction = "umap", group.by = "smc_phenotypes",
  pt.size = 1, shuffle = TRUE, label = TRUE, label.box = TRUE,
  colors_use = pal_celltypes_pig, repel = TRUE, label.size = 5) +
  NoLegend() +
  scale_color_manual(values = pal_celltypes_pig) +
  theme(
    plot.title = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank())

ggsave(file.path(output_dir, "dimplot_pig_smc_phenotypes.pdf"),
       plot = p_umap, device = "pdf", width = 4, height = 4, units = "cm", scale = 3)

message("Done. Outputs saved to: ", output_dir)