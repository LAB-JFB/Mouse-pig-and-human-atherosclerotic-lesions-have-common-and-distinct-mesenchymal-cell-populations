# =============================================================================
# 01 - Processing of Mouse Alencar (Owens) Dataset
# =============================================================================
# This script processes the mouse Alencar/Owens scRNA-seq dataset (GSE150644):
#   1. Load 10X data and create Seurat objects
#   2. QC metrics and visualisation
#   3. Doublet identification (scDblFinder)
#   4. QC filtering, cleaning, and feature filtering
#   5. RPCA integration across samples
#   6. Clustering and marker identification
#   7. SMC phenotype annotation and final visualisation
# =============================================================================

# --- Setup -------------------------------------------------------------------

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(sctransform)
library(tidyverse)
library(MatrixGenerics)
library(openxlsx)
library(RColorBrewer)
library(MAST)
library(psych)
library(MetBrewer)
library(data.table)
library(pheatmap)
library(corrplot)
library(ggpubr)
library(scCustomize)

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse")

# --- Load 10X data (Alencar GSE150644) ---------------------------------------

smc_plaque_85 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555585_KLF4_WTeYFPpos_filtered_gene_bc_matrices/mm10witheyfp',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
smc_plaque_86 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555586_KLF4WTeYFPpos_filtered_gene_bc_matrices/mm10witheyfp',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
smc_plaque_92 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555592_SMC_WT_eYFP_POS_filtered_gene_bc_matrices/mm10witheyfp',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
smc_plaque_93 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555593_01_green_filtered_gene_bc_matrices/mm10',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
smc_plaque_94 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555594_02_green_filtered_gene_bc_matrices/mm10',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
# smc_plaque_99 <- Read10X(
#   '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555599_03_bca_filtered_gene_bc_matrices/mm10',
#   gene.column = 2, cell.column = 1,
#   unique.features = TRUE, strip.suffix = FALSE
# )
smc_plaque_95 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555595_03_gfp_filtered_gene_bc_matrices/mm10',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
smc_plaque_96 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555596_01_tdtomato_filtered_gene_bc_matrices/mm10',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
smc_plaque_97 <- Read10X(
  '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555597_03_tomato_filtered_gene_bc_matrices/mm10',
  gene.column = 2, cell.column = 1,
  unique.features = TRUE, strip.suffix = FALSE
)
# smc_plaque_98 <- Read10X(
#   '/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/GSM4555598_01_bca_unsorted_01_filtered_gene_bc_matrices/mm10',
#   gene.column = 2, cell.column = 1,
#   unique.features = TRUE, strip.suffix = FALSE
# )

# --- Create Seurat objects ---------------------------------------------------

smc_plaque_85 <- CreateSeuratObject(counts = smc_plaque_85, project = "Klf4_WT_sorted_SMCs", min.cells = 3, min.features = 200) #94 cells
smc_plaque_86 <- CreateSeuratObject(counts = smc_plaque_86, project = "Klf4_WT_sorted_SMCs", min.cells = 3, min.features = 200)
smc_plaque_92 <- CreateSeuratObject(counts = smc_plaque_92, project = "SMC_WT_eYFP_sorted", min.cells = 3, min.features = 200)
smc_plaque_93 <- CreateSeuratObject(counts = smc_plaque_93, project = "gfp_SMC_sorted_cells", min.cells = 3, min.features = 200) #64
smc_plaque_94 <- CreateSeuratObject(counts = smc_plaque_94, project = "gfp_SMC_sorted_cells", min.cells = 3, min.features = 200) #72
smc_plaque_95 <- CreateSeuratObject(counts = smc_plaque_95, project = "gfp_SMC_sorted_cells", min.cells = 3, min.features = 200)
smc_plaque_96 <- CreateSeuratObject(counts = smc_plaque_96, project = "tdTomato_SMC_sorted_cells01", min.cells = 3, min.features = 200)
smc_plaque_97 <- CreateSeuratObject(counts = smc_plaque_97, project = "tdTomato_SMC_sorted_cells02", min.cells = 3, min.features = 200)
# smc_plaque_98 <- CreateSeuratObject(counts = smc_plaque_98, project = "bca_unsorted_01", min.cells = 3, min.features = 200)
# smc_plaque_99 <- CreateSeuratObject(counts = smc_plaque_99, project = "bca_unsorted_03", min.cells = 3, min.features = 200)

# --- Merge and QC metrics ----------------------------------------------------

smc_plaque <- merge(
  x = smc_plaque_85,
  y = c(smc_plaque_86, smc_plaque_92, smc_plaque_93, smc_plaque_94,
        smc_plaque_95, smc_plaque_96, smc_plaque_97
        # smc_plaque_98,
        # smc_plaque_99
  )
)

smc_plaque[["percent.mt"]] <- PercentageFeatureSet(smc_plaque, pattern = "^mt-")
smc_plaque[["percent.ribo"]] <- PercentageFeatureSet(smc_plaque, pattern = "^Rp[sl]")
smc_plaque[["percent.hb"]] <- PercentageFeatureSet(smc_plaque, pattern = "^Hb[^(p)]")
smc_plaque[["percent.Malat1"]] <- PercentageFeatureSet(smc_plaque, pattern = "^Malat1")

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt",
           "percent.ribo", "percent.hb", "percent.Malat1")

dim(smc_plaque)

# --- QC visualisation --------------------------------------------------------

VlnPlot(smc_plaque, features = feats, ncol = 3)

plot1 <- FeatureScatter(smc_plaque, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(smc_plaque, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# --- Doublet identification (scDblFinder) ------------------------------------


require(DoubletFinder)
require(scDblFinder)

seu_split <- SplitObject(smc_plaque, split.by = "orig.ident")
seu_split <- map(seu_split, function(x) {
  sce <- as.SingleCellExperiment(x, assay = "RNA")
  sce <- scDblFinder(sce, sce$seurat_clusters, dbr.sd = 1)
  scdbl_cols <- grep("scDblFinder", colnames(colData(sce)), value = TRUE)
  x@meta.data <- cbind(x@meta.data, colData(sce)[, scdbl_cols])
  return(x)
})

# --- Pre-filtering normalisation and embedding (with doublets) ---------------

smc_plaque <- Merge_Seurat_List(seu_split)
smc_plaque <- NormalizeData(smc_plaque,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
smc_plaque <- FindVariableFeatures(smc_plaque,
                                   selection.method = "vst",
                                   nfeatures = 2000)
all.genes <- rownames(smc_plaque)
smc_plaque <- ScaleData(smc_plaque, features = all.genes)



smc_plaque <-
  RunPCA(smc_plaque, npcs = 40, verbose = FALSE) %>%
  RunUMAP(reduction = "pca",
          dims = 1:22,
          verbose = FALSE) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:28,
                verbose = FALSE) %>%
  FindClusters(resolution = c(0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5, 0.7, 1),
               verbose = FALSE)

DimPlot_scCustom(smc_plaque, group.by = "RNA_snn_res.0.35") +
  DimPlot_scCustom(smc_plaque, group.by = "scDblFinder.class")

saveRDS(smc_plaque, file = "Alencar_SMCs_withDoublets.rds")

# --- QC filtering and cleaning -----------------------------------------------

smc_plaque <-
  smc_plaque[, smc_plaque$nFeature_RNA > 800 &
               smc_plaque$percent.mt < 15 &
               smc_plaque$nCount_RNA < 25000 &
               smc_plaque$percent.hb < 1 &
               smc_plaque$scDblFinder.class == "singlet"]

# Remove contaminating immune cells
smc_plaque <- subset(smc_plaque, Ptprc == 0 & Tyrobp == 0)

# Remove ribosomal, mitochondrial, and MALAT1 genes
smc_plaque <- smc_plaque[!grepl("^Rp[sl]", rownames(smc_plaque)), ]
smc_plaque <- smc_plaque[!grepl("^mt-", rownames(smc_plaque)), ]
smc_plaque <- smc_plaque[!grepl("Malat1", rownames(smc_plaque)), ]

# --- Gene filtering ----------------------------------------------------------

# Filter genes expressed in fewer than 5 cells
counts <- GetAssayData(object = smc_plaque, slot = "counts", assay = "RNA")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 5
filtered_counts <- counts[keep_genes, ]

smc_plaque <- CreateSeuratObject(filtered_counts, meta.data = smc_plaque@meta.data)

# --- Normalisation -----------------------------------------------------------

all.genes <- rownames(smc_plaque)
smc_plaque <- NormalizeData(smc_plaque,
                            normalization.method = "LogNormalize",
                            scale.factor = 10000)
smc_plaque <- FindVariableFeatures(smc_plaque,
                                   selection.method = "vst",
                                   nfeatures = 2000)
smc_plaque <- ScaleData(smc_plaque, features = all.genes)

# --- Integration (RPCA) ------------------------------------------------------

ifnb.list <- SplitObject(smc_plaque, split.by = "orig.ident")

# Normalize and identify variable features per sample
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# Select integration features and run PCA per sample
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Find anchors and integrate
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list,
                                         anchor.features = features,
                                         reduction = "rpca",
                                         k.anchor = 5)
smc_plaque <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(smc_plaque) <- "integrated"

# --- Dimensionality reduction and clustering ---------------------------------

smc_plaque <- ScaleData(smc_plaque, verbose = FALSE)
smc_plaque <- RunPCA(smc_plaque, npcs = 40, verbose = FALSE)

# Determine significant PCs
stdv <- smc_plaque[["pca"]]@stdev
sum.stdv <- sum(smc_plaque[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 85 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                     percent.stdv[2:length(percent.stdv)]) > 0.1),
            decreasing = T)[1] + 1

print(co1)
print(co2)

# UMAP, neighbours, and multi-resolution clustering
smc_plaque <-
  RunUMAP(object = smc_plaque, reduction = "pca",
          dims = 1:co1, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:co1, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5, 0.7, 1),
               verbose = FALSE)

# Fine-tuning resolutions
smc_plaque <- FindClusters(object = smc_plaque,
                           resolution = c(0.21, 0.22, 0.23, 0.24))

# --- Clustering visualisation ------------------------------------------------

DimPlot_scCustom(smc_plaque, group.by = "orig.ident") +
  DimPlot_scCustom(smc_plaque, group.by = "integrated_snn_res.0.24")

DimPlot_scCustom(smc_plaque, label = TRUE,
                 group.by = "integrated_snn_res.0.24", repel = TRUE,
                 reduction = "umap", pt.size = 1, label.size = 7, shuffle = TRUE) +
  DimPlot_scCustom(smc_plaque, label = TRUE,
                   group.by = "integrated_snn_res.0.3", repel = TRUE,
                   reduction = "umap", pt.size = 1, label.size = 7, shuffle = TRUE) +
  DimPlot_scCustom(smc_plaque, label = FALSE,
                   group.by = "orig.ident", repel = TRUE,
                   reduction = "umap", pt.size = 1, label.size = 7, shuffle = TRUE)

# --- Marker identification ---------------------------------------------------

smc_plaque$res.final.0.7 <- paste0("C", smc_plaque$integrated_snn_res.0.7)
smc_plaque$res.final.0.24 <- paste0("C", smc_plaque$integrated_snn_res.0.24)
smc_plaque$res.final.0.3 <- paste0("C", smc_plaque$integrated_snn_res.0.3)

smc_plaque$res.final.0.7 <- factor(smc_plaque$res.final.0.7)
smc_plaque$res.final.0.24 <- factor(smc_plaque$res.final.0.24)
smc_plaque$res.final.0.3 <- factor(smc_plaque$res.final.0.3)

# Markers at resolution 0.7
Idents(smc_plaque) <- "res.final.0.7"
markers_0.7 <- FindAllMarkers(smc_plaque, assay = "RNA",
                               logfc.threshold = 0.25, test.use = "MAST",
                               min.pct = 0.3, only.pos = TRUE,
                               max.cells.per.ident = 1000)
markers_0.7 <- markers_0.7[markers_0.7$p_val_adj < 0.05, ]

# Markers at resolution 0.24
Idents(smc_plaque) <- "res.final.0.24"
markers_0.3 <- FindAllMarkers(smc_plaque, assay = "RNA",
                               logfc.threshold = 0.25, test.use = "MAST",
                               min.pct = 0.3, only.pos = TRUE,
                               max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05, ]

# --- Dot plot (resolution 0.24) ----------------------------------------------

top_markers <- c()
for (i in unique(markers_0.3$cluster)) {
  cluster_markers <- subset(markers_0.3, subset = markers_0.3$cluster == i)
  cluster_markers <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = TRUE), ]
  top_markers <- c(top_markers, head(cluster_markers$gene, 10))
}

p1 <- DotPlot(smc_plaque,
              features = unique(top_markers),
              assay = "RNA", scale = TRUE,
              group.by = "res.final.0.24") +
  ggtitle("Top markers for SMC phenotypes res 0.24") +
  coord_flip() +
  RotatedAxis() +
  scale_colour_gradientn(name = "log2 (count + 1)",
                         colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  cowplot::theme_cowplot() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank()) +
  ylab("")
print(p1)

# --- SMC phenotype annotation ------------------------------------------------

Idents(smc_plaque) <- "res.final.0.24"
smc_phenotypes <- c("Transitional", "Terminal", "Contractile", "Contractile_Rgs5")
names(smc_phenotypes) <- levels(Idents(smc_plaque))
smc_plaque <- RenameIdents(smc_plaque, smc_phenotypes)
smc_plaque$smc_phenotypes <- Idents(smc_plaque)

# Short names
Idents(smc_plaque) <- "smc_phenotypes"
new_names <- c("Trans", "Term", "SMC1", "SMC2")
names(new_names) <- levels(Idents(smc_plaque))
smc_plaque <- Seurat::RenameIdents(object = smc_plaque, new_names)
smc_plaque$smc_phenotypes_short <- Idents(smc_plaque)

# --- Markers by phenotype ----------------------------------------------------

Idents(smc_plaque) <- "smc_phenotypes"
markers_0.3 <- FindAllMarkers(smc_plaque, assay = "RNA",
                               logfc.threshold = 0.25, test.use = "MAST",
                               min.pct = 0.3, only.pos = TRUE,
                               max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05, ]

top_markers <- c()
for (i in unique(markers_0.3$cluster)) {
  cluster_markers <- subset(markers_0.3, subset = markers_0.3$cluster == i)
  cluster_markers <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = TRUE), ]
  top_markers <- c(top_markers, head(cluster_markers$gene, 6))
}

# --- Final visualisations ----------------------------------------------------

# NOTE: pal_celltypes_mouse must be defined BEFORE use
pal_celltypes_mouse <- c("#E18727FF", "#BC3C29FF", "#0072B5FF", "#79AF97FF")

p1 <- DotPlot(smc_plaque,
              features = unique(top_markers),
              assay = "RNA", scale = TRUE,
              group.by = "smc_phenotypes") +
  ggtitle("Top markers for SMC phenotypes in mouse plaques") +
  scale_colour_gradientn(name = "log2 (count + 1)",
                         colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank()) +
  ylab("")
print(p1)

ggsave(filename = "dotplot_Alencar_0.24.pdf",
       plot = p1, device = "pdf", width = 8, height = 4,
       units = "cm", scale = 3)

p2 <- DimPlot_scCustom(seurat_object = smc_plaque,
                        reduction = "umap", group.by = "smc_phenotypes",
                        pt.size = 1, shuffle = TRUE, label = TRUE,
                        label.box = TRUE, colors_use = pal_celltypes_mouse,
                        repel = TRUE, label.size = 5) +
  NoLegend() +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes_mouse) +
  theme(axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
print(p2)

ggsave(filename = "dimplot_mouse_Alencar_smc_phenotypes_short.pdf",
       plot = p2, device = "pdf", width = 3, height = 3,
       units = "cm", scale = 3)


p2 <- DimPlot_scCustom(seurat_object = smc_plaque,
                        reduction = "umap", group.by = "orig.ident",
                        pt.size = 1, shuffle = TRUE,
                        colors_use = pal_celltypes, repel = TRUE) +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes) +
  theme(axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
print(p2)

ggsave(filename = "dimplot_mouse_Alencar_origin.pdf",
       plot = p2, device = "pdf", width = 5, height = 3,
       units = "cm", scale = 3)

# --- Add metadata and save ---------------------------------------------------

smc_plaque$species <- "mmusculus"
smc_plaque$Origin <- "Mouse"
smc_plaque$dataset <- smc_plaque$orig.ident
smc_plaque$rough_org <- paste("M", smc_plaque$smc_phenotypes, sep = "_")
smc_plaque$fine_org <- paste("M", smc_plaque$smc_phenotypes, sep = "_")

saveRDS(smc_plaque, file = "smc_plaque_Alencar_seurat.rds")
