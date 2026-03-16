#!/usr/bin/env Rscript
# =============================================================================
# 05_coronaries_preprocessing.R
#
# Compare coronary artery plaque SMCs from pigs and humans
# =============================================================================

# --- Setup -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(cluster)
  library(Seurat)
  library(cowplot)
  library(stats)
  library(data.table)
  library(tidyverse)
  library(ggsci)
  library(scCustomize)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(sctransform)
  library(MatrixGenerics)
  library(openxlsx)
  library(RColorBrewer)
  library(MAST)
  library(psych)
  library(MetBrewer)
  library(pheatmap)
  library(corrplot)
  library(ggpubr)
  library(DESeq2)
  library(SeuratObject)
  library(viridis)
  library(scales)
  library(SCpubr)
  library(ggplotify)
  library(DoubletFinder)
  library(slingshot)
  library(gridExtra)
  library(ggvenn)
  library(ggrepel)
  library(nebula)
  library(readxl)
})
set.seed(8899)
my_pal1 <- MetBrewer::met.brewer('Egypt', n = 4, type = 'discrete')
my_pal2 <- MetBrewer::met.brewer('Signac', n = 14, type = 'discrete')
options(width = 200)
options(scipen = 999)
options(repr.matrix.max.rows = 20)
options(repr.matrix.max.cols = 50)

setwd("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries")

# =============================================================================
# Preparation of human coronary arteries (Wirka)
# =============================================================================

athero_comb <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/humans/Carotids_Dani_human/coronaries_carotids_processed.rds")

athero_comb[["percent.Malat1"]] <- PercentageFeatureSet(athero_comb, pattern = "^MALAT1")
athero_comb$mito_ribo_ratio <- athero_comb$percent_mito/(athero_comb$percent_mito + athero_comb$percent_ribo)
athero_comb$novelty_score <- log10(athero_comb$nFeature_RNA) / log10(athero_comb$nCount_RNA)
athero_comb$is_NS_good <- ifelse(athero_comb$novelty_score < 0.8, 'Bad', 'Good')

athero_comb <- subset(athero_comb,  author == 'wirka')

# --- QC plots ----------------------------------------------------------------

FeatureScatter(athero_comb, feature1 = "novelty_score", feature2 = "mito_ribo_ratio",  plot.cor = FALSE) + theme(legend.position = "bottom") + geom_hline(yintercept = 0.75, lty =2) + geom_vline(xintercept = 0.8, lty = 2)
QC_Plot_UMIvsGene(seurat_object = athero_comb, low_cutoff_gene = 250, high_cutoff_gene = 4000, low_cutoff_UMI = 800, high_cutoff_UMI = 20000, group.by = 'is_NS_good')

# --- Clustering of full dataset ----------------------------------------------

DefaultAssay(athero_comb) <- 'integrated'
athero_comb <- RunPCA(athero_comb, npcs = 30, verbose = FALSE)
athero_comb <- RunUMAP(athero_comb, reduction = "pca", dims = 1:22)
athero_comb <- FindNeighbors(athero_comb, reduction = "pca", dims = 1:22)
athero_comb <- FindClusters(athero_comb, resolution = c(0.1, 0.15,  0.2, 0.3, 0.35, 0.4, 0.45, 0.46))

# --- Cluster annotation ------------------------------------------------------

smc <-
  c(
    'Endothelial cell 1',
    'Macrophages 1',
    'Mod.SMC 1',
    'Fibroblast',
    'SMC',
    'Macrophages 2',
    'Pericytes 1',
    'Mod.SMC 2',
    'Pericytes 2',
    'Endothelial cell 2' ,
    'T cell',
    'Neuron',
    'T cell',
    'Plasma cell 1',
    'Plasma cell 2',
     'Macrophages 3'

  )
Idents(athero_comb) <- 'integrated_snn_res.0.4'
names(smc) <- levels(Idents(athero_comb))
athero_comb <- RenameIdents(athero_comb, smc)
smc <- Idents(athero_comb)
athero_comb$new_name <- smc

# --- Whole dataset of coronary arteries from Wirka et al 2019: 4 patients, 8 datasets (1 F, 3 M) ---

Idents(athero_comb) <- 'superclusters'
DimPlot(athero_comb,  pt.size = 1)
DimPlot(athero_comb,  pt.size = 1, label.box = T, label = T, group.by = 'integrated_snn_res.0.4') + NoLegend()
DimPlot(athero_comb,  pt.size = 1, label.box = T, label = T, group.by = 'new_name', repel = T) & NoLegend()

# --- Find markers per cluster ------------------------------------------------

Idents(athero_comb) <- 'integrated_snn_res.0.4'
marker_list_tcell <- list()
for(i in levels(athero_comb$integrated_snn_res.0.4)){
  tmp <- FindMarkers(athero_comb, ident.1 =  i , assay = 'RNA',  max.cells.per.ident = 3000, only.pos = T)
  tmp <- tmp[tmp$p_val_adj < 0.05, ]
  marker_list_tcell[[i]] <- tmp
}
marker_list_tcell <- map(marker_list_tcell, ~ rownames_to_column(.x, var = "GeneName"))
openxlsx::write.xlsx(marker_list_tcell, file = "Wirka_all_celltypes.0.4.xlsx",
                     asTable = TRUE, overwrite = TRUE)

# =============================================================================
# Smooth muscle cell subset and reclusterisation
# =============================================================================

DefaultAssay(athero_comb) <- 'RNA'

athero_comb$keep_cells <-
  ifelse(
    athero_comb$nCount_RNA < 25000 &
      athero_comb$nCount_RNA > 1200 &
      athero_comb$nFeature_RNA < 4000 &
      athero_comb$nFeature_RNA > 1000 &
      athero_comb$percent_mito < 10 &
      athero_comb$percent_hb < 0.1  &
      athero_comb$novelty_score > 0.8 &
      athero_comb$mito_ribo_ratio < 0.75 &
      athero_comb$doublet_finder_sct_15pct == 'Singlet' &
      athero_comb$integrated_snn_res.0.4 %in% c(3,4,2,7,6,8),
    'Keep',
    'Discard'
  )
DimPlot_scCustom(athero_comb, group.by = 'keep_cells')

# --- QC feature plots --------------------------------------------------------

DimPlot_scCustom(athero_comb, group.by = 'integrated_snn_res.0.4', label = T)
FeaturePlot_scCustom(
  athero_comb,
  features = c(
    'percent_mito',
    'percent_ribo',
    'percent.Malat1',
    'mito_ribo_ratio',
    'novelty_score',
    'nFeature_RNA'
  )
)

# --- Subset and clean SMCs ---------------------------------------------------

my_data <- athero_comb[, athero_comb$keep_cells == 'Keep']
my_data <- subset(my_data, sample_id != 'GSM3819857')
FeaturePlot_scCustom(my_data, features = c('PTPRC', 'TYROBP', 'PECAM1'))
DefaultAssay(my_data) <- 'RNA'
my_data <- subset(my_data, subset = PTPRC < 1, slot = 'counts')
my_data <- subset(my_data, subset = TYROBP < 1, slot = 'counts')

DimPlot(my_data,  label = T, label.box = T) + NoLegend()+ plot_annotation(title = 'SMC subset')

# --- Filter genes and re-integrate -------------------------------------------

DefaultAssay(my_data) <- 'RNA'
my_data <- my_data[!grepl('^RP[SL]', rownames(my_data)),]
my_data <- my_data[!grepl('^MT-', rownames(my_data)),]
my_data <- my_data[!grepl('MALAT1', rownames(my_data)),]

# Extract counts
counts <- GetAssayData(object = my_data, slot = "counts", assay = 'RNA')

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
my_data <- CreateSeuratObject(filtered_counts, meta.data = my_data@meta.data)
ifnb.list <- SplitObject(my_data, split.by = "sample_id")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca", k.anchor = 20)
# this command creates an 'integrated' data assay
my_data <- IntegrateData(anchorset = immune.anchors, k.weight = 80)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(my_data) <- "integrated"


# Run the standard workflow for visualization and clustering
my_data <- ScaleData(my_data, verbose = FALSE)
my_data <- RunPCA(my_data, npcs = 30, verbose = FALSE)

# Find significant PCs
stdv <- my_data[["pca"]]@stdev
sum.stdv <- sum(my_data[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                     percent.stdv[2:length(percent.stdv)]) > 0.1),
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
min.pc



my_data <- RunUMAP(my_data, reduction = "pca", dims = 1:co1)
my_data <- FindNeighbors(my_data, reduction = "pca", dims = 1:co1)
my_data <- FindClusters(my_data, resolution = c(0.1, 0.15,  0.2, 0.3, 0.35, 0,4, 0.5, 0.7, 1))

# --- Reclusterisation. Resolution 0.2 works better --------------------------

DimPlot(my_data, group.by = 'integrated_snn_res.0.1', cols = my_pal2, label = T) + DimPlot(my_data, group.by = 'integrated_snn_res.0.2', cols = my_pal2, label = T) +  DimPlot(my_data, group.by = 'superclusters', cols = sample(my_pal2, 5))

# --- QC metrics --------------------------------------------------------------

Idents(my_data) <- 'integrated_snn_res.0.2'
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb","percent.Malat1","mito_ribo_ratio","novelty_score")
FeaturePlot_scCustom(my_data, features = feats)
VlnPlot_scCustom(my_data, features = feats)
DefaultAssay(my_data) <- 'RNA'

# Markers
FeaturePlot_scCustom(my_data, features = c('SERPINF1', 'LUM', 'ACTA2','VCAM1', 'TNFRSF11B', 'BNIP3', 'NDUFA4L2', 'GAPDH'))

# --- DE analysis at two resolutions -----------------------------------------

my_data$res.final.0.1 <- paste('C', my_data$integrated_snn_res.0.1, sep = '')
my_data$res.final.0.2 <- paste('C', my_data$integrated_snn_res.0.2, sep = '')


my_data$res.final.0.1 <- factor(my_data$res.final.0.1)
my_data$res.final.0.2 <- factor(my_data$res.final.0.2)


Idents(my_data) <- 'res.final.0.2'
markers_0.7 <- FindAllMarkers(my_data, assay = 'RNA', logfc.threshold = 0.25,  test.use = 'MAST', min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.7 <- markers_0.7[markers_0.7$p_val_adj < 0.05,]

Idents(my_data) <- 'res.final.0.1'
markers_0.3 <- FindAllMarkers(my_data, assay = 'RNA', logfc.threshold = 0.25,  test.use = 'MAST', min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

# low res
top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 10))
}

p1 <- DotPlot(
  my_data,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'res.final.0.1'
) + ggtitle('Top markers for SMC phenotypes res 0.1 Wirka') + coord_flip() + theme_bw() +
  RotatedAxis() + scale_colour_gradientn(name = 'log2 (count + 1)',colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank())
print(p1)

ggsave(filename = 'human_coronaries_smcs.res0.1.pdf', plot = p1, device = 'pdf', width = 5, height = 10)

# high res
top_markers <- c()
for (i in unique(markers_0.7$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.7, subset = markers_0.7$cluster == i)[order(subset(markers_0.7, subset = markers_0.7$cluster == i)$avg_log2FC, decreasing = T),])$gene, 10))
}

p1 <- DotPlot(
  my_data,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'res.final.0.2'
) + ggtitle('Top markers for SMC phenotypes res 0.2 Wirka') + coord_flip() + theme_bw() +
  RotatedAxis() + scale_colour_gradientn(name = 'log2 (count + 1)',colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank())
print(p1)
ggsave(filename = 'human_coronaries_smcs.res0.2.pdf', plot = p1, device = 'pdf', width = 6, height = 13)

saveRDS(my_data, file = 'human_coronaries_seurat.rds')

# --- Rename clusters to SMC phenotypes ---------------------------------------

Idents(my_data) <- 'res.final.0.2'
smc_phenotypes <- c('Fibroblasts','Pericytes','Transitional', "Contractile" , 'LBH_SMC',"Terminal"  )
names(smc_phenotypes)<- levels(Idents(my_data))
my_data <- RenameIdents(my_data, smc_phenotypes)
smc_phenotypes <- Idents(my_data)
my_data$smc_phenotypes <- smc_phenotypes ### make new clusterisation where all SMCs are in one group
###
Idents(my_data) <- 'smc_phenotypes'

new_names <-
  c( "FB",
    "PC",
    "Trans",
    'SMC_1',
    "SMC_2",
    'Term')
names(new_names) <- levels(Idents(my_data))

my_data <-
  Seurat::RenameIdents(object = my_data, new_names)
new_names <- Idents(my_data)
my_data$smc_phenotypes_short <- new_names


Idents(my_data) <- 'smc_phenotypes'
markers_0.3 <- FindAllMarkers(my_data, assay = 'RNA', logfc.threshold = 0.25,  test.use = 'MAST', min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 6))
}

pal_celltypes <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_celltypes.rds")

p1 <- DotPlot(
  my_data,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'smc_phenotypes'
) + ggtitle('Top markers for SMC phenotypes of human coronary plaques') +
  # coord_flip() +
  # RotatedAxis() +
  scale_colour_gradientn(name = 'log2 (count + 1)',colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme_bw() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank())
print(p1)
ggsave(filename = 'dorplot_human_coronaries_smc_phenotypes.pdf', plot = p1, device = 'pdf', width = 8, height = 3.2, units = 'cm', scale = 3)

pal_celltypes_human <- c( "#80796BFF","#20854EFF","#E18727FF","#0072B5FF",   "#374E55FF" ,"#BC3C29FF" )

p2 <- DimPlot_scCustom(seurat_object = my_data,
                 reduction = 'umap',
                 group.by = 'smc_phenotypes',
                 pt.size = 1,
                 shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes_human, repel = T, label.size = 5
) + NoLegend()+
  # guides(color = guide_legend(nrow = 6)) +
  # theme(legend.position = 'none', legend.text = element_text(size = 9)) +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes_human) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
print(p2)
ggsave(filename = 'dimplot_human_coronaries_smc_phenotypes.pdf', plot = p2, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 3)

my_data$species <- 'hsapiens'
my_data$Origin <- 'Human'
my_data$dataset <-  my_data$orig.ident
my_data$rough_org <- paste('H', my_data$smc_phenotypes, sep = '_')
my_data$fine_org <- paste('H', my_data$smc_phenotypes, sep = '_')

saveRDS(my_data, file = 'human_coronaries_seurat.rds')

# =============================================================================
# Pig RCA: all cell types and SMC subset (from SMCs RCA.Rmd)
# =============================================================================

pigs_2_os <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SC.Analysis.V4.ThirdPass.RNA.RCA.All.seuratSC.Final.rds")

# --- Visualise all cell types ------------------------------------------------

DimPlot(pigs_2_os, group.by = 'Condition') + DimPlot(pigs_2_os, group.by = 'Sample')
DimPlot(pigs_2_os, group.by = 'res.Final', split.by = 'Condition', label = T)

VlnPlot(pigs_2_os, features = c('nCount_RNA', 'nFeature_RNA'))
FeaturePlot(pigs_2_os, features = c('nCount_RNA', 'nFeature_RNA'))
FeaturePlot(subset(pigs_2_os, nCount_RNA < 20000 & nFeature_RNA < 4000), features = c( 'nFeature_RNA')) + ggtitle(label = 'Subset', subtitle =  'nCount_RNA < 20000 & nFeature_RNA < 4000')

describeBy(pigs_2_os$nCount_RNA, group = pigs_2_os$res.Final, mat = T)
describeBy(pigs_2_os$nFeature_RNA, group = pigs_2_os$res.Final, mat = T)

# --- DE for all cell types ----------------------------------------------------

DefaultAssay(pigs_2_os) <- 'RNA'
Idents(pigs_2_os) <- 'res.Final'
markers <- FindAllMarkers(pigs_2_os, assay = 'RNA', logfc.threshold = 0.25,  test.use = 'MAST', min.pct = 0.3,max.cells.per.ident = 800, only.pos = T)
markers <- markers[markers$p_val_adj < 0.05,]
top_markers <- c()
for (i in unique(markers$cluster)){
  top_markers <- c(top_markers, head((subset(markers, subset = markers$cluster == i)[order(subset(markers, subset = markers$cluster == i)$avg_log2FC, decreasing = T),])$gene, 10))
}

p1 <- DotPlot(
  pigs_2_os,
  features = c(unique(top_markers)),
  assay = 'RNA',
  scale = T,
  group.by = 'res.Final'
) + ggtitle('Top markers for res.Final') + coord_flip() +
  RotatedAxis() + scale_colour_gradientn(name = 'log2 (count + 1)',colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  cowplot::theme_cowplot() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank())
ggsave(filename = 'RCA_pig_all_celltypes_DOTPLOT.pdf', plot = p1, width = 12, height = 20)
print(p1)

# --- Annotate superclusters ---------------------------------------------------

new.cluster.ids <- c('Modulated.SMCs' ,'ECs','Macrophages', 'Contractile.SMCs', 'T.cells',  'ECs.2','Fibroblasts','Proliferating.1','Neurons', 'Proliferating.2')
Idents(pigs_2_os) <- 'res.Final'
pigs_2_os$res.Final <- as.factor(pigs_2_os$res.Final )
names(new.cluster.ids) <- levels(pigs_2_os$res.Final)
pigs_2_os <- RenameIdents(pigs_2_os, new.cluster.ids)
pigs_2_os$superclusters <- Idents(pigs_2_os)

colourCount <-
  length(unique(pigs_2_os$superclusters)) # number of levels
getPalette <- colorRampPalette(met.brewer(name = 'Klimt', n =6 , direction = 1))

p1 <- DimPlot_scCustom(
  seurat_object =
    pigs_2_os,
  label.size = 6,
  group.by = 'superclusters',
  aspect_ratio = 1,
  label.box = T,
  label = T,
  colors_use =  colorRampPalette(met.brewer(
    palette_name = 'Cross',
    n = 8 ,
    direction = 1
  ))(colourCount)
) & NoLegend()


p2 <- FeaturePlot_scCustom(pigs_2_os,
                           features = c('ACTA2', 'LUM'),
                           order = F, aspect_ratio = 1)

ggsave(filename = 'pig.coronary.pdf', plot = p1, width = 5.2, height = 5.2)
ggsave(filename = 'pig.coronary_feat.pdf', plot = p2, width = 8, height = 4)

# --- Subset SMCs (ACTA2+/MYH11+ clusters, remove PTPRC+/TYROBP+, filter QC) -

my_data <- subset(pigs_2_os,  res.Final %in% c('C0', 'C3', 'C6') &  nFeature_RNA > 1000  & nCount_RNA > 1200 & nFeature_RNA < 4500 & nCount_RNA < 20000)
my_data <- subset(my_data, PTPRC == 0 & TYROBP == 0)
my_data <- my_data[!grepl('^RP[SL]', rownames(my_data)),]

# Extract counts
counts <- GetAssayData(object = my_data, slot = "counts", assay = 'RNA')

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 5

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
my_data <- CreateSeuratObject(filtered_counts, meta.data = my_data@meta.data)
my_data$Condition_Sample <- paste(my_data$Condition, my_data$Sample, sep = '.')

# --- Normalise, cluster, and embed -------------------------------------------

all.genes <- rownames(my_data)
my_data <-
  NormalizeData(my_data,
                normalization.method = "LogNormalize",
                scale.factor = 10000)
my_data <-
  FindVariableFeatures(my_data,  selection.method = "vst", nfeatures = 2000)
my_data <- ScaleData(my_data,  features = all.genes)


my_data <-
  RunPCA(my_data, npcs = 40, verbose = FALSE) %>%
  RunUMAP(reduction = "pca",
          dims = 1:30,
          verbose = FALSE) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:30,
                verbose = FALSE) %>%
  FindClusters(resolution = c(0.1, 0.15,  0.2, 0.3, 0.35, 0,4, 0.5, 0.7, 1),
               verbose = FALSE)
my_data <- FindClusters(my_data, resolution = c(0.31,0.32,0.33,0.34),
               verbose = FALSE)

# --- SMC cluster overview -----------------------------------------------------

DimPlot(my_data, label = T,  group.by = 'RNA_snn_res.0.35', repel = T, reduction = 'umap', pt.size = 1, label.size = 7, shuffle = T)
FeaturePlot(my_data, features = c('nCount_RNA', 'nFeature_RNA', 'subsets_RB_percent', 'subsets_MT_percent'), pt.size = 1)

DimPlot(my_data, label = T, split.by = 'Sample', group.by = 'RNA_snn_res.0.35', repel = T, reduction = 'umap', pt.size = 1, label.size = 7)
FeaturePlot(my_data, features = c( 'nFeature_RNA'), pt.size = 1,split.by = 'Sample')
DimPlot(my_data,  group.by = 'Sample', split.by = 'RNA_snn_res.0.35', repel = T, reduction = 'umap', pt.size = 1, label.size = 7,shuffle = T )

saveRDS(my_data, 'SMC.RCA.processed.rds')

# =============================================================================
# Preparation of pig coronary arteries (our data)
# =============================================================================

SMC.RCA <- readRDS("SMC.RCA.processed.rds")
SMC.RCA[["percent.Malat1"]] <- PercentageFeatureSet(SMC.RCA, pattern = "^MALAT1")
SMC.RCA$mito_ribo_ratio <- SMC.RCA$subsets_MT_percent/(SMC.RCA$subsets_MT_percent + SMC.RCA$subsets_RB_percent)
SMC.RCA$novelty_score <- log10(SMC.RCA$nFeature_RNA) / log10(SMC.RCA$nCount_RNA)
SMC.RCA$is_NS_good <- ifelse(SMC.RCA$novelty_score < 0.8, 'Bad', 'Good')
SMC.RCA <- SMC.RCA[,SMC.RCA$Condition == 'HFD' & SMC.RCA$smc_0.3 != 'OGN_CLIC5']

dim(SMC.RCA)

# Find significant PCs
stdv <- SMC.RCA[["pca"]]@stdev
sum.stdv <- sum(SMC.RCA[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                     percent.stdv[2:length(percent.stdv)]) > 0.1),
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
min.pc


SMC.RCA <- RunUMAP(SMC.RCA, reduction = "pca", dims = 1:min.pc, verbose = FALSE)
SMC.RCA <- FindNeighbors(SMC.RCA, reduction = "pca", dims = 1:min.pc)
SMC.RCA <- FindClusters(SMC.RCA, resolution = c(0.1,0.2,0.3,0.4,0.5,0.7,1))

# --- QC plots ----------------------------------------------------------------

FeatureScatter(SMC.RCA, feature1 = "novelty_score", feature2 = "mito_ribo_ratio",  plot.cor = FALSE) + theme(legend.position = "bottom") + geom_hline(yintercept = 0.75, lty =2) + geom_vline(xintercept = 0.8, lty = 2)
QC_Plot_UMIvsGene(seurat_object = SMC.RCA, low_cutoff_gene = 250, high_cutoff_gene = 4000, low_cutoff_UMI = 800, high_cutoff_UMI = 20000, group.by = 'is_NS_good')

# --- Visualisation -----------------------------------------------------------

colourCount <-
  length(unique(athero_comb$superclusters)) # number of levels
getPalette <- colorRampPalette(met.brewer(name = 'Klimt', n =6 , direction = 1))

p1 <- DimPlot_scCustom(seurat_object =
                         pig_plaque,
                       group.by = 'superclusters', aspect_ratio = 1,label.box = T,
                       label = T,colors_use =  colorRampPalette(met.brewer(
                         palette_name = 'Cross',
                         n = 8 ,
                         direction = 1
                       ))(colourCount)
) & NoLegend()


p2 <- FeaturePlot_scCustom(athero_comb,
                           features = c('ACTA2', 'LUM'),
                           order = F, aspect_ratio = 1)

ggsave(filename = 'human.carotid.coronary.pdf', plot = p1, width = 6, height = 6)
ggsave(filename = 'human.carotid.coronary_feat.pdf', plot = p2, width = 8, height = 4)

DimPlot(SMC.RCA, group.by = 'RNA_snn_res.0.4')
DimPlot(SMC.RCA, group.by = 'smc_0.3')

FeaturePlot_scCustom(SMC.RCA, features = c('VCAM1', 'TNFRSF11B', 'LUM', 'MYH11'))

# --- Label transfer from human to pig ----------------------------------------

DefaultAssay(human_coronaries_seurat) <- 'RNA'
pancreas.anchors <- FindTransferAnchors(reference = human_coronaries_seurat, query = SMC.RCA,
    dims = 1:30, reference.reduction = "pca", reference.assay = 'RNA', query.assay = 'RNA', project.query = T)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = human_coronaries_seurat$smc_phenotypes,
    dims = 1:30)
SMC.RCA <- AddMetaData(SMC.RCA, metadata = predictions)

DimPlot(SMC.RCA, group.by = 'predicted.id')

# --- Rename pig clusters and find markers ------------------------------------

Idents(SMC.RCA) <- 'smc_0.3'
smc_phenotypes <- c( "Contractile",      "Transitional"    , "Transitional_2" ,  "Terminal_IGFBP2", "Terminal"   ,     "IFN_SMC", "Pericytes",  "Fibroblasts"  ,   "OGN_CLIC5")


names(smc_phenotypes)<- levels(Idents(SMC.RCA))
SMC.RCA <- RenameIdents(SMC.RCA, smc_phenotypes)
smc_phenotypes <- Idents(SMC.RCA)
SMC.RCA$smc_phenotypes <- smc_phenotypes ### make new clusterisation where all SMCs are in one group

Idents(SMC.RCA) <- 'smc_phenotypes'
markers_0.3 <- FindAllMarkers(SMC.RCA, assay = 'RNA', logfc.threshold = 0.25,  test.use = 'MAST', min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 6))
}

pal_celltypes <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_celltypes.rds")

p1 <- DotPlot(
  SMC.RCA,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'smc_phenotypes'
) + ggtitle('Top markers for SMC phenotypes of pig coronary plaques') +
  # coord_flip() +
  # RotatedAxis() +
  scale_colour_gradientn(name = 'log2 (count + 1)',colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme_bw() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank())
print(p1)
ggsave(filename = 'dorplot_pig_coronaries_smc_phenotypes.pdf', plot = p1, device = 'pdf', width = 8, height = 3.2, units = 'cm', scale = 3)

pal_celltypes_pig <- c( "#79AF97FF","#E18727FF", "#0072B5FF",  "#BC3C29FF", "#20854EFF", "#7876B1FF" ,"#FFDC91FF", "#EE4C97FF", "#80796BFF")

p2 <- DimPlot_scCustom(seurat_object = SMC.RCA,
                 reduction = 'umap',
                 group.by = 'smc_phenotypes',
                 pt.size = 1,
                 shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes_pig, repel = T, label.size = 5
) + NoLegend()+
  # guides(color = guide_legend(nrow = 6)) +
  # theme(legend.position = 'none', legend.text = element_text(size = 9)) +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes_pig) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
print(p2)
ggsave(filename = 'dimplot_pig_coronaries_smc_phenotypes.pdf', plot = p2, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 3)

SMC.RCA$species <- 'sscrofa'
SMC.RCA$Origin <- 'Pig'
SMC.RCA$dataset <-  SMC.RCA$Sample
# SMC.RCA$rough_org <- paste('H', SMC.RCA$smc_phenotypes, sep = '_')
# SMC.RCA$fine_org <- paste('H', SMC.RCA$smc_phenotypes, sep = '_')

saveRDS(SMC.RCA, file = 'pig_coronaries_seurat.rds')
