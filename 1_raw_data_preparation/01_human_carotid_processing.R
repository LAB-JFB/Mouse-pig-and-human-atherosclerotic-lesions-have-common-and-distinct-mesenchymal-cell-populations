# =============================================================================
# 01 - Processing of Human Carotid Dataset
# =============================================================================
# This script processes the human carotid plaque scRNA-seq dataset:
#   1. QC filtering and doublet removal
#   2. SMC subset extraction
#   3. Integration across samples (RPCA)
#   4. Clustering and marker identification
#   5. SMC phenotype annotation
# =============================================================================

# --- Setup -------------------------------------------------------------------

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/human")
set.seed(1984)

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(Matrix)
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
library(DESeq2)
library(gdata)
library(SeuratDisk)
library(SingleR)
library(celldex)
library(scCATCH)
library(DoubletFinder)

# =============================================================================
# Raw data loading and initial processing (from scrnaseq_2_seurat-processing.R)
# =============================================================================

# --- Load raw 10X counts -----------------------------------------------------

sc_paths <- list(
  alsaigh = list.dirs("~/THOR/scrna-seq/raw/2_counts/Alsaigh_2020_GSE159677", recursive = FALSE),
  pan = list.dirs("~/THOR/scrna-seq/raw/2_counts/Pan_2020_GSE155512", recursive = FALSE),
  wirka = list.dirs("~/THOR/scrna-seq/raw/2_counts/Wirka_2019_GSE131778", recursive = FALSE)
)
sc_samples <- sc_paths %>% map(~ gsub(".+/GSM", "GSM", .x))

sample_info <- read.delim("sample_info_basic.txt")
rownames(sample_info) <- sample_info$sample_id

# --- Create Seurat objects per dataset ----------------------------------------

sc_data <- list()
for (i in 1:length(sc_paths)) {
  my_ds <- sc_paths[[i]]
  sc_data[[i]] <- lapply(sc_paths[[i]], function(my_path) {
    my_sample <- gsub(".+/GSM", "GSM", my_path)
    my_data <- Read10X(my_path)
    my_seurat <- CreateSeuratObject(counts = my_data, min.cells = 3, min.features = 200, project = names(sc_paths)[i])
    my_seurat$sample_id <- my_sample
    my_seurat <- RenameCells(my_seurat, add.cell.id = my_sample)
    return(my_seurat)
  })
  names(sc_data[[i]]) <- sc_samples[[i]]
}
rm(i, my_ds)
names(sc_data) <- names(sc_paths)

# --- Merge datasets -----------------------------------------------------------

seurat_alsaigh <- merge(x = sc_data$alsaigh[[1]], y = sc_data$alsaigh[2:length(sc_data$alsaigh)])
seurat_pan <- merge(x = sc_data$pan[[1]], y = sc_data$pan[2:length(sc_data$pan)])
seurat_wirka <- merge(x = sc_data$wirka[[1]], y = sc_data$wirka[2:length(sc_data$wirka)])

save(seurat_alsaigh, seurat_pan, seurat_wirka, sample_info,
  file = "sc_athero_1.RData"
)

# --- Add QC metrics -----------------------------------------------------------

load("sc_athero_1.RData")

seurat_alsaigh <- PercentageFeatureSet(seurat_alsaigh, "^MT-", col.name = "percent_mito")
seurat_alsaigh <- PercentageFeatureSet(seurat_alsaigh, "^RP[SL]", col.name = "percent_ribo")
seurat_alsaigh <- PercentageFeatureSet(seurat_alsaigh, "^HB[^(P)]", col.name = "percent_hb")

seurat_pan <- PercentageFeatureSet(seurat_pan, "^MT-", col.name = "percent_mito")
seurat_pan <- PercentageFeatureSet(seurat_pan, "^RP[SL]", col.name = "percent_ribo")
seurat_pan <- PercentageFeatureSet(seurat_pan, "^HB[^(P)]", col.name = "percent_hb")

seurat_wirka <- PercentageFeatureSet(seurat_wirka, "^MT-", col.name = "percent_mito")
seurat_wirka <- PercentageFeatureSet(seurat_wirka, "^RP[SL]", col.name = "percent_ribo")
seurat_wirka <- PercentageFeatureSet(seurat_wirka, "^HB[^(P)]", col.name = "percent_hb")

# --- QC visualisation (raw) ---------------------------------------------------

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(seurat_alsaigh, group.by = "sample_id", features = feats, pt.size = 0.1, ncol = 3, combine = TRUE) + NoLegend()
p1 <- FeatureScatter(seurat_alsaigh, feature1 = "nCount_RNA", feature2 = "percent_mito")
p2 <- FeatureScatter(seurat_alsaigh, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2

VlnPlot(seurat_pan, group.by = "sample_id", features = feats, pt.size = 0.1, ncol = 5, combine = TRUE) + NoLegend()
p1 <- FeatureScatter(seurat_pan, feature1 = "nCount_RNA", feature2 = "percent_mito")
p2 <- FeatureScatter(seurat_pan, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 <- FeatureScatter(seurat_pan, feature1 = "percent_ribo", feature2 = "percent_mito")
p1 + p2 + p3

VlnPlot(seurat_wirka, group.by = "sample_id", features = feats, pt.size = 0.1, ncol = 5, combine = TRUE) + NoLegend()
p1 <- FeatureScatter(seurat_wirka, feature1 = "nCount_RNA", feature2 = "percent_mito")
p2 <- FeatureScatter(seurat_wirka, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 <- FeatureScatter(seurat_wirka, feature1 = "percent_ribo", feature2 = "percent_mito")
p1 + p2 + p3

# --- QC filtering -------------------------------------------------------------

seurat_alsaigh <- subset(seurat_alsaigh, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA < 4000 & percent_mito < 10 & percent_hb < 1)
seurat_pan <- subset(seurat_pan, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA < 4000 & percent_mito < 10 & percent_hb < 1)
seurat_wirka <- subset(seurat_wirka, subset = nCount_RNA > 300 & nCount_RNA < 30000 & nFeature_RNA < 4000 & percent_mito < 10 & percent_hb < 1)

athero <- list(alsaigh = seurat_alsaigh, pan = seurat_pan, wirka = seurat_wirka)
save(athero, sample_info, file = "sc_athero_2.RData")

# --- SCTransform, PCA, UMAP per dataset ---------------------------------------

load("sc_athero_2.RData")

athero <- lapply(athero, function(x) {
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent_mito", verbose = TRUE)
  x <- RunPCA(x)
  x <- RunUMAP(x, dims = 1:40)
})
lapply(athero, dim)

# --- Experiments with Alsaigh: automatic annotation ---------------------------

seurat_alsaigh <- athero$alsaigh
DefaultAssay(seurat_alsaigh) <- "RNA"
seurat_alsaigh <- FindNeighbors(seurat_alsaigh, dims = 1:40)
seurat_alsaigh <- FindClusters(seurat_alsaigh, resolution = 0.3)
table(seurat_alsaigh$seurat_clusters)

# SingleR
hpca_se <- readRDS("hpca_se.rds")
alsaigh_singler <- SingleR(
  test = seurat_alsaigh@assays$RNA@data,
  assay.type.test = 1,
  ref = hpca_se, labels = hpca_se$label.main
)
seurat_alsaigh$singler <- alsaigh_singler$labels

# scCATCH
alsaigh_scc_dat <- rev_gene(seurat_alsaigh@assays$RNA@data,
  data_type = "data",
  species = "Human", geneinfo = geneinfo
)
obj <- createscCATCH(data = alsaigh_scc_dat, cluster = formatC(seurat_alsaigh$seurat_clusters, flag = 0, width = 2))
obj <- findmarkergene(
  object = obj, species = "Human", marker = cellmatch,
  tissue = c("Blood vessel", "Heart", "Myocardium", "Serum")
)
obj <- findcelltype(object = obj) ## does not work :-(

# Seurat / Tabula sapiens
ts_vasculature <- LoadH5Seurat("~/THOR/scrna-seq/reference/tabula-sapiens/TS_Vasculature.h5seurat",
  assays = "RNA"
)
ts_vasculature <- FindVariableFeatures(ts_vasculature, selection.method = "vst")
ts_vasculature <- ScaleData(ts_vasculature)
ts_vasculature <- RunPCA(ts_vasculature, npcs = 40, verbose = TRUE)
DimPlot(ts_vasculature, reduction = "umap", group.by = "Annotation", label = TRUE, repel = TRUE)
alsaigh_anchors <- FindTransferAnchors(
  reference = ts_vasculature, query = seurat_alsaigh,
  normalization.method = "LogNormalize", dims = 1:40
)
predictions <- TransferData(anchorset = alsaigh_anchors, refdata = ts_vasculature$Annotation, dims = 1:40)
seurat_alsaigh <- AddMetaData(seurat_alsaigh, metadata = predictions)
table(seurat_alsaigh$predicted.id)

reduc_method <- "umap"
p1 <- DimPlot(seurat_alsaigh, reduction = reduc_method, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DimPlot(seurat_alsaigh, reduction = reduc_method, group.by = "singler", label = TRUE, repel = TRUE)
p3 <- DimPlot(seurat_alsaigh, reduction = reduc_method, group.by = "predicted.id", label = TRUE, repel = TRUE)
p1 + p2 + p3 # bad predictions for Seurat/TS

# --- DoubletFinder on Alsaigh (exploratory) -----------------------------------

ElbowPlot(athero[[1]], ndims = 40)
sweep_res_alsaigh <- paramSweep_v3(athero$alsaigh, PCs = 1:40, sct = TRUE)
sweep_stats_alsaigh <- summarizeSweep(sweep_res_alsaigh, GT = FALSE)
bcmvn_alsaigh <- find.pK(sweep_stats_alsaigh)
par(mar = c(4, 4, 3, 2))
plot(as.numeric(as.character(bcmvn_alsaigh$pK)), bcmvn_alsaigh$BCmetric,
  type = "b", col = 4, lty = 1, xlab = "pK", ylab = "BCmetric"
)
pK <- as.numeric(as.character(bcmvn_alsaigh$pK))[which.max(bcmvn_alsaigh$BCmetric)]
nExp_poi <- round(0.075 * ncol(athero$alsaigh))
seurat_alsaigh <- doubletFinder_v3(athero$alsaigh, PCs = 1:40, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

homotypic.prop <- modelHomotypic(seurat_alsaigh$singler)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
seurat_alsaigh <- doubletFinder_v3(seurat_alsaigh,
  PCs = 1:40, pN = 0.25, pK = pK, nExp = 10000,
  reuse.pANN = "pANN_0.25_0.14_3407", sct = TRUE
)
par(mfrow = c(1, 2))
FeatureScatter(seurat_alsaigh, feature1 = "pANN_0.25_0.14_3407", feature2 = "nCount_RNA") + NoLegend() +
  FeatureScatter(seurat_alsaigh, feature1 = "pANN_0.25_0.14_3407", feature2 = "nFeature_RNA") + NoLegend()
dev.off()
FeatureScatter(seurat_alsaigh, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()

reduc_method <- "umap"
p1 <- DimPlot(seurat_alsaigh, reduction = reduc_method, group.by = "singler", label = TRUE, repel = TRUE) + NoLegend()
p2 <- DimPlot(seurat_alsaigh, reduction = reduc_method, group.by = "DF.classifications_0.25_0.14_3407")
p3 <- FeaturePlot(seurat_alsaigh, reduction = reduc_method, features = c("CD3D"))
p1 + p2 + p3

# --- DoubletFinder + SingleR on all datasets ----------------------------------

hpca_se <- readRDS("hpca_se.rds")

athero_info <- athero_2 <- dbf_param <- list()
for (i in 1:length(athero)) {
  tmp_obj <- athero[[i]]
  DefaultAssay(tmp_obj) <- "RNA"
  tmp_obj <- FindNeighbors(tmp_obj, dims = 1:40)
  cat("! Automatic annotation using SingleR \n----------\n")
  tmp_singler <- SingleR(
    test = tmp_obj@assays$RNA@data,
    assay.type.test = 1,
    ref = hpca_se, labels = hpca_se$label.main
  )
  tmp_obj$singler <- tmp_singler$labels
  cat("! Identification of doublets using DoubletFinder \n-----------\n")
  tmp_sweep_res <- paramSweep_v3(tmp_obj, PCs = 1:40, sct = TRUE)
  tmp_sweep_stats <- summarizeSweep(tmp_sweep_res, GT = FALSE)
  tmp_bcmvn <- find.pK(tmp_sweep_stats)
  pK <- as.numeric(as.character(tmp_bcmvn$pK))[which.max(tmp_bcmvn$BCmetric)]
  nExp_poi <- round(0.075 * ncol(tmp_obj))
  homotypic.prop <- modelHomotypic(tmp_obj$singler)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  tmp_obj <- doubletFinder_v3(tmp_obj,
    PCs = 1:40, pN = 0.25, pK = pK, nExp = nExp_poi,
    reuse.pANN = FALSE, sct = TRUE
  )
  tmp_obj <- doubletFinder_v3(tmp_obj,
    PCs = 1:40, pN = 0.25, pK = pK, nExp = nExp_poi.adj,
    reuse.pANN = names(tmp_obj@meta.data)[11], sct = TRUE
  )
  DefaultAssay(tmp_obj) <- "SCT"
  tmp_obj <- FindNeighbors(tmp_obj, dims = 1:40)
  tmp_obj <- FindClusters(tmp_obj, resolution = 0.3)
  cat("! Saving results \n-----------\n")
  athero_2[[i]] <- tmp_obj
  athero_info[[i]] <- data.frame(
    cell.id = Cells(tmp_obj),
    prior = tmp_obj@meta.data[, 10:ncol(tmp_obj@meta.data)]
  )
  dbf_param[[i]] <- list(
    pK = pK, nExp_poi = nExp_poi, nExp_poi.homo.adj = nExp_poi.adj,
    bcmvn = tmp_bcmvn, sweep_stats = tmp_sweep_stats
  )
  cat("! Done! \n-----------\n\n")
}
names(athero_2) <- names(athero_info) <- names(dbf_param) <- names(athero)
save(athero_info, dblf_param, file = "athero_prior_info.RData")
rm(grep("tmp", ls(), value = TRUE))

# --- Integrate all three datasets ---------------------------------------------

features <- SelectIntegrationFeatures(object.list = athero)
athero <- PrepSCTIntegration(object.list = athero, anchor.features = features)
athero_anchors <- FindIntegrationAnchors(
  object.list = athero, normalization.method = "SCT",
  anchor.features = features
)
athero_comb <- IntegrateData(anchorset = athero_anchors, normalization.method = "SCT")
athero_comb <- RunPCA(athero_comb, verbose = TRUE)
athero_comb <- RunUMAP(athero_comb, reduction = "pca", dims = 1:40)
athero_comb <- RunTSNE(athero_comb, reduction = "pca", dims = 1:40)
athero_comb$dataset_author <- str_to_title(athero_comb$orig.ident)

# --- Annotate integrated object with SingleR ----------------------------------

hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
athero_singler <- SingleR(
  test = athero_comb@assays$RNA@data,
  assay.type.test = 1,
  ref = hpca.se, labels = hpca.se$label.main
)
athero_comb$singler <- athero_singler$labels

# --- Add PopV / Tabula Sapiens annotation -------------------------------------

ts_anno <- read.csv("annotated_query.csv", row.names = 1)
athero_comb <- AddMetaData(athero_comb, metadata = ts_anno[, 21:28])
table(athero_comb$consensus_prediction)
athero_comb$consensus_prediction[athero_comb$consensus_prediction == "antibody secreting cell"] <- "plasma cell"
athero_comb$consensus_prediction[athero_comb$consensus_prediction == "lymphocyte of b lineage"] <- "b cell"
athero_comb$consensus_prediction[athero_comb$consensus_prediction == "gamma-delta t cell"] <- "t cell"

# --- Add sample metadata ------------------------------------------------------

sc_anno <- sample_info[gsub("_.*", "", Cells(athero_comb)), ]
rownames(sc_anno) <- Cells(athero_comb)
athero_comb <- AddMetaData(athero_comb, metadata = sc_anno)

save(athero_comb, sample_info, file = "sc_athero_comb_new.rda")

# =============================================================================
# Carotid plaque SMC processing (original 01 pipeline)
# =============================================================================

# --- Load data ---------------------------------------------------------------
# athero_comb is now already in memory from the processing above

athero_comb$sample_id %>% table()
athero_comb$orig.ident %>% table()
athero_comb$group %>% table()

# --- QC metrics --------------------------------------------------------------

athero_comb[["percent.mt"]] <- PercentageFeatureSet(athero_comb, pattern = "^MT-")
athero_comb[["percent.ribo"]] <- PercentageFeatureSet(athero_comb, pattern = "^RP[SL]")
athero_comb[["percent.hb"]] <- PercentageFeatureSet(athero_comb, pattern = "^HB[^(P)]")
athero_comb[["percent.Malat1"]] <- PercentageFeatureSet(athero_comb, pattern = "^MALAT1")
athero_comb$mito_ribo_ratio <- athero_comb$percent.mt / (athero_comb$percent.mt + athero_comb$percent.ribo)
athero_comb$novelty_score <- log10(athero_comb$nFeature_RNA) / log10(athero_comb$nCount_RNA)
athero_comb$is_NS_good <- ifelse(athero_comb$novelty_score < 0.8, "Bad", "Good")

# Subset to carotid plaque samples only
athero_comb <- athero_comb[, athero_comb$group == "Carotid Plaque"]

# --- QC visualisation --------------------------------------------------------

FeatureScatter(athero_comb,
               feature1 = "novelty_score",
               feature2 = "mito_ribo_ratio",
               plot.cor = FALSE) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0.75, lty = 2) +
  geom_vline(xintercept = 0.8, lty = 2)

QC_Plot_UMIvsGene(seurat_object = athero_comb,
                  low_cutoff_gene = 250, high_cutoff_gene = 4000,
                  low_cutoff_UMI = 800, high_cutoff_UMI = 20000,
                  group.by = "is_NS_good")

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt",
           "percent.ribo", "percent.hb", "percent.Malat1")

dim(athero_comb)

VlnPlot(athero_comb, features = feats, ncol = 3)

plot1 <- FeatureScatter(athero_comb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(athero_comb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

summary(athero_comb$nCount_RNA)
summary(athero_comb$nFeature_RNA)
summary(athero_comb$percent.mt)
summary(athero_comb$percent.ribo)

DimPlot_scCustom(athero_comb)

# --- Doublet identification (scDblFinder) ------------------------------------

require(DoubletFinder)
require(scDblFinder)

DefaultAssay(athero_comb) <- "integrated"
athero_comb <- FindClusters(athero_comb, resolution = 1)
dr <- 0.075 # Assuming standard doublet rate (7.5% by default)

DimPlot_scCustom(athero_comb, group.by = "integrated_snn_res.1")

athero_comb$seurat_clusters <- athero_comb$integrated_snn_res.1
seu_split <- SplitObject(athero_comb, split.by = "sample_id")

seu_split <- map(seu_split, function(x) {
  sce <- as.SingleCellExperiment(x, assay = "RNA")
  sce <- scDblFinder(sce, sce$seurat_clusters, dbr.sd = 1)
  scdbl_cols <- grep("scDblFinder", colnames(colData(sce)), value = TRUE)
  x@meta.data <- cbind(x@meta.data, colData(sce)[, scdbl_cols])
  return(x)
})

athero_comb_2 <- Merge_Seurat_List(seu_split)
v_modified <- substring(colnames(athero_comb_2), 2)
athero_comb_2$cellID <- v_modified
colnames(athero_comb_2) <- v_modified

athero_comb$scDblFinder.class <- athero_comb_2$scDblFinder.class
rm(athero_comb_2)

DimPlot_scCustom(athero_comb, group.by = "integrated_snn_res.1") +
  DimPlot_scCustom(athero_comb, group.by = "scDblFinder.class")

FeaturePlot_scCustom(athero_comb,
                     features = c("ACTA2", "MYH11", "POSTN", "PTPRC", "TYROBP", "PECAM1"),
                     num_columns = 4)

saveRDS(athero_comb, file = "athero_comb_carotids_withDoublets.rds")

# --- SMC subset extraction ---------------------------------------------------

DimPlot_scCustom(athero_comb, group.by = "integrated_snn_res.0.3", label = TRUE)

FeaturePlot_scCustom(athero_comb,
                     features = c("percent.mt", "percent.ribo", "percent.Malat1",
                                  "mito_ribo_ratio", "novelty_score", "nFeature_RNA"))

# Define cells to keep based on QC thresholds and cluster identity
athero_comb$keep_cells <- ifelse(
  athero_comb$nCount_RNA < 25000 &
    athero_comb$nCount_RNA > 1200 &
    athero_comb$nFeature_RNA < 4000 &
    athero_comb$nFeature_RNA > 1000 &
    athero_comb$percent.mt < 10 &
    athero_comb$percent.hb < 0.1 &
    athero_comb$novelty_score > 0.8 &
    athero_comb$mito_ribo_ratio < 0.75 &
    athero_comb$scDblFinder.class == "singlet" &
    athero_comb$integrated_snn_res.0.3 %in% c(5, 7, 4, 14, 8),
  "Keep", "Discard"
)

DimPlot_scCustom(athero_comb, group.by = "keep_cells")

athero_comb_smc <- athero_comb[, athero_comb$keep_cells == "Keep"]
FeaturePlot_scCustom(athero_comb_smc, features = c("PTPRC", "TYROBP", "PECAM1"))

# Remove contaminating immune / endothelial cells
DefaultAssay(athero_comb_smc) <- "RNA"
athero_comb_smc <- subset(athero_comb_smc, subset = PTPRC < 1, slot = "counts")
athero_comb_smc <- subset(athero_comb_smc, subset = TYROBP < 1, slot = "counts")

# --- Feature cleaning --------------------------------------------------------

# Remove ribosomal, mitochondrial, and MALAT1 genes
athero_comb_smc <- athero_comb_smc[!grepl("^RP[SL]", rownames(athero_comb_smc)), ]
athero_comb_smc <- athero_comb_smc[!grepl("^MT-", rownames(athero_comb_smc)), ]
athero_comb_smc <- athero_comb_smc[!grepl("MALAT1", rownames(athero_comb_smc)), ]

# Filter genes expressed in fewer than 5 cells
counts <- GetAssayData(object = athero_comb_smc, slot = "counts", assay = "RNA")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 5
filtered_counts <- counts[keep_genes, ]

athero_comb_smc <- CreateSeuratObject(filtered_counts, meta.data = athero_comb_smc@meta.data)

# --- Integration (RPCA) ------------------------------------------------------

athero_comb_smc$orig.ident <- athero_comb_smc$sample_id
ifnb.list <- SplitObject(athero_comb_smc, split.by = "orig.ident")

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
athero_comb_smc <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(athero_comb_smc) <- "integrated"

# --- Dimensionality reduction and clustering ---------------------------------

athero_comb_smc <- ScaleData(athero_comb_smc, verbose = FALSE)
athero_comb_smc <- RunPCA(athero_comb_smc, npcs = 40, verbose = FALSE)

# Determine significant PCs
stdv <- athero_comb_smc[["pca"]]@stdev
sum.stdv <- sum(athero_comb_smc[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 80 & percent.stdv < 10)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                     percent.stdv[2:length(percent.stdv)]) > 0.1),
            decreasing = TRUE)[1] + 1

print(co1)
print(co2)
ElbowPlot(athero_comb_smc, ndims = 40)

# UMAP, neighbours, and multi-resolution clustering
athero_comb_smc <-
  RunUMAP(object = athero_comb_smc, reduction = "pca",
          dims = 1:co1, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca",
                dims = 1:co1, verbose = FALSE) %>%
  FindClusters(resolution = c(0.1, 0.15, 0.2, 0.3, 0.35, 0.4, 0.5, 0.7, 1),
               verbose = FALSE)

athero_comb_smc <- FindClusters(object = athero_comb_smc,
                                resolution = c(0.12, 0.13, 0.14, 0.15,
                                               0.31, 0.32, 0.33, 0.34))

# --- Clustering visualisation ------------------------------------------------

DimPlot_scCustom(athero_comb_smc, label = TRUE,
                 group.by = "integrated_snn_res.0.12", repel = TRUE,
                 reduction = "umap", pt.size = 1, label.size = 7, shuffle = TRUE) +
  DimPlot_scCustom(athero_comb_smc, label = TRUE,
                   group.by = "integrated_snn_res.0.3", repel = TRUE,
                   reduction = "umap", pt.size = 1, label.size = 7, shuffle = TRUE) +
  DimPlot_scCustom(athero_comb_smc, label = FALSE,
                   group.by = "orig.ident", repel = TRUE,
                   reduction = "umap", pt.size = 1, label.size = 7, shuffle = TRUE)

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo",
           "percent.hb", "percent.Malat1", "mito_ribo_ratio", "novelty_score")
FeaturePlot_scCustom(athero_comb_smc, features = feats)
VlnPlot_scCustom(athero_comb_smc, features = feats)

DefaultAssay(athero_comb_smc) <- "RNA"

# --- Marker identification ---------------------------------------------------

FeaturePlot_scCustom(athero_comb_smc,
                     features = c("SERPINF1", "LUM", "ACTA2", "VCAM1",
                                  "TNFRSF11B", "BNIP3", "NDUFA4L2", "GAPDH"))

athero_comb_smc$res.final.0.12 <- paste0("C", athero_comb_smc$integrated_snn_res.0.12)
athero_comb_smc$res.final.0.3 <- paste0("C", athero_comb_smc$integrated_snn_res.0.3)
athero_comb_smc$res.final.0.12 <- factor(athero_comb_smc$res.final.0.12)
athero_comb_smc$res.final.0.3 <- factor(athero_comb_smc$res.final.0.3)

# Markers at resolution 0.3
Idents(athero_comb_smc) <- "res.final.0.3"
markers_0.7 <- FindAllMarkers(athero_comb_smc, assay = "RNA",
                               logfc.threshold = 0.25, test.use = "MAST",
                               min.pct = 0.3, only.pos = TRUE,
                               max.cells.per.ident = 1000)
markers_0.7 <- markers_0.7[markers_0.7$p_val_adj < 0.05, ]

# Markers at resolution 0.12
Idents(athero_comb_smc) <- "res.final.0.12"
markers_0.3 <- FindAllMarkers(athero_comb_smc, assay = "RNA",
                               logfc.threshold = 0.25, test.use = "MAST",
                               min.pct = 0.3, only.pos = TRUE,
                               max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05, ]

# --- Dot plots (low resolution) ---------------------------------------------

top_markers <- c()
for (i in unique(markers_0.3$cluster)) {
  cluster_markers <- subset(markers_0.3, subset = markers_0.3$cluster == i)
  cluster_markers <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = TRUE), ]
  top_markers <- c(top_markers, head(cluster_markers$gene, 10))
}

p1 <- DotPlot(athero_comb_smc,
              features = unique(top_markers),
              assay = "RNA", scale = TRUE,
              group.by = "res.final.0.12") +
  ggtitle("Top markers for SMC phenotypes res 0.12") +
  coord_flip() +
  theme_bw() +
  RotatedAxis() +
  scale_colour_gradientn(name = "log2 (count + 1)",
                         colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank()) +
  ylab("")
print(p1)

ggsave(filename = "human_carotid_smcs.res0.12.pdf",
       plot = p1, device = "pdf", width = 5, height = 10)

# --- Dot plots (high resolution) --------------------------------------------

top_markers <- c()
for (i in unique(markers_0.7$cluster)) {
  cluster_markers <- subset(markers_0.7, subset = markers_0.7$cluster == i)
  cluster_markers <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = TRUE), ]
  top_markers <- c(top_markers, head(cluster_markers$gene, 10))
}

p1 <- DotPlot(athero_comb_smc,
              features = unique(top_markers),
              assay = "RNA", scale = TRUE,
              group.by = "res.final.0.3") +
  ggtitle("Top markers for SMC phenotypes res 0.3") +
  coord_flip() +
  theme_bw() +
  RotatedAxis() +
  scale_colour_gradientn(name = "log2 (count + 1)",
                         colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank()) +
  ylab("")
print(p1)

ggsave(filename = "human_carotid_smcs.res0.3.pdf",
       plot = p1, device = "pdf", width = 6, height = 13)

saveRDS(athero_comb_smc, file = "athero_comb_smc_seurat.rds")

# --- SMC phenotype annotation ------------------------------------------------

Idents(athero_comb_smc) <- "res.final.0.3"

smc_phenotypes <- c("Contractile", "Terminal", "Transitional",
                    "DLX5_SMC", "Pericytes", "CRTAC1_SMC")
names(smc_phenotypes) <- levels(Idents(athero_comb_smc))
athero_comb_smc <- RenameIdents(athero_comb_smc, smc_phenotypes)
athero_comb_smc$smc_phenotypes <- Idents(athero_comb_smc)

# Short names
Idents(athero_comb_smc) <- "smc_phenotypes"
new_names <- c("SMC1", "Term", "Trans", "FCh1", "PC", "FCh2")
names(new_names) <- levels(Idents(athero_comb_smc))
athero_comb_smc <- Seurat::RenameIdents(object = athero_comb_smc, new_names)
athero_comb_smc$smc_phenotypes_short <- Idents(athero_comb_smc)

# --- Markers by phenotype ----------------------------------------------------

Idents(athero_comb_smc) <- "smc_phenotypes"
markers_0.3 <- FindAllMarkers(athero_comb_smc, assay = "RNA",
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

pal_celltypes <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_celltypes.rds")
pal_celltypes_human <- c("#0072B5FF", "#BC3C29FF", "#E18727FF",
                         "#80796BFF", "#20854EFF", "#374E55FF")

p1 <- DotPlot(athero_comb_smc,
              features = unique(top_markers),
              assay = "RNA", scale = TRUE,
              group.by = "smc_phenotypes") +
  ggtitle("Top markers for SMC phenotypes human carotids") +
  scale_colour_gradientn(name = "log2 (count + 1)",
                         colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.ticks = element_blank()) +
  ylab("")
print(p1)

ggsave(filename = "dorplot_human_carotids_smc_phenotypes.pdf",
       plot = p1, device = "pdf", width = 8, height = 3.2,
       units = "cm", scale = 3)

p2 <- DimPlot_scCustom(seurat_object = athero_comb_smc,
                        reduction = "umap", group.by = "smc_phenotypes",
                        pt.size = 1, shuffle = TRUE, label = TRUE,
                        label.box = TRUE, colors_use = pal_celltypes_human,
                        repel = TRUE, label.size = 5) +
  NoLegend() +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes_human) +
  theme(axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
print(p2)

ggsave(filename = "dimplot_human_carotids_smc_phenotypes.pdf",
       plot = p2, device = "pdf", width = 4, height = 4,
       units = "cm", scale = 3)

# --- Add metadata and save ---------------------------------------------------

athero_comb_smc$species <- "hsapiens"
athero_comb_smc$Origin <- "Human"
athero_comb_smc$dataset <- athero_comb_smc$orig.ident
athero_comb_smc$rough_org <- paste0("H_", athero_comb_smc$smc_phenotypes)
athero_comb_smc$fine_org <- paste0("H_", athero_comb_smc$smc_phenotypes)

saveRDS(athero_comb_smc, file = "athero_comb_smc_seurat.rds")