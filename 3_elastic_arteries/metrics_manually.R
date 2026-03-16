# ============================================================================
# Integration quality metrics for elastic arteries (30 method-homology combos)
# Metrics: kBET, PC Regression, LISI (iLISI + cLISI), ARI, NMI, cASW, bASW
# Adapted from: github.com/MillerLab-CPHG/Human_athero_scRNA_meta
# ============================================================================

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration")

# --- Libraries ---
library(kBET)
library(FNN)
library(lisi)
library(cluster)
library(Seurat)
library(scRNAutils)
library(stats)
library(DescTools)
library(data.table)
library(tidyverse)
library(ggsci)
library(bluster)
library(aricode)
library(SingleCellExperiment)

# --- Load all 30 integrated datasets ---
all_datasets <- readRDS("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/all_datasets_ALL.rds")

# --- Map each dataset to its embedding reduction name ---
embs <- c('pca','pca','pca','pca','pca','pca',
          'X_scANVI', 'X_scVI' ,'X_scANVI', 'X_scVI' ,'X_scANVI', 'X_scVI' ,
          'inmf' , 'inmf' ,'inmf' , 'mnn', 'iNMF','mnn', 'iNMF','mnn', 'iNMF',
          'X_pca_harmony', 'X_scanorama' ,'X_pca_harmony', 'X_scanorama' ,'X_pca_harmony', 'X_scanorama',
          'pca','pca','pca')
names(embs) <- names(all_datasets)

# Annotate each dataset with method, homology, and embedding type
for(i in names(all_datasets)[1:27]){
  all_datasets[[i]]$method <- str_split(i, pattern = '_')[[1]][length(str_split(i, pattern = '_')[[1]])]
  all_datasets[[i]]$homology <- paste0(str_split(i, pattern = '_')[[1]][-length(str_split(i, pattern = '_')[[1]])], collapse = '_')
  all_datasets[[i]]$embedding <- as.character(embs[i])
}

for(i in names(all_datasets)[28:30]){
  all_datasets[[i]]$method <- str_split(i, pattern = '_')[[1]][1]
  all_datasets[[i]]$homology <- paste0(str_split(i, pattern = '_')[[1]][-1], collapse = '_')
  all_datasets[[i]]$embedding <- as.character(embs[i])
}

# --- Short display names for 30 combinations ---
new_names <-
  c(
    "HE_CCA",
    "HE_RPCA",
    "HHC_CCA",
    "HHC_RPCA",
    "O2O_CCA",
    "O2O_RPCA",
    "HE_scANVI",
    "HE_scVI",
    "HHC_scANVI" ,
    "HHC_scVI",
    "O2O_scANVI"  ,
    "O2O_scVI" ,

    "HE_LIGER_UINMF",
    "HHC_LIGER_UINMF",
    "O2O_LIGER_UINMF",

    "HE_fastMNN",
    "HE_LIGER",
    "HHC_fastMNN",
    "HHC_LIGER",
    "O2O_fastMNN",
    "O2O_LIGER",

    "HE_Harmony",
    "HE_Scanorama",
    "HHC_Harmony",
    "HHC_Scanorama",
    "O2O_Harmony",
    "O2O_Scanorama",

    "HE_Unintegrated",
    "HHC_Unintegrated",
    "O2O_Unintegrated"
  )
names(new_names) <- names(all_datasets)

new_names <- new_names[c(2,4,6,
                         1,3,5,
                         7,9,11,
                         8,10,12,
                         13, 14, 15,
                         16, 18, 20,
                         17, 19, 21,
                         22, 24, 26,
                         23, 25, 27,
                         28, 29, 30)]

# --- Color palette (grouped by method) ---
my_col <- c("#E64B35FF","#E64B35FF","#E64B35FF",
                       "#EE4C97FF",  "#EE4C97FF","#EE4C97FF",
                       "coral1","coral1","coral1",
                       "#E18727FF","#E18727FF","#E18727FF",
                       "#FFDC91FF","#FFDC91FF","#FFDC91FF",
                       "#00A087FF","#00A087FF","#00A087FF",
                       "#4DBBD5FF","#4DBBD5FF","#4DBBD5FF",
                       "#3C5488FF","#3C5488FF","#3C5488FF",
                       "#79AF97FF","#79AF97FF","#79AF97FF",
                       "#80796BFF","#80796BFF","#80796BFF"   )

names(my_col) <- new_names


# ============================================================================
# Extract and save embeddings + metadata per method
# ============================================================================

embeddings_dir = "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/embeddings/"
metadata_dir = "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/metadata/"

# Helper: save embedding and metadata for a set of datasets
save_embeddings_metadata <- function(dataset_names, reduction_name) {
  for(i in dataset_names) {
    my_data <- all_datasets[[i]]
    pca_list <- my_data@reductions[[reduction_name]]
    saveRDS(pca_list, paste0(embeddings_dir, i, ".rds"))
    my_metadata = my_data@meta.data
    saveRDS(my_metadata, paste0(metadata_dir, i, ".rds"))
  }
}

# --- Unintegrated ---
save_embeddings_metadata(
  c("unintegrated_many_higher_expr", "unintegrated_many_higher_homology_conf","unintegrated_one2one_only"),
  "pca")

# --- CCA ---
save_embeddings_metadata(
  c("many_higher_expr_seuratCCA", "many_higher_homology_conf_seuratCCA", "one2one_only_seuratCCA"),
  "pca")

# --- rPCA ---
save_embeddings_metadata(
  c("many_higher_expr_seuratRPCA", "many_higher_homology_conf_seuratRPCA", "one2one_only_seuratRPCA"),
  "pca")

# --- scANVI ---
save_embeddings_metadata(
  c("many_higher_homology_conf_scANVI", "many_higher_expr_scANVI", "one2one_only_scANVI"),
  "X_scANVI")

# --- scVI ---
save_embeddings_metadata(
  c("many_higher_homology_conf_scVI", "many_higher_expr_scVI", "one2one_only_scVI"),
  "X_scVI")

# --- rligerUINMF ---
save_embeddings_metadata(
  c("liger_many_higher_expr_rligerUINMF","liger_many_higher_homology_conf_rligerUINMF" ,"liger_one2one_only_rligerUINMF"),
  "inmf")

# --- LIGER ---
save_embeddings_metadata(
  c("many_higher_expr_LIGER","many_higher_homology_conf_LIGER" ,"one2one_only_LIGER"),
  "iNMF")

# --- fastMNN ---
save_embeddings_metadata(
  c("many_higher_expr_fastMNN","many_higher_homology_conf_fastMNN" ,"one2one_only_fastMNN"),
  "mnn")

# --- Scanorama ---
save_embeddings_metadata(
  c("many_higher_expr_scanorama","many_higher_homology_conf_scanorama" ,"one2one_only_scanorama"),
  "X_scanorama")

# --- Harmony ---
save_embeddings_metadata(
  c("many_higher_expr_harmony","many_higher_homology_conf_harmony" ,"one2one_only_harmony"),
  "X_pca_harmony")


# ============================================================================
# Load saved embeddings and metadata (first 30 PCs)
# ============================================================================

files <- list.files(path = embeddings_dir, pattern = '.rds')
my_embeddings <- list()
for(i in files){
  x <- readRDS(paste0(embeddings_dir, i))
  x <- x@cell.embeddings
  if (ncol(x) > 30) {
    x <- x[, 1:30]
  }
  name <-  sub("\\.rds$", "", i)
  my_embeddings[[name]] <- x
}
rm(x)

files <- list.files(path = metadata_dir, pattern = '.rds')
my_metadata <- list()
for(i in files){
  x <- readRDS(paste0(metadata_dir, i))
  name <-  sub("\\.rds$", "", i)
  my_metadata[[name]] <- x
}
rm(x)


# ============================================================================
# kBET — batch effect test
# ============================================================================

# Wrapper function for kBET across a range of k values
calc_kbet = function(dim_red_embeddings, batch_labels,
                     k_range, method_id) {
  kbet_list = list()
  for (i in seq_along(k_range)) {
    msg = paste("Calculating kBET for", k_range[i],
                "neighbors, please bear with us...")
    message(msg)
    batch_est = kBET(df=dim_red_embeddings,
                     batch=batch_labels,
                     k0 = k_range[i],
                     knn = NULL,
                     do.pca=FALSE,
                     heuristic = FALSE,
                     plot = FALSE,
                     adapt = FALSE,
                     verbose = TRUE)

    plotting_df = data.frame(class = rep(c("observed", "expected"),
                                         each=length(batch_est$stats$kBET.observed)),
                             data = c(batch_est$stats$kBET.observed,
                                      batch_est$stats$kBET.expected))
    plotting_df$k = paste("k_", k_range[i], sep="")
    plotting_df$method = method_id
    kbet_list[[i]] = plotting_df
    names(kbet_list)[i] = paste("kBET_", k_range[i], "_neighbors",
                                sep = "")

  }
  return(kbet_list)
}

# --- Calculate kBET for all methods ---
my_kbet_result <- list()
for(method in names(my_metadata)){
  k_range = c(10, 25, 50, 100, 500, 1000)
  my_batch_lebels <- my_metadata[[method]]$species
  k0 = floor(mean(table(my_batch_lebels)) / 4)
  k_range = c(k_range, k0)

  kbet_df = calc_kbet(my_embeddings[[method]],
                           batch_labels=my_batch_lebels,
                           k_range=k_range,
                           method_id=method)

  kbet_merged_df = rbindlist(kbet_df)
  my_kbet_result[[method]] <- kbet_merged_df
}

saveRDS(my_kbet_result,
        "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/kBET_rejection_rates_df.rds")

# --- Plot kBET rejection rate curves ---
my_kbet_result <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/kBET_rejection_rates_df.rds")
merged_res = rbindlist(my_kbet_result)
merged_res$method_2 <- merged_res$method
for(i in 1:length(new_names)) {
  merged_res$method_2[merged_res$method_2 == names(new_names)[i]] <- new_names[i]
}

methods_ordered <- sort(unique(merged_res$method_2))
methods_ordered <- methods_ordered[c(1,11,21,2,12,22,3,13,23,4,14,24,5,15,25,6,16,26,7,17,27,8,18,28,9,19,29,10,20,30)]
merged_res$method_2 <- factor(merged_res$method_2, levels = methods_ordered)

merged_res$method <- merged_res$method_2

kbet_rates_curve = merged_res %>%
  filter(class == "observed") %>%
  group_by(method, k) %>%
  summarize(mean_obs_rejection_rate = mean(data)) %>%
  ggplot(aes(x=reorder(k, mean_obs_rejection_rate),
             mean_obs_rejection_rate,
             color=method)) +
  geom_line(aes(x=reorder(k, mean_obs_rejection_rate),
                y=mean_obs_rejection_rate, group=method)) +
  geom_point(size=1.5) +
  labs(x="Neighborhood size (k)",
       y="Mean kBET observed rejection rate") +
  custom_theme() +
  theme(aspect.ratio = 1.3,
        legend.position = "bottom", legend.text = element_text(size=6),
        axis.text=element_text(size=10)) +
  scale_color_manual(values= my_col)

kbet_rates_curve
ggsave(filename = 'kBET_plot.pdf', plot = kbet_rates_curve, width = 10, height = 8, units = 'cm', scale = 2)

# --- kBET AUC (area under rejection rate curve) ---
methods = unique(as.character(merged_res$method))
auc_list = list()

for (i in seq_along(methods)) {
  sorted_rej_rates = merged_res %>%
    filter(class == "observed" & method == methods[i]) %>%
    group_by(method, k) %>%
    summarize(mean_obs_rejection_rate = mean(data)) %>%
    arrange(mean_obs_rejection_rate) %>%
    mutate(k_index = seq(1, 7, by=1))

  auc = DescTools::AUC(x=sorted_rej_rates$k_index,
                       y=sorted_rej_rates$mean_obs_rejection_rate,
                       method="spline")
  auc_list[[i]] = auc
  names(auc_list)[i] = methods[i]
}
saveRDS(object = auc_list, file = 'kBET_AUC.rds')


# ============================================================================
# PC Regression — variance explained by batch
# ============================================================================

# For non-Seurat methods (7-30): recompute PCA on RNA assay to get stdev
pcr_list <- list()

for(i in names(all_datasets)[-c(1:6)]) {
  my_data <- all_datasets[[i]]
  DefaultAssay(my_data) <-'RNA'
  my_data <- ScaleData(my_data)
  my_data <-
    FindVariableFeatures(my_data, assay = 'RNA', nfeatures = 2000)

  my_data <- RunPCA(my_data, assay = 'RNA', npcs = 40)

  method_reduc_data = list(x=my_data@reductions$pca@cell.embeddings,
                           sdev=my_data@reductions$pca@stdev)

  batch_labels = as.factor(my_data$species)

  method_pc_reg = kBET::pcRegression(pca.data=method_reduc_data,
                                     batch=batch_labels,
                                     n_top=30)
  df <- data.frame(R2Var = method_pc_reg[["R2Var"]],
                   method = i)

  pcr_list[[i]] <- df
  rm(my_data)
}

# For Seurat-based methods (1-6): PCA stdev already available
for(i in names(all_datasets)[c(1:6)]) {
  my_data <- all_datasets[[i]]
  DefaultAssay(my_data) <-'RNA'

  method_reduc_data = list(x=my_data@reductions$pca@cell.embeddings,
                           sdev=my_data@reductions$pca@stdev)

  batch_labels = as.factor(my_data$species)

  method_pc_reg = kBET::pcRegression(pca.data=method_reduc_data,
                                     batch=batch_labels,
                                     n_top=30)
  df <- data.frame(R2Var = method_pc_reg[["R2Var"]],
                   method = i)

  pcr_list[[i]] <- df
  rm(my_data)
}
saveRDS(pcr_list, 'pcr_list.rds')
pcr_list <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/pcr_list.rds")

# --- Plot PC regression results ---
r2var_df = rbindlist(pcr_list)
r2var_df$method <- factor(r2var_df$method)
levels(r2var_df$method) <- new_names[levels(r2var_df$method)]

r2var_plot = r2var_df %>%
  mutate(method=fct_reorder(method, -log10(R2Var), .desc=TRUE)) %>%
  ggplot(aes(x=method, y=-log10(R2Var),
             fill=method)) +
  geom_col(width = 0.6) +
  xlab("Method") +
  ylab("-log10(R^2 Variance)") +
  custom_theme() +
  theme(aspect.ratio = 0.6,
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust=0.5,
                                   hjust=1, size = 10)) +
  scale_fill_manual(values= my_col)
r2var_plot
ggsave("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/r2var_PC_regression_plot.pdf",
       plot = r2var_plot, width = 3, height = 2, scale = 2)


# ============================================================================
# LISI — Local Inverse Simpson Index (iLISI + cLISI)
# ============================================================================

# Reload embeddings (full, not truncated to 30 PCs)
files <- list.files(path = embeddings_dir, pattern = '.rds')
my_embeddings <- list()
for(i in files){
  x <- readRDS(paste0(embeddings_dir, i))
  name <-  sub("\\.rds$", "", i)
  my_embeddings[[name]] <- x
}
rm(x)

files <- list.files(path = metadata_dir, pattern = '.rds')
my_metadata <- list()
for(i in files){
  x <- readRDS(paste0(metadata_dir, i))
  name <-  sub("\\.rds$", "", i)
  my_metadata[[name]] <- x
}
rm(x)

# --- Calculate iLISI (batch mixing) and cLISI (cell type purity) ---
LISI_list <- list()
for (i in names(my_embeddings)) {
  batch_metadata <-
    data.frame(batch_label = my_metadata[[i]]$species,
               cluster = my_metadata[[i]]$smc_phenotypes)
  rownames(batch_metadata) = rownames(my_metadata[[i]])
  batch_metadata$method = rep(i,
                              length(rownames(batch_metadata)))

  # iLISI: higher = better batch mixing
  ilisi_res = lisi::compute_lisi(my_embeddings[[i]]@cell.embeddings,
                                 batch_metadata,
                                 c("batch_label"))
  batch_metadata$ilisi = ilisi_res$batch_label

  # cLISI: lower = better cell type separation
  clisi_res = lisi::compute_lisi(my_embeddings[[i]]@cell.embeddings,
                                 batch_metadata,
                                 c("cluster"))
  batch_metadata$clisi = clisi_res$cluster
  LISI_list[[i]] <- batch_metadata
}

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration")

merged_lisi_df = rbindlist(LISI_list, use.names = TRUE)
saveRDS(merged_lisi_df,
        "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/LISI_scores_per_cell.rds")

# --- Plot iLISI ---
merged_lisi_df <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/LISI_scores_per_cell.rds")
merged_lisi_df$method <- factor(merged_lisi_df$method)
levels(merged_lisi_df$method) <- new_names[levels(merged_lisi_df$method)]

mean_ilisi_scores = merged_lisi_df %>%
  group_by(method) %>%
  summarize(mean_iLISI = mean(ilisi)) %>%
  mutate(method = fct_reorder(as.factor(method), mean_iLISI, .desc=TRUE)) %>%
  ggplot(aes(x=method, y=mean_iLISI, fill=method)) +
  geom_col(width = 0.6) +
  xlab("Method") +
  ylab("Mean iLISI") +
  custom_theme() +
  scale_fill_manual(values= my_col) +  theme(aspect.ratio = 0.6,
                                                 legend.position = "none",
                                             axis.text.x = element_text(angle = 90, vjust=0.5,
                                                                        hjust=1, size = 10))
mean_ilisi_scores
ggsave(file="/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/integration_mean_iLISI_scores.pdf",
       plot = mean_ilisi_scores, width = 3, height = 3, scale = 3)

# --- Plot cLISI ---
mean_clisi_scores = merged_lisi_df %>%
  group_by(method) %>%
  summarize(mean_cLISI = mean(clisi)) %>%
  mutate(method = fct_reorder(as.factor(method), mean_cLISI)) %>%
  ggplot(aes(x=method, y=mean_cLISI, fill=method)) +
  geom_col(width = 0.6) +
  xlab("Method") +
  ylab("Mean cLISI") +
  custom_theme() +
  scale_fill_manual(values= my_col) +  theme(aspect.ratio = 0.6,
                                                 legend.position = "none",
                                             axis.text.x = element_text(angle = 90, vjust=0.5,
                                                                        hjust=1, size = 10))
mean_clisi_scores
ggsave(file="/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/integration_mean_cLISI_scores.pdf",
       plot = mean_clisi_scores, width = 3, height = 3, scale = 3)


# ============================================================================
# ARI — Adjusted Rand Index (cluster vs phenotype agreement)
# ============================================================================

ARI_list <- list()
for(i in names(all_datasets)){
  my_data <- all_datasets[[i]]
  my_data <- FindNeighbors(my_data, reduction = unique(my_data$embedding))
  my_data <- FindClusters(my_data, graph.name = names(my_data@graphs)[length(names(my_data@graphs))], resolution = 0.4)
  x <- bluster::pairwiseRand(ref = my_data$smc_phenotypes, alt = my_data@meta.data[,grep(pattern = '0.4', colnames(my_data@meta.data))]
               , mode="index")
  ARI_list[[i]] <- x
}

saveRDS(ARI_list, '~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/ARI_smc_phenotypes.rds')

# --- Plot ARI ---
ARI_list <- readRDS("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/ARI_smc_phenotypes.rds")
ARI_list <- data.frame(t(data.frame(ARI_list, use.names = TRUE)))
colnames(ARI_list) <- 'ARI'
ARI_list$method <- rownames(ARI_list)
ARI_list <- ARI_list[1:30,]

ARI_list$method <- factor(ARI_list$method)
levels(ARI_list$method) <- new_names[levels(ARI_list$method)]

ARI_scores = ARI_list %>% mutate(method = fct_reorder(as.factor(method), ARI, .desc=TRUE)) %>%
  ggplot(aes(x=method, y=ARI, fill=method)) +
  geom_col(width = 0.6) +
  xlab("Method") +
  ylab("ARI") +
  custom_theme() +
  scale_fill_manual(values= my_col) +  theme(aspect.ratio = 0.6,
                                             legend.position = "none",
                                             axis.text.x = element_text(angle = 90, vjust=0.5,
                                                                        hjust=1, size = 10))

ARI_scores

ggsave(file="/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/integration_ARI_scores.pdf",
       plot = ARI_scores, width = 3, height = 3, scale = 3)


# ============================================================================
# NMI — Normalized Mutual Information
# ============================================================================

NMI_list <- list()
for(i in names(all_datasets)){
  my_data <- all_datasets[[i]]
  my_data <- FindNeighbors(my_data, reduction = unique(my_data$embedding))
  my_data <- FindClusters(my_data, graph.name = names(my_data@graphs)[length(names(my_data@graphs))], resolution = 0.4)
  x <- NMI(c1 =   my_data$smc_phenotypes,c2 =  my_data@meta.data[,grep(pattern = '0.4', colnames(my_data@meta.data))], variant="sqrt")
  NMI_list[[i]] <- x
}

saveRDS(NMI_list, '~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/NMI_smc_phenotypes.rds')

# --- Plot NMI ---
NMI_list <- readRDS("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/NMI_smc_phenotypes.rds")
NMI_list <- data.frame(t(data.frame(NMI_list, use.names = TRUE)))
colnames(NMI_list) <- 'NMI'
NMI_list$method <- rownames(NMI_list)
NMI_list <- NMI_list[1:30,]

NMI_list$method <- factor(NMI_list$method)
levels(NMI_list$method) <- new_names[levels(NMI_list$method)]

NMI_scores = NMI_list %>% mutate(method = fct_reorder(as.factor(method), NMI, .desc=TRUE)) %>%
  ggplot(aes(x=method, y=NMI, fill=method)) +
  geom_col(width = 0.6) +
  xlab("Method") +
  ylab("NMI") +
  custom_theme() +
  scale_fill_manual(values= my_col) +  theme(aspect.ratio = 0.6,
                                             legend.position = "none",
                                             axis.text.x = element_text(angle = 90, vjust=0.5,
                                                                        hjust=1, size = 10))

NMI_scores

ggsave(file="/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/integration_NMI_scores.pdf",
       plot = NMI_scores, width = 3, height = 3, scale = 3)


# ============================================================================
# cASW — cell-type Average Silhouette Width (across clustering resolutions)
# ============================================================================

res <- seq(0.3, 1, by=0.1)

cASW_scores_list <- list()

for (i in names(all_datasets)[28:30]) {
  my_data <- all_datasets[[i]]
  embed_now <- unique(my_data$embedding)
  dist_matrix = dist(x = my_data@reductions[[embed_now]]@cell.embeddings)
  sil_perDataset_list <- list()

  for (j in seq_along(res)) {
    my_data <-
      FindNeighbors(my_data, reduction = unique(my_data$embedding))
    my_data <-
      FindClusters(my_data,
                   graph.name = names(my_data@graphs)[length(names(my_data@graphs))],
                   resolution = res[j])

    sil_scores <-
      cluster::silhouette(x = as.numeric(x = as.factor(x = my_data@meta.data[, grep(pattern = res[j], colnames(my_data@meta.data))[1]])),
                          dist = dist_matrix)

    sil_df <- data.frame(
      cluster = sil_scores[, 1],
      neighbor = sil_scores[, 2],
      sil_width = sil_scores[, 3],
      res = rep(res[j], length(sil_scores[, 3]))
    )

    sil_perDataset_list[[j]] <- sil_df

    names(sil_perDataset_list)[[j]] = paste("res",
                                            as.character(res[j]),
                                            sep = "_")
  }
  cASW_scores_list[[i]] <- sil_perDataset_list
}

saveRDS(cASW_scores_list, "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/cASW_scores_list.rds")

# --- Plot cASW ---
for(i in names(cASW_scores_list)){
  x <- cASW_scores_list[[i]]
  x <- rbindlist(x)
  x$method <- i
  cASW_scores_list[[i]] <- x
}

res_sil_scores_df <- rbindlist(cASW_scores_list)
res_sil_scores_df$method <- factor(res_sil_scores_df$method)
levels(res_sil_scores_df$method) <- new_names[levels(res_sil_scores_df$method)]

x <- res_sil_scores_df %>%
  group_by(res, method) %>%
  summarize(mean_sil = mean(sil_width))

res_sil_scores_plot = res_sil_scores_df %>%
  group_by(res, method) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x= factor(res), y=mean_sil, fill=method)) + geom_col(position = position_dodge(), width = 1) +
  xlab("Clustering resolution") +
  ylab("Mean Silhouette score") +
  custom_theme() +
  scale_fill_manual(values = my_col) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5,
                                                 hjust=1, size = 10), legend.position = "bottom",  legend.text = element_text(size = 8))

res_sil_scores_plot

ggsave(file="/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/cASW_plot.pdf",
       plot = res_sil_scores_plot, width = 3, height = 3, scale = 3)


# ============================================================================
# bASW — batch Average Silhouette Width (species as grouping)
# ============================================================================

bASW_scores_list <- list()
for (i in names(all_datasets)) {
  my_data <- all_datasets[[i]]
  embed_now <- unique(my_data$embedding)
  dist_matrix = dist(x = my_data@reductions[[embed_now]]@cell.embeddings)
  print(i)

  sil_scores <-
    cluster::silhouette(x = as.numeric(x = as.factor(x = my_data@meta.data$species)),
                        dist = dist_matrix)

  sil_df <- data.frame(method = rep(i, length(sil_scores[, 2])),
                       cluster = sil_scores[, 1],
                       neighbor = sil_scores[, 2],
                       sil_width = sil_scores[, 3])

  bASW_scores_list[[i]] <- sil_df
}

saveRDS(bASW_scores_list, "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/bASW_scores_list.rds")

# --- Plot bASW ---
res_sil_scores_df <- rbindlist(bASW_scores_list)
res_sil_scores_df$method <- factor(res_sil_scores_df$method)
levels(res_sil_scores_df$method) <- new_names[levels(res_sil_scores_df$method)]

x <- res_sil_scores_df %>%
  group_by( method) %>%
  summarize(mean_sil = mean(sil_width))

res_sil_scores_plot = res_sil_scores_df %>%
  group_by( method) %>%
  summarize(mean_sil = mean(sil_width)) %>%
  ggplot(aes(x= method, y=mean_sil, fill=method)) + geom_col(position = position_dodge(), width = 1) +
  xlab("Clustering resolution") +
  ylab("Mean Silhouette score") +
  custom_theme() +
  scale_fill_manual(values = my_col) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5,
                                   hjust=1, size = 10), legend.position = "bottom",  legend.text = element_text(size = 8))

res_sil_scores_plot

ggsave(file="/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/not_succesfull_integration/bASW_plot.pdf",
       plot = res_sil_scores_plot, width = 3, height = 3, scale = 3)
