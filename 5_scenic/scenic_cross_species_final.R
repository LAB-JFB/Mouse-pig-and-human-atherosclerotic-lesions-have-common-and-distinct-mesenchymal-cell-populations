# =============================================================================
# SCENIC Cross-Species Regulon Analysis (Human vs Mouse, Elastic Arteries)
# =============================================================================
# Compare SCENIC regulon activity between human and mouse SMCs in integrated
# clusters. Two analyses:
#   1) Binary regulon activity — species contribution ratio per cluster
#   2) Continuous AUC scores — differential regulon markers + dot plot
#
# Input:
#   - Binary regulon matrices: 4.2_binaryRegulonActivity_nonDupl.Rds (human, mouse)
#   - AUC regulon matrices:    3.4_regulonAUC.Rds (human, mouse)
#   - Integrated Seurat:       merged_for_heatmap.rds
#   - Species palette:         pal_species.rds
#   - Proliferating SMC Seurat object `prol` (expected pre-loaded in environment)
#   - Human SMC Seurat object `human` (expected pre-loaded for AUC cell filtering)
#
# Output:
#   - binar_regulons_contrib_of_species_to_clusters.rds / .xlsx
#   - Per-cluster PDF tile plots: <cluster> binar_regulon.pdf
#   - regulon_markers.xlsx (ROC-based markers per cluster)
#   - dotplot_merged_species-specific_regulons_vert.pdf
#   - merged_withTF.rds (Seurat with AUC + AUCBinary assays)
# =============================================================================

library(Seurat)
library(ggplot2)
library(stringr)    # word()
library(scales)     # rescale(), squish
library(openxlsx)
library(RColorBrewer)
library(AUCell)      # onlyNonDuplicatedExtended()

mouse <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/smc_plaque_Alencar_seurat.rds")
human <-  readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/human/athero_comb_smc_seurat.rds")

# TF comparison across species based on SCENIC results for humans and mice
vector_of_specificity <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/SCENIC_cross-species/vector_of_specificity.rds")
my_data <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/HE_RPCA_cell_categorization.rds")
prol <- my_data[,my_data$species == 'human' & my_data$shared_clusters == 'Proliferating']

mouse$specificity <-  vector_of_specificity$Mouse[colnames(mouse)]
human$specificity <-  vector_of_specificity$Human[colnames(human)]
# NOTE: pig object is not loaded in this script (SCENIC is human+mouse only).
# If pig specificity is needed downstream, load pig Seurat here first.
# pig$specificity <-  vector_of_specificity$Pig[colnames(pig)]

table(mouse$specificity)

human$smc_phenotypes_specificity <- paste(human$smc_phenotypes, human$specificity, sep = '_')
# pig$smc_phenotypes_specificity <- paste(pig$smc_phenotypes, pig$specificity, sep = '_')
mouse$smc_phenotypes_specificity <- paste(mouse$smc_phenotypes, mouse$specificity, sep = '_')

my_vec <- c('Contractile_Human-pig-specific', 'DLX5_SMC_Human-specific', 'Pericytes_Human-pig-specific', 'Terminal_Mouse-specific', 'CRTAC1_SMC_NA')

mouse$for_heatmap <- ifelse(mouse$smc_phenotypes_specificity %in% my_vec, mouse$smc_phenotypes_specificity, 'Shared')
human$for_heatmap <- ifelse(human$smc_phenotypes_specificity %in% my_vec, human$smc_phenotypes_specificity, 'Shared')
my_data <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/HE_scVI_shared.rds") # reclustered shared clusters


# We need to devide 'rest categories' to shared phenotypes
mouse$for_heatmap <- factor(mouse$for_heatmap)
human$for_heatmap <- factor(human$for_heatmap)

# For mice
mouse_shar <- my_data[,my_data$species == 'mmusculus']$shared_clusters 
names(mouse_shar) <-  substr(names(mouse_shar), start = 1, stop = 16)
mouse_shar <- mouse_shar[substr(colnames(mouse), start = 1, stop = 16)]
names(mouse_shar) <- colnames(mouse)
mouse$shared_clusters <- mouse_shar
mouse$shared_clusters <- as.character(mouse$shared_clusters)
mouse$shared_clusters[is.na(mouse$shared_clusters)] <- 'specific'
table(mouse$shared_clusters, mouse$for_heatmap)

mouse$xxx <- paste(mouse$shared_clusters, mouse$for_heatmap, sep = '_')

mouse$xxx <- ifelse(mouse$xxx %in% c('Contractile_Shared',  'specific_Shared'), 'Contractile_shared', mouse$xxx) # because mouse contractile cells contribute to pig-uman contractile cluster, but too little, technically theu are just shared contractile
mouse$xxx <- ifelse(mouse$xxx %in% c('Transitional_Shared'), 'Transitional_shared', mouse$xxx )
mouse$xxx <- ifelse(mouse$xxx %in% c('Terminal_Shared', 'Transitional_Terminal_Mouse-specific'), 'Terminal_shared', mouse$xxx )
table(mouse$xxx)

# names(merged$xxx[merged$xxx == 'Terminal_shared'])
mouse$for_heatmap_2 <- mouse$for_heatmap
mouse$for_heatmap_2 <- as.character(mouse$for_heatmap_2)
mouse$for_heatmap_2[names(mouse$xxx[mouse$xxx == 'Terminal_shared'])] <- 'Terminal_shared'
mouse$for_heatmap_2[names(mouse$xxx[mouse$xxx == 'Contractile_shared'])] <- 'Contractile_shared'
mouse$for_heatmap_2[names(mouse$xxx[mouse$xxx == 'Transitional_shared'])] <- 'Transitional_shared'
table(mouse$for_heatmap_2, mouse$smc_phenotypes) # 
table(mouse$for_heatmap_2)

mouse$for_heatmap_2 <-
  factor(mouse$for_heatmap_2, levels = 
           c(
             "Terminal_Mouse-specific",
             'Contractile_shared',
             'Transitional_shared',
             'Terminal_shared' ))

# For humans
human_shar <- my_data[,my_data$species == 'hsapiens']$shared_clusters
names(human_shar) <-  substr(names(human_shar), start = 12, stop = 27)
human_shar <- human_shar[substr(colnames(human), start = 12, stop = 27)]
names(human_shar) <- colnames(human)
human$shared_clusters <- human_shar
human$shared_clusters <- as.character(human$shared_clusters)
human$shared_clusters[is.na(human$shared_clusters)] <- 'specific'
table(human$shared_clusters, human$for_heatmap)

human$xxx <- paste(human$shared_clusters, human$for_heatmap, sep = '_')

human$xxx <- ifelse(human$xxx %in% c('Contractile_Shared'), 'Contractile_shared', human$xxx) # because mouse contractile cells contribute to pig-uman contractile cluster, but too little, technically theu are just shared contractile
human$xxx <- ifelse(human$xxx %in% c('Transitional_Shared',  'specific_Shared'), 'Transitional_shared', human$xxx )
human$xxx <- ifelse(human$xxx %in% c('Terminal_Shared', 'Transitional_Terminal_Mouse-specific'), 'Terminal_shared', human$xxx )
table(human$xxx)

# names(merged$xxx[merged$xxx == 'Terminal_shared'])
human$for_heatmap_2 <- human$for_heatmap
human$for_heatmap_2 <- as.character(human$for_heatmap_2)
human$for_heatmap_2[names(human$xxx[human$xxx == 'Terminal_shared'])] <- 'Terminal_shared'
human$for_heatmap_2[names(human$xxx[human$xxx == 'Contractile_shared'])] <- 'Contractile_shared'
human$for_heatmap_2[names(human$xxx[human$xxx == 'Transitional_shared'])] <- 'Transitional_shared'
table(human$for_heatmap_2, human$smc_phenotypes) # 
table(human$for_heatmap_2)


human$for_heatmap_2 <-
  factor(human$for_heatmap_2, levels = 
           c(
             "DLX5_SMC_Human-specific",
             "Pericytes_Human-pig-specific",
             'Contractile_Human-pig-specific',
             "CRTAC1_SMC_NA",
             'Contractile_shared',
             'Transitional_shared',
             'Terminal_shared' ))




# =============================================================================
# PART 1: Binary Regulon Activity — Species Contribution per Cluster
# =============================================================================

# --- 1.1 Load and harmonise binary regulon matrices --------------------------
# Standardise TF names: extract first word, uppercase, prefix with "tf_"

regulonAUC_human <- readRDS("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/SCENIC_human_carotids/int/4.2_binaryRegulonActivity_nonDupl.Rds")
b <- gsub('_', ' ', rownames(regulonAUC_human))
rownames(regulonAUC_human) <- word(b, 1)
rownames(regulonAUC_human) <- toupper(rownames(regulonAUC_human))
rownames(regulonAUC_human) <- paste0('tf_', rownames(regulonAUC_human))
rownames(regulonAUC_human) <- gsub("[(+)]", "", rownames(regulonAUC_human))

regulonAUC_mouse <- readRDS("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/SCENIC_mouse_Alencar/int/4.2_binaryRegulonActivity_nonDupl.Rds")
b <- gsub('_', ' ', rownames(regulonAUC_mouse))
rownames(regulonAUC_mouse) <- word(b, 1)
rownames(regulonAUC_mouse) <- toupper(rownames(regulonAUC_mouse))
rownames(regulonAUC_mouse) <- paste0('tf_', rownames(regulonAUC_mouse))
rownames(regulonAUC_mouse) <- gsub("[(+)]", "", rownames(regulonAUC_mouse))

# --- 1.2 Merge binary matrices (outer join, fill missing with 0) -------------

merged_matrix <- merge(
  as.data.frame(regulonAUC_human),
  as.data.frame(regulonAUC_mouse),
  by = "row.names",
  all = TRUE
)
rownames(merged_matrix) <- merged_matrix$Row.names
merged_matrix <- merged_matrix[, -1]
colnames(merged_matrix) <- paste0('_', colnames(merged_matrix))
merged_matrix[is.na(merged_matrix)] <- 0

# --- 1.3 Load integrated Seurat and add proliferating cells ------------------

merged <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/merged_for_heatmap.rds")
merged <- merged[, merged$species != 'sscrofa']  # Keep only human + mouse

# Prepare proliferating SMC object (rename barcodes to match merged_matrix format)
metadata_2 <- prol@meta.data
counts_2 <- GetAssay(object = prol, assay = 'RNA')
colnames(counts_2) <- paste0("_", colnames(counts_2))
colnames(counts_2) <- gsub(pattern = '-hsapiens', replacement = '', x = colnames(counts_2))
prol <- CreateSeuratObject(counts = counts_2, assay = 'RNA', meta.data = metadata_2)
prol$species <- 'hsapiens'
prol$for_heatmap_2 <- "Proliferating_shared"

merged <- merge(x = merged, y = prol)
table(merged$species, merged$for_heatmap_2)

# --- 1.4 Attach binary regulon matrix as assay -------------------------------

merged_matrix <- merged_matrix[, colnames(merged_matrix) %in% colnames(merged)]
merged@assays[["AUCBinary"]] <- CreateAssayObject(data = as.matrix(merged_matrix))

# --- 1.5 Compute species contribution ratio per regulon per cluster ----------
# For each integrated cluster:
#   - Keep regulons active in >30% of total cells (both species combined)
#   - Calculate fraction of active cells per species
#   - Compute normalised ratio: (human - mouse) / (human + mouse)
#     Range [-1, 1]: negative = mouse-enriched, positive = human-enriched

my_df_res <- data.frame(
  cluster = character(),
  tf = character(),
  mmusculus_active = numeric(),
  mmusculus_total = numeric(),
  mmusculus_ratio = numeric(),
  hsapiens_active = numeric(),
  hsapiens_total = numeric(),
  hsapiens_ratio = numeric(),
  ratio = numeric()
)

for (k in unique(merged$for_heatmap_2)) {
  metadata <- merged@meta.data
  subset_metadata <- metadata[metadata$for_heatmap_2 == k, ]
  subset_matrix <- merged_matrix[, colnames(merged_matrix) %in% rownames(subset_metadata)]

  # Filter regulons: keep only those active in >30% of cells in this cluster
  good_regulons <- c()
  for (i in rownames(subset_matrix)) {
    xxx <- subset_matrix[i, ]
    good_regulons <- c(good_regulons, ifelse(sum(xxx) / length(xxx) > 0.3, i, NA))
  }
  good_regulons <- na.omit(good_regulons)
  subset_matrix_2 <- subset_matrix[rownames(subset_matrix) %in% good_regulons, ]

  # Calculate per-species activity fractions for each passing regulon
  my_df <- data.frame(
    cluster = character(),
    tf = character(),
    mmusculus_active = numeric(),
    mmusculus_total = numeric(),
    mmusculus_ratio = numeric(),
    hsapiens_active = numeric(),
    hsapiens_total = numeric(),
    hsapiens_ratio = numeric(),
    ratio = numeric()
  )

  for (i in 1:nrow(subset_matrix_2)) {
    gene_expression <- subset_matrix_2[i, ]

    my_df[i, 'cluster'] <- k
    my_df[i, 'tf'] <- rownames(subset_matrix_2)[i]

    # Mouse: count active cells and total cells
    my_df[i, 'mmusculus_active'] <- sum(gene_expression[, colnames(gene_expression) %in% rownames(subset_metadata[subset_metadata$species == 'mmusculus', ])])
    my_df[i, 'mmusculus_total'] <- length(gene_expression[, colnames(gene_expression) %in% rownames(subset_metadata[subset_metadata$species == 'mmusculus', ])])
    my_df[i, 'mmusculus_ratio'] <- my_df[i, 'mmusculus_active'] / my_df[i, 'mmusculus_total']
    my_df[i, 'mmusculus_ratio'][is.na(my_df[i, 'mmusculus_ratio'])] <- 0

    # Human: count active cells and total cells
    my_df[i, 'hsapiens_active'] <- sum(gene_expression[, colnames(gene_expression) %in% rownames(subset_metadata[subset_metadata$species == 'hsapiens', ])])
    my_df[i, 'hsapiens_total'] <- length(gene_expression[, colnames(gene_expression) %in% rownames(subset_metadata[subset_metadata$species == 'hsapiens', ])])
    my_df[i, 'hsapiens_ratio'] <- my_df[i, 'hsapiens_active'] / my_df[i, 'hsapiens_total']
    my_df[i, 'hsapiens_ratio'][is.na(my_df[i, 'hsapiens_ratio'])] <- 0

    # Normalised species ratio: (human - mouse) / (human + mouse)
    my_df[i, 'ratio'] <- (my_df[i, 'hsapiens_ratio'] - my_df[i, 'mmusculus_ratio']) /
      (my_df[i, 'hsapiens_ratio'] + my_df[i, 'mmusculus_ratio'])
  }

  my_df_res <- rbind(my_df_res, my_df)
}

# --- 1.6 Save binary regulon results -----------------------------------------

saveRDS(my_df_res, file = 'binar_regulons_contrib_of_species_to_clusters.rds')
openxlsx::write.xlsx(x = my_df_res, file = 'binar_regulons_contrib_of_species_to_clusters.xlsx')

# --- 1.7 Tile plot: species contribution per regulon per cluster --------------
# Colour scale: teal = mouse-enriched, yellow = balanced, green = human-enriched

pal_species <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_species.rds")

for (i in unique(my_df_res$cluster)) {
  yyy <- my_df_res[my_df_res$cluster == i, c('tf', 'ratio')]
  yyy <- yyy[order(yyy$ratio, decreasing = FALSE), ]

  # Arrange TFs into a roughly square grid
  yyy$x <- rep(1:(trunc(sqrt(length(yyy$tf)))), length.out = nrow(yyy))
  yyy$y <- rep(1:round(length(yyy$tf) / max(yyy$x)),
               each = trunc(sqrt(length(yyy$tf))))[1:length(yyy$ratio)]
  yyy$tf_2 <- word(yyy$tf, start = 2, sep = '_')  # Strip "tf_" prefix for labels

  p <- ggplot(yyy, aes(x = x, y = y, fill = ratio)) +
    geom_tile() +
    geom_text(aes(label = tf_2), color = "black", size = 6) +
    scale_fill_gradientn(
      colours = c('#065465', '#06768d', '#009092',
                  "#D0AC1C", "#D0AC1C",
                  '#006400', '#306844', '#2c4c3b'),
      values = scales::rescale(c(-1, -0.9, -0.5, -0.3, 0.3, 0.5, 0.9, 1)),
      limits = c(-1, 1),
      na.value = "grey50",
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      aspect.ratio = 1
    ) +
    labs(title = paste("Species Contribution to Important Regulons in Integrated", i))

  print(p)
  ggsave(
    filename = paste(i, 'binar_regulon.pdf'),
    plot = p, device = 'pdf',
    units = 'cm', dpi = 400,
    width = 6, height = 6, scale = 3
  )
}

# =============================================================================
# PART 2: Continuous AUC Scores — Differential Regulon Analysis
# =============================================================================

# --- 2.1 Load and harmonise continuous AUC matrices --------------------------

regulonAUC_human <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/SCENIC_human_carotids/int/3.4_regulonAUC.Rds")
regulonAUC_human <- regulonAUC_human[onlyNonDuplicatedExtended(rownames(regulonAUC_human)), ]
regulonAUC_human <- regulonAUC_human@assays@data@listData[["AUC"]]
regulonAUC_human <- regulonAUC_human[, colnames(regulonAUC_human) %in% colnames(human)]
b <- gsub('_', ' ', rownames(regulonAUC_human))
rownames(regulonAUC_human) <- word(b, 1)
rownames(regulonAUC_human) <- toupper(rownames(regulonAUC_human))
rownames(regulonAUC_human) <- paste0('tf_', rownames(regulonAUC_human))
rownames(regulonAUC_human) <- gsub("[(+)]", "", rownames(regulonAUC_human))

regulonAUC_mouse <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/SCENIC_mouse_Alencar/int/3.4_regulonAUC.Rds")
regulonAUC_mouse <- regulonAUC_mouse[onlyNonDuplicatedExtended(rownames(regulonAUC_mouse)), ]
regulonAUC_mouse <- regulonAUC_mouse@assays@data@listData[["AUC"]]
b <- gsub('_', ' ', rownames(regulonAUC_mouse))
rownames(regulonAUC_mouse) <- word(b, 1)
rownames(regulonAUC_mouse) <- toupper(rownames(regulonAUC_mouse))
rownames(regulonAUC_mouse) <- paste0('tf_', rownames(regulonAUC_mouse))
rownames(regulonAUC_mouse) <- gsub("[(+)]", "", rownames(regulonAUC_mouse))

# --- 2.2 Merge AUC matrices and attach as assay ------------------------------

merged_matrix <- merge(
  as.data.frame(regulonAUC_human),
  as.data.frame(regulonAUC_mouse),
  by = "row.names",
  all = TRUE
)
rownames(merged_matrix) <- merged_matrix$Row.names
merged_matrix <- merged_matrix[, -1]
colnames(merged_matrix) <- paste0('_', colnames(merged_matrix))
merged_matrix[is.na(merged_matrix)] <- 0

merged_matrix <- merged_matrix[, colnames(merged_matrix) %in% colnames(merged)]
merged@assays[["AUC"]] <- CreateAssayObject(data = as.matrix(merged_matrix))

# --- 2.3 Create species-split cluster labels (for_heatmap_3) -----------------
# Split shared clusters by species for species-aware ordering

merged$for_heatmap_3 <- merged$for_heatmap_2

for (label in c('Contractile_shared', 'Transitional_shared', 'Terminal_shared', 'Proliferating_shared')) {
  mask <- merged$for_heatmap_3 == label
  merged$for_heatmap_3[mask] <- paste(merged$for_heatmap_3[mask],
                                       merged[, mask]$species, sep = '_')
}

merged$for_heatmap_3 <- factor(merged$for_heatmap_3, levels = c(
  'DLX5_SMC_Human-specific',
  'Pericytes_Human-pig-specific',
  'Terminal_Mouse-specific',
  'Contractile_Human-pig-specific',
  'Contractile_shared_hsapiens',
  'Contractile_shared_mmusculus',
  'Transitional_shared_hsapiens',
  'Transitional_shared_mmusculus',
  'Terminal_shared_hsapiens',
  'Terminal_shared_mmusculus',
  'Proliferating_shared_hsapiens',
  'Proliferating_shared_mmusculus'
))

merged$for_heatmap_2 <- factor(merged$for_heatmap_2, levels = c(
  'DLX5_SMC_Human-specific',
  'Pericytes_Human-pig-specific',
  'Terminal_Mouse-specific',
  'Contractile_Human-pig-specific',
  'Contractile_shared',
  'Transitional_shared',
  'Terminal_shared',
  'Proliferating_shared'
))

# --- 2.4 Find marker regulons per cluster (ROC test) -------------------------

Idents(merged) <- 'for_heatmap_2'
DefaultAssay(merged) <- 'AUC'

tmp <- FindAllMarkers(merged, assay = 'AUC', test.use = 'roc',
                      only.pos = TRUE, min.pct = 0.3)

openxlsx::write.xlsx(tmp, file = "regulon_markers.xlsx",
                     asTable = TRUE, overwrite = TRUE)

# --- 2.5 Select top 9 markers per cluster by avg_log2FC ----------------------

top_markers <- c()
for (i in unique(tmp$cluster)) {
  cluster_markers <- subset(tmp, cluster == i)
  cluster_markers <- cluster_markers[order(cluster_markers$avg_log2FC, decreasing = TRUE), ]
  top_markers <- c(top_markers, head(cluster_markers$gene, 9))
}

# --- 2.6 Dot plot of top regulon markers --------------------------------------

p1 <- DotPlot(
  merged,
  features = unique(top_markers),
  assay = 'AUC',
  scale = TRUE,
  group.by = 'for_heatmap_2'
) +
  coord_flip() +
  RotatedAxis() +
  scale_colour_gradientn(
    name = 'log2 (count + 1)',
    colours = rev(brewer.pal(n = 11, name = "BrBG"))
  ) +
  theme_bw() +
  theme(
    legend.position = 'left',
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks = element_blank()
  ) +
  ylab('')

print(p1)

ggsave(filename = 'dotplot_merged_species-specific_regulons_vert.pdf',
       plot = p1, device = 'pdf',
       height = 9.6, width = 4, units = 'cm', scale = 3)

# --- 2.7 Save final Seurat object with TF assays -----------------------------

saveRDS(merged, 'merged_withTF.rds')
