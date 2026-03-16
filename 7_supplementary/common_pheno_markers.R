# ============================================================================
# Common Phenotype Markers: ACTA2/LUM co-expression & marker validation
# ============================================================================
# Analyses ACTA2-LUM co-expression across species and artery types, then
# validates transitional markers (VCAM1, TNFRSF11B, SERPINF1) and generates
# phenotype-specific FeaturePlots for elastic and coronary datasets.
#
# Input:
#   - /faststorage/.../cluster_specificity/r_py_all_datasets_processed.rds
#   - /faststorage/.../coronaries/after_bengal/r_py_all_datasets_processed_with_specific_clusters.rds
#
# Output (saved to working directory):
#   - ACTA2/LUM co-expression heatmap & stacked bar plots
#   - mesenchymal_lum_acta2_composition.pdf
#   - common_pheno_markers_3_species.pdf / .tiff
#   - macroph_pheno_markers_3_species.pdf / .tiff
#   - myofibr_msc_pheno_markers_3_species.pdf / .tiff
#   - fibromyoc_pheno_markers_3_species.pdf / .tiff
#   - osteob_pheno_markers_3_species3.pdf / .tiff
#   - SEM_markers_3_species.pdf / .tiff
#   - common_pheno_markers_2_species_cor.pdf / .tiff
#   - macroph_pheno_markers_2_species_cor.pdf / .tiff
#   - myofibr_msc_pheno_markers_2_species_cor.pdf / .tiff
#   - fibromyoc_pheno_markers_2_species_cor.pdf / .tiff
#   - osteob_pheno_markers_2_species_cor.pdf / .tiff
# ============================================================================

library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(SCpubr)
library(scCustomize)
library(pheatmap)
library(scales)
library(RColorBrewer)

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/Common_phenotypes")

# ============================================================================
# 1. Load data
# ============================================================================

# Elastic arteries (3 species)
all_datasets <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/r_py_all_datasets_processed.rds")
my_data_el <- all_datasets$many_higher_expr_scVI
rm(all_datasets)

# Coronary arteries (pig + human)
all_datasets_cor <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/after_bengal/r_py_all_datasets_processed_with_specific_clusters.rds")
my_data_cor <- all_datasets_cor$many_higher_expr_seuratRPCA
rm(all_datasets_cor)
gc()

DefaultAssay(my_data_cor) <- 'RNA_genes'
DefaultAssay(my_data_el)  <- 'RNA_genes'

# ============================================================================
# 2. ACTA2 / LUM blend FeaturePlots across species and artery types
# ============================================================================

# Blend plots show ACTA2 (green) vs LUM (magenta) co-expression per species
p1 <- FeaturePlot(my_data_el[, my_data_el$species == 'hsapiens'],
                  features = c('ACTA2', 'LUM'), reduction = 'umap',
                  pt.size = 0.7, order = TRUE, blend = TRUE,
                  cols = c("darkblue", "green", "magenta"), blend.threshold = 0.5) &
  DarkTheme() & NoAxes() & plot_annotation(title = 'Carotid HS')
p1 <- wrap_elements(grid::textGrob('carotid HS')) + p1 + plot_layout(widths = c(1, 3))

p2 <- FeaturePlot(my_data_el[, my_data_el$species == 'sscrofa'],
                  features = c('ACTA2', 'LUM'), reduction = 'umap',
                  pt.size = 0.7, order = TRUE,
                  min.cutoff = "q1", max.cutoff = "q99",
                  blend = TRUE, cols = c("darkblue", "green", "magenta"), blend.threshold = 0) &
  DarkTheme() & NoAxes() & plot_annotation(title = 'Aorta SS')
p2 <- wrap_elements(grid::textGrob('aorta SS')) + p2 + plot_layout(widths = c(1, 3))

p3 <- FeaturePlot(my_data_el[, my_data_el$species == 'mmusculus'],
                  features = c('ACTA2', 'LUM'), reduction = 'umap',
                  pt.size = 0.7, order = TRUE,
                  min.cutoff = "q1", max.cutoff = "q99",
                  blend = TRUE, cols = c("darkblue", "green", "magenta"), blend.threshold = 0) &
  DarkTheme() & NoAxes() & plot_annotation(title = 'BCT MM')
p3 <- wrap_elements(grid::textGrob('BCT MM')) + p3 + plot_layout(widths = c(1, 3))

p4 <- FeaturePlot(my_data_cor[, my_data_cor$species == 'hsapiens'],
                  features = c('ACTA2', 'LUM'), reduction = 'umap',
                  pt.size = 0.7, order = TRUE,
                  min.cutoff = "q1", max.cutoff = "q99",
                  blend = TRUE, cols = c("darkblue", "green", "magenta"), blend.threshold = 0) &
  DarkTheme() & NoAxes() & plot_annotation(title = 'Coronary HS')
p4 <- wrap_elements(grid::textGrob('coronary HS')) + p4 + plot_layout(widths = c(1, 3))

p5 <- FeaturePlot(my_data_cor[, my_data_cor$species == 'sscrofa'],
                  features = c('ACTA2', 'LUM'), reduction = 'umap',
                  pt.size = 0.7, order = TRUE,
                  min.cutoff = "q1", max.cutoff = "q99",
                  blend = TRUE, cols = c("darkblue", "green", "magenta"), blend.threshold = 0) &
  DarkTheme() & NoAxes() & plot_annotation(title = 'Coronary SS')
p5 <- wrap_elements(grid::textGrob('coronary SS')) + p5 + plot_layout(widths = c(1, 3))

p6 <- p1 / p2 / p3 / p4 / p5
p6

# --- Transitional marker FeaturePlots (VCAM1, TNFRSF11B, ASPN) ---
p7 <- SCpubr::do_FeaturePlot(
  sample = my_data_el[, my_data_el$species == 'hsapiens'],
  features = c('VCAM1', 'TNFRSF11B', 'ASPN'), reduction = 'umap', ncol = 2,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.7, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H')

p8 <- SCpubr::do_FeaturePlot(
  sample = my_data_el[, my_data_el$species == 'sscrofa'],
  features = c('VCAM1', 'TNFRSF11B'), reduction = 'umap', ncol = 2,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.7, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H')

p9 <- SCpubr::do_FeaturePlot(
  sample = my_data_el[, my_data_el$species == 'mmusculus'],
  features = c('VCAM1', 'TNFRSF11B'), reduction = 'umap', ncol = 2,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.7, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H')

p10 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor[, my_data_cor$species == 'hsapiens'],
  features = c('VCAM1', 'TNFRSF11B'), reduction = 'umap', ncol = 2,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.7, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H')

p11 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor[, my_data_cor$species == 'sscrofa'],
  features = c('VCAM1', 'TNFRSF11B'), reduction = 'umap', ncol = 2,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.7, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H')

p12 <- p7 / p8 / p9
p13 <- p10 / p11

# ============================================================================
# 3. ACTA2 / LUM expression level heatmap (pig elastic arteries)
# ============================================================================

data_expr <- GetAssayData(my_data_el[, my_data_el$species == 'sscrofa'], slot = "data")[c("ACTA2", "LUM"), ]

expression_bins <- 0:6

bin_expression <- function(expression_data) {
  cut(expression_data,
      breaks = c(expression_bins, Inf),
      include.lowest = TRUE, right = FALSE,
      labels = paste0(expression_bins, "-", expression_bins + 1))
}

acta2_bins <- bin_expression(data_expr["ACTA2", ])
lum_bins   <- bin_expression(data_expr["LUM", ])
expr_df    <- data.frame(ACTA2 = acta2_bins, LUM = lum_bins)

contingency_table <- table(expr_df)

pheatmap(contingency_table,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE, show_colnames = TRUE)

# ============================================================================
# 4. ACTA2/LUM co-expression classification
# ============================================================================
# Classify cells as ACTA2+, LUM+, double-positive, or neither (threshold = 2)

# --- Elastic arteries ---
data_expr <- GetAssayData(my_data_el, slot = "data", assay = 'RNA_genes')

my_data_el$ACTA2_LUM <- ifelse(data_expr["ACTA2", ] > 2 & data_expr["LUM", ] > 2, "ACTA2_LUM", 'NA')
my_data_el$ACTA2_LUM <- ifelse(data_expr["ACTA2", ] > 2 & data_expr["LUM", ] < 2, "ACTA2", my_data_el$ACTA2_LUM)
my_data_el$ACTA2_LUM <- ifelse(data_expr["ACTA2", ] < 2 & data_expr["LUM", ] > 2, "LUM", my_data_el$ACTA2_LUM)

table(my_data_el$ACTA2_LUM)

# --- Coronary arteries ---
data_expr <- GetAssayData(my_data_cor, slot = "data", assay = 'RNA_genes')

my_data_cor$ACTA2_LUM <- ifelse(data_expr["ACTA2", ] > 2 & data_expr["LUM", ] > 2, "ACTA2_LUM", 'NA')
my_data_cor$ACTA2_LUM <- ifelse(data_expr["ACTA2", ] > 2 & data_expr["LUM", ] < 2, "ACTA2", my_data_cor$ACTA2_LUM)
my_data_cor$ACTA2_LUM <- ifelse(data_expr["ACTA2", ] < 2 & data_expr["LUM", ] > 2, "LUM", my_data_cor$ACTA2_LUM)

table(my_data_cor$ACTA2_LUM)

# ============================================================================
# 5. Stacked bar plots: ACTA2/LUM composition across phenotypes
# ============================================================================

# --- Per-phenotype composition (Contractile / Transitional / Terminal) ---
filtered_data <- subset(my_data_el, subset = smc_phenotypes %in% c("Contractile", "Transitional", "Terminal"))

data_for_plot <- filtered_data@meta.data %>%
  group_by(species, smc_phenotypes, ACTA2_LUM) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup() %>%
  mutate(smc_phenotypes = fct_relevel(smc_phenotypes, "Contractile", "Transitional", "Terminal"))

filtered_data <- subset(my_data_cor, subset = smc_phenotypes %in% c("Contractile", "Transitional", "Terminal"))

data_for_plot_cor <- filtered_data@meta.data %>%
  group_by(species, smc_phenotypes, ACTA2_LUM) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup() %>%
  mutate(smc_phenotypes = fct_relevel(smc_phenotypes, "Contractile", "Transitional", "Terminal"))

# --- All-SMC composition ---
data_for_plot_all <- my_data_el@meta.data %>%
  group_by(species, ACTA2_LUM) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

data_for_plot_all_cor <- my_data_cor@meta.data %>%
  group_by(species, ACTA2_LUM) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

bar_p1 <- ggplot(data_for_plot_all, aes(x = species, y = Percentage, fill = ACTA2_LUM)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Species", y = "Percentage of Expression", fill = "LUM and ACTA2 Expression") +
  theme_minimal() +
  scale_fill_manual(values = c("#79AF97FF", "#374E55FF", "#80796BFF", "#7788ac")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) +
  ggtitle(label = 'Elastic')

bar_p2 <- ggplot(data_for_plot_all_cor, aes(x = species, y = Percentage, fill = ACTA2_LUM)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Species", y = "Percentage of Expression", fill = "LUM and ACTA2 Expression") +
  theme_minimal() +
  scale_fill_manual(values = c("#79AF97FF", "#374E55FF", "#80796BFF", "#7788ac")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12)) +
  ggtitle(label = 'Coronary')

ggsave(filename = 'Acta2_Lum_composition_allSMC.pdf', plot = bar_p1 + bar_p2, scale = 2, width = 4, height = 1.7)

# --- Combined elastic + coronary, per-species barplot with percentages ---
filtered_data <- subset(my_data_el, subset = smc_phenotypes %in% c("Contractile", "Transitional", "Terminal"))
data_for_plot <- filtered_data@meta.data[, c('species', 'ACTA2_LUM')] %>%
  group_by(species, ACTA2_LUM) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(species) %>%
  mutate(Total_Count = sum(Count),
         Percentage = (Count / Total_Count) * 100) %>%
  ungroup()

filtered_data <- subset(my_data_cor, subset = smc_phenotypes %in% c("Contractile", "Transitional", "Terminal"))
data_for_plot_cor <- filtered_data@meta.data[, c('species', 'ACTA2_LUM')] %>%
  group_by(species, ACTA2_LUM) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(species) %>%
  mutate(Total_Count = sum(Count),
         Percentage = (Count / Total_Count) * 100) %>%
  ungroup()

# Label species with artery type and combine
data_for_plot_cor$species <- paste(data_for_plot_cor$species, 'Coronary', sep = '.')
data_for_plot$species     <- paste(data_for_plot$species, 'Elastic', sep = '.')
data_combined <- rbind(data_for_plot, data_for_plot_cor)

# NOTE: my_pal1 must be defined before running this section
# (e.g., loaded from a palette .rds or defined as a named vector)
bar_combined <- ggplot(data_combined, aes(x = species, y = Percentage, fill = ACTA2_LUM)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format()) +
  geom_text(aes(label = percent(Percentage / 100, accuracy = 1)),
            position = position_fill(vjust = 0.5), size = 3, color = "white") +
  labs(x = "SMC Phenotypes", y = "Percentage of Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        aspect.ratio = 2) +
  scale_fill_manual(values = my_pal1[c(3, 4, 1, 2)])

ggsave(filename = 'mesenchymal_lum_acta2_composition.pdf',
       plot = bar_combined, width = 6, height = 3, scale = 4, units = 'cm')

# ============================================================================
# 6. VCAM1 / TNFRSF11B / SERPINF1 marker status across phenotypes
# ============================================================================
# Add binary marker status columns to elastic artery data

DefaultAssay(my_data_el) <- 'RNA_genes'

my_data_el <- my_data_el %>%
  AddMetaData(metadata = ifelse(GetAssayData(., slot = "data")["VCAM1", ] > 1,
                                "VCAM1", "noVCAM1"), col.name = "VCAM1_Status")

my_data_el <- my_data_el %>%
  AddMetaData(metadata = ifelse(GetAssayData(., slot = "data")["TNFRSF11B", ] > 1,
                                "TNFRSF11B", "noTNFRSF11B"), col.name = "TNFRSF11B_Status")

my_data_el <- my_data_el %>%
  AddMetaData(metadata = ifelse(GetAssayData(., slot = "data")["SERPINF1", ] > 1,
                                "SERPINF1", "noSERPINF1"), col.name = "SERPINF1_Status")

# Subset to main phenotypes
filtered_data <- subset(my_data_el, subset = smc_phenotypes %in% c("Contractile", "Transitional", "Terminal"))
filtered_data$smc_phenotypes <- factor(as.character(filtered_data$smc_phenotypes),
                                       levels = c("Contractile", "Transitional", "Terminal"))

# --- VCAM1 ---
data_for_plot <- filtered_data@meta.data %>%
  group_by(smc_phenotypes, ACTA2_LUM, species, VCAM1_Status) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(smc_phenotypes, ACTA2_LUM, species) %>%
  mutate(Percentage = Count / sum(Count) * 100)

p_vcam1 <- ggplot(data_for_plot, aes(x = smc_phenotypes, y = Percentage, fill = VCAM1_Status)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(species ~ ACTA2_LUM, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "SMC Phenotypes", y = "Percentage of Marker-Positive Cells", fill = "Marker Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

# --- TNFRSF11B ---
data_for_plot <- filtered_data@meta.data %>%
  group_by(smc_phenotypes, ACTA2_LUM, species, TNFRSF11B_Status) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(smc_phenotypes, ACTA2_LUM, species) %>%
  mutate(Percentage = Count / sum(Count) * 100)

p_tnfrsf <- ggplot(data_for_plot, aes(x = smc_phenotypes, y = Percentage, fill = TNFRSF11B_Status)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(species ~ ACTA2_LUM, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "SMC Phenotypes", y = "Percentage of Marker-Positive Cells", fill = "Marker Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

# --- SERPINF1 ---
data_for_plot <- filtered_data@meta.data %>%
  group_by(smc_phenotypes, ACTA2_LUM, species, SERPINF1_Status) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(smc_phenotypes, ACTA2_LUM, species) %>%
  mutate(Percentage = Count / sum(Count) * 100)

p_serpinf <- ggplot(data_for_plot, aes(x = smc_phenotypes, y = Percentage, fill = SERPINF1_Status)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(species ~ ACTA2_LUM, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "SMC Phenotypes", y = "Percentage of Marker-Positive Cells", fill = "Marker Status") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10))

p_vcam1 / p_tnfrsf / p_serpinf

# ============================================================================
# 7. CNN1/TNFRSF11B/LUM triple classification (elastic arteries)
# ============================================================================
# Classify cells into 7 categories based on CNN1, TNFRSF11B, LUM (threshold = 1)

data_expr <- GetAssayData(my_data_el, slot = "data", assay = 'RNA_genes')

cnn1_pos     <- data_expr["CNN1", ] > 1
tnfrsf11b_pos <- data_expr["TNFRSF11B", ] > 1
lum_pos      <- data_expr["LUM", ] > 1

my_data_el$CNN1_LUM <- case_when(
  cnn1_pos & !tnfrsf11b_pos & !lum_pos  ~ "CNN1",
  cnn1_pos &  tnfrsf11b_pos & !lum_pos  ~ "CNN1_TNFRSF11B",
  cnn1_pos & !tnfrsf11b_pos &  lum_pos  ~ "CNN1_LUM",
  !cnn1_pos &  tnfrsf11b_pos & !lum_pos ~ "TNFRSF11B",
  !cnn1_pos &  tnfrsf11b_pos &  lum_pos ~ "TNFRSF11B_LUM",
  !cnn1_pos & !tnfrsf11b_pos &  lum_pos ~ "LUM",
  cnn1_pos &  tnfrsf11b_pos &  lum_pos  ~ "Triple_pos",
  TRUE                                   ~ "NA"
)

table(my_data_el$CNN1_LUM)
DimPlot_scCustom(my_data_el, group.by = 'CNN1_LUM', split.by = 'species',
                 pt.size = 1, reduction = 'umap')

# ============================================================================
# 8. Phenotype marker FeaturePlots — Elastic arteries (3 species)
# ============================================================================

# --- Common phenotype markers (CNN1, MYH11, VCAM1, TNFRSF11B, LUM, COL6A3) ---
P1 <- SCpubr::do_FeaturePlot(
  sample = my_data_el,
  features = c('CNN1', 'MYH11', 'VCAM1', 'TNFRSF11B', 'LUM', 'COL6A3'),
  split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'common_pheno_markers_3_species.tiff', plot = P1, device = 'tiff',
       dpi = 600, width = 20, height = 10, unit = "cm", bg = 'white', scale = 3.5)
ggsave(filename = 'common_pheno_markers_3_species.pdf', plot = P1, device = 'pdf',
       dpi = 600, width = 20, height = 10, unit = "cm", bg = 'white', scale = 3.5)

# --- Macrophage-like markers (LGALS3, LAMP2) ---
P3 <- SCpubr::do_FeaturePlot(
  sample = my_data_el,
  features = c('LGALS3', 'LAMP2'), split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'macroph_pheno_markers_3_species.tiff', plot = P3, device = 'tiff',
       dpi = 600, width = 10, height = 10, unit = "cm", bg = 'white', scale = 1.8)
ggsave(filename = 'macroph_pheno_markers_3_species.pdf', plot = P3, device = 'pdf',
       dpi = 600, width = 10, height = 10, unit = "cm", bg = 'white', scale = 1.8)

# --- Myofibroblast / MSC markers (ENG, PDGFRB, S100A4) ---
P4 <- SCpubr::do_FeaturePlot(
  sample = my_data_el,
  features = c('ENG', 'PDGFRB', 'S100A4'), split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'myofibr_msc_pheno_markers_3_species.tiff', plot = P4, device = 'tiff',
       dpi = 600, width = 15, height = 10, unit = "cm", bg = 'white', scale = 3)
ggsave(filename = 'myofibr_msc_pheno_markers_3_species.pdf', plot = P4, device = 'pdf',
       dpi = 600, width = 15, height = 10, unit = "cm", bg = 'white', scale = 3)

# --- Fibromyocyte markers (FN1, TNFRSF11B, COL1A1, BGN, DCN) ---
P5 <- SCpubr::do_FeaturePlot(
  sample = my_data_el,
  features = c('FN1', 'TNFRSF11B', 'COL1A1', 'BGN', 'DCN'),
  split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'fibromyoc_pheno_markers_3_species.tiff', plot = P5, device = 'tiff',
       dpi = 600, width = 20, height = 10, unit = "cm", bg = 'white', scale = 2)
ggsave(filename = 'fibromyoc_pheno_markers_3_species.pdf', plot = P5, device = 'pdf',
       dpi = 600, width = 20, height = 10, unit = "cm", bg = 'white', scale = 2)

# --- Osteoblast-like markers (SOX9, TNFRSF11B, RUNX2) ---
P6_el <- SCpubr::do_FeaturePlot(
  sample = my_data_el,
  features = c('SOX9', 'TNFRSF11B', 'RUNX2'), split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  sequential.palette = 'YlOrRd') & theme(aspect.ratio = 1)

P6_cor <- SCpubr::do_FeaturePlot(
  sample = my_data_cor,
  features = c('SOX9', 'TNFRSF11B', 'RUNX2'), split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  sequential.palette = 'YlOrRd') & theme(aspect.ratio = 1)

ggsave(filename = 'osteob_pheno_markers_3_species3.tiff', plot = P6_el, device = 'tiff',
       dpi = 600, width = 13, height = 10, unit = "cm", bg = 'white', scale = 3)
ggsave(filename = 'osteob_pheno_markers_3_species3.pdf', plot = P6_el, device = 'pdf',
       dpi = 600, width = 13, height = 10, unit = "cm", bg = 'white', scale = 3)

# --- SEM markers (LGALS3, VCAM1, FBLN5, TIMP1, FN1) ---
P6_sem <- SCpubr::do_FeaturePlot(
  sample = my_data_el,
  features = c('LGALS3', 'VCAM1', 'FBLN5', 'TIMP1', 'FN1'),
  split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'SEM_markers_3_species.tiff', plot = P6_sem, device = 'tiff',
       dpi = 600, width = 20, height = 10, unit = "cm", bg = 'white', scale = 3)
ggsave(filename = 'SEM_markers_3_species.pdf', plot = P6_sem, device = 'pdf',
       dpi = 600, width = 20, height = 10, unit = "cm", bg = 'white', scale = 3)

# ============================================================================
# 9. Phenotype marker FeaturePlots — Coronary arteries (2 species)
# ============================================================================

# --- Common phenotype markers ---
P7 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor,
  features = c('CNN1', 'MYH11', 'VCAM1', 'TNFRSF11B', 'LUM', 'COL6A3'),
  split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'common_pheno_markers_2_species_cor.tiff', plot = P7, device = 'tiff',
       dpi = 600, width = 20, height = 8, unit = "cm", bg = 'white', scale = 3.5)
ggsave(filename = 'common_pheno_markers_2_species_cor.pdf', plot = P7, device = 'pdf',
       dpi = 600, width = 20, height = 8, unit = "cm", bg = 'white', scale = 3.5)

# --- Macrophage-like markers ---
P9 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor,
  features = c('LGALS3', 'LAMP2'), split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'macroph_pheno_markers_2_species_cor.tiff', plot = P9, device = 'tiff',
       dpi = 600, width = 15, height = 8, unit = "cm", bg = 'white', scale = 2.8)
ggsave(filename = 'macroph_pheno_markers_2_species_cor.pdf', plot = P9, device = 'pdf',
       dpi = 600, width = 15, height = 8, unit = "cm", bg = 'white', scale = 2.8)

# --- Myofibroblast / MSC markers ---
P10 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor,
  features = c('ENG', 'PDGFRB', 'S100A4'), split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'myofibr_msc_pheno_markers_2_species_cor.tiff', plot = P10, device = 'tiff',
       dpi = 600, width = 15, height = 8, unit = "cm", bg = 'white', scale = 2.8)
ggsave(filename = 'myofibr_msc_pheno_markers_2_species_cor.pdf', plot = P10, device = 'pdf',
       dpi = 600, width = 15, height = 8, unit = "cm", bg = 'white', scale = 2.8)

# --- Fibromyocyte markers ---
P11 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor,
  features = c('FN1', 'TNFRSF11B', 'COL1A1', 'BGN', 'DCN'),
  split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  use_viridis = TRUE, viridis.palette = 'H') & theme(aspect.ratio = 1)

ggsave(filename = 'fibromyoc_pheno_markers_2_species_cor.tiff', plot = P11, device = 'tiff',
       dpi = 600, width = 15, height = 8, unit = "cm", bg = 'white', scale = 3.5)
ggsave(filename = 'fibromyoc_pheno_markers_2_species_cor.pdf', plot = P11, device = 'pdf',
       dpi = 600, width = 15, height = 8, unit = "cm", bg = 'white', scale = 3.5)

# --- Osteoblast-like markers ---
P12 <- SCpubr::do_FeaturePlot(
  sample = my_data_cor,
  features = c('RUNX2', 'MSX2', 'CYTL1', 'TRPV4', 'SOX9'),
  split.by = 'species', reduction = 'umap', ncol = 1,
  order = FALSE, assay = 'RNA_genes', pt.size = 0.6, plot.axes = FALSE,
  sequential.palette = 'YlOrRd') & theme(aspect.ratio = 1)

ggsave(filename = 'osteob_pheno_markers_2_species_cor.tiff', plot = P12, device = 'tiff',
       dpi = 600, width = 12, height = 8, unit = "cm", bg = 'white', scale = 2.8)
ggsave(filename = 'osteob_pheno_markers_2_species_cor.pdf', plot = P12, device = 'pdf',
       dpi = 600, width = 12, height = 8, unit = "cm", bg = 'white', scale = 2.8)
