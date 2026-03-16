# ============================================================================
# Figure 2: Integration overview across 12 method-homology combinations
# Composition of 4-panel plots per dataset: UMAP smc_phenotypes, barplot
# species contribution, UMAP stable clusters, UMAP specificity categories
# ============================================================================

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity")

# --- Libraries ---
library(cluster)
library(Seurat)
library(SeuratObject)
library(cowplot)
library(tidyverse)
library(patchwork)
library(ggsci)
library(scCustomize)
library(sctransform)
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
library(extrafont)
library(ggrastr)
library(ggside)
library(viridis)
library(scales)
library(SCpubr)
library(ggplotify)

# --- Font setup ---
font_import(paths = "/faststorage/project/THOR/diana/Fonts") # import all your fonts
fonts() #get a list of fonts
fonttable()
Sys.setenv(R_GSCMD = "~/miniforge3/envs/r422/bin/gs")

options(scipen = 999)
options(repr.matrix.max.rows = 20)
options(repr.matrix.max.cols = 50)
set.seed(8899)

# --- Load data ---
all_datasets <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/r_py_all_datasets_processed_with_specific_clusters.rds")

# --- Common plot themes ---
my_pal1 <- MetBrewer::met.brewer('Egypt', n = 4, type = 'discrete')
my_pal2 <- MetBrewer::met.brewer('Signac', n = 14, type = 'discrete')

text_sizes <- theme(axis.text.x=element_text(size=9,colour="black"),
                    axis.text.y=element_text(size=9,colour="black"),
                    axis.title.y=element_text(size=9,colour="black", margin = margin(t = 2, l = 1, r = 1, b =2, unit = "pt")),
                    axis.title.x=element_text(size=9,colour="black", margin = margin(t = 1, l = 2, r = 2, b = 1, unit = "pt")),
                    legend.text = element_text(size=9,colour="black"),
                    legend.title = element_text(size=9,colour="black", margin = margin(t = 5, l = 0, r = 0, b = 5, unit = "pt")),
                    legend.key = element_rect(colour="transparent", fill = "transparent"),
                    strip.text.x = element_text(size=6,color = 'black',face="bold", angle=0),
                    strip.text.y = element_text(size=6,color = 'black', face="bold", angle=0, vjust=0.5, hjust=0),
                    axis.ticks= element_line(color = 'black', size=0.2),
                    axis.line = element_line(colour = "gray50", size = 0.2, linetype = "solid"),
                    plot.margin=unit(c(0,0,0,0),"pt"),
                    plot.title=element_text(size=7, face='plain', colour="black"))

common_minimal <- text_sizes + theme(
  plot.background = element_rect(fill = NA,colour = NA),
  strip.background = element_rect(fill = NA,colour = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank()) + theme(
    legend.spacing = unit(0.15, 'cm'),
    legend.key.size = unit(0.2, "cm"))

common_0x <- common_minimal + theme(axis.text.x = element_text(angle=0))

# --- Color palettes ---
pal_celltypes <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_celltypes.rds")
pal_species <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_species.rds")
names(pal_species) <- c('pig', 'human', 'mouse')

pal_cat <- c("#868686FF",   "#8F7700FF",  "#EFC000FF",  "#0073C2FF", "#4ebdab", "#7788ac", "#E64B35B2")
names(pal_cat) <-  names(table(all_datasets$many_higher_homology_conf_scANVI$specificity_category_original))

pal_celltypes_all <- c("#0072B5FF",  "#79AF97FF",  "#374E55FF","#80796BFF", "#FFDC91FF","#EE4C97FF", "#7876B1FF","#20854EFF", "#BC3C29FF",  "#E18727FF")
names(pal_celltypes_all) <- names(table(all_datasets$many_higher_expr_scVI$smc_phenotypes))[c(1,3,4,5,6,7,2,9,10,8)]

pal_celltypes <- brewer.pal(n = 12, name = 'Paired')

# --- Method-homology short names (12 combinations) ---
new_names <- c("HE_CCA","HE_RPCA","HHC_CCA",  "HHC_RPCA","O2O_CCA", "O2O_RPCA", "HE_scANVI","HE_scVI",   "HHC_scANVI" ,    "HHC_scVI", "O2O_scANVI"  ,"O2O_scVI" )
names(new_names) <- names(all_datasets)
new_names <- new_names[c(2,4,6,
                         1,3,5,
                         7,9,11,
                         8,10,12)]

# --- Pre-computed species contribution barplots ---
cluster_prop <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/cluster_proportions_list_plots.rds")

# ============================================================================
# 4-panel composite plots — selected methods (HE_RPCA + HHC_scVI)
# Each panel: SMC phenotypes UMAP, stable clusters UMAP,
#             specificity category UMAP, species contribution barplot
# ============================================================================

dev.new(width = 90, height = 180, unit = "mm",noRStudioGD=TRUE)

myplots <- list()

for(name in names(new_names)[c(1,  10)]) {
  # Subsample human cells to 2500 for balanced visualisation
  my_vec <- sample(names(all_datasets[[name]]$species[all_datasets[[name]]$species == 'human']), size = 2500, replace = F)

  # Panel 1: UMAP colored by SMC phenotypes
  my_smc_pheno = DimPlot_scCustom(
    seurat_object = all_datasets[[name]][, all_datasets[[name]]$species %in% c('pig', 'mouse') |
                                           colnames(all_datasets[[name]]) %in% my_vec],
    reduction = 'umap',
    group.by = 'smc_phenotypes',
    pt.size = 1,
    shuffle = TRUE,
    colors_use = pal_celltypes_all,
  ) + NoLegend() + plot_annotation(title = 'SMC phenotypes')+
    theme(plot.title = element_blank()) +
    scale_color_manual(values = pal_celltypes_all) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text =  element_text(family = "Helvetica")
    )& theme(legend.background = element_rect(fill = "transparent"),
             legend.box.background = element_rect(fill = "transparent"),
             panel.background = element_rect(fill = "transparent"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.background = element_rect(fill = "transparent",
                                            color = NA))

  # Panel 2: UMAP colored by stable clusters (bootstrap-derived)
  my_cell_type = DimPlot_scCustom(seurat_object = all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
                                  reduction = 'umap',
                                  group.by = 'specific_clusters',
                                  pt.size = 1,
                                  shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes, repel = T, label.size = 5
  ) + NoLegend()+ plot_annotation(title = 'Stable clusters') +
    theme(plot.title = element_blank()) +
    scale_color_manual(values = pal_celltypes) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text =  element_text(family = "Helvetica")& theme(legend.background = element_rect(fill = "transparent"),
                                                        legend.box.background = element_rect(fill = "transparent"),
                                                        panel.background = element_rect(fill = "transparent"),
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        plot.background = element_rect(fill = "transparent",
                                                                                       color = NA))
    )

  # Panel 3: species contribution barplot
  my_prop <-
    cluster_prop[[name]] +
    theme(plot.title = element_blank()) + plot_annotation(title = 'Species contribution') +
    theme(legend.position = 'none', text =  element_text(family = "Helvetica", size = 12))& theme(legend.background = element_rect(fill = "transparent"),
                                                                                                  legend.box.background = element_rect(fill = "transparent"),
                                                                                                  panel.background = element_rect(fill = "transparent"),
                                                                                                  panel.grid.major = element_blank(),
                                                                                                  panel.grid.minor = element_blank(),
                                                                                                  plot.background = element_rect(fill = "transparent",
                                                                                                                                 color = NA))

  # Panel 4: UMAP colored by species-specificity category
  my_specific_cl <- DimPlot_scCustom(
    all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
    reduction = 'umap',
    group.by = c('specificity_category_original'),
    pt.size = 1,
    shuffle = TRUE, label = T, label.box = T,  colors_use = pal_cat, repel = T, label.size = 5
  ) +  plot_annotation(title = 'Species-specific SMCs') +
    scale_color_manual(values = pal_cat) +
    theme(legend.position = 'none') +
    theme(plot.title = element_blank(),  text =  element_text(family = "Helvetica")) +
    scale_color_manual(values = pal_cat) + theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) & theme(legend.background = element_rect(fill = "transparent"),
              legend.box.background = element_rect(fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent",
                                             color = NA))

  # Assemble 4-panel composite with cowplot
  fig = ggdraw() +
    draw_plot(my_smc_pheno, x = 0.0, y = .5, width = .5, height = .5) +
    draw_plot(my_cell_type, x = 0.5, y = .5, width = .5, height = .5) +
    draw_plot(my_specific_cl, x = 0.0, y = 0.0, width = .5, height = .5)  +
    draw_plot(my_prop , x = 0.5, y = 0.0, width = .5, height = .5)  +
    draw_plot_label(label = new_names[name], size = 14, family = 'Helvetica',
                    x = c(0), y = c(1.01))

  myplots[[name]] <- fig

}

# saveRDS(myplots, file = 'plot_of_4_specificity_selected.rds')

p <- plot_grid(plotlist = myplots, nrow = 1, ncol = 2, align = 'hv',  axis = 'b', greedy = T,  label_x = -0.1, label_fontfamily = 'Helvetica', label_size = 16)
ggsave(filename = 'plot_of_4_specificity_selected_only2.pdf',  plot = p, device = 'pdf', dpi = 600, width = 10, height = 5, unit = "cm", bg = 'transparent', scale = 3)

# --- Save standalone legends ---
my_smc_pheno = DimPlot_scCustom(seurat_object = all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
                                reduction = 'umap',
                                group.by = 'smc_phenotypes',
                                pt.size = 1,
                                shuffle = TRUE,  colors_use = pal_celltypes_all, repel = T
) +  plot_annotation(title = 'SMC phenotypes') +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes_all) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text =  element_text(family = "Helvetica")
  )& theme(legend.background = element_rect(fill = "transparent"),
           legend.box.background = element_rect(fill = "transparent"),
           panel.background = element_rect(fill = "transparent"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.background = element_rect(fill = "transparent",
                                          color = NA))
ggsave(filename = 'plot_of_4_specificity_selected_only2_LEGEND.pdf',  plot = my_smc_pheno, device = 'pdf', dpi = 600, width = 3.5, height = 2, unit = "cm", bg = 'transparent', scale = 4)

my_prop <-
  cluster_prop[[name]] +
  theme(plot.title = element_blank()) + plot_annotation(title = 'Species contribution') +
  theme(legend.position = 'right', text =  element_text(family = "Helvetica", size = 12))& theme(legend.background = element_rect(fill = "transparent"),
                                                                                                legend.box.background = element_rect(fill = "transparent"),
                                                                                                panel.background = element_rect(fill = "transparent"),
                                                                                                panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),
                                                                                                plot.background = element_rect(fill = "transparent",
                                                                                                                               color = NA))
ggsave(filename = 'plot_of_4_specificity_selected_only2_LEGEND_PROP.pdf',  plot = my_prop, device = 'pdf', dpi = 600, width = 3.5, height = 2, unit = "cm", bg = 'transparent', scale = 4)

# --- ggrastr test ---
ggplot() +
  rasterise(geom_point(aes(carat, price, colour = cut), data=diamonds), dpi=30) +
  geom_point(aes(x=runif(20, 0, 5), y=runif(20, 0, 20000)), size=10, color="black", shape=8)

# ============================================================================
# Rename phenotypes for consistency across figures
# ============================================================================
for(name in names(new_names)) {
  all_datasets[[name]]$smc_phenotypes[all_datasets[[name]]$smc_phenotypes == "Medial"] <- 'Contractile_2'
  all_datasets[[name]]$smc_phenotypes[all_datasets[[name]]$smc_phenotypes == "Transitional"] <- 'Intermediate'
}

# ============================================================================
# 4-panel composites — remaining 10 methods (rasterised for file size)
# ============================================================================
myplots <- list()

for(name in names(new_names)[-c(1,   10)]) {
  my_vec <- sample(names(all_datasets[[name]]$species[all_datasets[[name]]$species == 'human']), size = 2500, replace = F)

  # Panel 1: UMAP SMC phenotypes (rasterised)
  my_smc_pheno = DimPlot_scCustom(seurat_object = all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
                                  reduction = 'umap',
                                  group.by = 'smc_phenotypes',
                                  pt.size = 1,
                                  shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes_all, repel = T, label.size = 5, raster = T, raster.dpi = c(110,110)
  ) + NoLegend()+
    theme(plot.title = element_blank())  +
    scale_color_manual(values = pal_celltypes_all) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text =  element_text(family = "Helvetica")
    )& theme(legend.background = element_rect(fill = "transparent"),
             legend.box.background = element_rect(fill = "transparent"),
             panel.background = element_rect(fill = "transparent"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.background = element_rect(fill = "transparent",
                                            color = NA))

  # Panel 2: UMAP stable clusters (rasterised)
  my_cell_type = DimPlot_scCustom(seurat_object = all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
                                  reduction = 'umap',
                                  group.by = 'specific_clusters',
                                  pt.size = 1,
                                  shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes, repel = T, label.size = 5, raster = T, raster.dpi = c(110,110)
  ) + NoLegend()+
    theme(plot.title = element_blank()) +
    scale_color_manual(values = pal_celltypes) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      theme(legend.background = element_rect(fill = "transparent"),
                                                        legend.box.background = element_rect(fill = "transparent"),
                                                        panel.background = element_rect(fill = "transparent"),
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        plot.background = element_rect(fill = "transparent",
                                                                                       color = NA))
    )

  # Panel 3: species contribution barplot
  my_prop <-
    cluster_prop[[name]] +
    theme(plot.title = element_blank()) +
    theme(legend.position = 'none', text =  element_text(family = "Helvetica", size = 12))& theme(legend.background = element_rect(fill = "transparent"),
                                                                                                  legend.box.background = element_rect(fill = "transparent"),
                                                                                                  panel.background = element_rect(fill = "transparent"),
                                                                                                  panel.grid.major = element_blank(),
                                                                                                  panel.grid.minor = element_blank(),
                                                                                                  plot.background = element_rect(fill = "transparent",
                                                                                                                                 color = NA))

  # Panel 4: UMAP specificity category (rasterised)
  my_specific_cl <- DimPlot_scCustom(
    all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
    reduction = 'umap',
    group.by = c('specificity_category_original'),
    pt.size = 1,
    shuffle = TRUE, label = T, label.box = T,  colors_use = pal_cat, repel = T, label.size = 5, raster = T, raster.dpi = c(110,110)
  ) +
    scale_color_manual(values = pal_cat) +
    theme(legend.position = 'none') +
    theme(plot.title = element_blank(),  text =  element_text(family = "Helvetica")) +
    scale_color_manual(values = pal_cat) + theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) & theme(legend.background = element_rect(fill = "transparent"),
              legend.box.background = element_rect(fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent",
                                             color = NA))

  # Assemble 4-panel composite
  fig = ggdraw() +
    draw_plot(my_smc_pheno, x = 0.0, y = .5, width = .5, height = .5) +
    draw_plot(my_cell_type, x = 0.5, y = .5, width = .5, height = .5) +
    draw_plot(my_specific_cl, x = 0.0, y = 0.0, width = .5, height = .5)  +
    draw_plot(my_prop , x = 0.5, y = 0.0, width = .5, height = .5)  +
    draw_plot_label(label = new_names[name], size = 14, family = 'Helvetica',
                    x = c(0), y = c(1.005))

  myplots[[name]] <- fig

}

saveRDS(myplots, file = 'plot_of_4_specificity_selected_rest.rastr.rds')

p <- plot_grid(plotlist = myplots , nrow = 4, ncol = 3, align = 'hv',  axis = 'b', greedy = T,  label_x = -0.1, label_fontfamily = 'Helvetica', label_size = 16)
ggsave(filename = 'plot_of_4_specificity_selected_rest.rastr.pdf',  plot = p, device = 'pdf', dpi = 600, width = 16, height = 18, unit = "cm", bg = 'transparent', scale = 3)

# ============================================================================
# Individual panel grids — remaining 10 methods saved separately
# ============================================================================

# --- SMC phenotypes UMAPs (no labels, for supplementary grid) ---
myplots <- list()

for(name in names(new_names)[-c(1,   10)]) {
  my_vec <- sample(names(all_datasets[[name]]$species[all_datasets[[name]]$species == 'human']), size = 2500, replace = F)

  my_smc_pheno = DimPlot_scCustom(seurat_object = all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
                                  reduction = 'umap',
                                  group.by = 'smc_phenotypes',
                                  pt.size = 1,
                                  shuffle = TRUE, label = F, label.box = F, colors_use = pal_celltypes_all, repel = T, label.size = 6
  ) + NoLegend()+
    theme(plot.title = element_blank())  +
    scale_color_manual(values = pal_celltypes_all) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text =  element_text(family = "Helvetica")
    )& theme(legend.background = element_rect(fill = "transparent"),
             legend.box.background = element_rect(fill = "transparent"),
             panel.background = element_rect(fill = "transparent"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             plot.background = element_rect(fill = "transparent",
                                            color = NA),
             aspect.ratio = 1)

  fig = ggdraw() +
    draw_plot(my_smc_pheno, x = 0.0, y = .0, width = .95, height = .95) +
    draw_plot_label(label = new_names[name], size = 14, family = 'Helvetica',
                    x = c(0), y = c(1))

  myplots[[name]] <- fig

}

p <- plot_grid(plotlist = myplots, nrow = 2, ncol = 5, align = 'hv',  axis = 'b', greedy = T,  label_x = -0.1)
p
ggsave(filename = 'plot_of_SMC_pheno_specificity_selected_rest.tiff',  plot = p, device = 'tiff', dpi = 600, width = 13, height = 5, unit = "cm", bg = 'white', scale = 3)

# --- Stable clusters UMAPs ---
myplots <- list()

for(name in names(new_names)[-c(1,   10)]) {
  my_vec <- sample(names(all_datasets[[name]]$species[all_datasets[[name]]$species == 'human']), size = 2500, replace = F)

  my_cell_type = DimPlot_scCustom(seurat_object = all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
                                  reduction = 'umap',
                                  group.by = 'specific_clusters',
                                  pt.size = 1,
                                  shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes, repel = T, label.size = 5
  ) + NoLegend()+
    theme(plot.title = element_blank()) +
    scale_color_manual(values = pal_celltypes) +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      text =  element_text(family = "Helvetica")& theme(legend.background = element_rect(fill = "transparent"),
                                                        legend.box.background = element_rect(fill = "transparent"),
                                                        panel.background = element_rect(fill = "transparent"),
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(),
                                                        plot.background = element_rect(fill = "transparent",
                                                                                       color = NA),
                                                        aspect.ratio = 1)
    )

  fig = ggdraw() +
    draw_plot(my_cell_type, x = 0, y = 0, width = 0.95, height = .95) +
    draw_plot_label(label = new_names[name], size = 14, family = 'Helvetica',
                    x = c(0), y = c(1))

  myplots[[name]] <- fig

}
p <- plot_grid(plotlist = myplots, nrow = 2, ncol = 5, align = 'hv',  axis = 'b', greedy = T,  label_x = -0.1)
p
ggsave(filename = 'plot_of_stable_cl_selected_rest.tiff',  plot = p, device = 'tiff', dpi = 600, width = 13, height = 5, unit = "cm", bg = 'white', scale = 3)

# --- Species-specificity category UMAPs ---
# for(name in names(new_names)){
#   all_datasets[[name]]$specificity_category <- gsub(pattern = 'human', replacement = 'hs', x = all_datasets[[name]]$specificity_category)
#   all_datasets[[name]]$specificity_category <- gsub(pattern = 'mouse', replacement = 'mm', x = all_datasets[[name]]$specificity_category)
#   all_datasets[[name]]$specificity_category <- gsub(pattern = 'pig', replacement = 'ss', x = all_datasets[[name]]$specificity_category)
#
# }
#
# for(i in 1:length(pal_cat)){
#   names(pal_cat)[i] <- gsub(pattern = 'human', replacement = 'hs', x = names(pal_cat)[i])
#   names(pal_cat)[i] <- gsub(pattern = 'mouse', replacement = 'mm', x = names(pal_cat)[i])
#   names(pal_cat)[i] <- gsub(pattern = 'pig', replacement = 'ss', x = names(pal_cat)[i])
# }

myplots <- list()

for(name in names(new_names)[-c(1,   10)]) {
  my_vec <- sample(names(all_datasets[[name]]$species[all_datasets[[name]]$species == 'human']), size = 2500, replace = F)

  my_specific_cl <- DimPlot_scCustom(
    all_datasets[[name]][,all_datasets[[name]]$species %in% c('pig', 'mouse')| colnames(all_datasets[[name]]) %in% my_vec],
    reduction = 'umap',
    group.by = c('specificity_category_original'),
    pt.size = 1,
    shuffle = TRUE, label = T, label.box = T,  colors_use = pal_cat, repel = T, label.size = 5
  ) +
    scale_color_manual(values = pal_cat) +
    theme(legend.position = 'none') +
    theme(plot.title = element_blank(),  text =  element_text(family = "Helvetica")) +
    scale_color_manual(values = pal_cat) + theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) & theme(legend.background = element_rect(fill = "transparent"),
              legend.box.background = element_rect(fill = "transparent"),
              panel.background = element_rect(fill = "transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_rect(fill = "transparent",
                                             color = NA),
              aspect.ratio = 1)

  fig = ggdraw() +
    draw_plot(my_specific_cl, x = 0.0, y = 0.0, width = .95, height = .95)  +
    draw_plot_label(label = new_names[name], size = 14, family = 'Helvetica',
                    x = c(0), y = c(1))

  myplots[[name]] <- fig

}
p <- plot_grid(plotlist = myplots, nrow = 2, ncol = 5, align = 'hv',  axis = 'b', greedy = T,  label_x = -0.1)
p
ggsave(filename = 'plot_of_spec_cat_cl_selected_rest.tiff',  plot = p, device = 'tiff', dpi = 600, width = 13, height = 5, unit = "cm", bg = 'white', scale = 3)
