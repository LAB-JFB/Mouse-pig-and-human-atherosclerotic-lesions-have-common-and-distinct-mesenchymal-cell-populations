# ============================================================================
# Figure 1: Individual species dotplots and UMAP dimplots for elastic arteries
# SMC phenotype marker expression in human carotid, pig aorta, and mouse aortic root
# ============================================================================

setwd("~/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC")
set.seed(1984)
library(Seurat)
library(tidyverse)
library(patchwork)
library(scCustomize)
library(Matrix)
library(sctransform)
library(MatrixGenerics)
library(openxlsx)
library(RColorBrewer)
library(MetBrewer)
library(data.table)
library(pheatmap)
library(corrplot)
library(DESeq2)
library(ggh4x)
library(scales)
library(purrr)
library(SCpubr)

# --- Common theme for publication-quality plots ---
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

# --- Load individual species Seurat objects ---
athero_comb_smc <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/human/athero_comb_smc_seurat.rds")
pig_plaque <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SMC/bengal/pig_plaque_seurat.rds")
smc_plaque <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/smc_plaque_Alencar_seurat.rds")

pal_celltypes <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_celltypes.rds")

# ============================================================================
# Human carotid SMCs — dotplot and UMAP
# ============================================================================
setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/human")

# Reorder phenotype levels for plotting
Idents(athero_comb_smc) <- 'smc_phenotypes'
levels(Idents(athero_comb_smc))
athero_comb_smc$smc_phenotypes <- factor(athero_comb_smc$smc_phenotypes, levels = levels(athero_comb_smc$smc_phenotypes)[c(1,3,2,4,5,6)])
Idents(athero_comb_smc) <- 'smc_phenotypes'

# Find top markers per phenotype
markers_0.3 <- FindAllMarkers(athero_comb_smc, assay = 'RNA', logfc.threshold = 0.25,  min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 5))
}

# Manual marker replacements for display
top_markers[3] <- 'ACTA2'
top_markers[5] <- 'LMOD1'
top_markers[16] <- 'DLX5'
top_markers[21] <- 'TINAGL1'
top_markers[22] <- 'SOCS3'
top_markers[23] <- 'NOTCH3'
top_markers[24] <- 'GJA4'
top_markers[25] <- 'RGS5'
# top_markers[25] <- 'ABCC9'

# Dotplot
p1 <- DotPlot(
  athero_comb_smc,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'smc_phenotypes'
) +
  scale_colour_gradientn(name = 'log2 (count + 1)',
                         colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme_bw() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank()) + coord_fixed(ratio = 1)+ text_sizes
print(p1)
ggsave(filename = 'dorplot_human_carotids_smc_phenotypes.pdf', plot = p1, device = 'pdf', width = 8, height = 3, units = 'cm', scale = 3.3)

# UMAP dimplot
pal_celltypes_human <- c("#0072B5FF", "#BC3C29FF" ,"#E18727FF", "#80796BFF", "#20854EFF",  "#374E55FF")[c(1,3,2,4,5,6)]

p2 <- DimPlot_scCustom(seurat_object = athero_comb_smc,
                       reduction = 'umap',
                       group.by = 'smc_phenotypes',
                       pt.size = 1.5,
                       shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes_human, repel = T, label.size = 6
) + NoLegend()+
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
ggsave(filename = 'dimplot_human_carotids_smc_phenotypes.pdf', plot = p2, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 2.5)


# ============================================================================
# Pig aorta SMCs — dotplot and UMAP
# ============================================================================
setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/pig")

# Reorder phenotype levels for plotting
levels(Idents(pig_plaque))
pig_plaque$smc_phenotypes <- factor(pig_plaque$smc_phenotypes, levels = levels(pig_plaque$smc_phenotypes)[c(2,1,3,5,4,6,7)])
Idents(pig_plaque) <- 'smc_phenotypes'
pal_celltypes <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/pal_celltypes.rds")

pal_celltypes_pig <- c("#E18727FF", "#0072B5FF",  "#BC3C29FF", "#20854EFF", "#7876B1FF" ,"#FFDC91FF", "#EE4C97FF")[c(2,1,3,5,4,6,7)]

# Find top markers per phenotype
markers_0.3 <- FindAllMarkers(pig_plaque, assay = 'RNA', logfc.threshold = 0.25, min.pct = 0.4, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 4))
}

# Manual marker replacements for display
top_markers[18] <- 'RGS16'
top_markers[19] <- 'TINAGL1'
top_markers[17] <- 'SOCS3'

# Dotplot
p3 <- DotPlot(
  pig_plaque,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'smc_phenotypes'
) +
  scale_colour_gradientn(name = 'log2 (count + 1)',colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme_bw() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank()) + coord_fixed(ratio = 1)+ text_sizes
print(p3)

# Rename Ensembl ID to gene symbol for display
levels(p3@data$features.plot)[levels(p3@data$features.plot) == "ENSSSCG00000032436"] <- "IFITM1"

ggsave(filename = 'dorplot_pig_smc_phenotypes.pdf', plot = p1, device = 'pdf', width = 8, height = 4.5, units = 'cm', scale = 3.2)

# UMAP dimplot
p4 <- DimPlot_scCustom(seurat_object = pig_plaque,
                       reduction = 'umap',
                       group.by = 'smc_phenotypes',
                       pt.size = 1.5,
                       shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes_pig, repel = T, label.size = 6
) + NoLegend()+
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
print(p4)
ggsave(filename = 'dimplot_pig_smc_phenotypes.pdf', plot = p4, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 2.5)

# ============================================================================
# Mouse aortic root SMCs (Alencar et al.) — dotplot and UMAP
# ============================================================================
setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse")

# Reorder phenotype levels for plotting
levels(Idents(smc_plaque))
smc_plaque$smc_phenotypes <- factor(smc_plaque$smc_phenotypes, levels = levels(smc_plaque$smc_phenotypes)[c(3,4,1,2)])

Idents(smc_plaque) <- 'smc_phenotypes'

# Find top markers per phenotype
markers_0.3 <- FindAllMarkers(smc_plaque, assay = 'RNA', logfc.threshold = 0.25,   min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 8))
}
top_markers <- top_markers[-30]
top_markers[19] <- 'Lum'

# Dotplot
p5 <- DotPlot(
  smc_plaque,
  features = unique(top_markers),
  assay = 'RNA',
  scale = T,
  group.by = 'smc_phenotypes'
) +
  scale_colour_gradientn(name = 'log2 (count + 1)', colours = rev(brewer.pal(n = 11, name = "Spectral")))  +
  theme_bw() +
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab('') +
  theme(axis.ticks = element_blank()) + coord_fixed(ratio = 1) + text_sizes
print(p5)
ggsave(filename = 'dotplot_Alencar_0.24.pdf', plot = p1, device = 'pdf', width = 11, height = 6, units = 'cm', scale = 2.5)

# UMAP dimplot — SMC phenotypes
pal_celltypes_mouse <- c("#E18727FF","#BC3C29FF",  "#0072B5FF", "#79AF97FF")[c(3,4,1,2)]
p6 <- DimPlot_scCustom(seurat_object = smc_plaque,
                       reduction = 'umap',
                       group.by = 'smc_phenotypes',
                       pt.size = 1.5,
                       shuffle = TRUE, label = T, label.box = T, colors_use = pal_celltypes_mouse, repel = T, label.size = 6
) + NoLegend()+
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes_mouse) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
print(p6)
ggsave(filename = 'dimplot_mouse_Alencar_smc_phenotypes.pdf', plot = p6, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 2.5)

# UMAP dimplot — colored by sample origin
pal_celltypes_mouse <- c("#E18727FF","#BC3C29FF",  "#0072B5FF", "#79AF97FF")

p2 <- DimPlot_scCustom(seurat_object = smc_plaque,
                       reduction = 'umap',
                       group.by = 'orig.ident',
                       pt.size = 1,
                       shuffle = TRUE, colors_use = pal_celltypes, repel = T
) +
  theme(plot.title = element_blank()) +
  scale_color_manual(values = pal_celltypes) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), aspect.ratio = 1
  )
print(p2)
ggsave(filename = 'dimplot_mouse_Alencar_origin.pdf', plot = p2, device = 'pdf', width = 5, height = 3, units = 'cm', scale = 3)

# ============================================================================
# Combined cross-species dotplot (human + pig + mouse)
# ============================================================================

# Merge dotplot data from all three species
df135 <- rbind(mutate(p1$data, organism = "human"),
              mutate(p3$data, organism = "pig"),
              mutate(p5$data, organism = "mouse"))

# Rename Ensembl IDs / aliases to gene symbols
levels(df135$features.plot)[levels(df135$features.plot) == "ENSSSCG00000032436"] <- "IFITM1"
levels(df135$features.plot)[levels(df135$features.plot) == "3110079O15Rik"] <- "Snorc"

# Get number of unique clusters per organism (for panel sizing)
nids <- df135 %>%
  group_by(organism) %>%
  distinct(id) %>%
  dplyr::count() %>%
  pull(n)

# Combined faceted dotplot with ggh4x::force_panelsizes
p135 <- df135 %>%
  ggplot(aes(x = features.plot, y = id, fill = avg.exp.scaled, size = pct.exp)) +
  geom_point(shape = 21, stroke = 0.1, color = "grey30") +
  scale_fill_gradientn(colors = rev(brewer.pal(9, "RdBu"))) +
                       #limits = c(-1.5, 1.5), oob = squish) +
  facet_wrap(~ organism, ncol = 1, scales = "free") +
  force_panelsizes(rows = nids/2) +
  theme_bw() +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, face = "italic"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  labs(x = NULL, y = NULL, title = NULL)
p135
ggsave(filename = 'dotplot_all_species.pdf',
       plot = p135, width = 9, height = 9, units = 'cm', scale = 2.5)

# --- Diverging color scale centered at zero for combined dotplots ---
# mypal <- met.brewer("OKeeffe1", n = 11)
mypal <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
grad101 <- colorRampPalette(rev(mypal))(101)

# Helper: apply zero-centered diverging gradient to a dotplot
apply_grad <- function(p) {
  rng <- range(p$data$avg.exp.scaled, na.rm = TRUE)
  vals <- rescale(c(seq(rng[1], 0, length.out = 51),
                    seq(0,      rng[2], length.out = 51)[-1]))
  p + scale_color_gradientn(colors = grad101, values = vals) +
    labs(x = NULL, y = NULL, title = NULL)
}

# Attempt 1: patchwork with apply_grad
pw135 <- (apply_grad(p1) + apply_grad(p3) + apply_grad(p5)) +
  patchwork::plot_layout(ncol = 1, guides = "collect", heights = c(1,1,1)) &
  theme(legend.position = "bottom")

# Attempt 2: manual per-panel rescaling
pw135 <- (p1 + scale_color_gradientn(
  colors = colorRampPalette(rev(mypal))(101),
  values = rescale(c(seq(min(p1$data$avg.exp.scaled), 0, length.out = 51),
                     seq(0, max(p1$data$avg.exp.scaled), length.out = 51)[-1]))) +
    labs(x = NULL, y = NULL, title = NULL)) +
  (p5 + scale_color_gradientn(
    colors = colorRampPalette(rev(mypal))(101),
    values = rescale(c(seq(min(p3$data$avg.exp.scaled), 0, length.out = 51),
                       seq(0, max(p3$data$avg.exp.scaled), length.out = 51)[-1]))) +
     labs(x = NULL, y = NULL, title = NULL)) +
  (p3 + scale_color_gradientn(
    colors = colorRampPalette(rev(mypal))(101),
    values = rescale(c(seq(min(p5$data$avg.exp.scaled), 0, length.out = 51),
                       seq(0, max(p5$data$avg.exp.scaled), length.out = 51)[-1]))) +
     labs(x = NULL, y = NULL, title = NULL)) +
  plot_layout(ncol = 1, heights = nids/2, widths = rep(1, 3))
pw135
ggsave(filename = 'dotplot_all_species_pw_4.pdf',
       plot = pw135, width = 9, height = 12, units = 'cm', scale = 2.5)

# Attempt 3: map + wrap_plots
list(p1, p3, p5) |>
  map(apply_grad) |>
  wrap_plots(ncol = 1,
             heights = if (length(nids) >= 3) nids[1:3] / 2 else rep(1, 3),
             widths = 1)

pw135 <- (apply_grad(p1) + apply_grad(p3) + apply_grad(p5)) +
  patchwork::plot_layout(ncol = 1, heights = c(1,1,1))

# ============================================================================
# Blended ACTA2/LUM FeaturePlots — contractile vs terminal gradient
# ============================================================================

SCpubr::do_FeaturePlot(sample = PCSK9_aorta_pigs,
                       features =  c(  'PECAM1',  'GJA4', 'RGS5'),   ncol = 3,
                       order = F,  pt.size = 0.7, plot.axes = F, use_viridis = T, viridis.palette = 'H')

FeaturePlot(athero_comb_smc, features = c('ACTA2', 'LUM'),reduction = 'umap',pt.size = 0.7, order = T,  blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0.01) &DarkTheme() &NoAxes() & plot_annotation(title = 'Carotid HS')
FeaturePlot(pig_plaque, features = c('ACTA2', 'LUM'),reduction = 'umap',pt.size = 0.7, order = T,  blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0.01) &DarkTheme() &NoAxes() & plot_annotation(title = 'Aorta SS')
FeaturePlot(smc_plaque, features = c('Acta2', 'Lum'),reduction = 'umap',pt.size = 0.7, order = T,  blend = T, cols = c("darkblue", "green", "magenta"), blend.threshold = 0.01) &DarkTheme() &NoAxes() & plot_annotation(title = 'AoRoot MM')

# ============================================================================
# Additional marker identification (exploratory)
# ============================================================================

Idents(athero_comb_smc) <- 'smc_phenotypes'
tmp <- FindAllMarkers(athero_comb_smc, assay = 'RNA', only.pos = T)

Idents(pig_plaque) <- 'smc_phenotypes'
tmp2 <- FindAllMarkers(pig_plaque, assay = 'RNA', only.pos = T)

# --- Pericyte marker overlap across datasets ---
# Top 100 pericyte markers per dataset, then intersect
h.markers.filtered <- h.markers[h.markers$cluster == "Pericytes" & h.markers$p_val_adj < 0.05, ][order(h.markers$avg_log2FC[h.markers$cluster == "Pericytes" & h.markers$p_val_adj < 0.05], decreasing = TRUE), ][1:100,]
p.markers.filtered  <- p.markers[p.markers$cluster == "Pericytes" & p.markers$p_val_adj < 0.05, ][order(p.markers$avg_log2FC[p.markers$cluster == "Pericytes" & p.markers$p_val_adj < 0.05], decreasing = TRUE), ][1:100,]
hc.markers.filtered  <- hc.markers[hc.markers$cluster == "Pericytes" & hc.markers$p_val_adj < 0.05, ][order(hc.markers$avg_log2FC[hc.markers$cluster == "Pericytes" & hc.markers$p_val_adj < 0.05], decreasing = TRUE), ][1:100,]
pc.markers.filtered  <- pc.markers[pc.markers$cluster == "Pericytes" & pc.markers$p_val_adj < 0.05, ][order(pc.markers$avg_log2FC[pc.markers$cluster == "Pericytes" & pc.markers$p_val_adj < 0.05], decreasing = TRUE), ][1:100,]

xxx <- intersect(intersect(intersect(pc.markers.filtered$gene, h.markers.filtered$gene) , hc.markers.filtered$gene), p.markers.filtered$gene)
