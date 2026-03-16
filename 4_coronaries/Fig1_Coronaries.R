# ============================================================================
# Figure: Coronary artery SMC phenotype dotplots and UMAP dimplots
# Individual species (human Wirka + pig RCA) and combined cross-species dotplot
# ============================================================================

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries")
set.seed(1984)

# --- Libraries ---
library(Seurat)
library(tidyverse)
library(patchwork)
library(scCustomize)
library(Matrix)
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
library(ggh4x)
library(scales)

# --- Common theme ---
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

# --- Load coronary SMC Seurat objects ---
human <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/human_coronaries_seurat.rds")
pig <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/pig_coronaries_seurat.rds")


# ============================================================================
# Human coronary SMCs — dotplot and UMAP
# ============================================================================

# Reorder phenotype levels for plotting
Idents(human) <- 'smc_phenotypes'
levels(Idents(human))
human$smc_phenotypes <- factor(human$smc_phenotypes, levels = levels(human$smc_phenotypes)[c(4,3,6,2,5,1)])
Idents(human) <- 'smc_phenotypes'

# Find top markers per phenotype
markers_0.3 <- FindAllMarkers(human, assay = 'RNA', logfc.threshold = 0.25,  min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 5))
}

# Manual marker replacements for display
top_markers[2] <- 'ACTA2'
top_markers[14] <- 'LUM'
top_markers[16] <- 'TINAGL1'
top_markers[17] <- 'NOTCH3'
# top_markers[18] <- 'RGS5'
top_markers[19] <- 'RGS16'
top_markers[20] <- 'GJA4'

# Dotplot
p1 <- DotPlot(
  human,
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
ggsave(filename = 'dotplot_human_coronaries_smc_phenotypes.pdf', plot = p1, device = 'pdf', width = 8, height = 3, units = 'cm', scale = 3.3)

# UMAP dimplot
pal_celltypes_human <- c( "#FFDC91FF","#20854EFF","#E18727FF","#0072B5FF",   "#374E55FF" ,"#BC3C29FF" )[c(4,3,6,2,5,1)]

p2 <- DimPlot_scCustom(seurat_object = human,
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
ggsave(filename = 'dimplot_human_coronaries_smc_phenotypes.pdf', plot = p2, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 2.5)


# ============================================================================
# Pig coronary (RCA) SMCs — dotplot and UMAP
# ============================================================================

# Reorder phenotype levels for plotting
levels(Idents(pig))
pig$smc_phenotypes <- factor(pig$smc_phenotypes, levels = levels(pig$smc_phenotypes)[c(3,4,2,1,5,8,7,6)])
Idents(pig) <- 'smc_phenotypes'

# Find top markers per phenotype
markers_0.3 <- FindAllMarkers(pig, assay = 'RNA', logfc.threshold = 0.25, min.pct = 0.3, only.pos = T, max.cells.per.ident = 1000)
markers_0.3 <- markers_0.3[markers_0.3$p_val_adj < 0.05,]

top_markers <- c()
for (i in unique(markers_0.3$cluster)){
  top_markers <- c(top_markers, head((subset(markers_0.3, subset = markers_0.3$cluster == i)[order(subset(markers_0.3, subset = markers_0.3$cluster == i)$avg_log2FC, decreasing = T),])$gene, 4))
}

# Manual marker replacements for display
top_markers[3] <- 'ACTA2'
top_markers[16] <- 'LUM'
top_markers[24] <- 'RGS16'
top_markers[23] <- 'GJA4'
top_markers[21] <- 'TINAGL1'

# Dotplot
p3 <- DotPlot(
  pig,
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
ggsave(filename = 'dotplot_pig_smc_phenotypes.pdf', plot = p1, device = 'pdf', width = 8, height = 4.5, units = 'cm', scale = 3.2)

# UMAP dimplot
pal_celltypes_pig <- c("#0072B5FF","#E18727FF", "coral1",  "#BC3C29FF","#7876B1FF" ,"#20854EFF", "#EE4C97FF" ,"#FFDC91FF")

p4 <- DimPlot_scCustom(seurat_object = pig,
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
ggsave(filename = 'dimplot_pig_smc_phenotypes_coronary.pdf', plot = p4, device = 'pdf', width = 4, height = 4, units = 'cm', scale = 2.5)


# ============================================================================
# Combined cross-species dotplot (human + pig coronaries)
# ============================================================================

# Merge dotplot data from both species
df135 <- rbind(mutate(p1$data, organism = "human"),
              mutate(p3$data, organism = "pig"))

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
  scale_fill_gradientn(colors = rev(brewer.pal(11, "Spectral"))) +
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

# --- Diverging color scale centered at zero ---
# mypal <- met.brewer("OKeeffe1", n = 11)
mypal <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")

# Attempt with apply_grad (defined in elastic arteries script)
pw135 <- (apply_grad(p1) + apply_grad(p3) + apply_grad(p5)) +
  patchwork::plot_layout(ncol = 1, heights = c(1,1,1))

# Manual per-panel rescaling
pw135 <- (p3 + scale_color_gradientn(
  colors = colorRampPalette(rev(mypal))(101),
  values = scales::rescale(c(seq(min(p3$data$avg.exp.scaled), 0, length.out = 51),
                     seq(0, max(p3$data$avg.exp.scaled), length.out = 51)[-1]))) +
    labs(x = NULL, y = NULL, title = NULL)) +
  (p1 + scale_color_gradientn(
    colors = colorRampPalette(rev(mypal))(101),
    values = scales::rescale(c(seq(min(p1$data$avg.exp.scaled), 0, length.out = 51),
                       seq(0, max(p1$data$avg.exp.scaled), length.out = 51)[-1]))) +
     labs(x = NULL, y = NULL, title = NULL)) +
  plot_layout(ncol = 1, heights = nids/2, widths = rep(1, 2))
pw135
ggsave(filename = 'dotplot_all_species_coronaries.2.pdf',
       plot = pw135, width = 9, height = 12, units = 'cm', scale = 2.5)
