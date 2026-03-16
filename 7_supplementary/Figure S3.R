# ============================================================================
# Figure S3: Carramolino et al. SMC phenotype validation
# ============================================================================
# Reproduces FeaturePlots and DimPlot from Carramolino et al. dataset
# to validate DLX5-associated SMC markers.
#
# Input:
#   - /faststorage/project/THOR/diana/integration_pigs_with_WT/hdfc/
#     SC.Analysis.V4.FirstPass.RNA.W4.MT15NbW034.W4.SMCs.Final.rds
#
# Output (saved to working directory):
#   - FS3_Carramolinoetal_pheno.pdf   — DimPlot of clusters
#   - FS3_Carramolinoetal_markers.pdf — FeaturePlot of DLX5-associated markers
# ============================================================================

library(Seurat)
library(ggplot2)
library(RColorBrewer)

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity")

# --- Load Carramolino et al. Seurat object ---
carr <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/hdfc/SC.Analysis.V4.FirstPass.RNA.W4.MT15NbW034.W4.SMCs.Final.rds")

# --- FeaturePlot: DLX5-associated markers ---
p1 <- FeaturePlot(carr,
                  features = c('DLX5', 'SOST', 'PTN', 'SUCNR1', 'LMO2'),
                  reduction = 'umap',
                  ncol = 4, order = FALSE, pt.size = 0.5) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) &
  theme(aspect.ratio = 1) &
  NoAxes()

# --- DimPlot: cluster identity ---
p3 <- DimPlot(carr,
              group.by = 'RNA_snn_res.0.34',
              reduction = 'umap',
              cols = brewer.pal(n = 12, name = 'Paired')[1:5],
              label = TRUE, pt.size = 1, label.box = TRUE) +
  theme(plot.title = element_blank()) +
  NoLegend() + NoAxes() + theme(aspect.ratio = 1)

# --- Save ---
ggsave(filename = 'FS3_Carramolinoetal_pheno.pdf', plot = p3)
ggsave(filename = 'FS3_Carramolinoetal_markers.pdf', plot = p1, scale = 2)
