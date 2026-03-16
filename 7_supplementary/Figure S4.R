# ============================================================================
# Figure S4: Species-specific SMC phenotype marker FeaturePlots
# ============================================================================
# Generates FeaturePlots for phenotype-defining markers across all 5 datasets
# (mouse elastic, human elastic, human coronary, pig elastic, pig coronary).
# Mouse uses lowercase gene symbols; pig/human use uppercase.
#
# Input:
#   - /faststorage/.../mouse/smc_plaque_Alencar_seurat.rds
#   - /faststorage/.../human/athero_comb_smc_seurat.rds
#   - /faststorage/.../SMC/coronaries/human_coronaries_seurat.rds
#   - /faststorage/.../Pigs/.../pig_plaque_seurat.rds
#   - /faststorage/.../SMC/coronaries/pig_coronaries_seurat.rds
#
# Output:
#   - 1–30 species-specific_markers.pdf in Figure_S4/ subfolder
# ============================================================================

library(Seurat)
library(ggplot2)
library(RColorBrewer)

# --- Load Seurat objects ---
mouse     <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/mouse/smc_plaque_Alencar_seurat.rds")
human     <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/human/athero_comb_smc_seurat.rds")
pig       <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/Pigs/Pigs_regression/Carlos_new_analysis/SC.Analysis.V4.ThirdPass.RNA.New/SMC/bengal/pig_plaque_seurat.rds")
human_cor <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/human_coronaries_seurat.rds")
pig_cor   <- readRDS("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/coronaries/pig_coronaries_seurat.rds")

DefaultAssay(mouse)     <- 'RNA'
DefaultAssay(human)     <- 'RNA'
DefaultAssay(pig)       <- 'RNA'
DefaultAssay(human_cor) <- 'RNA'
DefaultAssay(pig_cor)   <- 'RNA'

# --- Define marker gene sets per phenotype ---
# Mouse genes are lowercase; pig/human are uppercase
marker_sets <- list(
  Chondromyocyte   = list(mouse = c('Alpl', 'Col2a1', 'Pth1r'),
                          upper = c('ALPL', 'COL2A1', 'PTH1R')),
  DLX5_SMC         = list(mouse = c('Dlx5', 'Sost', 'Ptn'),
                          upper = c('DLX5', 'SOST', 'PTN')),
  Pericyte         = list(mouse = c('Gja4', 'Pgf', 'Rgs16'),
                          upper = c('GJA4', 'PGF', 'RGS16')),
  Contractile      = list(mouse = c('Pln', 'Filip1l', 'Sparcl1'),
                          upper = c('PLN', 'FILIP1L', 'SPARCL1')),
  LBH_Pericyte     = list(mouse = c('Rergl', 'Fabp4', 'Ptp4a3'),
                          upper = c('RERGL', 'FABP4', 'PTP4A3')),
  VEGFA_Terminal   = list(mouse = c('Loxl2', 'Vegfa', 'Egln3'),
                          upper = c('LOXL2', 'VEGFA', 'EGLN3'))
)

# Ordered list of datasets: mouse uses lowercase, rest use uppercase
datasets <- list(
  list(obj = mouse,     use_mouse = TRUE),
  list(obj = human,     use_mouse = FALSE),
  list(obj = human_cor, use_mouse = FALSE),
  list(obj = pig,       use_mouse = FALSE),
  list(obj = pig_cor,   use_mouse = FALSE)
)

# --- Helper: create a single FeaturePlot with consistent styling ---
make_feature_plot <- function(seurat_obj, genes) {
  FeaturePlot(seurat_obj,
              features = genes, reduction = 'umap',
              ncol = 3, order = FALSE, pt.size = 0.5) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) &
    theme(aspect.ratio = 1) &
    NoAxes()
}

# --- Generate all 30 plots (6 phenotypes x 5 datasets) ---
output_dir <- "/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/after_bengal/cluster_specificity/Figure_S4/"

plot_idx <- 1
for (pheno in names(marker_sets)) {
  for (ds in datasets) {
    genes <- if (ds$use_mouse) marker_sets[[pheno]]$mouse else marker_sets[[pheno]]$upper
    p <- make_feature_plot(ds$obj, genes)
    ggsave(filename = paste0(plot_idx, ' species-specific_markers.pdf'),
           plot = p, path = output_dir,
           units = 'cm', width = 14, height = 14, bg = 'white', scale = 1.5)
    plot_idx <- plot_idx + 1
  }
}
