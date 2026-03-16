# ============================================================================
# Bulk RNA-seq: Visualize marker gene expression across species
# ============================================================================
# Input:  xspecies_bulk_rnaseq.rds (from bulk-rnaseq_2_deseq-analysis.R)
# Output: species_specific_gene_expression_bulk.pdf
# ============================================================================

setwd("/faststorage/project/THOR/diana/integration_pigs_with_WT/integrated_pig_mice_hum/closer/more_closer/SMC/cross-species-analysis/6_bulk")

library(DESeq2)
library(tidyverse)
library(ggsci)

# --- Load DESeq2 results ---
my_data <- readRDS("xspecies_bulk_rnaseq.rds")

# Uppercase mouse gene names to match human/pig nomenclature
rownames(my_data$mouse$vst_norm_count) <- toupper(rownames(my_data$mouse$vst_norm_count))
rownames(my_data$mouse$norm_count) <- toupper(rownames(my_data$mouse$norm_count))

# --- Exploratory plot: immune marker genes (VST-normalized) ---
my_genes <- c("CD68", "TYROBP", "CD8A", "SUCNR1")

gene_exprs <- imap_dfr(my_data, ~ {
  t(.x$vst_norm_count[toupper(rownames(.x$vst_norm_count)) %in% my_genes, ]) %>%
    as.data.frame() %>%
    mutate(species = .y) %>%
    cbind(colData(my_data[[.y]]$deseq_object)[, c("sample", "group")])
})

# Harmonize group labels across species
gene_exprs$group <- fct_recode(
  .f = gene_exprs$group,
  Unstable = "unstable",
  Stable = "stable",
  Plaque = "hfd3m",
  Healthy_artery = "chow",
  Healthy_artery = "Normal_chow_diet",
  Plaque = "High_fat_high_cholesterol_and_high_fructose_diet"
)

gene_exprs %>%
  pivot_longer(cols = -c(sample, group, species), names_to = "gene", values_to = "exprs") %>%
  ggplot(aes(x = group, y = exprs, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  facet_grid(rows = vars(gene), cols = vars(species), drop = TRUE, scales = "free") +
  theme_bw()

# --- Final figure: SMC marker genes in plaque samples (normalized counts) ---
my_genes <- factor(c(
  "UBC", "ALPL", "COL2A1", "DLX5", "SUCNR1",
  "GJA4", "PGF", "SPARCL1", "ACTA2"
))

my_colors <- c(Plaque = "#E64B35FF", Healthy_artery = "#4DBBD5FF")

gene_exprs <- imap_dfr(my_data, ~ {
  t(.x$norm_count[toupper(rownames(.x$norm_count)) %in% my_genes, ]) %>%
    as.data.frame() %>%
    mutate(species = .y) %>%
    cbind(colData(my_data[[.y]]$deseq_object)[, c("sample", "group")])
})

# Harmonize group labels and keep only plaque samples
gene_exprs$group <- fct_recode(
  .f = gene_exprs$group,
  Plaque = "unstable",
  Plaque = "stable",
  Plaque = "hfd3m",
  Healthy_artery = "chow",
  Healthy_artery = "Normal_chow_diet",
  Plaque = "High_fat_high_cholesterol_and_high_fructose_diet"
)
gene_exprs <- gene_exprs[gene_exprs$group == "Plaque", ]

p1 <- gene_exprs %>%
  pivot_longer(
    cols = -c(sample, group, species),
    names_to = "gene",
    values_to = "exprs"
  ) %>%
  mutate(gene = factor(gene, levels = my_genes)) %>%
  ggplot(aes(x = group, y = exprs, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) +
  facet_grid(
    rows = vars(gene),
    cols = vars(species),
    drop = TRUE,
    scales = "free"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "Normalized counts") +
  scale_fill_manual(values = my_colors) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

print(p1)

ggsave(
  filename = "species_specific_gene_expression_bulk.pdf",
  plot = p1,
  device = "pdf",
  width = 2.7,
  height = 8,
  units = "cm",
  scale = 3.2
)
