# Cross-Species SMC Analysis in Atherosclerotic Plaques

This repository contains the analysis code for cross-species comparison of smooth muscle cells (SMCs) in atherosclerotic plaques across pig, mouse, and human. The analysis covers both elastic (aorta/carotid) and muscular (coronary) arteries.

All scripts were originally run on an HPC cluster with absolute paths. Input data files are not included in this repository due to size. See the [Input Data](#input-data) section for details on required files.

## Folder Structure

```
cross-species-analysis/
|
|-- 1_raw_data_preparation/       # Preprocessing of individual species datasets
|   |-- scrnaseq_1_cellranger-pipeline.sh     # Cell Ranger pipeline for raw FASTQ processing
|   |-- 01_human_carotid_processing.R         # Process human carotid plaque scRNA-seq (includes raw data loading, QC, DoubletFinder, SCT integration, SMC subsetting, phenotype annotation)
|   |-- 01_mouse_Alencar_processing.R         # Process mouse SMC lineage-tracing data (Alencar et al.) from raw 10X counts
|   |-- 01_pig_plaque_processing.R            # Process pig aorta plaque scRNA-seq (SMC subsetting, pseudotime, phenotype annotation)
|   |-- 01_coronaries_preprocessing.R         # Preprocess coronary artery data (pig RCA + human Wirka): SMC subsetting, clustering, label transfer, phenotype annotation
|   |-- 02_seurat_to_anndata_conversion.Rmd   # Convert Seurat objects to .h5ad for BENGAL input (Python/scanpy)
|   |-- sample_info.txt                       # Sample metadata for human carotid datasets
|   '-- annotated_query.csv                   # Cell type annotation reference
|
|-- 2_bengal/                     # Cross-species integration with BENGAL
|   |-- BENGAL/                   # Only our run-specific files; pipeline code is external
|   |   |-- config/               # Run configuration files
|   |   |   |-- mph_lessp_mice_2.config       # Config for elastic arteries (3 species)
|   |   |   '-- coronaries.config             # Config for coronary arteries (2 species)
|   |   |-- mph_metadata_nf.tsv              # Input metadata for elastic artery run
|   |   '-- mph_metadata_coronaries.tsv      # Input metadata for coronary artery run
|   |-- elastic/results/          # BENGAL output for elastic arteries (not in repo)
|   '-- coronary/results/         # BENGAL output for coronary arteries (not in repo)
|
|-- 3_elastic_arteries/           # Post-integration analysis of elastic arteries (aorta/carotid)
|   |-- cluster_specificity.Rmd   # Cluster specificity across integration methods and gene homology
|   |-- Fig1_elastic_arteries.R   # Figure 1: individual species dotplots and dimplots
|   |-- Fig2_integration_overview.R  # Figure 2: integrated dimplots, cluster specificity barplots
|   '-- metrics_manually.R        # Integration quality metrics (kBET, LISI, ARI, NMI, PCR)
|
|-- 4_coronaries/                 # Post-integration analysis of coronary arteries
|   |-- cluster_specificity.Rmd   # Cluster specificity for coronary datasets
|   |-- Fig1_Coronaries.R         # Figure: coronary artery dimplots and dotplots
|   '-- metrics_manually_coronaries.R  # Integration quality metrics for coronaries
|
|-- 5_scenic/                     # SCENIC transcription factor analysis
|   |-- mouse_scenic_smc.Rmd      # Run SCENIC on mouse SMC data
|   |-- SCENIC_human_carotids_part.Rmd  # Run SCENIC on human carotid SMC data
|   '-- scenic_cross_species_final.R    # Cross-species TF regulon comparison (human vs mouse)
|
|-- 6_bulk/                       # Bulk RNA-seq validation across species
|   |-- bulk-rnaseq_1_nfcore-pipeline.sh   # nf-core/rnaseq pipeline for human, mouse, pig
|   |-- bulk-rnaseq_2_deseq-analysis.R     # DESeq2 analysis per species (GEO datasets)
|   |-- bulk-rnaseq_3_visualization.R      # Marker gene expression visualization
|   |-- human_samplesheet.csv              # Human bulk RNA-seq sample metadata
|   |-- mouse_samplesheet.csv              # Mouse bulk RNA-seq sample metadata
|   '-- pig_samplesheet.csv               # Pig bulk RNA-seq sample metadata
|
|-- 7_supplementary/              # Supplementary figure scripts
|   |-- common_pheno_markers.R    # Shared SMC phenotype markers across elastic + coronary
|   |-- Figure S3.R               # Supp. Fig S3: Carramolino et al. external validation
|   '-- Figure S4.R               # Supp. Fig S4: species-specific marker features
|
'-- 8_revision/                   # Revision analyses
    |-- ROC_analysis_cross_species.R   # ROC analysis of ACTA2-LUM marker pair across tissues
    |-- MarkerPairEval_v3.R            # Paired-gene ROC evaluation function (sourced by ROC script)
    |-- SingleMarkerEval_v3.R          # Single-gene ROC evaluation function (sourced by ROC script)
    '-- sex_stratification.R           # Sex-stratified analysis of SMC phenotypes
```

## Execution Order

1. **Raw data preparation** (`1_raw_data_preparation/`): The four `01_` scripts can be run independently (in any order). Then run `02_seurat_to_anndata_conversion.Rmd` to produce the `.h5ad` files needed for BENGAL.
2. **BENGAL integration** (`2_bengal/`): Run the BENGAL Nextflow pipeline using the provided config and metadata files. Produces integrated datasets under `elastic/results/` and `coronary/results/`.
3. **Downstream analysis** (`3_elastic_arteries/`, `4_coronaries/`, `5_scenic/`, `6_bulk/`, `7_supplementary/`, `8_revision/`): These scripts consume the BENGAL output and individual species objects. They can be run in any order within each folder. `6_bulk/` uses independent bulk RNA-seq data processed with nf-core/rnaseq and DESeq2.

## BENGAL Pipeline

The BENGAL cross-species integration tool is available at: **[BENGAL GitHub repository URL]**

Only our run-specific configuration and metadata files are included in this repository. The BENGAL pipeline code, Singularity containers, Nextflow workflows, and all pipeline outputs are excluded via `.gitignore`.

### Files tracked in `2_bengal/BENGAL/`

| File | Description |
|------|-------------|
| `config/mph_lessp_mice_2.config` | Nextflow config for elastic arteries run (pig, mouse, human) |
| `config/coronaries.config` | Nextflow config for coronary arteries run (pig, human) |
| `mph_metadata_nf.tsv` | Input metadata mapping species to `.h5ad` files (elastic arteries) |
| `mph_metadata_coronaries.tsv` | Input metadata mapping species to `.h5ad` files (coronary arteries) |
| `homology_tbl_bengal_sscrofa.rds` | Homology table for pig (Sus scrofa) gene ortholog mapping |
| `mart_bengal_sscrofa.rds` | BioMart query results for pig gene annotations |
| `avail_attributes_bengal_sscrofa.rds` | Available BioMart attributes for pig |

### To reproduce

1. Clone the [BENGAL repository] https://github.com/Papatheodorou-Group/BENGAL
2. Copy the config files into the BENGAL `config/` directory
3. Copy the metadata `.tsv` files into the BENGAL root directory
4. Adjust absolute paths in the config files to match your environment


## Input Data

All input files are Seurat `.rds` objects or AnnData `.h5ad` files. Scripts reference these via absolute paths on our HPC system. Users will need to adjust paths accordingly.

### Step 1 input data

| File | Used by | Description |
|------|---------|-------------|
| Raw 10X directories (Alsaigh, Pan, Wirka) | `01_human_carotid_processing.R` | Human carotid plaque raw scRNA-seq counts |
| `hpca_se.rds` | `01_human_carotid_processing.R` | Human Primary Cell Atlas (SingleR reference) |
| Raw 10X directories (8 samples: 85–97) | `01_mouse_Alencar_processing.R` | Mouse Alencar lineage-tracing raw scRNA-seq counts |
| `SC.Analysis.V4.ThirdPass.RNA.Aorta.All.seuratSC.Final.rds` | `01_pig_plaque_processing.R` | Full pig aorta Seurat object |
| `coronaries_carotids_processed.rds` | `01_coronaries_preprocessing.R` | Human coronary + carotid combined Seurat object |
| `SC.Analysis.V4.ThirdPass.RNA.RCA.All.seuratSC.Final.rds` | `01_coronaries_preprocessing.R` | Full pig RCA (coronary) Seurat object |
| `pal_celltypes.rds` | `01_human_carotid_processing.R`, `01_coronaries_preprocessing.R` | Cell type color palette |

### Step 1 output data

| File | Produced by | Description |
|------|-------------|-------------|
| `athero_comb_smc_seurat.rds` | `01_human_carotid_processing.R` | Human carotid SMCs (final, with phenotype annotation) |
| `smc_plaque_Alencar_seurat.rds` | `01_mouse_Alencar_processing.R` | Mouse aorta SMCs (final, with phenotype annotation) |
| `pig_plaque_seurat.rds` | `01_pig_plaque_processing.R` | Pig aorta SMCs (final, with phenotype annotation) |
| `SMC.aorta.processed.RNA.norm.rds` | `01_pig_plaque_processing.R` | Pig aorta SMC normalized Seurat object |
| `human_coronaries_seurat.rds` | `01_coronaries_preprocessing.R` | Human coronary SMCs (final, with phenotype annotation) |
| `pig_coronaries_seurat.rds` | `01_coronaries_preprocessing.R` | Pig coronary SMCs (final, with phenotype annotation) |
| `SMC.RCA.processed.rds` | `01_coronaries_preprocessing.R` | Pig RCA SMC processed Seurat object |

### Step 2 input data (AnnData conversion)

| File | Used by | Description |
|------|---------|-------------|
| `sscrofa_i.h5ad` | `02_seurat_to_anndata_conversion.Rmd` | Pig SMC intermediate .h5ad (converted from Seurat) |
| `athero_comb_i.h5ad` | `02_seurat_to_anndata_conversion.Rmd` | Human SMC intermediate .h5ad (converted from Seurat) |
| `mice_hfdc_i.h5ad` | `02_seurat_to_anndata_conversion.Rmd` | Mouse SMC intermediate .h5ad (converted from Seurat) |

### Step 2 output data (BENGAL-ready)

| File | Produced by | Description |
|------|-------------|-------------|
| `sscrofa.h5ad` | `02_seurat_to_anndata_conversion.Rmd` | Pig SMC .h5ad with QC metrics (BENGAL input) |
| `hsapiens.h5ad` | `02_seurat_to_anndata_conversion.Rmd` | Human SMC .h5ad with QC metrics (BENGAL input) |
| `mmusculus.h5ad` | `02_seurat_to_anndata_conversion.Rmd` | Mouse SMC .h5ad with QC metrics (BENGAL input) |

### Processed species-level objects (produced by Step 1, used in Steps 3-8)

The following `.rds` files are outputs of the `01_*` scripts and serve as inputs for downstream analysis:

| File | Produced by | Used by | Description |
|------|-------------|---------|-------------|
| `athero_comb_smc_seurat.rds` | `01_human_carotid_processing.R` | `Fig1_elastic_arteries.R`, `Figure S4.R` | Human carotid SMCs (processed) |
| `smc_plaque_Alencar_seurat.rds` | `01_mouse_Alencar_processing.R` | `Fig1_elastic_arteries.R`, `Figure S4.R` | Mouse aorta SMCs (processed) |
| `pig_plaque_seurat.rds` | `01_pig_plaque_processing.R` | `Fig1_elastic_arteries.R`, `Figure S4.R` | Pig aorta SMCs (processed) |
| `human_coronaries_seurat.rds` | `01_coronaries_preprocessing.R` | `Fig1_Coronaries.R`, `Figure S4.R` | Human coronary SMCs |
| `pig_coronaries_seurat.rds` | `01_coronaries_preprocessing.R` | `Fig1_Coronaries.R`, `Figure S4.R` | Pig coronary SMCs |



### Post-BENGAL integrated objects (used in Steps 3-8)

| File | Used by | Description |
|------|---------|-------------|
| `r_py_all_datasets_processed_with_specific_clusters.rds` | `Fig2_integration_overview.R`, `common_pheno_markers.R` | All elastic artery integration results with cluster specificity |
| `r_py_all_datasets_processed_with_specific_clusters_coronary_2.rds` | `Fig1_Coronaries.R` | All coronary integration results with cluster specificity |
| `all_datasets_ALL.rds` (coronaries) | `metrics_manually_coronaries.R` | All coronary integration datasets |
| `merged_pig_human_smcs_new.rds` | `Fig1_Coronaries.R` | Merged pig-human coronary SMCs |
| `merged_pig_human_mouse_smcs.rds` | `Fig2_integration_overview.R` | Merged 3-species elastic SMCs |
| `HE_RPCA_cell_categorization.rds` | `Fig2_integration_overview.R` | Elastic artery cell categorization |
| `shared_Cor_HE_RPCA_cell_categorization.rds` | `Fig1_Coronaries.R` | Coronary cell categorization |
| `specific_clusters_df_10_80.rds` | `Fig2_integration_overview.R` | Cluster specificity scores (elastic) |
| `specific_clusters_df_90_10.rds` | `Fig1_Coronaries.R` | Cluster specificity scores (coronary) |
| `cluster_proportions_list_plots.rds` | `Fig2_integration_overview.R` | Pre-computed cluster proportion plots |
| `correlation_shared_SMC_phenotypes_LIST.rds` | `Fig2_integration_overview.R` | Correlation of shared SMC phenotypes |
| `pal_cat.rds`, `pal_celltypes.rds`, `pal_species.rds` | `Fig2_integration_overview.R` | Color palettes |

### Integration metrics (used in Steps 3-4 metrics scripts)

| File | Used by | Description |
|------|---------|-------------|
| `kBET_rejection_rates_df.rds` | `metrics_manually*.R` | kBET batch effect scores |
| `LISI_scores_per_cell.rds` | `metrics_manually*.R` | LISI integration scores |
| `pcr_list.rds` | `metrics_manually*.R` | Principal component regression results |
| `ARI_smc_phenotypes.rds` | `metrics_manually*.R` | Adjusted Rand Index |
| `NMI_smc_phenotypes.rds` | `metrics_manually*.R` | Normalized Mutual Information |

### SCENIC input (Step 5)

| File | Used by | Description |
|------|---------|-------------|
| `4.2_binaryRegulonActivity_nonDupl.Rds` (human) | `scenic_cross_species_final.R` | Human SCENIC binary regulon activity |
| `4.2_binaryRegulonActivity_nonDupl.Rds` (mouse) | `scenic_cross_species_final.R` | Mouse SCENIC binary regulon activity |
| `3.4_regulonAUC.Rds` (human) | `mouse_scenic_smc.Rmd` | Human SCENIC regulon AUC scores |
| `3.4_regulonAUC.Rds` (mouse) | `mouse_scenic_smc.Rmd` | Mouse SCENIC regulon AUC scores |

### Bulk RNA-seq data (Step 6)

| File | Used by | Description |
|------|---------|-------------|
| `xspecies_bulk_rnaseq.rds` | `bulk-rnaseq_3_visualization.R` | Combined DESeq2 results for all three species |
| `human_samplesheet.csv` | `bulk-rnaseq_1_nfcore-pipeline.sh` | Human bulk RNA-seq sample sheet (GEO: GSE226790) |
| `mouse_samplesheet.csv` | `bulk-rnaseq_1_nfcore-pipeline.sh` | Mouse bulk RNA-seq sample sheet (GEO: GSE205929) |
| `pig_samplesheet.csv` | `bulk-rnaseq_1_nfcore-pipeline.sh` | Pig bulk RNA-seq sample sheet (GEO: GSE199592) |

### Revision data (Step 8)

| File | Used by | Description |
|------|---------|-------------|
| `coronaries_carotids_processed_2.rds` | `ROC_analysis_cross_species.R` | Human coronary + carotid (updated) |
| `biomart_109_orth_new.rds` | `sex_stratification.R` | BioMart orthology table (Ensembl 109) |

## Software Requirements

### R packages
Seurat (v4+), sctransform, scCustomize, patchwork, ggplot2, tidyverse, cowplot, MetBrewer, RColorBrewer, pheatmap, corrplot, ggpubr, MAST, DESeq2, openxlsx, psych, data.table, ggsci, pROC, extrafont, ggrastr, speckle, limma, mclust, tximport, PCAtools, khroma, GEOquery

### Python packages
scanpy, numpy (used in `02_seurat_to_anndata_conversion.Rmd`)

### External tools
- [BENGAL](https://github.com/TODO/BENGAL) — cross-species scRNA-seq integration pipeline (Nextflow)
- [SCENIC](https://scenic.aertslab.org/) — single-cell regulatory network inference
- [nf-core/rnaseq](https://nf-co.re/rnaseq) v3.14.0 — bulk RNA-seq processing pipeline (used in `6_bulk/`)
