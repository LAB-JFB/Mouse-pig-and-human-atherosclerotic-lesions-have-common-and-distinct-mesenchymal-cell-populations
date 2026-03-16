### Check the information about samples

library(GEOquery)
gse_list <- c("GSE226790", "GSE205929", "GSE199592")
tmp <- sapply(gse_list, \(x) getGEO(x)[[1]])
View(pData(tmp$GSE226790))


### Analysis using DESeq2

library(tidyverse)
library(tximport)
library(DESeq2)
library(khroma)
library(PCAtools)

xsp <- list("human", "mouse", "pig")
deseq_xsp <- list()

#x = "pig"
for(x in xsp) {
  # make annotation
  tx_info <- read_tsv(
    file.path(x, "star_salmon/tx2gene.tsv"), col_names = FALSE) %>% 
    as.data.frame()
  colnames(tx_info) <- c("tx_id", "gene_id", "gene_name")
  rownames(tx_info) <- tx_info$tx_id
  gene_info <- tx_info %>% 
    distinct(gene_id, gene_name)
  rownames(gene_info) <- gene_info$gene_id 
  # get sample info
  sample_info <- read_csv(file.path(x, "samplesheet.csv"))
  files <- file.path(x,"star_salmon", sample_info$sample, "quant.sf")
  # read quantified counts
  txi <- tximport(files, type = "salmon", tx2gene = tx_info)
  # DESeq2 analysis
  dds <- DESeqDataSetFromTximport(txi, sample_info, ~ group)
  colnames(dds) <- colData(dds)$sample
  rowData(dds) <- DataFrame(gene_info[rownames(dds), ])
  rownames(dds) <- make.unique(rowData(dds)$gene_name)
  dds <- DESeq(dds)
  dds_vsd <- assay(vst(dds))
  dds_nc <- counts(dds, normalized = TRUE)
  res_shrink <- lfcShrink(dds, res = results(dds), type = "ashr")
  deseq_xsp[[x]] <- list(norm_count = dds_nc,
                         vst_norm_count = dds_vsd,
                         dge_results = as.data.frame(res_shrink),
                         deseq_object = dds)
  rm(tx_info, gene_info, sample_info, files, txi, dds, dds_vsd, dds_nc, res_shrink)
}

saveRDS(deseq_xsp, file = "xspecies_bulk_rnaseq.rds")
