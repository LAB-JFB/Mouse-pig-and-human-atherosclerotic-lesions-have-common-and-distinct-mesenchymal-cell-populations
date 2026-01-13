MarkerPairEval <- function(genes, labels, expression_matrix, 
                           combination = c("AND", "OR", "sum", "product"),
                           expression_threshold = 0) {
  
  combination <- match.arg(combination)
  
  # Check for NA's in labels
  if(any(is.na(labels))) {
    sel_na <- which(is.na(labels))
    labels <- labels[-sel_na]
    expression_matrix <- expression_matrix[, -sel_na]
  }
  labels <- as.logical(labels)
  
  # Limit the number of genes
  if(length(genes) > 100) { genes <- genes[1:100] }
  
  # Filter to genes present in matrix
  genes <- genes[genes %in% rownames(expression_matrix)]
  
  # Generate all gene pairs
  gp <- combn(genes, 2)
  n_pairs <- ncol(gp)
  
  # Pre-allocate results
  metrics <- data.frame(
    Gene.1 = gp[1, ],
    Gene.2 = gp[2, ],
    AUC = numeric(n_pairs),
    Sensitivity = numeric(n_pairs),
    Specificity = numeric(n_pairs),
    Accuracy = numeric(n_pairs),
    Precision = numeric(n_pairs),
    F1.score = numeric(n_pairs),
    stringsAsFactors = FALSE
  )
  
  # Calculate metrics for each pair
  for(i in seq_len(n_pairs)) {
    expr1 <- expression_matrix[gp[1, i], ]
    expr2 <- expression_matrix[gp[2, i], ]
    
    # Create combined score based on combination method
    combined_score <- switch(combination,
                             "AND" = pmin(as.numeric(expr1), as.numeric(expr2)), # AND (pmin): both genes need high expression — good for co-expression markers; OR (pmax): either gene high is sufficient — good for redundant markers; sum: additive effect — treats genes as contributing equally; product: multiplicative — emphasizes cells where both are expressed
                             "OR" = pmax(as.numeric(expr1), as.numeric(expr2)),
                             "sum" = as.numeric(expr1) + as.numeric(expr2),
                             "product" = as.numeric(expr1) * as.numeric(expr2)
    )
    
    # Calculate true ROC AUC
    roc_obj <- pROC::roc(labels, combined_score, quiet = TRUE, direction = "<")
    metrics$AUC[i] <- as.numeric(pROC::auc(roc_obj))
    
    # Get optimal threshold metrics (Youden's J)
    best_coords <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "all")
    
    # If multiple optimal points, take the first
    if(is.data.frame(best_coords) && nrow(best_coords) > 1) {
      best_coords <- best_coords[1, ]
    }
    
    metrics$Sensitivity[i] <- best_coords$sensitivity
    metrics$Specificity[i] <- best_coords$specificity
    
    # Calculate binary metrics at optimal threshold
    threshold <- best_coords$threshold
    predicted <- combined_score >= threshold
    
    tp <- sum(predicted & labels)
    tn <- sum(!predicted & !labels)
    fp <- sum(predicted & !labels)
    fn <- sum(!predicted & labels)
    
    metrics$Accuracy[i] <- (tp + tn) / (tp + tn + fp + fn)
    metrics$Precision[i] <- ifelse((tp + fp) > 0, tp / (tp + fp), 0)
    metrics$F1.score[i] <- ifelse((2*tp + fp + fn) > 0, (2*tp) / (2*tp + fp + fn), 0)
  }
  
  # Add predictive power
  metrics$Predictive.power <- abs(metrics$AUC - 0.5) * 2
  
  # Sort by AUC descending
  metrics <- metrics[order(metrics$AUC, decreasing = TRUE), ]
  
  return(metrics)
}

my_res <- list()

my_ds <- list(h_car = human_carotid,
              h_cor = human_coronary,
              p_ao = pig_aorta, 
              p_cor = pig_coronary,
              mouse = mouse)


for(k in names(my_ds)){
  seurat_obj <- my_ds[[k]]
  
  tmp <- FindMarkers(seurat_obj, ident.1 = 'Mesenchymal', only.pos = T, assay = 'RNA',group.by = 'superduperclusters' )
  tmp <- tmp[tmp$avg_log2FC > 0.6 & tmp$pct.1 > 0.3,]
  tmp <- tmp[order(tmp$avg_log2FC, decreasing = T),]
  tmp$gene <- rownames(tmp)
  my_data <-
    seurat_obj[, seurat_obj$superduperclusters == 'Mesenchymal']
  
  tmp[, 'exp'] <-
    matrixStats::rowMeans2(as.matrix(expm1(my_data[["RNA"]]@data[tmp$gene,])))
  tmp <- tmp[tmp$exp > 2, ]
  
  
  Idents(seurat_obj) <- 'superduperclusters'
  
  group1_cells <- WhichCells(seurat_obj, idents = "Mesenchymal")
  group2_cells <- colnames(seurat_obj[,!colnames(seurat_obj) %in% group1_cells])
  
  # Fetch data
  data_use <- GetAssayData(seurat_obj, slot = "data")
  genes <- unique(c("LUM",rownames(tmp)))
  
  best_genes <- MarkerPairEval(
    genes = genes,
    labels = seurat_obj$superduperclusters == "Mesenchymal",
    expression_matrix = seurat_obj@assays$RNA@data, expression_threshold = 0, combination = "sum"
  )
  
  
  
  for(i in 1:nrow(best_genes)){
    
    # Fetch data for the two genes of interest
    data_target <- FetchData(seurat_obj, vars = c(best_genes[i, 'Gene.1'], best_genes[i, 'Gene.2']), cells = group1_cells)
    data_nontarget <- FetchData(seurat_obj, vars = c(best_genes[i, 'Gene.1'], best_genes[i, 'Gene.2']), cells = group2_cells)
    
    # Determine cells where both genes are expressed
    # Assuming expression means non-zero counts
    coexpressed_cells_t <- which(data_target[[best_genes[i, 'Gene.1']]] > 0 & data_target[[best_genes[i, 'Gene.2']]] > 0)
    coexpressed_cells_nt <- which(data_nontarget[[best_genes[i, 'Gene.1']]] > 0 & data_nontarget[[best_genes[i, 'Gene.2']]] > 0)
    
    # Calculate percentage
    best_genes[i, 'pct.Target'] <- length(coexpressed_cells_t) / length(group1_cells) * 100
    best_genes[i, 'pct.NonTarget'] <- length(coexpressed_cells_nt) / length(group2_cells) * 100
  }
  
  tmp <-
    FindMarkers(
      seurat_obj,
      ident.1 = 'Mesenchymal',
      only.pos = T,
      assay = 'RNA',
      group.by = 'superduperclusters',
      test.use = 'roc'
    )
  
  
  best_genes[, 'pct_gene1'] <- tmp[best_genes[, 'Gene.1'], 'pct.1']
  best_genes[, 'pct_gene2'] <- tmp[best_genes[, 'Gene.2'], 'pct.1']
  
  
  
  my_res[[k]] <- best_genes
  
}
