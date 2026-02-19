MarkerPairEval <- function(genes, labels, expression_matrix, 
                           combination = c("AND", "OR", "sum", "product"),
                           hard_threshold = 0) {
  
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
    
    # ROC-based metrics (optimal threshold via Youden's J)
    AUC = numeric(n_pairs),
    ROC.Threshold = numeric(n_pairs),
    ROC.Sensitivity = numeric(n_pairs),
    ROC.Specificity = numeric(n_pairs),
    ROC.Accuracy = numeric(n_pairs),
    ROC.Precision = numeric(n_pairs),
    ROC.F1.score = numeric(n_pairs),
    
    # Hard threshold metrics (based on expression > hard_threshold)
    Hard.Sensitivity = numeric(n_pairs),
    Hard.Specificity = numeric(n_pairs),
    Hard.Accuracy = numeric(n_pairs),
    Hard.Precision = numeric(n_pairs),
    Hard.F1.score = numeric(n_pairs),
    
    # Co-expression percentages (both genes > hard_threshold)
    pct.Target = numeric(n_pairs),
    pct.NonTarget = numeric(n_pairs),
    
    # Single gene percentages (each gene > hard_threshold in target)
    pct.Gene1.Target = numeric(n_pairs),
    pct.Gene2.Target = numeric(n_pairs),
    pct.Gene1.NonTarget = numeric(n_pairs),
    pct.Gene2.NonTarget = numeric(n_pairs),
    
    stringsAsFactors = FALSE
  )
  
  # Calculate metrics for each pair
  for(i in seq_len(n_pairs)) {
    expr1 <- as.numeric(expression_matrix[gp[1, i], ])
    expr2 <- as.numeric(expression_matrix[gp[2, i], ])
    
    # Create combined score based on combination method
    combined_score <- switch(combination,
                             "AND" = pmin(expr1, expr2),
                             "OR" = pmax(expr1, expr2),
                             "sum" = expr1 + expr2,
                             "product" = expr1 * expr2
    )
    
    # =========================================
    # ROC-BASED METRICS (Youden's J threshold)
    # =========================================
    
    roc_obj <- pROC::roc(labels, combined_score, quiet = TRUE, direction = "<")
    metrics$AUC[i] <- as.numeric(pROC::auc(roc_obj))
    
    # Get optimal threshold metrics (Youden's J)
    best_coords <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "all")
    
    # If multiple optimal points, take the first
    if(is.data.frame(best_coords) && nrow(best_coords) > 1) {
      best_coords <- best_coords[1, ]
    }
    
    metrics$ROC.Threshold[i] <- best_coords$threshold
    metrics$ROC.Sensitivity[i] <- best_coords$sensitivity
    metrics$ROC.Specificity[i] <- best_coords$specificity
    
    # Binary metrics at ROC optimal threshold
    roc_threshold <- best_coords$threshold
    roc_predicted <- combined_score >= roc_threshold
    
    tp_roc <- sum(roc_predicted & labels)
    tn_roc <- sum(!roc_predicted & !labels)
    fp_roc <- sum(roc_predicted & !labels)
    fn_roc <- sum(!roc_predicted & labels)
    
    metrics$ROC.Accuracy[i] <- (tp_roc + tn_roc) / (tp_roc + tn_roc + fp_roc + fn_roc)
    metrics$ROC.Precision[i] <- ifelse((tp_roc + fp_roc) > 0, tp_roc / (tp_roc + fp_roc), 0)
    metrics$ROC.F1.score[i] <- ifelse((2*tp_roc + fp_roc + fn_roc) > 0, 
                                      (2*tp_roc) / (2*tp_roc + fp_roc + fn_roc), 0)
    
    # =========================================
    # HARD THRESHOLD METRICS (practical/staining-relevant)
    # =========================================
    
    # For "sum" combination: cell is "positive" if EITHER gene > hard_threshold
    # This matches the biological interpretation of sum method
    if(combination == "sum") {
      hard_positive <- (expr1 > hard_threshold) | (expr2 > hard_threshold)
    } else if(combination == "AND") {
      hard_positive <- (expr1 > hard_threshold) & (expr2 > hard_threshold)
    } else if(combination == "OR") {
      hard_positive <- (expr1 > hard_threshold) | (expr2 > hard_threshold)
    } else if(combination == "product") {
      hard_positive <- (expr1 > hard_threshold) & (expr2 > hard_threshold)
    }
    
    tp_hard <- sum(hard_positive & labels)
    tn_hard <- sum(!hard_positive & !labels)
    fp_hard <- sum(hard_positive & !labels)
    fn_hard <- sum(!hard_positive & labels)
    
    metrics$Hard.Sensitivity[i] <- ifelse((tp_hard + fn_hard) > 0, 
                                          tp_hard / (tp_hard + fn_hard), 0)
    metrics$Hard.Specificity[i] <- ifelse((tn_hard + fp_hard) > 0, 
                                          tn_hard / (tn_hard + fp_hard), 0)
    metrics$Hard.Accuracy[i] <- (tp_hard + tn_hard) / (tp_hard + tn_hard + fp_hard + fn_hard)
    metrics$Hard.Precision[i] <- ifelse((tp_hard + fp_hard) > 0, 
                                        tp_hard / (tp_hard + fp_hard), 0)
    metrics$Hard.F1.score[i] <- ifelse((2*tp_hard + fp_hard + fn_hard) > 0, 
                                       (2*tp_hard) / (2*tp_hard + fp_hard + fn_hard), 0)
    
    # =========================================
    # CO-EXPRESSION PERCENTAGES
    # =========================================
    
    # Cells where BOTH genes are expressed (> hard_threshold)
    coexpr <- (expr1 > hard_threshold) & (expr2 > hard_threshold)
    
    n_target <- sum(labels)
    n_nontarget <- sum(!labels)
    
    metrics$pct.Target[i] <- sum(coexpr & labels) / n_target * 100
    metrics$pct.NonTarget[i] <- sum(coexpr & !labels) / n_nontarget * 100
    
    # Single gene percentages
    gene1_pos <- expr1 > hard_threshold
    gene2_pos <- expr2 > hard_threshold
    
    metrics$pct.Gene1.Target[i] <- sum(gene1_pos & labels) / n_target * 100
    metrics$pct.Gene2.Target[i] <- sum(gene2_pos & labels) / n_target * 100
    metrics$pct.Gene1.NonTarget[i] <- sum(gene1_pos & !labels) / n_nontarget * 100
    metrics$pct.Gene2.NonTarget[i] <- sum(gene2_pos & !labels) / n_nontarget * 100
  }
  
  # Add derived metrics
  metrics$Predictive.power <- abs(metrics$AUC - 0.5) * 2
  metrics$pct.Difference <- metrics$pct.Target - metrics$pct.NonTarget
  
  # Sort by AUC descending
  metrics <- metrics[order(metrics$AUC, decreasing = TRUE), ]
  
  return(metrics)
}
