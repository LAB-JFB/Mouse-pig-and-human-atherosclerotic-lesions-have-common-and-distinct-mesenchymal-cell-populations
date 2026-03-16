# =============================================================================
# SingleMarkerEval v3: pct values at both ROC and hard thresholds
# =============================================================================

SingleMarkerEval <- function(genes, labels, expression_matrix, hard_threshold = 0) {
  
  # Check for NA's in labels
  if(any(is.na(labels))) {
    sel_na <- which(is.na(labels))
    labels <- labels[-sel_na]
    expression_matrix <- expression_matrix[, -sel_na]
  }
  labels <- as.logical(labels)
  
  # Filter to genes present in matrix
  genes <- genes[genes %in% rownames(expression_matrix)]
  n_genes <- length(genes)
  
  if(n_genes == 0) {
    stop("No genes found in expression matrix")
  }
  
  n_target <- sum(labels)
  n_nontarget <- sum(!labels)
  
  # Pre-allocate results
  metrics <- data.frame(
    Gene = genes,
    
    # ---------------------------------------------------------
    # ROC-BASED METRICS (optimal threshold via Youden's J)
    # ---------------------------------------------------------
    AUC = numeric(n_genes),
    ROC.Threshold = numeric(n_genes),
    ROC.Sensitivity = numeric(n_genes),
    ROC.Specificity = numeric(n_genes),
    ROC.Accuracy = numeric(n_genes),
    ROC.Precision = numeric(n_genes),
    ROC.F1.score = numeric(n_genes),
    
    # ---------------------------------------------------------
    # HARD THRESHOLD METRICS (expression > hard_threshold)
    # ---------------------------------------------------------
    Hard.Sensitivity = numeric(n_genes),
    Hard.Specificity = numeric(n_genes),
    Hard.Accuracy = numeric(n_genes),
    Hard.Precision = numeric(n_genes),
    Hard.F1.score = numeric(n_genes),
    
    # ---------------------------------------------------------
    # PERCENTAGES AT HARD THRESHOLD
    # ---------------------------------------------------------
    Hard.pct.Target = numeric(n_genes),
    Hard.pct.NonTarget = numeric(n_genes),
    
    # ---------------------------------------------------------
    # PERCENTAGES AT ROC THRESHOLD
    # ---------------------------------------------------------
    ROC.pct.Target = numeric(n_genes),
    ROC.pct.NonTarget = numeric(n_genes),
    
    # ---------------------------------------------------------
    # MEAN EXPRESSION IN RAW COUNTS: by population
    # ---------------------------------------------------------
    mean.Target = numeric(n_genes),
    mean.NonTarget = numeric(n_genes),
    log2FC = numeric(n_genes),
    
    # ---------------------------------------------------------
    # MEAN EXPRESSION IN RAW COUNTS: by ROC threshold
    # ---------------------------------------------------------
    mean.aboveROC = numeric(n_genes),
    mean.belowROC = numeric(n_genes),
    n.aboveROC = numeric(n_genes),
    n.belowROC = numeric(n_genes),
    
    # ---------------------------------------------------------
    # MEAN EXPRESSION IN RAW COUNTS: by hard threshold
    # ---------------------------------------------------------
    mean.aboveHard = numeric(n_genes),
    mean.belowHard = numeric(n_genes),
    n.aboveHard = numeric(n_genes),
    n.belowHard = numeric(n_genes),
    
    stringsAsFactors = FALSE
  )
  
  # Calculate metrics for each gene
  for(i in seq_len(n_genes)) {
    
    expr <- as.numeric(expression_matrix[genes[i], ])
    
    # Convert to raw counts for mean calculations
    expr_raw <- expm1(expr)
    
    # =========================================
    # ROC-BASED METRICS (Youden's J threshold)
    # =========================================
    
    roc_obj <- pROC::roc(labels, expr, quiet = TRUE, direction = "<")
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
    roc_predicted <- expr >= roc_threshold
    
    tp_roc <- sum(roc_predicted & labels)
    tn_roc <- sum(!roc_predicted & !labels)
    fp_roc <- sum(roc_predicted & !labels)
    fn_roc <- sum(!roc_predicted & labels)
    
    metrics$ROC.Accuracy[i] <- (tp_roc + tn_roc) / (tp_roc + tn_roc + fp_roc + fn_roc)
    metrics$ROC.Precision[i] <- ifelse((tp_roc + fp_roc) > 0, tp_roc / (tp_roc + fp_roc), 0)
    metrics$ROC.F1.score[i] <- ifelse((2*tp_roc + fp_roc + fn_roc) > 0, 
                                       (2*tp_roc) / (2*tp_roc + fp_roc + fn_roc), 0)
    
    # =========================================
    # HARD THRESHOLD METRICS
    # =========================================
    
    hard_positive <- expr > hard_threshold
    
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
    # PERCENTAGES AT HARD THRESHOLD
    # =========================================
    
    metrics$Hard.pct.Target[i] <- sum(hard_positive & labels) / n_target * 100
    metrics$Hard.pct.NonTarget[i] <- sum(hard_positive & !labels) / n_nontarget * 100
    
    # =========================================
    # PERCENTAGES AT ROC THRESHOLD
    # =========================================
    
    metrics$ROC.pct.Target[i] <- sum(roc_predicted & labels) / n_target * 100
    metrics$ROC.pct.NonTarget[i] <- sum(roc_predicted & !labels) / n_nontarget * 100
    
    # =========================================
    # MEAN EXPRESSION (RAW COUNTS): by population
    # =========================================
    
    metrics$mean.Target[i] <- mean(expr_raw[labels])
    metrics$mean.NonTarget[i] <- mean(expr_raw[!labels])
    metrics$log2FC[i] <- log2((metrics$mean.Target[i] + 1) / (metrics$mean.NonTarget[i] + 1))
    
    # =========================================
    # MEAN EXPRESSION (RAW COUNTS): by ROC threshold
    # =========================================
    
    metrics$n.aboveROC[i] <- sum(roc_predicted)
    metrics$n.belowROC[i] <- sum(!roc_predicted)
    
    metrics$mean.aboveROC[i] <- ifelse(sum(roc_predicted) > 0, 
                                        mean(expr_raw[roc_predicted]), NA)
    metrics$mean.belowROC[i] <- ifelse(sum(!roc_predicted) > 0, 
                                        mean(expr_raw[!roc_predicted]), NA)
    
    # =========================================
    # MEAN EXPRESSION (RAW COUNTS): by hard threshold
    # =========================================
    
    metrics$n.aboveHard[i] <- sum(hard_positive)
    metrics$n.belowHard[i] <- sum(!hard_positive)
    
    metrics$mean.aboveHard[i] <- ifelse(sum(hard_positive) > 0, 
                                         mean(expr_raw[hard_positive]), NA)
    metrics$mean.belowHard[i] <- ifelse(sum(!hard_positive) > 0, 
                                         mean(expr_raw[!hard_positive]), NA)
  }
  
  # Add derived metrics
  metrics$Predictive.power <- abs(metrics$AUC - 0.5) * 2
  metrics$Hard.pct.Difference <- metrics$Hard.pct.Target - metrics$Hard.pct.NonTarget
  metrics$ROC.pct.Difference <- metrics$ROC.pct.Target - metrics$ROC.pct.NonTarget
  
  # Sort by AUC descending
  metrics <- metrics[order(metrics$AUC, decreasing = TRUE), ]
  
  return(metrics)
}
