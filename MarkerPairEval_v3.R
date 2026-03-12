# =============================================================================
# MarkerPairEval v3
#
# Evaluate all pairwise combinations of candidate marker genes for their ability
# to distinguish a target cell population (e.g. mesenchymal) from non-target
# cells using ROC analysis on a combined two-gene score.
#
# For each gene pair the function reports:
#   - ROC-based metrics   : AUC, optimal threshold (Youden's J), sensitivity,
#                           specificity, accuracy, precision, F1
#   - Hard-threshold metrics: same classification metrics using a fixed
#                             expression cutoff (default 0 = any detection)
#   - Detection percentages: co-expression and per-gene % of target vs
#                            non-target cells, at both ROC and hard thresholds
#   - Mean expression      : per-gene raw-scale means split by population,
#                            by ROC-positive/negative, and by hard-positive/negative
#   - Derived metrics      : predictive power, pct difference, log2 fold-change
#
# Arguments:
#   genes             – character vector of candidate gene names
#   labels            – logical vector (TRUE = target cell); same length as
#                       ncol(expression_matrix)
#   expression_matrix – log-normalised expression matrix (genes x cells),
#                       e.g. seurat@assays$RNA@data
#   combination       – how to combine the two genes' expression into a single
#                       score: "AND" (min), "OR" (max), "sum", or "product"
#   hard_threshold    – fixed expression cutoff for the hard-threshold branch
#                       (default 0, i.e. any non-zero expression counts as positive)
#
# Returns:
#   A data.frame with one row per gene pair, sorted by AUC (descending).
# =============================================================================

MarkerPairEval <- function(genes, labels, expression_matrix,
                           combination = c("AND", "OR", "sum", "product"),
                           hard_threshold = 0) {

  combination <- match.arg(combination)

  # --- Input cleaning --------------------------------------------------------
  # Remove cells with NA labels and coerce to logical
  if (any(is.na(labels))) {
    sel_na <- which(is.na(labels))
    labels <- labels[-sel_na]
    expression_matrix <- expression_matrix[, -sel_na]
  }
  labels <- as.logical(labels)

  # Cap at 100 genes to keep runtime manageable (C(100,2) = 4950 pairs)
  if (length(genes) > 100) genes <- genes[1:100]

  # Keep only genes present in the expression matrix
  genes <- genes[genes %in% rownames(expression_matrix)]

  # --- Enumerate all gene pairs ----------------------------------------------
  gp      <- combn(genes, 2)
  n_pairs <- ncol(gp)

  n_target    <- sum(labels)
  n_nontarget <- sum(!labels)

  # --- Pre-allocate output data.frame ----------------------------------------
  metrics <- data.frame(
    Gene.1 = gp[1, ],
    Gene.2 = gp[2, ],

    # ROC-based classification metrics (optimal threshold via Youden's J)
    AUC             = numeric(n_pairs),
    ROC.Threshold   = numeric(n_pairs),
    ROC.Sensitivity = numeric(n_pairs),
    ROC.Specificity = numeric(n_pairs),
    ROC.Accuracy    = numeric(n_pairs),
    ROC.Precision   = numeric(n_pairs),
    ROC.F1.score    = numeric(n_pairs),

    # Hard-threshold classification metrics (expression > hard_threshold)
    Hard.Sensitivity = numeric(n_pairs),
    Hard.Specificity = numeric(n_pairs),
    Hard.Accuracy    = numeric(n_pairs),
    Hard.Precision   = numeric(n_pairs),
    Hard.F1.score    = numeric(n_pairs),

    # Co-expression % at hard threshold (both genes > hard_threshold)
    Hard.pct.Target         = numeric(n_pairs),
    Hard.pct.NonTarget      = numeric(n_pairs),
    Hard.pct.Gene1.Target   = numeric(n_pairs),
    Hard.pct.Gene2.Target   = numeric(n_pairs),
    Hard.pct.Gene1.NonTarget = numeric(n_pairs),
    Hard.pct.Gene2.NonTarget = numeric(n_pairs),

    # Co-expression % at ROC threshold (both genes >= ROC threshold)
    ROC.pct.Target          = numeric(n_pairs),
    ROC.pct.NonTarget       = numeric(n_pairs),
    ROC.pct.Gene1.Target    = numeric(n_pairs),
    ROC.pct.Gene2.Target    = numeric(n_pairs),
    ROC.pct.Gene1.NonTarget = numeric(n_pairs),
    ROC.pct.Gene2.NonTarget = numeric(n_pairs),

    # Mean raw-scale expression by population (target vs non-target)
    mean.Gene1.Target    = numeric(n_pairs),
    mean.Gene1.NonTarget = numeric(n_pairs),
    mean.Gene2.Target    = numeric(n_pairs),
    mean.Gene2.NonTarget = numeric(n_pairs),

    # Mean raw-scale expression by ROC classification
    mean.Gene1.aboveROC = numeric(n_pairs),
    mean.Gene1.belowROC = numeric(n_pairs),
    mean.Gene2.aboveROC = numeric(n_pairs),
    mean.Gene2.belowROC = numeric(n_pairs),

    # Mean raw-scale expression by hard-threshold classification
    mean.Gene1.aboveHard = numeric(n_pairs),
    mean.Gene1.belowHard = numeric(n_pairs),
    mean.Gene2.aboveHard = numeric(n_pairs),
    mean.Gene2.belowHard = numeric(n_pairs),

    stringsAsFactors = FALSE
  )

  # --- Main loop: evaluate each gene pair ------------------------------------
  for (i in seq_len(n_pairs)) {

    # Extract log-normalised expression for the two genes
    expr1 <- as.numeric(expression_matrix[gp[1, i], ])
    expr2 <- as.numeric(expression_matrix[gp[2, i], ])

    # Back-transform to raw counts for mean-expression summaries
    expr1_raw <- expm1(expr1)
    expr2_raw <- expm1(expr2)

    # Combine the two genes into a single score
    #   AND     → min (both genes must be high)
    #   OR      → max (either gene suffices)
    #   sum     → additive
    #   product → multiplicative
    combined_score <- switch(combination,
      "AND"     = pmin(expr1, expr2),
      "OR"      = pmax(expr1, expr2),
      "sum"     = expr1 + expr2,
      "product" = expr1 * expr2
    )

    # ----- ROC analysis & optimal threshold (Youden's J) ---------------------
    roc_obj <- pROC::roc(labels, combined_score, quiet = TRUE, direction = "<")
    metrics$AUC[i] <- as.numeric(pROC::auc(roc_obj))

    best_coords <- pROC::coords(roc_obj, "best", best.method = "youden", ret = "all")
    # If multiple optimal points exist, take the first
    if (is.data.frame(best_coords) && nrow(best_coords) > 1) {
      best_coords <- best_coords[1, ]
    }

    roc_threshold <- best_coords$threshold
    metrics$ROC.Threshold[i]   <- roc_threshold
    metrics$ROC.Sensitivity[i] <- best_coords$sensitivity
    metrics$ROC.Specificity[i] <- best_coords$specificity

    # Confusion-matrix entries at ROC threshold
    roc_predicted <- combined_score >= roc_threshold
    tp_roc <- sum(roc_predicted & labels)
    tn_roc <- sum(!roc_predicted & !labels)
    fp_roc <- sum(roc_predicted & !labels)
    fn_roc <- sum(!roc_predicted & labels)

    metrics$ROC.Accuracy[i]  <- (tp_roc + tn_roc) / length(labels)
    metrics$ROC.Precision[i] <- ifelse((tp_roc + fp_roc) > 0,
                                       tp_roc / (tp_roc + fp_roc), 0)
    metrics$ROC.F1.score[i]  <- ifelse((2 * tp_roc + fp_roc + fn_roc) > 0,
                                       (2 * tp_roc) / (2 * tp_roc + fp_roc + fn_roc), 0)

    # ----- Hard-threshold classification -------------------------------------
    # For AND/product both genes must exceed; for OR/sum either gene suffices
    hard_positive <- switch(combination,
      "AND"     = (expr1 > hard_threshold) & (expr2 > hard_threshold),
      "OR"      = (expr1 > hard_threshold) | (expr2 > hard_threshold),
      "sum"     = (expr1 > hard_threshold) | (expr2 > hard_threshold),
      "product" = (expr1 > hard_threshold) & (expr2 > hard_threshold)
    )

    tp_hard <- sum(hard_positive & labels)
    tn_hard <- sum(!hard_positive & !labels)
    fp_hard <- sum(hard_positive & !labels)
    fn_hard <- sum(!hard_positive & labels)

    metrics$Hard.Sensitivity[i] <- ifelse((tp_hard + fn_hard) > 0,
                                          tp_hard / (tp_hard + fn_hard), 0)
    metrics$Hard.Specificity[i] <- ifelse((tn_hard + fp_hard) > 0,
                                          tn_hard / (tn_hard + fp_hard), 0)
    metrics$Hard.Accuracy[i]    <- (tp_hard + tn_hard) / length(labels)
    metrics$Hard.Precision[i]   <- ifelse((tp_hard + fp_hard) > 0,
                                          tp_hard / (tp_hard + fp_hard), 0)
    metrics$Hard.F1.score[i]    <- ifelse((2 * tp_hard + fp_hard + fn_hard) > 0,
                                          (2 * tp_hard) / (2 * tp_hard + fp_hard + fn_hard), 0)

    # ----- Detection percentages at hard threshold ---------------------------
    # Co-expression: both genes above hard_threshold
    hard_coexpr <- (expr1 > hard_threshold) & (expr2 > hard_threshold)
    metrics$Hard.pct.Target[i]    <- sum(hard_coexpr & labels) / n_target * 100
    metrics$Hard.pct.NonTarget[i] <- sum(hard_coexpr & !labels) / n_nontarget * 100

    # Per-gene detection at hard threshold
    gene1_hard_pos <- expr1 > hard_threshold
    gene2_hard_pos <- expr2 > hard_threshold
    metrics$Hard.pct.Gene1.Target[i]    <- sum(gene1_hard_pos & labels)  / n_target    * 100
    metrics$Hard.pct.Gene2.Target[i]    <- sum(gene2_hard_pos & labels)  / n_target    * 100
    metrics$Hard.pct.Gene1.NonTarget[i] <- sum(gene1_hard_pos & !labels) / n_nontarget * 100
    metrics$Hard.pct.Gene2.NonTarget[i] <- sum(gene2_hard_pos & !labels) / n_nontarget * 100

    # ----- Detection percentages at ROC threshold ----------------------------
    # Co-expression: both genes at or above ROC threshold
    roc_coexpr <- (expr1 >= roc_threshold) & (expr2 >= roc_threshold)
    metrics$ROC.pct.Target[i]    <- sum(roc_coexpr & labels) / n_target * 100
    metrics$ROC.pct.NonTarget[i] <- sum(roc_coexpr & !labels) / n_nontarget * 100

    # Per-gene detection at ROC threshold
    gene1_roc_pos <- expr1 >= roc_threshold
    gene2_roc_pos <- expr2 >= roc_threshold
    metrics$ROC.pct.Gene1.Target[i]    <- sum(gene1_roc_pos & labels)  / n_target    * 100
    metrics$ROC.pct.Gene2.Target[i]    <- sum(gene2_roc_pos & labels)  / n_target    * 100
    metrics$ROC.pct.Gene1.NonTarget[i] <- sum(gene1_roc_pos & !labels) / n_nontarget * 100
    metrics$ROC.pct.Gene2.NonTarget[i] <- sum(gene2_roc_pos & !labels) / n_nontarget * 100

    # ----- Mean raw-scale expression by population ---------------------------
    metrics$mean.Gene1.Target[i]    <- mean(expr1_raw[labels])
    metrics$mean.Gene1.NonTarget[i] <- mean(expr1_raw[!labels])
    metrics$mean.Gene2.Target[i]    <- mean(expr2_raw[labels])
    metrics$mean.Gene2.NonTarget[i] <- mean(expr2_raw[!labels])

    # ----- Mean raw-scale expression by ROC classification -------------------
    metrics$mean.Gene1.aboveROC[i] <- ifelse(sum(roc_predicted) > 0,
                                             mean(expr1_raw[roc_predicted]), NA)
    metrics$mean.Gene1.belowROC[i] <- ifelse(sum(!roc_predicted) > 0,
                                             mean(expr1_raw[!roc_predicted]), NA)
    metrics$mean.Gene2.aboveROC[i] <- ifelse(sum(roc_predicted) > 0,
                                             mean(expr2_raw[roc_predicted]), NA)
    metrics$mean.Gene2.belowROC[i] <- ifelse(sum(!roc_predicted) > 0,
                                             mean(expr2_raw[!roc_predicted]), NA)

    # ----- Mean raw-scale expression by hard classification ------------------
    metrics$mean.Gene1.aboveHard[i] <- ifelse(sum(hard_positive) > 0,
                                              mean(expr1_raw[hard_positive]), NA)
    metrics$mean.Gene1.belowHard[i] <- ifelse(sum(!hard_positive) > 0,
                                              mean(expr1_raw[!hard_positive]), NA)
    metrics$mean.Gene2.aboveHard[i] <- ifelse(sum(hard_positive) > 0,
                                              mean(expr2_raw[hard_positive]), NA)
    metrics$mean.Gene2.belowHard[i] <- ifelse(sum(!hard_positive) > 0,
                                              mean(expr2_raw[!hard_positive]), NA)
  }

  # --- Derived summary metrics -----------------------------------------------
  # Predictive power: rescaled AUC from [0.5, 1] → [0, 1]
  metrics$Predictive.power <- abs(metrics$AUC - 0.5) * 2

  # Difference in co-expression % between target and non-target
  metrics$Hard.pct.Difference <- metrics$Hard.pct.Target - metrics$Hard.pct.NonTarget
  metrics$ROC.pct.Difference  <- metrics$ROC.pct.Target  - metrics$ROC.pct.NonTarget

  # Per-gene log2 fold-change (target vs non-target), pseudocount = 1
  metrics$log2FC.Gene1 <- log2((metrics$mean.Gene1.Target + 1) /
                               (metrics$mean.Gene1.NonTarget + 1))
  metrics$log2FC.Gene2 <- log2((metrics$mean.Gene2.Target + 1) /
                               (metrics$mean.Gene2.NonTarget + 1))

  # Sort by AUC descending
  metrics <- metrics[order(metrics$AUC, decreasing = TRUE), ]

  return(metrics)
}
