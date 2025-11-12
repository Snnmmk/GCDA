#' CSEGA: Cell-type Specific Expression-Guided Analysis 

#' @description
#' This function analyzes how the expression of a given gene is associated with the composition of each cell type
#' in a single-cell dataset. It extends the conventional thresholding methods (mean, median, quantile) by
#' introducing two advanced modes:
#' (1) **k-means-based thresholding**, which automatically partitions expression values into two clusters (High/Low);
#' (2) **"learn" mode**, which searches across quantile thresholds (from 1% to 99%) to identify the optimal cutoff
#' that maximizes the cell-type compositional divergence.
#'
#' For interpretability, the function further fits a logistic regression model for each cell type:
#' \deqn{I(celltype == c_t) ~ expression}
#' estimating the β coefficient and its statistical significance, thereby quantifying how strongly the gene’s
#' expression level predicts the probability of belonging to that specific cell type.
#'
#' @param seu A Seurat object containing single-cell transcriptomic data.
#' @param group_var Metadata column indicating the grouping variable (e.g., "celltype").
#' @param gene Target gene symbol (e.g., "DDIT3").
#' @param threshold_method Thresholding method, one of:
#'   - "median", "mean", "q0.xx" (quantile-based);
#'   - "kmeans" (two-cluster k-means);
#'   - "learn" (automatically learns the best cutoff from quantiles 0.01–0.99);
#'   - or a numeric value.
#'
#' @return A list containing:
#'   \item{summary}{Data frame of proportions, fold changes, odds ratios, logistic β values, and adjusted p-values.}
#'   \item{heatmap}{A ComplexHeatmap object visualizing compositional differences.}
#'   \item{bubbleplot}{A ggplot2 bubble plot showing High/Low proportions and significance.}
#'   \item{adaptive_cutoff}{The learned or applied cutoff threshold.}
#'   \item{models}{List of per-celltype logistic regression model objects.}
#'
#' @details
#' The "learn" mode iteratively evaluates candidate thresholds and selects the one that maximizes
#' the chi-square statistic between High/Low groups and cell-type composition, thus making the threshold
#' data-driven rather than fixed. Logistic regression enhances interpretability by providing
#' effect size (β) and significance for each cell type.
#'
#' @export


CSEGA <- function(seu,
                  group_var,
                  gene,
                  threshold_method = "kmeans",
                  learn_grid = seq(0.01, 0.99, by = 0.01),
                  min_cells_per_type = 5,
                  verbose = TRUE) {
  # ---- dependencies ----
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Please install Seurat.")
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Please install ComplexHeatmap.")
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Please install circlize.")
  # load libs (quietly)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tidyr)
  # ---- checks ----
  if (!is(seu, "Seurat")) stop("`seu` must be a Seurat object.")
  if (!gene %in% rownames(seu)) stop(paste0("Gene '", gene, "' not found in Seurat object rownames."))
  if (!group_var %in% colnames(seu@meta.data)) stop(paste0("Metadata column '", group_var, "' not found in Seurat object."))
  # ---- extract expression vector ----
  expr_df <- tryCatch(FetchData(seu, vars = gene), error = function(e) stop("FetchData failed: ", e$message))
  expr <- as.numeric(expr_df[[gene]])
  # handle NA
  if (all(is.na(expr))) stop("Expression vector is all NA.")
  # if constant expression -> warn and fallback to median
  if (var(expr, na.rm = TRUE) == 0) {
    warning("Expression is constant; falling back to median split.")
    threshold_method <- "median"
  }
  # ---- determine threshold ----
  adaptive_cutoff <- NA_real_
  if (is.numeric(threshold_method)) {
    adaptive_cutoff <- as.numeric(threshold_method)
    if (verbose) message("Using numeric cutoff: ", round(adaptive_cutoff, 6))
  } else if (threshold_method == "median") {
    adaptive_cutoff <- median(expr, na.rm = TRUE)
    if (verbose) message("Using median cutoff: ", round(adaptive_cutoff, 6))
  } else if (threshold_method == "mean") {
    adaptive_cutoff <- mean(expr, na.rm = TRUE)
    if (verbose) message("Using mean cutoff: ", round(adaptive_cutoff, 6))
  } else if (grepl("^q[0-9.]+$", threshold_method)) {
    qval <- as.numeric(sub("q", "", threshold_method))
    if (is.na(qval) || qval <= 0 || qval >= 1) stop("quantile must be between 0 and 1, e.g., 'q0.25'")
    adaptive_cutoff <- quantile(expr, probs = qval, na.rm = TRUE)
    if (verbose) message("Using quantile cutoff ", threshold_method, ": ", round(adaptive_cutoff, 6))
  } else if (threshold_method == "kmeans") {
    # optional log1p handling for count-like data: do not force, leave user to preprocess
    set.seed(123)
    km <- kmeans(expr, centers = 2, nstart = 10)
    high_cluster <- which.max(km$centers)
    # cluster centers are numeric -> choose midpoint between centers as cutoff for stability
    centers <- sort(km$centers)
    adaptive_cutoff <- mean(centers)
    if (verbose) message("Using kmeans cutoff (midpoint of centers): ", round(adaptive_cutoff, 6))
  } else if (threshold_method == "learn") {
    # grid search thresholds by quantiles; choose threshold maximizing chi-square statistic between High/Low and group_var
    if (verbose) message("Learning cutoff by grid search (max chi-square association)...")
    qs <- quantile(expr, probs = learn_grid, na.rm = TRUE)
    qs <- unique(as.numeric(qs))
    best_stat <- -Inf
    best_q <- NA_real_
    for (t in qs) {
      lab <- ifelse(expr > t, "High", "Low")
      tb <- table(as.character(seu@meta.data[[group_var]]), lab)
      # require at least 2x2
      if (nrow(tb) < 2 || ncol(tb) < 2) next
      # use chisq.test; catch warnings for low counts
      ct <- tryCatch(chisq.test(tb), error = function(e) NULL, warning = function(w) suppressWarnings(chisq.test(tb)))
      if (is.null(ct)) next
      stat <- as.numeric(ct$statistic)
      if (!is.na(stat) && stat > best_stat) {
        best_stat <- stat; best_q <- t
      }
    }
    if (is.na(best_q)) {
      warning("Learn failed to find discriminative threshold; fallback to median.")
      adaptive_cutoff <- median(expr, na.rm = TRUE)
    } else {
      adaptive_cutoff <- best_q
      if (verbose) message("Learned cutoff (best chi-square) = ", round(adaptive_cutoff, 6))
    }
  } else {
    stop("threshold_method must be numeric, 'median','mean','q0.xx','kmeans', or 'learn'.")
  }
  # ---- set expr_group in Seurat metadata ----
  seu$expr_group <- ifelse(expr > adaptive_cutoff, "High", "Low")
  # ---- build contingency / proportions ----
  tab <- table(seu@meta.data[[group_var]], seu$expr_group)
  prop_by_group <- prop.table(tab, margin = 2) # proportions within each column (High/Low)
  mat <- as.matrix(prop_by_group)
  storage.mode(mat) <- "numeric"
  # ---- per-celltype logistic models: model P(celltype == ct | expr) = logit^{-1}(α + β * expr) ----
  celltypes <- rownames(tab)
  models <- list()
  summ_list <- vector("list", length(celltypes))
  names(summ_list) <- celltypes
  all_meta <- seu@meta.data
  all_meta$expr <- expr
  for (i in seq_along(celltypes)) {
    ct <- celltypes[i]
    # construct binary label across all cells
    y <- as.integer(as.character(all_meta[[group_var]]) == ct)
    df <- data.frame(y = y, expr = all_meta$expr)
    # require variation in y and expr not constant
    if (length(unique(y)) < 2 || var(df$expr, na.rm = TRUE) == 0) {
      models[[ct]] <- NULL
      summ_list[[ct]] <- tibble::tibble(celltype = ct, beta = NA_real_, p = NA_real_)
      next
    }
    glm_fit <- tryCatch(glm(y ~ expr, data = df, family = binomial), error = function(e) NULL, warning = function(w) tryCatch(glm(y ~ expr, data = df, family = binomial), error = function(e) NULL))
    if (is.null(glm_fit)) {
      models[[ct]] <- NULL
      summ_list[[ct]] <- tibble::tibble(celltype = ct, beta = NA_real_, p = NA_real_)
      next
    }
    models[[ct]] <- glm_fit
    coefs <- summary(glm_fit)$coefficients
    if ("expr" %in% rownames(coefs)) {
      beta <- coefs["expr", "Estimate"]
      pval <- coefs["expr", "Pr(>|z|)"]
    } else {
      beta <- NA_real_; pval <- NA_real_
    }
    summ_list[[ct]] <- tibble::tibble(celltype = ct, beta = beta, p = pval)
  }
  summary_df <- dplyr::bind_rows(summ_list) %>%
    # join proportion metrics from res_df logic (x_high etc.)
    left_join(
      (as.data.frame(tab) %>%
         tidyr::pivot_wider(names_from = Var2, values_from = Freq) %>%
         dplyr::rename(celltype = Var1) %>%
         mutate(n_high = ifelse(is.na(High), 0, High),
                n_low = ifelse(is.na(Low), 0, Low),
                n_total = n_high + n_low,
                p_high = ifelse(n_total==0, NA, n_high / sum(n_high, n_low)),
                p_low  = ifelse(n_total==0, NA, n_low / sum(n_high, n_low))
         ) %>%
         select(celltype, n_high, n_low)),
      by = "celltype"
    )
  # compute fold_change and odds_ratio and p_value by prop.test per celltype (keep for compatibility)
  prop_stats <- lapply(celltypes, function(ct) {
    if (!ct %in% rownames(tab)) {
      return(tibble::tibble(celltype = ct, x_high = NA_integer_, n_high = NA_integer_, x_low = NA_integer_, n_low = NA_integer_, fold_change = NA_real_, odds_ratio = NA_real_, p_value = NA_real_))
    }
    x_high <- tab[ct, "High"]
    n_high <- sum(tab[, "High"])
    x_low  <- tab[ct, "Low"]
    n_low  <- sum(tab[, "Low"])
    # safe prop.test
    safe_pt <- tryCatch(prop.test(x = c(x_high, x_low), n = c(n_high, n_low)), error = function(e) NULL)
    pval <- if (is.null(safe_pt)) NA_real_ else safe_pt$p.value
    p_high <- if (n_high == 0) NA_real_ else x_high / n_high
    p_low  <- if (n_low == 0) NA_real_ else x_low / n_low
    fold_change <- if (is.na(p_low) || p_low == 0) NA_real_ else p_high / p_low
    or <- ((x_high + 0.5)/(n_high - x_high + 0.5)) / ((x_low + 0.5)/(n_low - x_low + 0.5))
    tibble::tibble(celltype = ct, x_high = x_high, n_high = n_high, x_low = x_low, n_low = n_low, fold_change = fold_change, odds_ratio = or, p_value = pval)
  })
  prop_stats_df <- dplyr::bind_rows(prop_stats)
  # combine
  final_df <- dplyr::left_join(summary_df, prop_stats_df, by = "celltype") %>%
    dplyr::mutate(p_adj = ifelse(is.na(p), NA_real_, p.adjust(p, method = "BH")),
                  sig_label = dplyr::case_when(
                    !is.na(p_adj) & p_adj < 0.001 ~ "***",
                    !is.na(p_adj) & p_adj < 0.01  ~ "**",
                    !is.na(p_adj) & p_adj < 0.05  ~ "*",
                    TRUE ~ ""
                  ))
  # ---- heatmap (proportions) ----
  # ensure columns "High" and "Low" exist
  if (!"High" %in% colnames(mat)) mat <- cbind(mat, High = rep(0, nrow(mat)))
  if (!"Low" %in% colnames(mat))  mat <- cbind(mat, Low  = rep(0, nrow(mat)))
  order_by <- order(mat[, "High"], decreasing = TRUE)
  mat_ordered <- mat[order_by, , drop = FALSE]
  col_fun <- colorRamp2(c(min(mat_ordered), mean(mat_ordered), max(mat_ordered)), c("#057dcd", "#FFFFBF", "#e50000"))
  sig_cells <- final_df %>% select(celltype, sig_label) %>% tibble::column_to_rownames("celltype")
  # heatmap cell_fun must handle missing sig_cells rows robustly
  ht <- Heatmap(mat_ordered,
                name = "Proportion",
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_gp = gpar(fontsize = 12),
                column_title = paste0(gene, " expression"),
                column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                heatmap_legend_param = list(title = "Proportion\nwithin group", title_gp = gpar(fontsize = 12, fontface = "bold")),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.rect(x, y, width, height, gp = gpar(col = "grey", fill = NA, lwd = 0.5))
                  ct <- rownames(mat_ordered)[i]
                  lab <- if (ct %in% rownames(sig_cells)) sig_cells[ct, "sig_label"] else ""
                  if (!is.null(lab) && nzchar(lab)) grid.text(lab, x, y, gp = gpar(col = "black", fontsize = 10, fontface = "bold"))
                })
  # ---- bubble plot (uses proportion per celltype in High/Low and colors by beta) ----
  bubble_df <- prop_stats_df %>%
    tidyr::pivot_longer(cols = c("x_high", "x_low"), names_to = "count_name", values_to = "count") %>%
    mutate(expr_group = ifelse(count_name == "x_high", "High", "Low")) %>%
    group_by(expr_group) %>%
    mutate(proportion = ifelse(sum(count, na.rm = TRUE) == 0, 0, count / sum(count, na.rm = TRUE))) %>%
    ungroup() %>%
    rename(celltype = celltype)
  # join beta
  bubble_df <- bubble_df %>% left_join(final_df %>% select(celltype, beta, sig_label), by = "celltype")
  # replace NA beta with 0 for coloring but keep NA record
  bubble_df$beta_for_plot <- ifelse(is.na(bubble_df$beta), 0, bubble_df$beta)
  x_max <- mean(as.numeric(factor(bubble_df$expr_group)))
  bp <- ggplot(bubble_df, aes(x = expr_group, y = celltype)) +
    geom_point(aes(size = proportion, fill = beta_for_plot), shape = 21, color = "grey30") +
    scale_size_continuous(range = c(3, 10), name = "Proportion") +
    scale_fill_gradient2(low = "#057dcd", mid = "#FFFFBF", high = "#e50000", midpoint = 0, name = "β (expr → P(celltype))") +
    geom_text(data = bubble_df %>% distinct(celltype, sig_label), aes(x = x_max, y = celltype, label = sig_label), inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 4, fontface = "bold") +
    theme_minimal(base_family = "Arial") +
    theme(axis.text.y = element_text(face = "bold", size = 10), axis.text.x = element_text(face = "bold", size = 10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color = "grey30", fill = NA, size = 0.5)) +
    labs(x = "Expression group", y = "Cell type")
  # ---- return ----
  return(list(
    summary = final_df %>% arrange(p_adj),
    adaptive_cutoff = adaptive_cutoff,
    heatmap = ht,
    bubbleplot = bp,
    models = models
  ))
}
