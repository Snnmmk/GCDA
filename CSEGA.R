#' Gene-based cell-type composition analysis
#'
#' @description
#' Given a Seurat object, a metadata column (cell type or cluster), and a gene name,
#' this function compares the composition of each cell type between High and Low
#' expression groups of the specified gene.
#'
#' @param seu A Seurat object.
#' @param group_var Metadata column name (e.g., "celltype_minor" or "orig.ident").
#' @param gene Target gene symbol (e.g., "DDIT3").
#' @param threshold_method Method to split expression into High/Low ("median" or a numeric cutoff).
#' @return A list with:
#'   \item{summary}{Data frame of proportions, fold change, odds ratio, and adjusted p-values.}
#'   \item{heatmap}{A ComplexHeatmap object visualizing proportions.}
#' @export

CSEGA <- function(seu, group_var, gene, threshold_method = "median") {
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  
  # 检查输入
  if (!gene %in% rownames(seu)) stop(paste("Gene", gene, "not found in object."))
  if (!group_var %in% colnames(seu@meta.data)) stop(paste("Column", group_var, "not found."))
  
  # 提取表达量
  expr <- FetchData(seu, vars = gene)
  
  # 阈值分组
  if (threshold_method == "median") {
    threshold <- median(expr[[gene]], na.rm = TRUE)
  } else if (is.numeric(threshold_method)) {
    threshold <- threshold_method
  } else {
    stop("threshold_method must be 'median' or a numeric cutoff.")
  }
  
  seu$expr_group <- ifelse(expr[[gene]] > threshold, "High", "Low")
  
  # 构造表格
  tab <- table(seu@meta.data[[group_var]], seu$expr_group)
  prop_by_group <- prop.table(tab, margin = 2)
  mat <- as.matrix(prop_by_group)
  storage.mode(mat) <- "numeric"
  
  # 比例检验（每个cell type）
  celltypes <- rownames(tab)
  res_list <- lapply(celltypes, function(ct){
    x_high <- tab[ct, "High"]
    n_high <- sum(tab[, "High"])
    x_low  <- tab[ct, "Low"]
    n_low  <- sum(tab[, "Low"])
    prop_res <- prop.test(x = c(x_high, x_low), n = c(n_high, n_low))
    p_high <- x_high / n_high
    p_low  <- x_low / n_low
    fold_change <- ifelse(p_low == 0, NA, p_high / p_low)
    or <- ((x_high + 0.5)/(n_high - x_high + 0.5)) /
      ((x_low + 0.5)/(n_low - x_low + 0.5))
    data.frame(celltype = ct,
               x_high = x_high, n_high = n_high, p_high = p_high,
               x_low = x_low, n_low = n_low, p_low = p_low,
               fold_change = fold_change,
               odds_ratio = or,
               p_value = prop_res$p.value)
  })
  res_df <- bind_rows(res_list) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_adj)
  
  # 显著性标注
  sig_cells <- res_df %>%
    mutate(sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )) %>%
    select(celltype, sig_label)
  rownames(sig_cells) <- sig_cells$celltype
  
  # 配色方案
  col_fun <- colorRamp2(
    c(min(mat), mean(mat), max(mat)),
    c("#4575B4", "#FFFFBF", "#D73027")
  )
  
  # 热图
  order_by <- order(mat[, "High"], decreasing = TRUE)
  mat <- mat[order_by, ]
  
  ht <- Heatmap(
    mat,
    name = "Proportion",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_title = paste0(gene, " expression group"),
    row_title = group_var,
    heatmap_legend_param = list(title = "Proportion\nwithin group"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      ct <- rownames(mat)[i]
      if (sig_cells[ct, "sig_label"] != "") {
        grid.text(sig_cells[ct, "sig_label"], x, y)
      }
    }
  )
  
  list(summary = res_df, heatmap = ht)
}