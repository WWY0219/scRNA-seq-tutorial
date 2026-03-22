#' Downsample Seurat Object by Group with Memory Optimization
#'
#' This function reduces the size of a Seurat object by randomly sampling cells from specified groups.
#' It is optimized for speed using vectorized operations and manages memory efficiently by subsetting
#' before running DietSeurat.
#'
#' @title Downsample Seurat Object by Group
#' @description
#' Performs stratified sampling on a Seurat object. The function can either:
#' 1. Automatically calculate the number of cells per group based on a target total (`sp.total`).
#' 2. Use a fixed number of cells per group (`sp.size`).
#'
#' It also includes an optional `DietSeurat` step to reduce the file size of the resulting object by removing unused assays and reductions.
#'
#' @param obj A Seurat object to be downsampled.
#' @param group.by Character string. The column name in `meta.data` to group cells by (e.g., "seurat_clusters", "cell_type"). Default is "seurat_clusters".
#' @param sp.size Integer or NULL. The specific number of cells to sample per group. If NULL (default), the sample size is calculated automatically as `ceiling(sp.total / n_groups)`.
#' @param sp.total Integer. The approximate total number of cells desired in the final object. This parameter is ignored if `sp.size` is provided. Default is 1000.
#' @param diet Logical. Whether to run `DietSeurat` on the subsetted object to remove unused assays, slots, and graphs. Default is TRUE.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is 1234. Set to NULL to disable seeding.
#' @param keep.reductions Character vector. List of dimensional reductions to keep when `diet = TRUE`. Non-existent reductions will be ignored with a warning. Default is c('pca', 'umap', 'harmony').
#'
#' @return A subsetted (and optionally dieted) Seurat object containing the sampled cells.
#' @author WWY
#' @export
#' @import Seurat
#' @examples
#' \dontrun{
#' # Example 1: Auto-calculate sample size to get approximately 1000 cells total
#' # If you have 10 clusters, it will sample 100 cells per cluster.
#' small_obj <- sc_sample(seurat_obj, sp.total = 1000)
#'
#' # Example 2: Fixed size - Sample exactly 200 cells per cluster
#' fixed_obj <- sc_sample(seurat_obj, group.by = "cell_type", sp.size = 200)
#'
#' # Example 3: Sample without removing other slots (keeping raw counts, etc.)
#' full_slot_obj <- sc_sample(seurat_obj, diet = FALSE)
#' }

Sample_seob_opt <- function(obj, 
                            group.by = "seurat_clusters", 
                            sp.size = NULL, 
                            sp.total = 1000, 
                            diet = TRUE, 
                            seed = 1234,
                            keep.reductions = c('pca', 'umap', 'harmony')) {
  
  # 1. 依赖与参数校验
  if (!inherits(obj, "Seurat")) stop("错误: obj 必须是 Seurat 对象")
  if (!group.by %in% colnames(obj@meta.data)) stop(paste0("错误: meta.data 中找不到列: ", group.by))
  
  # 设置随机种子
  if (!is.null(seed)) set.seed(seed)
  
  meta_data <- obj@meta.data
  groups <- meta_data[[group.by]]
  unique_groups <- unique(groups)
  n_groups <- length(unique_groups)
  
  # 2. 计算每个分组的抽样数
  if (is.null(sp.size)) {
    if (n_groups == 0) stop("错误: 分组列无有效值")
    sp.size <- ceiling(sp.total / n_groups)
    message(sprintf("策略: 自动计算每组抽样数 %d (目标总数 %d / 分组数 %d)", sp.size, sp.total, n_groups))
  } else {
    message(sprintf("策略: 使用固定每组抽样数 %d", sp.size))
  }
  
  # 3. 高效抽样 (向量化操作替代 for 循环)
  # 获取所有细胞ID并按组分割
  cells_by_group <- split(colnames(obj), groups)
  
  # 定义抽样函数
  sample_func <- function(cells, size, group_name) {
    n <- length(cells)
    if (n > size) {
      keep <- sample(cells, size)
      # message(sprintf("  - 分组 %s: %d -> %d", group_name, n, size)) # 调试可开，大数据建议关闭减少IO
      return(keep)
    } else {
      # message(sprintf("  - 分组 %s: %d (保留全部)", group_name, n))
      return(cells)
    }
  }
  
  # 执行抽样
  message("正在执行分组抽样...")
  # mapply 可以同时传递细胞列表和组名，方便(如果需要)打印日志
  cells_to_keep <- unlist(mapply(sample_func, cells_by_group, names(cells_by_group), 
                                 MoreArgs = list(size = sp.size), SIMPLIFY = FALSE), use.names = FALSE)
  
  # 4. 子集化 (Subset)
  # 先 Subset 后 Diet，效率更高
  obj_sub <- subset(obj, cells = cells_to_keep)
  message(sprintf("抽样完成: 原始细胞 %d -> 剩余细胞 %d", ncol(obj), ncol(obj_sub)))
  
  # 5. 对象精简 (DietSeurat)
  if (isTRUE(diet)) {
    # 检查哪些降维结果实际存在于对象中，避免报错
    avail_reducs <- names(obj_sub@reductions)
    valid_reducs <- intersect(keep.reductions, avail_reducs)
    
    if (length(valid_reducs) < length(keep.reductions)) {
      missing <- setdiff(keep.reductions, avail_reducs)
      warning(paste("注意: 以下降维结果不存在，将被忽略:", paste(missing, collapse = ", ")))
    }
    
    message(paste("正在精简对象 (保留降维:", paste(valid_reducs, collapse=","), ")..."))
    obj_sub <- DietSeurat(obj_sub, dimreducs = valid_reducs)
  }
  
  return(obj_sub)
}
