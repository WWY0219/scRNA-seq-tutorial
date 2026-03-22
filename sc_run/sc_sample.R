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
