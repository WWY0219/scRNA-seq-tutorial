# =========================================LoadingPackage===========================================================
library(Seurat)
library(SCP)
library(qs)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(circlize)
library(GeneNMF)
library(Seurat)
library(UCell)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
library(qs)
library(msigdbr)
library(fgsea)



# =============================Prepare for GeneNMF=============================================
seurat_obj <- qread("seurat_obj_allcelltypes.qs")
seurat_obj <- SCTransform(seurat_obj,
                          verbose = T,
                          vars.to.regress = c("percent.mt", "percent.HB"),
                          return.only.var.genes = TRUE)
ndim <-15
seurat_obj <-FindVariableFeatures(seurat_obj, nfeatures = 1000)
seurat_obj <- runNMF(seurat_obj, 
                     k = ndim,    # 手动指定降维维度
                     assay="SCT")
seurat_obj@reductions$NMF
## 看NMF分群
seurat_obj <-RunUMAP(seurat_obj, reduction ="NMF", dims=1:ndim, reduction.name ="NMF_UMAP", reduction.key ="nmfUMAP_")
DimPlot(seurat_obj, reduction ="NMF_UMAP", label=T,group.by="subcelltype") + ggtitle("NMF UMAP")


# =======================================Run GeneNMF============================================
seu.list <- SplitObject(seurat_obj, split.by = "group")
geneNMF.programs <- multiNMF(seu.list, 
                             assay="SCT", 
                             slot="data", 
                             k=4:9, 
                             nfeatures = 1000)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nMP=10,
                                        weight.explained = 0.7,
                                        max.genes=100)

# =====================================Visulaization HeatMap=====================================
## 修改元程序的颜色
anno_colors <- list(
  Metaprogram = c(
    "MP1" = "#AED0DF",  # 匹配函数硬编码的MP1
    "MP2" = "#EBAEA9"   
  )
)

## 绘制热图（封装在heatmap）
ph <- plotMetaPrograms(
  mp.res = geneNMF.metaprograms,
  similarity.cutoff = c(0, 1),           # 相似性阈值
  scale = "none",                        # 不标准化
  downsample = Inf,                      # 关闭下采样，保留聚类树
  showtree = TRUE,                       # 显示聚类树
  annotation_colors = anno_colors,       # 注释色块颜色
  main = "Metaprogram Similarity Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE
)
geneNMF.metaprograms$metaprograms.metrics


# ====================================GSEA富集分析================================================
## 将每个MP中的基因构建成一个列表
lapply(geneNMF.metaprograms$metaprograms.genes, head)
## !!需要联网
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, 
          universe=rownames(seurat_obj), 
          category = "H")   # GSEA的database
})
head(top_p$MP2)

# =====================================基因程序的特征打分============================================
mp.genes <- geneNMF.metaprograms$metaprograms.genes
seurat_obj <- AddModuleScore_UCell(seurat_obj, 
                                   features = mp.genes, 
                                   assay="SCT", 
                                   ncores=4, 
                                   name = "")

VlnPlot(seurat_obj, 
        features=names(mp.genes), 
        group.by = "celltype",
        pt.size = 0, ncol=5)

# ==================================基于基因程序的特征对细胞进行分组===================================
## --------------- 提取MP1-MP6的评分列 -----------------
mp_cols <- paste0("MP", 1:6)  
mp_score_mat <- seurat_obj@meta.data[, mp_cols]

## ----------------按最高评分分配细胞群 ----------------
### 对每个细胞，找到评分最高的MP编号
mp_max <- apply(mp_score_mat, 1, function(x) {
  which.max(x)
})
### 转换为分组标签（MP1-MP6）
mp_group <- paste0("MP", mp_max)
seurat_obj$MP_group <- factor(mp_group, levels = paste0("MP", 1:6))  # 设定因子水平，保证顺序

## ----------------可选：处理“评分全0”的特殊细胞-----------------
zero_cells <- which(apply(mp_score_mat, 1, sum) == 0)
if (length(zero_cells) > 0) {
  seurat_obj$MP_group[zero_cells] <- "Unassigned"
  seu$MP_group <- factor(seu$MP_group, levels = c(paste0("MP", 1:6), "Unassigned"))
  message(paste("发现", length(zero_cells), "个无活性细胞，标记为Unassigned"))
}

## ------------------ 验证分群结果 ---------------------
cat("MP分组细胞数统计：\n")
print(table(seu$MP_group))
for (i in 1:6) {
  mp_name <- paste0("MP", i)
  p <- VlnPlot(
    seurat_obj,
    features = mp_name,
    group.by = "MP_group",
    pt.size = 0.1,
    cols = RColorBrewer::brewer.pal(7, "Set1")[1:7]
  ) +
    ggtitle(paste(mp_name, "评分在各MP_group中的分布")) +
    theme(plot.title = element_text(size = 14))
  print(p)
}

DimPlot(seurat_obj,group.by = "MP_group",reduction = 'tsne')
table(seurat_obj$MP_group,seurat_obj$subcelltype)


# ==================================整合空间分数的signature===========================================
matrix <- seurat_obj@meta.data[,names(mp.genes)]

#dimred <- scale(matrix)
dimred <- as.matrix(matrix)

colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
#New dim reduction
seurat_obj@reductions[["MPsignatures"]] <- new("DimReduc",
                                               cell.embeddings = dimred,
                                               assay.used = "RNA",
                                               key = "MP_",
                                               global = FALSE)
set.seed(1234)
seurat_obj <- RunUMAP(seurat_obj, reduction="MPsignatures", 
                      dims=1:length(seu@reductions[["MPsignatures"]]),
                      metric = "euclidean", 
                      reduction.name = "umap_MP")

FeaturePlot(seurat_obj, features = names(mp.genes), reduction = "umap_MP", ncol=5) &
  scale_color_viridis(option="B") &
   theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

a <- DimPlot(seurat_obj, reduction = "umap_MP", group.by = "celltype.l1", label=T) + theme(aspect.ratio = 1,
                                                            axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()) + ggtitle("Original cell types") + NoLegend()

b <- DimPlot(seurat_obj, reduction = "umap_MP", group.by = "donor", label=T) + theme(aspect.ratio = 1,
                                                            axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()) + ggtitle("Donor") + NoLegend()
a | b














