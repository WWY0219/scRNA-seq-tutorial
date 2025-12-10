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














