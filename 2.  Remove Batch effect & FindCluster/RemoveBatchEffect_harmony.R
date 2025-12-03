# ============================ 加载配置和工作环境 ============================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
list.files()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(harmony)
library(RColorBrewer)
set.seed(1234)


# ===============================================Load scData after QC-tutorial===================================================
seurat_obj <- readRDS("../03.Output/seurat_obj_merge_QC.rds")
seurat_obj[["RNA"]]=JoinLayers(seurat_obj[["RNA"]])
print(seurat_obj)
table(seurat_obj@meta.data$orig.ident)

# ===============================================Load scData after QC-tutorial===============================================
## 用细胞总 UMI 计数的中位数作为缩放因子消除细胞间测序差异
seurat_obj <- NormalizeData(seurat_obj, normalization.method ="LogNormalize", 
                            scale.factor = median(seurat_obj@meta.data$nCount_RNA))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000) 

## 计算细胞周期评分
cc.genes.updated.2019 <- cc.genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes                        
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes) 

## 回归掉不感兴趣的变量
seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo1",
                                                        "percent_ribo2","percent_mt","percent_RBC"))
## 使用HVG去跑PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj,ndims = 50) 

## Perform Harmony batch correction
seurat_obj <- seurat_obj %>% RunHarmony(
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony",    
  plot_convergence = TRUE,
  max.iter = 30,
  verbose = FALSE
)
## Save Harmony convergence plot
ElbowPlot(seurat_obj, reduction = "harmony", ndims = 50) + 
  theme_minimal() + 
  ggtitle("Elbow Plot for Harmony-Corrected Dimensions")
ggsave(file.path(out_dir, paste0(obj_name, "_Harmony_convergence.pdf")), width = 8, height = 6, dpi = 300)


# ===============================================Run sc_resolutionfinder.R===============================================
## Loading sc_resolutionfinder.R
source("sc_resolutionfinder.R")
res <- seq(0.1, 1, by = 0.2)
res
seurat_obj <- sc_resolutionfinder(seurat_obj = seurat_obj,
                                  max_dim_pca = 30,
                                  max_dim_harmony = 30,
                                  resolutions = 0.5, 
                                  top_n = 20,
                                  out_dir = "../03.Output/")
##Output profiles-seurat_obj only contains resolution data, not containing Reduction Data(UMAP/tSNE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30, 
                                reduction.name = "umap", verbose = FALSE)
seurat_obj <- RunTSNE(seurat_obj, reduction = "harmony", dims = 1:30, 
                                reduction.name = "tsne", verbose = FALSE)
qsave(seurat_obj, "../03.Output/seurat_obj_harmony.qs")


# ===============================================Find Best Resolution===============================================
table(seurat_obj@meta.data$seurat_clusters)
## clustree method
p_clustree <- clustree(seurat_obj, prefix = "cluster_res") + coord_flip()
p_clustree 
ggsave("../03.Output/clustree.pdf",width =30 ,height =30,dpi =300)

## other methods (待补充）


## Determined your resolution
levels(Idents(seurat_obj))
Idents(seurat_obj) <- "cluster_res0.5"
plot1=DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T)
plot2=DimPlot(seurat_obj,reduction = "umap",group.by = "cluster_res0.5",label = T)+NoLegend()
plot = plot1+plot2
ggsave("../03.Output/USOO_umap.pdf",width = 20, height = 15,dpi =300)




