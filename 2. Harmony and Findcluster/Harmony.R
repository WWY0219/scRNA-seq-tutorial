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
library(RColorBrewer)
library(DoubletFinder)
set.seed(1234)


# ============================Load scData after QC-tutorial=================
seurat_obj <- readRDS("../03.Output/seurat_obj_merge_QC.rds")
seurat_obj[["RNA"]]=JoinLayers(seurat_obj[["RNA"]])

# ============================Load sc_run===================================
source("sc_harmony.R")

# ============================Load scData after QC-tutorial=================
##用细胞总 UMI 计数的中位数作为缩放因子消除细胞间测序差异
seurat_obj <- NormalizeData(seurat_obj, normalization.method ="LogNormalize", 
                            scale.factor = median(seurat_obj@meta.data$nCount_RNA ))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000) 
                            
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = USOO))

USOO <- FindNeighbors(USOO, dims = 1:40)
USOO <- FindClusters(USOO, resolution = 1)
# 计算细胞周期评分
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
                            
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
# 回归掉不感兴趣的变量
seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo","percent_mt","percent_RBC"))
# 使用HVG去跑PCA
seurat_obj=RunPCA(seurat_obj)
ElbowPlot(seurat_obj,ndims = 50)
seurat_obj=RunHarmony(seurat_obj,group.by.vars="orig.ident",max_iter=20)
seurat_obj=RunUMAP(seurat_obj,dims=1:6,verbose = T,reduction = "harmony")
seurat_obj=FindNeighbors(seurat_obj,dims = 1:6,reduction = "harmony")
seurat_obj=FindClusters(seurat_obj,resolution = 0.1)
table(seurat_obj@meta.data$seurat_clusters)
plot1=DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T)
plot2=DimPlot(seurat_obj,reduction = "umap",group.by = "seurat_clusters",label = T)+NoLegend()
plot1|plot2
# ==== 细胞注释 ====
known_markers=list("T/NK"=c('CD3E','CD4','CD8A'),
                    B=c("CD79A","MS4A1"),
                    Plasma=c("MZB1"),
                    Epithelium=c('EPCAM',"KRT8","KRT18"),
                    Tumor=c("MKI67","FOXJ1","SOX2","SOX9"),
                    "Mono/Macro"=c("CD14","CD68","CD163"),
                    Mast=c("KIT"),
                    Fibroblast=c("COL5A2",'PDGFRB'),
                    Pericyte=c("CSPG4","RGS5"),
                    Endothelium=c("PECAM1")) 
seq <- seq(0.5, 2, by = 0.5)
                            seq
                            USOO <- sc_harmony(seurat_object = USOO,
                                               max.dim =40,
                                               max.dim_harmony =40, 
                                               resolutions = seq,   
                                               max.iter = 10, 
                                               top_n = 10,
                                               out_dir = "../03.Output/") 
                            top10 <- USOO.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
                            DoHeatmap(USOO, features = top10$gene) + NoLegend()                           # Doheatmap图
                            VlnPlot(USOO, features = top10$gene[1:5], pt.size=0)                         #小提琴图观察基因分布
                            marker2 <- FindMarkers(object =USOO, ident.1 = 2)
                            marker_DEG <- FindMarkers(object = USOO, ident.1 = 10, ident.2 = 11,           # 第二个cluster（对比组）
                                                      logfc.threshold = 0.25,  # 最小log2倍数变化（过滤微小差异）
                                                      min.pct = 0.1,           # 基因在至少10%的细胞中表达（过滤低表达基因）
                                                      test.use = "wilcox",     # 统计检验方法（默认wilcox，适用于单细胞数据）
                                                      assay = "RNA"            # 使用的assay（默认"RNA"）
                            )
                            ####寻找最适Resolution
                            ##clustree法
                            p1<- clustree(USOO,prefix = "clusters_res")+coord_flip()
                            p1
                            ggsave("../03.Output/clustree.pdf",width =30 ,height =30,dpi =300)
                            saveRDS(USOO,"../03.Output/USOO_harmony/USOO_harmony.rds")







                        
  
