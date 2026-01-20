# =================================== Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")  # replace your workspace
getwd()
## Loading Packages
library(qs)
library(tidyverse)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
set.seed(1234)
dir.create("../03.Output/Major-CellAnnotation/")

# =================================== Load scData with celltype_major  ===================================================
seurat_obj <- readRDS("../03.Output/seurat_obj_harmony/seurat_obj_harmony.rds")
Idents(seurat_obj) <- "clusters_res0.5"
DimPlot(seurat_obj,  reduction = 'umap', 
        group.by="clusters_res0.5", label = TRUE, pt.size = 0.5) + NoLegend()

# ===============================================Cellmarkers Figures======================================================
## AllcellMarkers-1
cell_markers <- list(MSC =c("CD34","PDGFRA"),                                 
                     immunocell= c("PTPRC"),                                  
                     T_cell = c("CD3E", "CD3D", "CD4", "CD8A"),
                     B_cell = c("MS4A1", "CD19", "CD79A"),
                     Monocyte = c("CD14", "CD68", "CSF1R", "FCGR3A", "LYZ","ITGAM"),
                     Dcs = c("CD1C", "CLEC9A", "HLA-DRA", "CD83"),
                     NK = c("NKG7", "GNLY", "KLRD1", "NCAM1"),                      
                     Mastcell = c("CPA3", "TPSAB1", "KIT"),
                     Fibroblast = c("COL1A1", "COL3A1", "THY1", "DCN", "FAP"),      
                     Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),     
                     SMC = c("ACTA2", "MYH11", "TAGLN"),                            
                     SMC_cancercell =c("ACTG2","PRLR","SFRP4"),
                     ESC =c("PLN","RGS5","SUSD2"),
                     Proliferating_cell = c("MKI67", "PCNA", "TOP2A", "CCNB1"),
                     Epithelial = c("EPCAM", "KRT8", "KRT18", "CDH1")
                    )

p_cellmarkers <- DotPlot(object = seurat_obj,
                         features = cell_markers,  
                         cols = c("grey", "red"),
                         cluster.idents = TRUE) +
RotatedAxis() +
theme( panel.border = element_rect(color = "black", fill = NA),  
       panel.spacing = unit(1, "mm"),
       strip.text = element_text(margin = margin(b = 3, unit = "mm")),
       strip.placement = 'outlet',
       axis.line = element_blank()
     ) +
labs(x = "", y = "")
print(p1)
ggsave("../03.Output/Major-CellAnnotation/seurat_obj_cellmarkers.pdf", plot = p_cellmarkers, width = 20, height = 12, dpi = 300)

## Cell-Markers2
major_markers <- list("T/NK" =c('CD3E','CD4','CD8A'),
                      B=c("CD79A","MS4A1"),
                      Plasma=c("MZB1"),
                      Epithelial = c('EPCAM',"KRT8","KRT18"),
                      "Mono/Macro"=c("CD14","CD68","CD163"),
                      Mast = c("KIT"),
                      Fibroblast = c("COL5A2",'PDGFRB'),
                      Pericyte = c("CSPG4","RGS5"),
                      Endothelial =c("PECAM1"),
                      Tumor=c("MKI67","FOXJ1","SOX2","SOX9")
                     ) 
plot3=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()

p_know <- plot2|plot3
ggsave("../03.Output/Major-CellAnnotation/seurat_obj_cellmarkers2.pdf", plot = p_know, width = 20, height = 12, dpi = 300)

# ===============================================Major Celltype Annotation===============================================
top10 <- seurat_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
## Doheatmap图
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()                            
## 小提琴图观察基因分布
VlnPlot(seurat_obj, features = top10$gene[1:5], pt.size=0)                           
marker2 <- FindMarkers(object = , ident.1 = 2)
marker_DEG <- FindMarkers(object = USOO, ident.1 = 10, ident.2 = 11,                 # 第二个cluster（对比组）
                                                      logfc.threshold = 0.25,        # 最小log2倍数变化（过滤微小差异）
                                                      min.pct = 0.1,                 # 基因在至少10%的细胞中表达（过滤低表达基因）
                                                      test.use = "wilcox",           # 统计检验方法（默认wilcox，适用于单细胞数据）
                                                      assay = "RNA" )                # 使用的assay（默认"RNA"）


## CellAnnotationsMethods--1
meta_supp = data.frame(seurat_cluster = 0:(length(unique(seurat_obj$cluster_res0.5)) - 1), celltype = NA)
meta_supp[meta_supp$seurat_cluster %in% c(0), 'celltype'] = 'B'
meta_supp[meta_supp$seurat_cluster %in% c(4), 'celltype'] = 'Plasma'
meta_supp[meta_supp$seurat_cluster %in% c(1), 'celltype'] = 'T/NK'
meta_supp[meta_supp$seurat_cluster %in% c(5), 'celltype'] = 'Stromal'
meta_supp[meta_supp$seurat_cluster %in% c(3), 'celltype'] = 'Monocyte'
meta_supp[meta_supp$seurat_cluster %in% c(2), 'celltype'] = 'Tumor'
meta_supp[meta_supp$seurat_cluster %in% c(2), 'celltype'] = 'db-like'

for (i in 1:nrow(meta_supp)) {
  seurat_obj@meta.data[which(seurat_obj$cluster_res0.5 == meta_supp$seurat_cluster[i]), 'celltype_major'] = meta_supp$celltype[i]
}
Idents(seurat_obj) <- 'celltype_major'   #replace your metadata@celltype-name
table(seurat_obj@meta.data$celltype_major,useNA = "always")                  #If NA, plot5 will error
table(seurat_obj@meta.data$celltype_major,
      seurat_obj@meta.data$celltype_major)

### Head annotation umap
plot4 <- DimPlot(seurat_obj, group.by = "celltype_major", label = T) & NoLegend()
plot5 <- DotPlot(object = seurat_obj,
                 features = known_markers,
                 scale=T,
                 group.by = "celltype_major") +
scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
theme_pubr()+
theme(axis.text.x = element_text(angle=90)) & NoLegend()

plot4 | plot5

qsave(seurat_obj, "../03.Output/Major-CellAnnotation/seurat_obj_beforeremovedblike.qs")

### remove db-like cell
seurat_obj <- subset(seurat_obj, subset = celltype_major != "db-like")


## ============!!!去除掉双细胞后需重新进行降维聚类！！！请重复上述操作！！！========================
### Auto Methods for re-Reduction
source("sc_harmony.R")
seurat_obj <- sc_harmony(seurat_obj= seurat_obj,
                         out_dir="")

### Manual Methods for re-Reduction
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

seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30, 
                      reduction.name = "umap", verbose = FALSE)
seurat_obj <- RunTSNE(seurat_obj, reduction = "harmony", dims = 1:30, 
                      reduction.name = "tsne", verbose = FALSE)
seurat_obj=FindNeighbors(seurat_obj,dims = 1:15,reduction = "harmony")
seurat_obj=FindClusters(seurat_obj,resolution = 0.5)
table(seurat_obj@meta.data$seurat_clusters)
table(seurat_obj@meta.data$orig.ident)
plot1 <- DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T,pt.size = 0.25)
plot2 <- DimPlot(seurat_obj,reduction = "umap",group.by = "seurat_clusters",label = T,pt.size = 0.25)+NoLegend()
plot1|plot2

### Re-Annotation

## save annotation
qsave(seurat_obj, file = 'seurat_obj_annotation.qs')


## CellAnnotationsMethods--2 (for choose)
new.cluster.ids <- c("0"= "Tumor cell",
                    "1"=" Endothelial_cell",
                    "2"= "Tumor cell",
                    "3"= "MSC",
                    ...)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
p_umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5, repel = TRUE,        
                  cols = NULL          #自定义：cols = c("细胞类型1"="red", ...)
) + 
  ggtitle("UMAP") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

print(p_umap)

# 保存图片
ggsave("../03.Output/USOO_celltype_umap.pdf", plot = p_umap,
  width = 10,
  height = 8,
  dpi = 300
)
seurat_obj@meta.data$celltype <- Idents(seurat_obj)



















