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

# ===============================================Find Best Resolution===============================================
table(seurat_obj@meta.data$seurat_clusters)
## clustree method
p_clustree <- clustree(seurat_obj, prefix = "cluster_res") + coord_flip()
p_clustree 

## other methods (待补充）


## Determined your resolution
levels(Idents(seurat_obj))
Idents(seurat_obj) <- "cluster_res0.5"

# ===============================================Cellmarkers Figures===============================================
dir.create("../03.Output/cellmarkers_fig/")

## AllcellMarkers-1
cell_markers <- list(
  #间充质细胞
  MSC =c("CD34","PDGFRA"),
  #免疫细胞
  immunocell= c("PTPRC"),
  T_cell = c("CD3E", "CD3D", "CD4", "CD8A"),
  B_cell = c("MS4A1", "CD19", "CD79A"),
  Monocyte_Macrophage = c("CD14", "CD68", "CSF1R", "FCGR3A", "LYZ","ITGAM"),
  Dendritic_cell = c("CD1C", "CLEC9A", "HLA-DRA", "CD83"),
  NK_cell = c("NKG7", "GNLY", "KLRD1", "NCAM1"),  # NCAM1即CD56
  Mast_cell = c("CPA3", "TPSAB1", "KIT"),
  # 基质细胞
  Fibroblast = c("COL1A1", "COL3A1", "THY1", "DCN", "FAP"),  # THY1即CD90
  Endothelial_cell = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"),  # PECAM1即CD31
  Smooth_muscle_cell = c("ACTA2", "MYH11", "TAGLN"),
  SMC_cancercell =c("ACTG2","PRLR","SFRP4"),
  ESC =c("PLN","RGS5","SUSD2"),
  #增值细胞 
  Proliferating_cell = c("MKI67", "PCNA", "TOP2A", "CCNB1"),
  # 其他常见细胞
  Epithelial_cell = c("EPCAM", "KRT8", "KRT18", "CDH1")
)

p1 <- DotPlot(object = seurat_obj,
              features = cell_markers,  
              cols = c("grey", "red"),
              luster.idents = TRUE) +
RotatedAxis() +
theme( panel.border = element_rect(color = "black", fill = NA),  
       panel.spacing = unit(1, "mm"),
       strip.text = element_text(margin = margin(b = 3, unit = "mm")),
       strip.placement = 'outlet',
       axis.line = element_blank()
     ) +
labs(x = "", y = "")
print(p)
ggsave("../03.Output/usoooo_cellmarker_0.5.pdf", plot = p, width = 20, height = 12, dpi = 300)


top10 <- seurat_obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()                            # Doheatmap图
VlnPlot(seurat_obj, features = top10$gene[1:5], pt.size=0)                           #小提琴图观察基因分布
marker2 <- FindMarkers(object =USOO, ident.1 = 2)
marker_DEG <- FindMarkers(object = USOO, ident.1 = 10, ident.2 = 11,                 # 第二个cluster（对比组）
                                                      logfc.threshold = 0.25,        # 最小log2倍数变化（过滤微小差异）
                                                      min.pct = 0.1,                 # 基因在至少10%的细胞中表达（过滤低表达基因）
                                                      test.use = "wilcox",           # 统计检验方法（默认wilcox，适用于单细胞数据）
                                                      assay = "RNA"                   # 使用的assay（默认"RNA"）)



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


ggsave("../03.Output/clustree.pdf",width =30 ,height =30,dpi =300)
saveRDS(USOO,"../03.Output/USOO_harmony/USOO_harmony.rds")







                        
  
