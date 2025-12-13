# ============================================ Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(qs)
library(CellChat)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(clustree)
library(cowplot)
library(stringr)
library(ggsci)
library(SCP)
library(pheatmap)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
set.seed(1234)
list.files()
dir.create("../03.Output/")


# ============================================ Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")
Idents(seurat_obj) <- "celltype"
DimPlot(seurat_obj,reduction = 'umap',
        label = TRUE,pt.size = 0.5) +NoLegend()

# ============================================ 为WGCNA设置Seurat对象 ========================================
## WGCNA分析的时候会把信息储存在seurat对象的@misc槽中
## variable: 使用存储在Seurat对象的VariableFeatures中的基因
## fraction: 使用在整个数据集或每组细胞中表达的基因，由 group.by 指定
## custom: 使用在Custom 列表中指定的基因
## 一个seurat对象可以包含多个hdWGCNA实验对象

## V5版本需要这行代码，V4不需要
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",        # fraction(自动覆盖适合筛选）;variable(seurat_HVG);custom(自定义)
  fraction = 0.05,                 # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "celltype_1"             # the name of the hdWGCNA experiment
)
## !!!!手动指定要纳入 WGCNA 分析的基因列表!!!!
custom_genes <- c("CD3D", "CD3E", "CD4", "IL2", "IFNG", "TNF", "FOXP3")
seurat_obj <- SetupForWGCNA(
  seurat_obj = seurat_obj,
  gene_select = "custom",        # 启用自定义基因筛选模式
  custom_genes = custom_genes,   # 传入自定义基因列表（关键参数！）
  wgcna_name = "CD4+T"           # hdWGCNA实验名称
)

















