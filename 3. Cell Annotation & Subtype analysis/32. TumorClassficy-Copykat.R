# =================================== Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")  #replace your workspace
getwd()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(infercnv)
library(CopyKAT)
library(dplyr)
library(ggsignif)          
library(circlize)
library(ComplexHeatmap)
set.seed(1234)

# =================================== Load scData with celltype_major  ===================================================
surat_obj <- qread("../seurat_obj.qs")
Idents(seurat_obj) <- "clusters_res0.5"
##绘制样本和cluster的umap图
