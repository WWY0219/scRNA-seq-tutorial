#============= Load environment =================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
gc()
setwd("workplace")     #replace your workplace
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
library(DoubletFinder)
library(tidyverse)
library(ggsci)
library(pheatmap)
library(scCustomize)
set.seed(1234)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


#==================================Load Running script=============================
source("sc_qc.R")
source("sc_doublefinder.R")

#==================================Read 10X data==================================
##10X .h5 
seurat_obj_1 <- Read10X_h5("seurat_obj_1.h5",use.names = T,unique.features = T)                                    
seurat_obj_1 <- CreateSeuratObject(counts = seurat_obj_1, project = "seurat_obj_1", min.cells = 3, min.features = 200)
##10X 标准读取方式
seurat_obj_2 <- Read10X("../01.RawData/seurat_obj_2/")                                    
seurat_obj_2 <- CreateSeuratObject(counts = seurat_obj_2, project = "seurat_obj_2", min.cells = 3, min.features = 200)
ncol(seurat_obj_1)
ncol(seurat_obj_2)

#==================================Run sc_run.R==================================
seurat_obj_1_qc <- sc_qc(seurat_obj = seurat_obj_1)
seurat_obj_2_qc <- sc_qc(seurat_obj = seurat_obj_2)

##Process QC based on above result
###Eg.seurat_obj_1
seurat_obj_1_filtered <- subset(seurat_obj_1_qc, 
                                subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA < 40000 & 
                                percent.mt < 10 & 
                                percent.ribo1 < 12 & percent.ribo2 <18 &
                                percent.RBC < 1)
ncol(seurat_obj_1_filtered)
p_seurat_obj_1_filtered <- VlnPlot(seurat_obj_1_filtered, 
                                   features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo1","percent.ribo2","percent.RBC"), 
                                   ncol = 3,pt.size = 0)
p_seurat_obj_1_filtered
ggsave("../03.Output/seurat_obj_1/seurat_obj_1_qc_mt_violin_filtered.pdf", plot = p_seurat_obj_1_filtered, width = 10, height = 8, dpi = 300)   #Head QC effect

#==================================Run sc_run.R==================================
seurat_obj_1_db <- sc_doublefinder(seurat_obj= seurat_obj_1_filtered ,max.dim = 40 , max.pcs=40, pN = 0.25)


#==================================Merge singledata==============================
seurat_obj_merge <- merge(seurat_obj_1_db, y = c(seurat_obj_2_db, ...), add.cell.ids = c("seurat_obj_1", "seurat_obj_2",...), project = "seurat_obj_merge")
ncol(seurat_obj_merge)

#==================================Output==============================
##JoinLayers
seurat_obj_merge[["RNA"]] <- JoinLayers(seurat_obj_merge[["RNA"]])
## rds
saveRDS(seurat_obj_merge, "../03.Output/seurat_obj_merge_QC.rds")
## qs
qsave(seurat_obj_merge, "../03.Output/seurat_obj_merge_QC.qs")
