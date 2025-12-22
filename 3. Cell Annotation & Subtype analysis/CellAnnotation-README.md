# Cell Annotation Tutorial
单细胞分析中最重要的一步-细胞注释。<br>
同时我们会分析亚群分组。<br>

## Prepare Environment
```R
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
```
## Load scData with celltype_major
```R
seurat_obj <- readRDS("../03.Output/seurat_obj_harmony/seurat_obj_harmony.rds")
Idents(seurat_obj) <- "clusters_res0.5"
DimPlot(seurat_obj,  reduction = 'umap', 
        group.by="clusters_res0.5", label = TRUE, pt.size = 0.5) + NoLegend()
```
