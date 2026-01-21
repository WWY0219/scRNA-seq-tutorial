# =========================================================== Prepare Environment ===================================================
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
library(patchwork)
library(RColorBrewer)
library(RcppML)
library(readr)
library(irGSEA)
library(ggsci)
library(doMC)
library(ggunchull)
set.seed(1234)
dir.create("./irGSVA/")

# ================================================ Load scData with celltype_major  ==================================================
seurat_obj <- qread("./03.Output/seurat_obj_harmony/seurat_obj_harmony.qs")
DimPlot(seurat_obj,  reduction = 'umap', 
        group.by="clusters_res0.5", label = TRUE, pt.size = 0.5) + NoLegend()

# ===============================================irGSEA通过内置的msigdbr包进行打分======================================================
## Hallmark基因集打分
seurat_obj <- irGSEA.score(object = seurat_obj,
                           assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 1,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T,                                # 指定是否使用自定义基因集
                           species = "Homo sapiens", category = "H",  
                           subcategory = NULL, geneid = "symbol",
                           method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL,                           # 用于设置 AUCell 和 UCell 方法的最大排名，这里设置为 NULL，表示使用默认值
                           kcdf = 'Gaussian')                                                     # 设置核密度估计的类型，这里使用高斯（Gaussian）分布
Seurat::Assays(seurat_obj)

## KEGG基因集打分 
seurat_obj <- irGSEA.score(object = seurat_obj,
                           assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 1,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T, 
                           species = "Homo sapiens", category = "C2",  
                           subcategory = "CP:KEGG", geneid = "symbol",
                           method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')
Seurat::Assays(seurat_obj)

## GO-BP基因集打分 
seurat_obj <- irGSEA.score(object = seurat_obj,
                           assay = "RNA", slot = "data", 
                           seeds = 123, ncores = 1,
                           min.cells = 3, min.feature = 0,
                           custom = F, geneset = NULL, msigdb = T, 
                           species = "Homo sapiens", category = "C5",  
                           subcategory = "GO:BP", geneid = "symbol",
                           method = c("AUCell","UCell","singscore","ssgsea", "JASMINE", "viper"),
                           aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                           kcdf = 'Gaussian')
Seurat::Assays(seurat_obj)


# ===============================================================整合差异基因集================================================================
## 如果报错，考虑加句代码options(future.globals.maxSize = 100000 * 1024^5)
result.dge <- irGSEA.integrate(object = seurat_obj,
                               group.by = "celltype",
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea", "JASMINE", "viper"))

## 查看RRA识别的在多种打分方法中都普遍认可的差异基因集
geneset.show <- result.dge$RRA %>% 
  dplyr::filter(pvalue <= 0.05) %>% 
  dplyr::pull(Name) %>% unique(.)


# =============================================================== Vislization =============================================================
## -----------------------Heatmap-----------------------------
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA", 
                                      top = 10,
                                      show.geneset = NULL)
heatmap.plot
### 从'RRA"换成其他的单独分析方法，这里尝试使用了“ssgsea”
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "ssgsea”
                                      top = 10,
                                      show.geneset = NULL)
heatmap.plot
### geneset.show。这种方式是可以指定在当前method下展示具有统计学意义的通路
heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA", #从'RRA"换成“ssgsea”
                                      top = 10,
                                      show.geneset = geneset.show)
heatmap.plot

## -----------------------BubblePlot-------------------------
bubble.plot <- irGSEA.bubble(object = result.dge,
                                    method = "RRA",
                                    top = 10,
                                    show.geneset = geneset.show)
bubble.plot

## -----------------------Upset图-------------------------
upset.plot <- irGSEA.upset(object = result.dge,
                                  method = "RRA")
upset.plot

## -----------------------堆叠条形图-------------------------
barplot.plot <- irGSEA.barplot(object = result.dge,
                                      method = c("AUCell","UCell","singscore",
                                                 "ssgsea", "JASMINE", "viper"))
barplot.plot

## -----------------------小提琴图-------------------------
Idents(seurat_obj) <- seurat_obj$celltype
halfvlnplot <- irGSEA.halfvlnplot(object = seurat_obj,
                                  method = "AUCell",
                                  show.geneset = "HALLMARK-NOTCH-SIGNALING")
halfvlnplot

vlnplot <- irGSEA.vlnplot(object = seurat_obj,
                          method = c("AUCell", "UCell", 
                                     "singscore", "ssgsea", 
                                     "JASMINE", "viper"),
                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE")
vlnplot

#山峦图
ridgeplot <- irGSEA.ridgeplot(object = sc_dataset2,
                              method = "AUCell",
                              show.geneset = "HALLMARK-NOTCH-SIGNALING")
ridgeplot

#密度热图
densityheatmap <- irGSEA.densityheatmap(object = sc_dataset2,
                                        method = "AUCell",
                                        show.geneset = "HALLMARK-NOTCH-SIGNALING")
densityheatmap
