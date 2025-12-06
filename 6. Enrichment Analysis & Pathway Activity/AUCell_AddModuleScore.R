
rm(list = ls())
library(Seurat)
library(tidyverse)
library(openxlsx)
load("step1.final.Rdata") #pbmc数据
sce <- step1.final
#check一下
DimPlot(sce,group.by = "celltype",label = T)+ NoLegend() +ggsci::scale_color_d3()
autophagy_genes <- read.xlsx("Autophagy.xlsx",colNames = T) 
g <- autophagy_genes$Symbol
head(autophagy_genes)
