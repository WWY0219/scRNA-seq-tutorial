# CellChat
为了预测显著的通讯，CellChat识别出每个细胞组中差异过表达的配体和受体。为了量化由这些信号基因介导的两个细胞组之间的通讯，CellChat将每个相互作用与一个概率值相关联。 后者是基于配体在一个细胞组中的平均表达值和受体在另一个细胞组中的平均表达值，以及它们的协同因子<br>
以下内容参考[知乎文章](https://zhuanlan.zhihu.com/p/717734779)<br>
## DESCRIPTION
### Version
CellChat V2
### Github Learning
* Github
## Usage
### Prepare Environment
```R
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
library(SCP)
set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("./CellChat/")
```
### Load scData for cellchat
```R
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)
DimPlot(seurat_obj, pt.size = 0.8,group.by = "celltype_major",label = T)
table(seurat_obj@meta.data$celltype_major)
```
### Prepare scData for cellchat
> CellChat需要两个用户输入:一个是细胞的基因表达数据，另一个是用户分配的细胞标签。<br>
> 对于基因表达数据矩阵，行为基因，列为细胞。<br>
> CellChat需要数据归一化数据。如果是count数据，使用normalizeData进行归一化。<br>
```R
data.input <- GetAssayData(seurat_obj, layer = 'data')           # normalized data matrix
meta <- seurat_obj@meta.data[,c("orig.ident","celltype_major")]  # your cellType
colnames(meta) <-  c("samples","labels")
table(meta$labels)
identical(rownames(meta),colnames(data.input))                   # 严格判断两个向量是否一致
celltype_order <- c(
  "Tumor", "Fibroblast", "SMC", "T/NK", 
  "LEC", "Monocyte", "Mastcell", "Neutrophil", 
  "Endothelial","Epithelial")
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)
ordered_indices <- order(meta$labels)
## 对 meta 和 data.input 进行排序
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))
```
