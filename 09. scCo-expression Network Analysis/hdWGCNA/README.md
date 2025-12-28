# hdWGCNA
WGCNA(Weighted Gene Co-expression Network Analysis)即加权基因共表达网络分析，能帮助我们找到各方面相似的基因模块(module)，探索基因网络与关注表型之间的关系，寻找核心基因。适用于组别较多样本量大的数据，推荐5组(或者15个样品)以上的数据.<br>
## DESCRIPTION
### Version
hdWGCNA 
### Learning
* Github :
* Document:<https://smorabit.github.io/hdWGCNA/>
* Zhihu-2:
## Usage
### 01. 加载包并准备环境
```R
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)
library(qs)
library(tidyverse)
library(cowplot)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("./hdWGCNA/)
enableWGCNAThreads(nThreads = 4)
```
将 ggplot2 的全局默认绘图主题设置为 cowplot 包提供的 theme_cowplot() 主题
```R
theme_set(theme_cowplot())
```
### 02. 加载每个数据集的cellChat对象，然后合并在一起
```R
seurat_obj <- readRDS('data/Zhou_2020.rds')

# 这个 Seurat 对象最初是使用 Seurat v4 创建的。如果你使用的是 Seurat v5, 可以运行下面的代码。
seurat_obj <- UpdateSeuratObject(seurat_obj)
seurat_obj
```
