# CytoTRACE2
CytoTRACE2是一种基于单细胞RNA测序数据预测细胞潜能类别及绝对发育潜力的计算方法。<br>
* 其定义的潜能类别将细胞分为多个层次：1. 根据细胞的分化潜能对细胞进行分类：totipotent(多谱系分化的全能)—pluripotent(多谱系分化的多能)—multipotent(谱系限制分化的多能)—oligopoten(谱系限制分化的寡能)—unipotent(谱系限制分化的单能)—differentiated(分化); 
* 该方法通过预测的潜能分数(范围0到1)提供连续发育潜力评估，其中0(differentiated)代表分化，1(totipotent)代表全能。<br>
* 其核心是一个新型可解释的深度学习框架，该框架基于覆盖发育全谱系的34个人类和小鼠单细胞数据集(涵盖24种组织类型)训练验证，可学习每个潜能类别的基因表达程序，并通过校准不同细胞发育起源的输出结果，实现跨数据集比较发育潜力.<br>
## DESCRIPTION
### Version
CytoTRACE2
### Learning
* Github :<https://github.com/jinworks/CellChat>
* ZhiHu-1:

## Usage
### 01. 加载包并准备环境
```R
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(CytoTRACE2)
library(patchwork)
library(ggplot2)
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
set.seed(1234)
# 查看工作路径下的文件
list.files()
dir.create("./CytoTRACE2/)
```

### 02. 加载每个数据集的cellChat对象，然后合并在一起
