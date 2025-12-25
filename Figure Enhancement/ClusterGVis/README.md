# ClusterGVis
Plink 是用于遗传学数据分析的开源工具软件，特别适用于大规模的基因组关联研究（Genome-Wide Association Studies, GWAS）和基因型数据分析。<br>
其主要功能包括质量控制、关联测试、群体结构分析、基因组数据格式转换以及基因型数据的基本统计分析。<br>
以下内容参考[知乎文章](https://zhuanlan.zhihu.com/p/717734779)<br>
## DESCRIPTION
### Version
ClusterGVis v0.1.4
### Document
下载地址：<https://www.cog-genomics.org/plink/>
## Usage
### 读入
```R
library(Seurat)
library(qs)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(circlize)
library(ClusterGVis)
```
