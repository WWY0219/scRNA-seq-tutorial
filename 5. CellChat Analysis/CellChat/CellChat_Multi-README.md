# CellChat
Cellchat能够从scRNA-seq数据中定量推断和分析细胞间通讯网络，其预测细胞的主要信号的输入及输出，以及这些细胞和信号是如何通过网络分析和模式识别方法进行协调从而执行功能。Cellchat中配体-受体分析主要基于质量作用定律模型量化两个细胞群体之间的信号通讯概率，其中的质量作用定律模型是指形成复合物的速率（即信号强度）与配体和受体的浓度有关.
以下内容参考[知乎文章]([https://zhuanlan.zhihu.com/p/717734779](https://zhuanlan.zhihu.com/p/554129661))<br>
以下内容参考[知乎文章2](https://zhuanlan.zhihu.com/p/1894789522887250489)
## DESCRIPTION
### Environment
Conda::R432
### Download
下载地址：<https://www.cog-genomics.org/plink/>
## Usage
### 加载库
```R
library(Seurat)
library(dplyr)
library(SeuratData)
library(patchwork) #最强大的拼图包
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
options(stringsAsFactors = FALSE)
rm(list=ls()) #清空所有变量
options(stringsAsFactors = F) #输入数据不自动转换成因子（防止数据格式错误）
setwd("/home/rstudio/project/scar_data_analysis/")
```
### 加载每个数据集的cellChat对象，然后合并在一起
> 需要在每个数据集中单独运行CellChat，然后将不同的CellChat对象合并在一起<br>

```R
cellchat.NL <- readRDS("/Users/jinsuoqin/Documents/CellChat/tutorial/cellchat_humanSkin_NL.rds")
cellchat.LS <- readRDS("/Users/jinsuoqin/Documents/CellChat/tutorial/cellchat_humanSkin_LS.rds")
object.list <- list(NL = cellchat.NL, LS = cellchat.LS)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
```
### Part Ⅰ：预测细胞间通讯的一般原理
CellChat 从全局出发，预测细胞间通信的一般原理。在比较多种生物学条件下的细胞间通讯时，它可以回答以下生物学问题：<br>
* 细胞间通讯是否增强；
* 哪些细胞类型之间的相互作用发生了显着变化；
* 主要来源和目标如何从一种情况变为另一种情况；
#### 比较交互总数和交互强度
```R
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
```
#### 比较不同细胞群之间的相互作用数量和相互作用强度
```R
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
```
```R
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
```
差异网络分析仅适用于成对数据集。如果有更多的数据集进行比较，我们可以直接显示每个数据集中任意两个细胞群之间的相互作用数量或相互作用强度。 为了更好地控制跨不同数据集的推断网络的节点大小和边缘权重，我们计算每个细胞组的最大细胞数和所有数据集的最大交互数（或交互权重）
```R
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
```
