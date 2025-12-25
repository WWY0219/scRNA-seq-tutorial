# CellChat
Cellchat能够从scRNA-seq数据中定量推断和分析细胞间通讯网络，其预测细胞的主要信号的输入及输出，以及这些细胞和信号是如何通过网络分析和模式识别方法进行协调从而执行功能。Cellchat中配体-受体分析主要基于质量作用定律模型量化两个细胞群体之间的信号通讯概率，其中的质量作用定律模型是指形成复合物的速率（即信号强度）与配体和受体的浓度有关.
以下内容参考[知乎文章][https://zhuanlan.zhihu.com/p/717734779](https://zhuanlan.zhihu.com/p/554129661)<br>
以下内容参考[知乎文章2](https://zhuanlan.zhihu.com/p/1894789522887250489)
## DESCRIPTION
### Environment
Conda::R432
### Download
下载地址：<https://www.cog-genomics.org/plink/>
Github: <https://github.com/jinworks/CellChat>
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
library(ComplexHeatmap)
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
#### 不同细胞类型之间的相互作用数量或相互作用强度不同
为了简化复杂的网络并深入了解细胞类型级别的细胞间通信，我们可以根据定义的细胞组聚合细胞间通信。这里我们将细胞群分为三种细胞类型，然后重新合并 CellChat 对象的列表。
```R
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
```
然后，我们可以显示每个数据集中任意两种细胞类型之间的相互作用数量或相互作用强度
```R
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
```
同样，我们还可以使用圆形图显示任何两种细胞类型之间相互作用的差异数量或相互作用强度。与第一个数据集相比，第二个数据集中的红色（或蓝色）边缘表示增加（或减少）的信号
```R
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)
```
#### 比较二维空间中的主要源和目标
比较 2D 空间中的传出和传入交互强度可以很容易地识别在不同数据集之间发送或接收信号发生显着变化的细胞群
```R
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
```
从散点图中，我们可以看到与 NL 相比，Inflam.DC 和 cDC1 成为 LS 中的主要来源和目标之一。成纤维细胞群也成为 LS 中发出信号的主要来源。 此外，我们可以识别 NL 和 LS 之间 Inflam.DC 和 cDC1 的特定信号变化。## 识别与一个细胞组相关的信号变化。
```R
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Inflam. DC", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "cDC1", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))
```
### Part Ⅱ：识别保守和特定环境的信号通路
CellChat 可以根据它们在多种生物条件下的细胞间通信网络识别具有较大（或较小）差异的信号网络、信号组以及保守和特定于上下文的信号通路
#### 根据功能/结构相似性识别具有较大（或较小）差异的信号网络以及信号组
CellChat 基于其功能和拓扑相似性对推断的通信网络执行联合流形学习和分类。注意：这种分析适用于两个以上的数据集。 功能相似性：高度的功能相似性表明主要的发送者和接收者相似，可以解释为两条信号通路或两个配体-受体对表现出相似和/或冗余的作用。注意：功能相似性分析不适用于具有不同细胞类型组成的多个数据集。 结构相似性：结构相似性用于比较它们的信令网络结构，不考虑发送者和接收者的相似性。注意：结构相似性分析适用于具有相同细胞类型组成或细胞类型组成截然不同的多个数据集。 在这里，我们可以基于功能相似性运行流形和分类学习分析，因为两个数据集具有相同的细胞类型组成
#### 根据功能相似性识别信号组
```R
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
```
```R
# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
```
#### 根据结构相似性识别信号组
```R
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
```
```R
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2
```
#### 计算并可视化学习的联合流形中的路径距离
我们可以根据共享二维空间中的欧几里德距离来识别差异较大（或较小）的信令网络。较大的距离意味着两个数据集之间的通信网络在功能或结构相似性方面的差异较大。注意：我们只计算两个数据集之间重叠信号通路的距离。那些仅在一个数据集中识别的信号通路在此不予考虑。如果有超过三个数据集，可以通过comparison在函数中定义来进行成对比较rankSimilarity。
```R
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2
```
#### 识别和可视化保守和特定环境的信号通路
通过比较每个信号通路的信息流/相互作用强度，我们可以识别信号通路，（i）关闭，（ii）减少，（iii）打开或（iv）增加，通过在一个条件下改变它们的信息流为与另一种情况相比。
#### 比较每个信号通路的整体信号流
我们可以通过简单地比较每个信号通路的信息流来识别保守和特定于上下文的信号通路，该信息流由推断网络中所有细胞组对之间的通信概率之和（即网络中的总权重）定义）。 此条形图可以以堆叠模式绘制，也可以不绘制。根据 NL 和 LS 皮肤之间推断网络内的整体信息流的差异，对重要的信号通路进行排序。红色的顶部信号通路在 NL 皮肤中富集，而这些绿色在 LS 皮肤中富集。
```R
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
```
#### 比较每个细胞群相关的传出（或传入）信号
上述分析将来自传出和传入信令的信息汇总在一起。我们还可以比较两个数据集之间的传出（或传入）信号模式，从而识别表现出不同信号模式的信号通路/配体受体。 我们可以组合来自不同数据集的所有已识别信号通路，从而通过将传出和传入信号聚合在一起，将它们并排比较，包括传出信号、传入信号和整体信号。注意：rankNet也显示了整体信号的比较，但它没有显示特定细胞群中的信号强度
```R
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
```
```R
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
```
```R
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
```
### Part Ⅲ：识别上调和下调的信号配体-受体对
通过比较通信概率来识别功能失调的信号
> 可以比较一些细胞群与其他细胞群的配体-受体对介导的通信概率。这可以通过`comparison`在函数中设置来完成`netVisual_bubble`<br>

```R
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object
```
### Part Ⅳ：使用层次图、圈图或者弦图直观地比较细胞间通信
与单个数据集的 CellChat 分析类似，我们可以使用层次图、圆图或弦图来可视化细胞间通信网络。 边缘颜色/权重，节点颜色/大小/形状：在所有可视化图中，边缘颜色与作为发送者的源一致，边缘权重与交互强度成正比。较粗的边缘线表示较强的信号。在Hierarchy plot 和 Circle plot中，圆圈大小与每个单元组中的单元数成正比。在层次图中，实心圆和空心圆分别代表源和目标。在和弦图中，内部较细的条颜色表示从相应外部条接收信号的目标。内部条的大小与目标接收到的信号强度成正比。
```R
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
```
```R
# Chord diagram
pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
```
### Part Ⅴ：比较不同数据集之间地信号基因表达分布
我们可以使用 Seurat 包装函数绘制与 LR 对或信号通路相关的信号基因的基因表达分布`plotGeneExpression`
```R
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
```
