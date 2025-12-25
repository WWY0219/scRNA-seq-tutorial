# CellChat
为了预测显著的通讯，CellChat识别出每个细胞组中差异过表达的配体和受体。为了量化由这些信号基因介导的两个细胞组之间的通讯，CellChat将每个相互作用与一个概率值相关联。 后者是基于配体在一个细胞组中的平均表达值和受体在另一个细胞组中的平均表达值，以及它们的协同因子<br>
以下内容参考[知乎文章](https://zhuanlan.zhihu.com/p/717734779)<br>
## DESCRIPTION
### Version
CellChat V2
### Github Learning
* Github:<https://github.com/jinworks/CellChat>
* ZhiHu :<https://zhuanlan.zhihu.com/p/1894789522887250489>
## Usage
### 01. 加载包并准备环境
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
### 02. 准备CellChat输入数据
> CellChat需要两个输入文件：1.细胞的基因表达数据;2.细胞标签。<br>
> 基因表达数据矩阵，行为基因，列为细胞。<br>
> CellChat需要数据归一化数据。如果是count数据，使用normalizeData进行归一化。<br>
#### 输入需要分析的单细胞数据
```R
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)
DimPlot(seurat_obj, pt.size = 0.8,group.by = "celltype_major",label = T)
table(seurat_obj@meta.data$celltype_major)
```
#### 准备cellchat需要的文件类型
```R
data.input <- GetAssayData(seurat_obj, layer = 'data')           # normalized data matrix
meta <- seurat_obj@meta.data[,c("orig.ident","celltype_major")]  # your cellType
colnames(meta) <-  c("samples","labels")
table(meta$labels)
identical(rownames(meta),colnames(data.input))                   # 严格判断两个向量是否一致
```
#### 对细胞类型进行排序
> 这一步是有用的，后续亚群用的到，尽量同一大类的亚群分在一起
```R
celltype_order <- c(
  "Tumor", "Fibroblast", "SMC", "T/NK", 
  "LEC", "Monocyte", "Mastcell", "Neutrophil", 
  "Endothelial","Epithelial")
meta$labels <- factor(meta$labels ,levels = celltype_order)
table(meta$labels)
ordered_indices <- order(meta$labels)
meta <- meta[ordered_indices, ]
data.input <- data.input[, ordered_indices]
identical(rownames(meta),colnames(data.input))
```
#### 构建cellchat Object
```R
cellchat <- createCellChat(object = data.input, 
                           meta = meta, 
                           group.by = "labels")
levels(cellchat@idents)
```
### 03. 设置配体-受体相互作用数据库
展示CellChat数据库
```R
CellChatDB <- CellChatDB.human          # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)        # Display database category-All
dplyr::glimpse(CellChatDB$interaction)
```
> use CellChatDB.mouse if running on mouse data<br>
> 不建议把非蛋白质信号纳入分析，可能是由于在单细胞转录组数据中无法准确检测或量化、对应的信号分子不直接由基因编码，而是代谢产物、离子等、生物学机制复杂，且缺乏统一可靠的注释和数据库支持等原因<br>

#### 选取需要的数据库
* 除“非蛋白信号”外，使用所有CellChatDB数据进行细胞-细胞通信分析
```R
CellChatDB.use <- subsetDB(CellChatDB) 
cellchat@DB <- CellChatDB.use # 在构建的cellchat中设定需要使用的数据库
```
* 使用CellChatDB的中特定的数据库进行细胞-细胞通信分析
```R
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
```
* Only uses the Secreted Signaling from CellChatDB v1
```R
CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))
```
* 使用所有CellChatDB数据进行细胞-细胞通信分析
```R
CellChatDB.use <- CellChatDB 
```
> 不建议以这种方式使用它，因为CellChatDB v2包含“非蛋白信号”(即代谢和突触信号)。<br>

### 04. 预处理细胞-细胞通讯分析的表达数据
```R
cellchat <- subsetData(cellchat,features = NULL)            # This step is necessary even if using the whole database
future::plan("multisession", workers = 1)                   # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
```
> 默认情况下,cellchat使用object@data.signaling进行网络推断<br>
> 同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。<br>
> 原来是projectData，新版是smoothData函数

### 05. 细胞-细胞通信网络的推理
参数设定：‘triMean’会产生更少但更强的相互作用；而‘truncatedMean’方法中，当‘trim’参数值较小时（例如 ‘trim = 0.1或0.05’），会输出更多的相互作用，从而能够检测到较弱的信号传导活动
```R
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
```
> 该分析的关键参数是类型，即计算每个细胞组的平均基因表达的方法。默认情况下，type = “triMean”，产生较少但更强的交互。当设置 type = “truncatedMean” 时，应为trim分配一个值，从而产生更多交互。请详细检查上述计算每个细胞组平均基因表达的方法。<br>
> 使用的是投射到PPI网络的模式时候需要用FALSE。如果使用了raw data就需要设置为TRUE<br>
> min.cells是设置阈值，最小是需要10个细胞参与通讯推断(可以自定义)

* 如果所研究的信号没有被测到，可以采用如下函数进行探查，trim设为0.1或者0.05
``` R
computeAveExpr(cellchat, features = c("CXCL12","CXCR4"),type =  "truncatedMean",trim = 0.1)
```
* 如果发现修改参数之后所研究的信号被测到了，那就修改代码如下
```R
cellchat <- computeCommunProb(cellchat, type =  "truncatedMean",trim = 0.1,raw.use = FALSE) 
```

### 06. 在信号通路水平上推断细胞间通讯
CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平上的通信概率
```R
cellchat <- computeCommunProbPathway(cellchat)
```
> NB:推断的每个配体-受体对的细胞间通信网络和每个信号通路分别存储在槽'net'和'netP'中。<br>

数据提取，subsetCommunication函数，一般全部提取并保存
```R
df.net <- subsetCommunication(cellchat)
qsave(cellchat,"cellchat.qs")
save(df.net,file = "df.net.Rdata")
write.csv(df.net,"df.net.csv")
```
* 表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯
```R
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
```
```R
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
```
### 07. CellChat Visulization
CellChat可以可视化聚合的蜂窝间通信网络。circle展示互作数目<br>
* 计算聚合细胞-细胞通信网络
```R
cellchat <- aggregateNet(cellchat)
```
> 互作网络整合,可以设置soure和target，不设置就是默认全部<br>

可视化
```R
groupSize <- as.numeric(table(cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
```
#### Heatmap 展示互作数据
```R
pheatmap::pheatmap(cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")
```
#### 贝克图
指定顺序和指定颜色,生成颜色向量（例如使用彩虹色）,将颜色向量命名为矩阵的行名
```R
celltype_order <- c(...)
color.use <- rainbow(nrow(mat))
names(color.use) <- rownames(mat)
```
```R
mat <- as.data.frame(cellchat@net$weight)
mat <- mat[celltype_order,]                #行排序
mat <- mat[,celltype_order] %>% as.matrix()
```
```R
par(mfrow = c(5,4), xpd=TRUE,mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i],
                   color.use = color.use)
}
```
> 如果图片显示不全,需要考虑是不是重新设置mfrow

#### 层次结构图
展示pathways
```R
cellchat@netP$pathways
levels(cellchat@idents) 
pathways.show <- "CXCL"
```
* Hierarchy plot
```R
vertex.receiver = seq(1:9) # a numeric vector
netVisual_aggregate(cellchat, signaling = pathways.show,
                    vertex.receiver = vertex.receiver,layout= "hierarchy")
                  # vertex.size = groupSize)
```
> vertex.receiver定义层次图的左边细胞<br>

* circle plot
```R
netVisual_aggregate(cellchat, signaling = pathways.show,layout = "circle")
```
* Chord diagram
```R
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
```
* 分组弦图
```R
levels(cellchat@idents)
group.cellType <- c(rep("T/NK", 9), "B","VSMCs","endothelial","epithelial/cancer",
                    "fibroblasts","mast","myeloid","plasma","proliferative" )
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show,
                     group = group.cellType,
                     title.name = paste0(pathways.show, " signaling network"))
```
* heatmap
```R
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
```
> 行： 热图的行代表了不同的细胞类型，这些细胞作为信号的发送者。<br>
> 列： 热图的列代表了不同的细胞类型，这些细胞作为信号的接收者。<br>
> 颜色深浅： 热图中的颜色深浅表示了通讯概率的大小。颜色越深，表示通讯概率越高，这意味着发送方细胞和接收方细胞之间的信号传递越强。<br>
> 通信概率（Communication Prob.）： 右侧的颜色条是颜色映射的参考。图中的深红色表示较高的通讯概率（靠近0.05），浅色表示较低的通讯概率（靠近0或更低）。<br>
> 顶部数值范围 (0 - 0.3)： 显示不同细胞类型接受到总传入信号强度(incoming）<br>
> 右侧部分的数值范围 (0 - 0.4)：显示不同细胞类型发出的总传出信号强度(outcoming）<br>

### 计算配-受体对信号通路的贡献并可视化
* 计算配-受体对信号通路的贡献并可视化
```R
p1 <- netAnalysis_contribution(cellchat, signaling = pathways.show,
                               title =  pathways.show)                 # 展现对特定通路的贡献程度
p2 <- netAnalysis_contribution(cellchat, signaling = df.net$pathway_name)
cowplot::plot_grid(p1, p2, align = "h",ncol=2)
```
* 可视化由单个配体-受体对介导的细胞间通讯
```R
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show,
                                 geneLR.return = FALSE)
pairLR
LR.show <- pairLR[2,] # show one ligand-receptor pair
```
* Hierarchy plot
```R
vertex.receiver = seq(1,9) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  
                     pairLR.use = LR.show, 
                     vertex.receiver = vertex.receiver,
                     layout = "hierarchy")
```
* Circle plot
```
netVisual_individual(cellchat, signaling = pathways.show, 
                     pairLR.use = LR.show, layout = "circle")
```
### 基于配-受体结果进一步可视化
```R
levels(cellchat@idents) 
netVisual_bubble(cellchat, sources.use = seq(1:9), 
                 targets.use = c(13), remove.isolate = FALSE)
ggsave("bubbleplot_nont.pdf",width = 7,height = 20)
```
> **sources.use**是发出信号的细胞系,**target.use**是接受信号的细胞系<br>

* 还可以增加signaling参数用于展示特定的配受体
```R
cellchat@netP$pathways
netVisual_bubble(cellchat, sources.use = seq(1:9), 
                 targets.use = c(13), 
                 signaling = c("CXCL"),
                 remove.isolate = FALSE)
ggsave("bubbleplot2.pdf",width = 5,height = 10)
````
* 自定义signaling输入展示-所有通路汇总之后
```R
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CXCL","COLLAGEN","MK","CD99",
                                                      "LAMININ","APP","DESMOSOME"))
netVisual_bubble(cellchat, sources.use = c(1:9),
                 targets.use = c(13,14), 
                 pairLR.use = pairLR.use,
                 remove.isolate = TRUE)
ggsave("bubbleplot-LR.pdf",width = 5,height = 10)
```

* 可以通过增加下面的参数去设置X轴上的顺序
```
# sort.by.target = T
# sort.by.source = T
# sort.by.source = T, sort.by.target = T
# sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE
```
### 使用小提琴/点图绘制信号转导基因表达分布
CellChat可以使用Seurat包装函数plotGeneExpression绘制与L-R对或信号通路相关的信号转导基因的基因表达分布。
```R
plotGeneExpression(cellchat, signaling = "VEGF", 
                   enriched.only = TRUE, 
                   type = "violin")
```
> 该功能提供 “violin”、“dot”、“bar” 三种类型的可视化。<br>
> 或用户可以使用 extractEnrichedLR 提取与推断的 L-R 对或信号通路相关的信号转导基因，然后使用Seurat或其他软件包绘制基因表达。

### 计算并可视化网络中心性得分
```R
cellchat@netP$pathways
pathways.show <- "VEGF"
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, 
                                          slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, 
                                  signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)                                  
```
* 行（Sender, Receiver, Mediator, Influencer）：行表示在信号通路网络中，不同细胞类型扮演的角色：
> Sender：信号的发送者，即哪些细胞类型是主要的信号发出者。<br>
> Receiver：信号的接收者，即哪些细胞类型是主要的信号接收者。<br>
> Mediator：中介者，用于识别在信号传播路径中起中介作用的细胞群体。<br>
> Influencer：影响者，用于识别在整个网络中对信息传播影响最大的细胞群体。 列（不同的细胞类型）：列表示具体的细胞类型。颜色深浅（Importance）：颜色的深浅表示每个细胞类型在特定角色中的重要性。颜色越深，表示该细胞类型在这个角色中的重要性越高（例如信号传递的强度或频率越大）；颜色越浅，表示该细胞类型在这个角色中的重要性较低。<br>

####  在二维空间中可视化占优势的发送者(源)和接收者(目标)
```R
# 从所有信号通路对聚合细胞-细胞通信网络的信号作用分析
gg1 <- netAnalysis_signalingRole_scatter(cellchat);gg1
# 对特定细胞间通讯网络的信号作用分析
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL"));gg2
gg1 + gg2
```
#### 识别对某些细胞群的输出或输入信号贡献最大的信号
```R
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
# ht1
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
# ht2
# ht1 + ht2
# class(ht1)

# 特定的signaling
cellchat@netP$pathways
htout <- netAnalysis_signalingRole_heatmap(cellchat, 
                                        pattern = "outgoing",
                                        signaling = c("ICAM","TGFb"))
htout

htcome <- netAnalysis_signalingRole_heatmap(cellchat, 
                                        pattern = "incoming",
                                        signaling = c("ICAM","TGFb"))
htcome
```
#### 识别整体通信模式/以探索多种细胞类型和信号通路如何协调运作
```R
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

# 当输出模式的数量为5时，Cophenetic值和Silhouette值都开始突然下降。
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "outgoing", 
                                          k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")

# 当输出模式的数量为6时，Cophenetic值和Silhouette值都开始突然下降。
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, 
                                          pattern = "incoming", 
                                          k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")
```




