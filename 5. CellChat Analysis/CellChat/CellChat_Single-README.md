# CellChat
为了预测显著的通讯，CellChat识别出每个细胞组中差异过表达的配体和受体。为了量化由这些信号基因介导的两个细胞组之间的通讯，CellChat将每个相互作用与一个概率值相关联。 后者是基于配体在一个细胞组中的平均表达值和受体在另一个细胞组中的平均表达值，以及它们的协同因子<br>
以下内容参考[知乎文章](https://zhuanlan.zhihu.com/p/717734779)<br>
## DESCRIPTION
### Version
CellChat V2
### Github Learning
* Github
* ZhiHu:<https://zhuanlan.zhihu.com/p/1894789522887250489>
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
### Prepare scData for cellchat
> CellChat需要两个用户输入:一个是细胞的基因表达数据，另一个是用户分配的细胞标签。<br>
> 对于基因表达数据矩阵，行为基因，列为细胞。<br>
> CellChat需要数据归一化数据。如果是count数据，使用normalizeData进行归一化。<br>
### Load scData for cellchat
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
### 设置配体-受体相互作用数据库
研究者不建议把非蛋白质信号纳入分析，这可能是由于在单细胞转录组数据中无法准确检测或量化、对应的信号分子不直接由基因编码，而是代谢产物、离子等、生物学机制复杂，且缺乏统一可靠的注释和数据库支持等原因
```R
CellChatDB <- CellChatDB.human          # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)        # Display database category-All
dplyr::glimpse(CellChatDB$interaction)
```
> use CellChatDB.mouse if running on mouse data<br>

#### 选取需要的数据库
```R
# 使用CellChatDB的中特定的数据库进行细胞-细胞通信分析
# 示例中使用了Secreted Signaling
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# 除“非蛋白信号”外，使用所有CellChatDB数据进行细胞-细胞通信分析
CellChatDB.use <- subsetDB(CellChatDB)

# 使用所有CellChatDB数据进行细胞-细胞通信分析
# 研究者不建议以这种方式使用它，因为CellChatDB v2包含“非蛋白信号”(即代谢和突触信号)。
# CellChatDB.use <- CellChatDB 

# 在构建的cellchat中设定需要使用的数据库
cellchat@DB <- CellChatDB.use
```
### 预处理细胞-细胞通讯分析的表达数据
```R
cellchat <- subsetData(cellchat)                            # This step is necessary even if using the whole databa
future::plan("multisession", workers = 1) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#默认情况下,cellchat使用object@data.signaling进行网络推断
#同时也提供了projectData函数,通过扩散过程基于高置信度实验验证的蛋白质互作网络中的邻近节点对基因表达值进行平滑处理。该功能在处理测序深度较浅的单细胞数据时尤为有用，因其能减少信号基因（特别是配体/受体亚基可能存在的零表达）的dropout效应。不担心其可能在扩散过程引入伪影，因其仅会引发极微弱的通讯信号。
# 原来是projectData，新版是smoothData函数
cellchat <- smoothData(cellchat, adj = PPI.human)
```
### 细胞-细胞通信网络的推理
参数设定：‘triMean’会产生更少但更强的相互作用；而‘truncatedMean’方法中，当‘trim’参数值较小时（例如 ‘trim = 0.1或0.05’），会输出更多的相互作用，从而能够检测到较弱的信号传导活动
```R
# 该分析的关键参数是类型，即计算每个细胞组的平均基因表达的方法。默认情况下，type = “triMean”，产生较少但更强的交互。当设置 type = “truncatedMean” 时，应为trim分配一个值，从而产生更多交互。请详细检查上述计算每个细胞组平均基因表达的方法。
# 使用的是投射到PPI网络的模式时候需要用FALSE。如果使用了raw data就需要设置为TRUE
cellchat <- computeCommunProb(cellchat, type = "triMean",raw.use = FALSE) 

# 如果所研究的信号没有被测到，可以采用如下函数进行探查，trim设为0.1或者0.05
# computeAveExpr(cellchat, features = c("CXCL12","CXCR4"),type =  "truncatedMean",trim = 0.1)
# 如果发现修改参数之后所研究的信号被测到了，那就修改代码如下
# cellchat <- computeCommunProb(cellchat, type =  "truncatedMean",trim = 0.1,raw.use = FALSE) 

# min.cells是设置阈值，最小是需要10个细胞参与通讯推断(可以自定义)
cellchat <- filterCommunication(cellchat, min.cells = 10)
```

### 在信号通路水平上推断细胞间通讯
```R
# CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平上的通信概率。 
# NB:推断的每个配体-受体对的细胞间通信网络和每个信号通路分别存储在槽'net'和'netP'中。
cellchat <- computeCommunProbPathway(cellchat)

#数据提取，subsetCommunication函数，一般全部提取并保存
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) #表示从细胞群 1 和 2 向细胞群 4 和 5 推断出的细胞间通讯。
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
df.net <- subsetCommunication(cellchat)

qsave(cellchat,"cellchat.qs")
save(df.net,file = "df.net.Rdata")
```

write.csv(df.net,"df.net.csv")
