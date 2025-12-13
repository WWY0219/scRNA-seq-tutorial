# ============================================ Prepare Environment ===================================================
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
library(clustree)
library(cowplot)
library(stringr)
library(ggsci)
library(SCP)
library(pheatmap)
library(ggrepel)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
set.seed(1234)
list.files()
dir.create("../03.Output/")


# ============================================ Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")
Idents(seurat_obj) <- "celltype"
DimPlot(seurat_obj,reduction = 'umap',
        label = TRUE,pt.size = 0.5) +NoLegend()

# ============================================ 为WGCNA设置Seurat对象 ========================================
## WGCNA分析的时候会把信息储存在seurat对象的@misc槽中
## variable: 使用存储在Seurat对象的VariableFeatures中的基因
## fraction: 使用在整个数据集或每组细胞中表达的基因，由 group.by 指定
## custom: 使用在Custom 列表中指定的基因
## 一个seurat对象可以包含多个hdWGCNA实验对象

## V5版本需要这行代码，V4不需要
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",        # fraction(自动覆盖适合筛选）;variable(seurat_HVG);custom(自定义)
  fraction = 0.05,                 # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "celltype_1"             # the name of the hdWGCNA experiment
)
## !!!!手动指定要纳入 WGCNA 分析的基因列表!!!!
custom_genes <- c("CD3D", "CD3E", "CD4", "IL2", "IFNG", "TNF", "FOXP3")  #Eg
seurat_obj <- SetupForWGCNA(
  seurat_obj = seurat_obj,
  gene_select = "custom",        # 启用自定义基因筛选模式
  custom_genes = custom_genes,   # 传入自定义基因列表（关键参数！）
  wgcna_name = "celltype_1"      
)


# ============================================ 各组构建metacell ============================================
## metacells是由来自同一个生物样本的、相似细胞组成的小群体聚合而成的
## 该过程使用k最近邻(KNN)算法来识别相似细胞的群体，然后计算这些细胞的平均或总表达量，从而生成一个metacell基因表达矩阵
## <1 万细胞：k=10~20；>5 万细胞：k=25~50
seurat_obj <- MetacellsByGroups(
        seurat_obj = seurat_obj,
        group.by = c("celltype", "orig.ident"), # 指定seurat_obj@meta.data中要分组的列
        reduction = 'harmony',                  # 选择要执行KNN的降维
        k = 25,                                 # KNN：k值越大，元细胞数量越少，聚合程度越高
        max_shared = 10,                        # 两个metacell之间共享细胞的最大数目
        ident.group = 'celltype',               # 等价于设置metacell的active.ident
        min_cells = 100                         # 排除数量小于100的细胞亚群
)

## normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


# ============================================ 共表达网络分析 ============================================
## 设置表达式矩阵，使用hdWGCNA对目标细胞亚群进行共表达网络分析
seurat_obj <- SetDatExpr(
        seurat_obj,
        group_name = C("celltype_1",...),                     # the name of the group of interest in the group.by column
        group.by='celltype',                                  # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
        assay = 'RNA', 
        slot = 'data'                                         # using normalized data
)

## -------------------选择软阈值------------------------
### "unsigned" ：不考虑相关性的正负，仅用相关性的绝对值（适用于研究基因共表达的强弱，不关注调控方向，示例中使用此类型）
### "signed"   ：考虑相关性的正负（正相关为激活，负相关为抑制，适用于研究基因调控的方向性）
### "signed hybrid"：强调正相关，弱化负相关（常用作折中方案）
seurat_obj <- TestSoftPowers(
  seurat_obj,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  networkType = 'unsigned'                                    
)

### plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

### assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)

### WGCNA和hdWGCNA的一般原则是选择使尺度自由拓扑模型拟合度(Scale Free Topology Model Fit)大于或等于0.8的最低软阈值(soft power threshold)
### 在构建网络时，如果用户没有提供软阈值，程序会自动选择一个软阈值

##-------------------构建共表达网络------------------------
### construct co-expression network
seurat_obj <- ConstructNetwork(
        seurat_obj,
        soft_power = 4,           # 自定义软阈值
        min_power = 3,            # 自动选择软阈值时的最小阈值
        tom_outdir = "TOM",       # TOM矩阵的输出目录
        tom_name = 'Treg',        # TOM矩阵的文件名
        overwrite_tom = TRUE              # 允许覆盖已存在的同名文件
)
seurat_obj <- ConstructNetwork(
 
  
       
         
  consensus = FALSE,        # 是否构建共识网络（多数据集整合）
  overwrite_tom = FALSE,    # 是否覆盖已存在的TOM文件
  wgcna_name = NULL,        # hdWGCNA实验名称
  blocks = NULL,            # 基因分块（处理大量基因时）
  maxBlockSize = 30000,     # 每个块的最大基因数
  randomSeed = 12345,       # 随机种子（保证结果可重复）
  corType = "pearson",      # 相关性计算方法
  consensusQuantile = 0.3,  # 共识网络的分位数
  networkType = "signed",   # 网络类型
  TOMType = "signed",       # TOM矩阵类型
  TOMDenom = "min",         # TOM分母的计算方式
  scaleTOMs = TRUE,         # 是否缩放TOM矩阵
  calibrationQuantile = 0.95,# 校准分位数
  sampleForCalibration = TRUE,# 是否抽样校准TOM
  sampleForCalibrationFactor = 1000,# 校准抽样的因子
  useDiskCache = TRUE,      # 是否使用磁盘缓存存储中间结果
  chunkSize = NULL,         # 分块处理的块大小
  deepSplit = 4,            # 模块检测的深度分割参数
  pamStage = FALSE,         # 是否使用PAM优化模块
  detectCutHeight = 0.995,  # 模块检测的剪切高度
  minModuleSize = 50,       # 最小模块基因数
  mergeCutHeight = 0.2,     # 模块合并的剪切高度
  saveConsensusTOMs = TRUE, # 是否保存共识TOM矩阵
  ...                       # 其他传递给WGCNA的参数
)
### 可视化WGCNA树状图
### ！！！“灰色”模块由那些未被归入任何共表达模块的基因组成。对于所有下游分析和解释，应忽略灰色模块！！！
PlotDendrogram(seurat_obj, main='Treg hdWGCNA Dendrogram')

### 检查拓扑重叠矩阵(topoligcal overlap matrix，TOM)
TOM <- GetTOM(seurat_obj)
TOM













