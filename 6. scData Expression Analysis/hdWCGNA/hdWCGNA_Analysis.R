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
  wgcna_name = "celltype_1"        # the name of the hdWGCNA experiment
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
### "unsigned" ：不考虑相关性的正负，仅用相关性的绝对值（适用于研究基因共表达的强弱，不关注调控方向）
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
        overwrite_tom = TRUE,     # 允许覆盖已存在的同名文件
        consensus = FALSE,        # 是否构建共识网络（多数据集整合）
        overwrite_tom = FALSE,    # 是否覆盖已存在的TOM文件
        blocks = NULL,            # 基因分块（处理大量基因时）
        maxBlockSize = 30000,     # 每个块的最大基因数
        randomSeed = 12345,       # 随机种子（保证结果可重复）
        corType = "pearson",      # 相关性计算方法（"pearson"/"spearman"）
        consensusQuantile = 0.3,  # 共识网络的分位数
        networkType = "signed",   # 网络类型（"signed"/"unsigned"/"signed hybrid"）
        TOMType = "signed",       # TOM矩阵类型
        TOMDenom = "min",         # TOM分母的计算方式（"min" or "mean")
        scaleTOMs = TRUE,         # 是否缩放TOM矩阵
        calibrationQuantile =0.95,# 校准分位数
        sampleForCalibration=TRUE,# 是否抽样校准TOM
  sampleForCalibrationFactor=1000,# 校准抽样的因子
        chunkSize = NULL,         # 分块处理的块大小
        deepSplit = 4,            # 模块检测的深度分割参数
        pamStage = FALSE,         # 是否使用PAM优化模块
        detectCutHeight = 0.995,  # 模块检测的剪切高度
        minModuleSize = 50,       # 最小模块基因数
        mergeCutHeight = 0.2,     # 模块合并的剪切高度
        saveConsensusTOMs = TRUE, # 是否保存共识TOM矩阵
)

### 可视化WGCNA树状图
### ！！！“灰色”模块由那些未被归入任何共表达模块的基因组成。对于所有下游分析和解释，应忽略灰色模块！！！
PlotDendrogram(seurat_obj, main='Celltye_1 hdWGCNA Dendrogram')

### 检查拓扑重叠矩阵(topoligcal overlap matrix，TOM)
TOM <- GetTOM(seurat_obj)
TOM


##-------------------计算模块特征基因------------------------
### 模块特征基因(Module Eigengenes，MEs)是用于总结整个共表达模块基因表达谱的常用指标
### 模块特征基因是通过在每个模块的基因表达矩阵子集上执行主成分(PCA)分析来计算的
### 这些PCA矩阵的第一个主成分就是MEs。此外对MEs应用Harmony批量校正，从而得到harmony后的模块特征基因(hMEs）
### 需要先运行ScaleData，否则harmony会报错
# seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

### 计算完整单细胞数据集中的所有MEs
seurat_obj <- ModuleEigengenes(
        seurat_obj,
        group.by.vars="orig.ident"
)

### 协调模特征基因:允许用户对MEs应用Harmony批量校正，生成协调模块特征基因(hMEs)
hMEs <- GetMEs(seurat_obj)

### 提取未经过批次校正的模块特征基因（ME）
#MEs <- GetMEs(seurat_obj, harmonized=FALSE)


##-------------------计算模块连接性------------------------
### 在共表达网络分析中，通常希望关注“枢纽基因”，即在每个模块内高度连接的基因。因此，希望确定每个基因的基于特征基因(eigengene)的连接性，也称为kME
### hdWGCNA提供了ModuleConnectivity 函数，用于在完整的单细胞数据集(而不是metacell数据集)中计算基因的kME值。这个函数本质上是计算基因与模块特征基因之间的成对相关性
### 虽然可以在整个数据集中计算所有细胞的kME，但建议在之前用于运行ConstructNetwork的细胞类型或分组中计算kME

### 计算基于特征基因的连接性(kME)：关注枢纽基因
seurat_obj <- ModuleConnectivity(
        seurat_obj,
        group.by = 'celltype', 
        group_name = 'celltype_1')

### 模块重命名
seurat_obj <- ResetModuleNames(
        seurat_obj,
        new_name = "celltype_NEW"
)

### 绘制每个模块按kME排序的基因
p <- PlotKMEs(seurat_obj, ncol=4)
p


##-------------------获取模块内部信息------------------------
### 获取模块内部信息:这一步去除了不需要的灰色模块基因
modules <- GetModules(seurat_obj) %>% 
  subset(module != 'grey')
head(modules[,1:6])

### 得到枢纽基因(可以提取按kME排序的前N个枢纽基因的表格,这里选择了10)
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
head(hub_df)

### 保存数据
qsave(seurat_obj, 'hdWGCNA_object.qs')


##-------------------计算hub基因siganture得分------------------------
### 计算每个模块前25个枢纽基因的kME得分
library(UCell)
seurat_obj <- ModuleExprScore(
        seurat_obj,
        n_genes = 25,
        method='UCell'                # "Seurat" or "UCell"
        #wgcna_name = NULL,
)



# ============================================ Visulization ============================================
## -------------------模块特征图------------------
### "hMEs"（默认）：批次校正后的模块特征基因（harmonized MEs），适用于有批次差异的数据，展示模块的生物学表达模式
### "MEs"：原始模块特征基因（未校正批次），适用于无批次差异的数据
### "kME"：基因的模块内连通性（intramodular connectivity），需额外指定基因（较少用）
### "module_score"：模块的平均表达分数（类似 Seurat 的AddModuleScore结果）

### 每个模块制作hMEs的特征图
plot_list <- ModuleFeaturePlot(
        seurat_obj,
        reduction = "umap",
        features='hMEs',      # 要可视化的模块特征类型
        order=TRUE            # order so the points with highest hMEs are on top
        seurat_obj,           # 输入的Seurat对象
        module_names = NULL,  # 要可视化的模块名称列表
        wgcna_name = NULL,    # hdWGCNA实验名称
        order_points = TRUE,  # 是否按特征值大小排序点（使高值点在上）
        restrict_range =TRUE, # 是否限制颜色范围（减少极端值影响）
        point_size = 0.5,     # 点的大小
        alpha = 1,            # 点的透明度
        label_legend = FALSE, # 是否为图例添加标签
        ucell = FALSE,        # 是否使用UCell分数可视化
        raster = FALSE,       # 是否光栅化绘图（减少内存占用）
        raster_dpi = 500,     # 光栅化的DPI
        raster_scale = 1,     # 光栅化的缩放比例
        plot_ratio = 1,       # 多模块绘图时的行列比例
        title = TRUE          # 是否显示标题
)
wrap_plots(plot_list, ncol=4)


## -------------------相同函数绘制hub基因特征得分------------------
### 每个模块制作hub scores的特征图
plot_list <- ModuleFeaturePlot(
        seurat_obj,
        features='scores', # plot the hub gene scores
        order='shuffle',   # order so cells are shuffled
        ucell = TRUE       # depending on Seurat vs UCell for gene scoring
)
wrap_plots(plot_list, ncol=4)


## -------------------每个模块在不同细胞亚群中的情况每个模块在不同细胞亚群中的情况------------------
### 每个模块在不同样本中的情况
seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$orig.ident), ' '))[,1]
ModuleRadarPlot(
        seurat_obj,
        group.by = 'cluster',
        barcodes = seurat_obj@meta.data %>% 
        subset(celltype == 'celltype_1') %>% 
        rownames(),
        axis.label.size=4,
        grid.label.size=4
)


## -------------------查看模块相关图------------------
ModuleCorrelogram(seurat_obj)


## -------------------气泡图------------------
### get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

### add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

### plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'celltype')
### flip the x/y axes, rotate the axis labels, and change color scheme
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
p


## -------------------单模块的网络图------------------
# 使用ModuleNetworkPlot可视化每个模块前50(数值可自定)的hub gene
ModuleNetworkPlot(
    seurat_obj, 
    outdir='ModuleNetworks',       # new folder name
    n_inner = 20,                  # number of genes in inner ring
    n_outer = 30,                  # number of genes in outer ring
    n_conns = Inf,                 # show all of the connections
    plot_size=c(10,10),            # larger plotting area
    vertex.label.cex=1             # font size
)


## -------------------结合hub基因的网络图------------------
options(future.globals.maxSize = 5 * 1024^3)  # 5GB
graphics.off()                                # 关闭绘图设备
### hubgene network(基因数可自定)
HubGeneNetworkPlot(
        seurat_obj,
        n_hubs = 2, 
        n_other=2,
        edge_prop = 0.75,
        mods = 'all'
)

### 可以选择模块数
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)
### get the list of modules
modules <- GetModules(seurat_obj)
mods <- levels(modules$module)
mods <- mods[mods != 'grey']

### hubgene network
HubGeneNetworkPlot(
        seurat_obj,
        n_hubs = 2, 
        n_other= 2,
        edge_prop = 0.75,
        mods = mods[1:5]    # only select 5 modules
)


## -------------------UMAP共表达网络------------------
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,       # number of hub genes to include for the UMAP embedding
  n_neighbors=15,    # neighbors parameter for UMAP
  min_dist=0.1       # min distance between points in UMAP space
)

### get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,  # color each point by WGCNA module
   size=umap_df$kME*2    # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1,        # proportion of edges to sample (20% here)
  label_hubs=2 ,        # how many hub genes to plot per module?
  keep_grey_edges=FALSE
)
