# ============================================ Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(qs)
library(tidyverse)
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
library(magrittr)
library(WGCNA)
library(hdWGCNA)
library(igraph)
library(JASPAR2020)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(JASPAR2024)
library(RSQLite)
library(EnsDb.Hsapiens.v86)
library(pheatmap)
library(ggrepel)
set.seed(1234)
list.files()
dir.create("14-hdWGCNA")

# ============================================ Load Data ===================================================
seurat_obj <- qread("hdWGCNA_object.qs")

# ============================================ 转录因子分析 ===================================================
## JASPAR 2020
pfm_core <- TFBSTools::getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", 
              tax_group = 'vertebrates', 
              all_versions = FALSE)
)

# JASPAR 2024 (not used for this tutorial)
# JASPAR2024 <- JASPAR2024()
# sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
# pfm_core <- TFBSTools::getMatrixSet(
#   x = sq24,
#   opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
# )

# 进行motif分析
seurat_obj <- MotifScan(
  seurat_obj,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)

# 获取motif df:
motif_df <- GetMotifs(seurat_obj)

# 保留所有TFs, 并去除灰色模块基因
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(seurat_obj)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
seurat_obj <- SetWGCNAGenes(seurat_obj, genes_use)
seurat_obj <- SetDatExpr(seurat_obj, group.by = 'celltype', group_name='Treg')

# define model params:
model_params <- list(
    objective = 'reg:squarederror',
    max_depth = 1,
    eta = 0.1,
    nthread=16,
    alpha=0.5
)

# 构建转录因子网络
seurat_obj <- ConstructTFNetwork(seurat_obj, model_params)
results <- GetTFNetwork(seurat_obj)
head(results)


# 策略“A”为每个基因选择前10个TF
# 策略“B”选择每个转录因子的顶级基因
# 策略“C”保留所有高于一定调控分数的TF-基因对
seurat_obj <- AssignTFRegulons(
    seurat_obj,
    strategy = "A", # 还有B和C策略
    reg_thresh = 0.01,
    n_tfs = 10
)

# 可视化
# 根据基因表达分为与TF正相关(右侧)或负相关(左侧)的目标基因
p1 <- RegulonBarPlot(seurat_obj, selected_tf='ZNF652')
p2 <- RegulonBarPlot(seurat_obj, selected_tf='NFKB2', cutoff=0.15)
p1 | p2



##============================================获取正/负共表达的TF-基因对======================================
# 正向regulons
seurat_obj <- RegulonScores(
    seurat_obj,
    target_type = 'positive',
    ncores=8
)

# 负向regulons
seurat_obj <- RegulonScores(
    seurat_obj,
    target_type = 'negative',
    cor_thresh = -0.05,
    ncores=8
)

# 获取数据结果:
pos_regulon_scores <- GetRegulonScores(seurat_obj, target_type='positive')
neg_regulon_scores <- GetRegulonScores(seurat_obj, target_type='negative')

# 选择感兴趣的TF
cur_tf <- 'FOXP3'

# 把regulon分数添加到Seurat metadata
seurat_obj$pos_regulon_score <- pos_regulon_scores[,cur_tf]
seurat_obj$neg_regulon_score <- neg_regulon_scores[,cur_tf]

# plot using FeaturePlot
p1 <- FeaturePlot(seurat_obj, feature=cur_tf) + umap_theme()
p2 <- FeaturePlot(seurat_obj, feature='pos_regulon_score', cols=c('lightgrey', 'red')) + umap_theme()
p3 <- FeaturePlot(seurat_obj, feature='neg_regulon_score', cols=c('lightgrey', 'seagreen')) + umap_theme()

p1 | p2 | p3
# select TF of interest
cur_tf <- 'FOXP3'

# plot with default settings
p <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf)
p

# plot the FOXP3 network with primary, secondary, and tertiary targets
p1 <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf, depth=1, no_labels=TRUE)
p2 <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf, depth=2, no_labels=TRUE)
p3 <- TFNetworkPlot(seurat_obj, selected_tfs=cur_tf, depth=3, no_labels=TRUE)

p1 | p2 | p3
