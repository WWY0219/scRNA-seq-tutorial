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
library(WGCNA)
library(hdWGCNA)
library(enrichR)
library(GeneOverlap)
library(pheatmap)
library(ggrepel)
set.seed(1234)
list.files()
dir.create("../03.Output/")

# ============================================ Load Data ===================================================
seurat_obj <- qread("hdWGCNA_object.qs")

# ============================================ GeneEnrichment ============================================
## 定义enrichr databases
dbs <- c('GO_Biological_Process_2023',
         'GO_Cellular_Component_2023',
         'GO_Molecular_Function_2023')

### 富集分析
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs,
  max_genes = 100           # use max_genes = Inf to choose all genes
)

### 检索输出表
enrich_df <- GetEnrichrTable(seurat_obj)

### 查看结果
head(enrich_df)

### make GO term plots:
EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots",                           # name of output directory
  n_terms = 10,                                       # number of enriched terms to show (sometimes more are shown if there are ties)
  plot_size = c(5,7),                                 # width, height of the output .pdfs
  logscale=TRUE                                       # do you want to show the enrichment as a log scale?
)

### enrichr dotplot
EnrichrDotPlot(
  seurat_obj,
  mods = "all", # use all modules (default)
  database = "GO_Biological_Process_2023",            # this must match one of the dbs used previously
  n_terms=2,                                          # number of terms per module
  term_size=8,                                        # font size for the terms
  p_adj = FALSE                                       # show the p-val or adjusted p-val?
)  + scale_color_stepsn(colors=rev(viridis::magma(256)))

# ============================================ GSEA ============================================
library(fgsea)
### load the GO Biological Pathways file (downloaded from EnrichR website)
pathways <- fgsea::gmtPathways('GO_Biological_Process_2023.txt')

### optionally, remove the GO term ID from the pathway names to make the downstream plots look cleaner
names(pathways) <- stringr::str_replace(names(pathways), " \\s*\\([^\\)]+\\)", "")


# get the modules table and remove grey genes
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# rank by Treg_NEW1 genes only by kME
cur_mod <- 'Treg_NEW1'
modules <- GetModules(seurat_obj) %>% subset(module == cur_mod)
cur_genes <- modules[,(c('gene_name', 'module', paste0('kME_', cur_mod)))]
ranks <- cur_genes$kME; names(ranks) <- cur_genes$gene_name
ranks <- ranks[order(ranks)]

# run fgsea to compute enrichments
gsea_df2 <- fgsea::fgsea(
  pathways = pathways, 
  stats = ranks,
  minSize = 3,
  maxSize = 500
)

# 可视化
top_pathways <- gsea_df2 %>% 
    subset(pval < 0.05) %>% 
    slice_max(order_by=NES, n=25) %>% 
    .$pathway

plotGseaTable(
    pathways[top_pathways], 
    ranks, 
    gsea_df2, 
    gseaParam=0.5,
    colwidths = c(10, 4, 1, 1, 1)
)

# name of the pathway to plot 
selected_pathway <- 'Cellular Respiration'
plotEnrichment(
    pathways[[selected_pathway]],
    ranks
) + labs(title=selected_pathway)
