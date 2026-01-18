# ============================================ Prepare Environment ===================================================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(qs)
library(tibble)
library(DESeq2)
library(msigdbr) 
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(cowplot)
library(stringr)
library(ggsci)
library(SCOP)
library(pheatmap)
library(ggrepel)
library(gplots)
library(clusterProfiler)
library(org.Hs.eg.db)
library(singleseqgset)
set.seed(1234)
list.files()
dir.create("./07.GE/")


# ============================================ Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")

# ============================================ ALL Celltypes ===================================================
ids = bitr(seurat_obj.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') 
seurat_obj.markers = merge(seurat_obj.markers,ids,by.x='gene',by.y='SYMBOL')
View(seurat_obj.markers)

## 函数split()可以按照分组因子，把向量，矩阵和数据框进行适当的分组。
## 它的返回值是一个列表，代表分组变量每个水平的观测。
gcSample=split(sce.markers$ENTREZID, sce.markers$cluster) 

## -------------------------------KEGG----------------------------------
kegg <- compareCluster(gcSample,
  fun = "enrichKEGG",
  organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(kegg)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))

## --------------------------------GO------------------------------------
go <- compareCluster(gcSample,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)
p <- dotplot(go)
p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5, hjust = 0.5
))


# ============================================ 差异分析后的GO以及KEGG分析 ===================================================
exp_m <- GetAssayData(seurat_obj, layer = 'counts')
gene.all <- rownames(exp_m)
mt.gene  = grep('^MT-', gene.all, value = TRUE)                       # 线粒体基因
ribosome.gene = grep('^RPL|^RPS|^MRPS|^MRPL', gene.all, value = TRUE) # 核糖体基因
features = setdiff(gene.all, c(mt.gene, ribosome.gene))               # 去除不需要的基因
unique(seurat_obj$group)

## ---------------------------------------全部细胞的差异分析 ----------------------------------------
Idents(seurat_obj) <- seurat_obj$group
### 循环处理
process_cell_type <- function(cell_type, seurat_obj) {
  # 子集化为特定细胞类型
  cell_subset <- subset(seurat_obj, subset = celltype == cell_type)
  # 查找差异表达基因
  degs <- FindMarkers(cell_subset, 
                      logfc.threshold = 0.25,
                      min.pct = 0.1, 
                      only.pos = FALSE, 
                      features = features,
                      ident.1 = "group1", 
                      ident.2 = "group2") %>%
    mutate(gene = rownames(.))
  # 过滤符合条件的差异表达基因
  degs_filtered <- degs %>%
    filter(pct.1 > 0.1 & p_val_adj < 0.05) %>%
    filter(abs(avg_log2FC) > 0.5)
  # 返回过滤后的结果，并标记是哪个细胞类型
  return(list(cell_type = cell_type, degs_filtered = degs_filtered))
}

### 获取所有独特的细胞类型
cell_types <- unique(surat_obj$celltype)
cell_types
### 使用 lapply 对每个细胞类型进行处理
results <- lapply(cell_types, 
                  function(cell_type) process_cell_type(cell_type, seurat_obj))

### 组合所有细胞类型的差异表达基因数据
combined_results <- do.call(rbind, lapply(results, function(x) {
  x$degs_filtered %>% mutate(cell_type = x$cell_type)
}))
### 将最后一列改成cluster
colnames(combined_results)[colnames(combined_results) == "cell_type"] <- "cluster"
combined_results$cluster <- as.factor(combined_results$cluster)
## 保存差异结果
head(combined_results)
table(combined_results$cluster)
write.csv(combined_results,file = "./07.GE/combined_results.csv")
combined_results <- read.csv(file = "./07.GE/combined_results.csv")
### 火山图
jjVolcano(diffData = combined_results,
          tile.col = corrplot::COL2('PuOr', 15)[4:12],
          size  = 4,
          fontface = 'italic',
          base_size = 16,
          legend.position = c(0.9, 0.9),
          cluster.order = rev(unique(combined_results$cluster)))


## 富集分析需要挂上代理
Sys.setenv("http_proxy"="http://10.16.46.126:7890")
Sys.setenv("https_proxy"="http://10.16.46.126:7890") 
## 读取差异分析之后的结果
tex_cd8_degs_fil <- read.csv(file = "./data/tex_cd8_差异基因.csv")
head(tex_cd8_degs_fil)
ids_tex <- bitr(tex_cd8_degs_fil$gene, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)
head(ids_tex)
tex_cd8_degs_fil <- merge(tex_cd8_degs_fil, ids_tex, by.x = 'gene', by.y = 'SYMBOL')
tex_cd8_kegg <- enrichKEGG(gene = tex_cd8_degs_fil$ENTREZID.x, organism = "hsa", pvalueCutoff = 1)
write.csv(as.data.frame(tex_cd8_kegg@result), file = "./data/tex_cd8_kegg_result.csv", row.names = FALSE)
dotplot(tex_cd8_kegg, showCategory = 10, title = "KEGG Enrichment for Tex_CD8")

## GSEA
tex_cd8_degs_fil <- tex_cd8_degs_fil[order(tex_cd8_degs_fil$avg_log2FC, decreasing = TRUE),]
tex_markers_list <- setNames(as.numeric(tex_cd8_degs_fil$avg_log2FC), tex_cd8_degs_fil$ENTREZID.x)
tex_cd8_gsea_go <- gseGO(geneList = tex_markers_list, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05)
tex_cd8_gsea_go_arrange <- arrange(as.data.frame(tex_cd8_gsea_go@result), desc(abs(NES)))
write.csv(tex_cd8_gsea_go_arrange, file = "./data/tex_cd8_gsea_go_results.csv", row.names = FALSE)
# 定义配色
color <- c("#f7ca64", "#43a5bf", "#86c697", "#a670d6", "#ef998a")
# 绘制 tex_cd8 的 GSEA-GO 图
## 上调
# 提取上调通路
upregulated_pathways <- tex_cd8_gsea_go@result %>%
  filter(NES > 0) %>%                       # 筛选上调通路
  arrange(desc(NES)) %>%                    # 按NES降序排列
  slice(1:5)                                # 选择前5条通路
# 绘制前5条上调通路的GSEA曲线
gsekp1_tex <- gseaplot2(
  tex_cd8_gsea_go,                          # GSEA结果对象
  geneSetID = upregulated_pathways$ID,      # 上调通路的ID向量
  pvalue_table = F,                         # 是否显示p值表
  base_size = 12,                           # 字体大小
  color = color # 可选颜色
)
upregulated_pathways$Description
# 显示图形
gsekp1_tex








sce = sce[, Idents(sce) %in% c("FCGR3A+ Mono", "CD14+ Mono")] # 挑选细胞
deg=FindMarkers(object = sce, 
                ident.1 = 'FCGR3A+ Mono',
                ident.2 = 'CD14+ Mono', 
                test.use='MAST' )  ## MAST在单细胞领域较为常用
head(deg)
save(deg,file = 'deg-by-MAST-for-mono-2-cluster.Rdata')
##火山图
degdf <- deg
degdf$symbol <- rownames(deg)
logFC_t=0
P.Value_t = 1e-28
degdf$change = ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC < 0,"down",
                      ifelse(degdf$p_val_adj < P.Value_t & degdf$avg_log2FC > 0,"up","stable"))
ggplot(degdf, aes(avg_log2FC,  -log10(p_val_adj))) +
  geom_point(alpha=0.4, size=3.5, aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()

##功能注释
## 获取上下调基因
gene_up=rownames(deg[deg$avg_log2FC > 0,])
gene_down=rownames(deg[deg$avg_log2FC < 0,])
## 把SYMBOL改为ENTREZID
library(org.Hs.eg.db)
gene_up=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                   keys = gene_up,
                                                   columns = 'ENTREZID',
                                                   keytype = 'SYMBOL')[,2]))
gene_down=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                                     keys = gene_down,
                                                     columns = 'ENTREZID',
                                                     keytype = 'SYMBOL')[,2]))
library(clusterProfiler)
## 以上调基因为例，下调基因同理
## KEGG
gene_up <- unique(gene_up)
kk.up <- enrichKEGG(gene = gene_up,
                    organism = "hsa",
                    pvalueCutoff = 0.9,
                    qvalueCutoff = 0.9)
dotplot(kk.up)

## GO
go.up <- enrichGO(gene = gene_up,
                OrgDb = org.Hs.eg.db,
                ont = "BP" ,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.99,
                qvalueCutoff = 0.99,
                readabl = TRUE)
dotplot(go.up)

##差异分析后的GSEA
## 上一步差异分析得到差异基因列表deg后取出，p值和log2FC
nrDEG = deg[,c('avg_log2FC', 'p_val')]
colnames(nrDEG)=c('log2FoldChange','pvalue') ##更改列名
library(org.Hs.eg.db)
library(clusterProfiler)
## 把SYMBOL转换为ENTREZID，可能有部分丢失
gene <- bitr(rownames(nrDEG),     
             fromType = "SYMBOL",     
             toType =  "ENTREZID",    
             OrgDb = org.Hs.eg.db)
## 基因名、ENTREZID、logFC一一对应起来
gene$logFC <- nrDEG$log2FoldChange[match(gene$SYMBOL,rownames(nrDEG))]
## 构建genelist
geneList=gene$logFC
names(geneList)=gene$ENTREZID 
geneList=sort(geneList,decreasing = T) # 降序，按照logFC的值来排序
## GSEA分析
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.9,
                  verbose      = FALSE)
kk_gse=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
sortkk<-kk_gse[order(kk_gse$enrichmentScore, decreasing = T),]
library(enrichplot)
gseaplot2(kk_gse, 
          "hsa04510", 
          color = "firebrick",
          rel_heights=c(1, .2, .6))

## 展示排名前四的通路
gseaplot2(kk_gse, row.names(sortkk)[1:4])
## 把p值标上去
gseaplot2(kk_gse, 
          "hsa04510", 
          color = "firebrick",
          rel_heights=c(1, .2, .6),
          pvalue_table = TRUE)
















































