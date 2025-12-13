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
set.seed(1234)
list.files()
dir.create("../03.Output/")


# ============================================ Load Data ===================================================
seurat_obj <- qread("seurat_obj.qs")


# ============================================ 首先观察各个细胞亚群的相关性（验证性工作） ===================================================
av <-AggregateExpression(seurat_obj,
                         group.by = c("orig.ident","celltype"),
                         assays = "RNA",
                         slot = "counts",                              # counts适用于DEseq2；data适用于可视化、相对表达量比较                       
                         return.seurat = FALSE)                        # 返回总的计数 

## colnames为orig.ident_celltype, rownames为gene symbol
av=as.data.frame(av[[1]])
head(av)[1:3,1:3] 
write.csv(av,file = 'AverageExpression_seurat_obj.csv')

## 找到sd最显著的排名前1000基因
cg <- names(tail(sort(apply(log(av+1), 
                            1,              # 1为按行计算；2为按列计算
                            sd)),
                 1000)) 

## 提前这些基因矩阵并进行log处理后计算不同样本对应的细胞的相关性值
df =cor(as.matrix(log(av[cg,]+1)))
colnames(df)
ac=as.data.frame(str_split(colnames(df),'_',simplify = T))
rownames(ac)=colnames(df)
colnames(ac)[1:2] <- c("orig.ident","celltype")

## seurat_obj 进行分组
ac$group = ifelse(grepl('GSM4942397|GSM4942396|GSM4942398|GSM4942399|GSM5023319',ac$orig.ident),"control","tumor")
table(ac$group)
head(ac)
## -------------------------------------------------visulazition---------------------------------------------------------
#### 未对行列按照group进行聚类
pheatmap::pheatmap(df ,
                   show_colnames = F,
                   show_rownames = F,
                   annotation_col = ac,
                   annotation_colors = list(
                     group = c(LM = "red", MM = "green"),
                     celltype = c(`FC-C0` = "orange", `SMC-C1` = "purple", `FC-C2` = "cyan")
                   )) 
dev.off()

#### 对行列按照group进行聚类
my_group_order <- c("LM", "MM")  
ac$group <- factor(ac$group, levels = my_group_order)
sample_order <- order(ac$group) 
df_ordered <- df[sample_order, sample_order] # 表达矩阵按样本顺序调整
ac_ordered <- ac[sample_order, ] # 注释数据按样本顺序调整

pheatmap(
  mat = df_ordered,
  show_colnames = FALSE,               # 不显示列名
  show_rownames = FALSE,               # 不显示行名
  annotation_col = ac_ordered,         # 列注释
  color = colorRampPalette(brewer.pal(9, "OrRd"))(50), 
  cluster_cols = FALSE,                # 关闭列聚类（已按group排序）
  cluster_rows = T,                    # 保留行聚类
  treeheight_col = 0,                  # 隐藏列聚类树
  treeheight_row = 20,                 # 调整行聚类树高度
  border_color = T                     # “NA” 去掉色块边框
)
dev.off()


# ============================================ 对于每个细胞亚群进行主成分分析 ===================================================
av <-AggregateExpression(seurat_obj,
                         group.by = c("orig.ident","celltype"),
                         assays = "RNA",
                         slot = "counts",                                                  
                         return.seurat = FALSE)                        

av=as.data.frame(av[[1]])
df=log(av +1) 
head(ac)
celltp = unique(ac$celltype);celltp

## -------------------------------------------------visulazition---------------------------------------------------------
library("FactoMineR") 
library("factoextra")  
library(ggstatsplot)
pca_list <- lapply(celltp, function(x){
  x
  exp <- df[,rownames(ac[ac$celltype==x,])]  
  cg=names(tail(sort(apply(exp, 1, sd)),1000)) 
  exp=exp[cg,]
  dat.pca <- PCA(as.data.frame(t(exp)) , graph = FALSE)
  group_list=ac[ac$celltype==x,'group']
  this_title <- paste0(x,'_PCA')
  p2 <- fviz_pca_ind(dat.pca,
                     geom.ind = "point",   # show points only (nbut not "text")
                     col.ind = group_list, # color by groups
                     palette = "Dark2",
                     addEllipses = TRUE, # Concentration ellipses
                     legend.title = "Groups")+
    ggtitle(this_title)+
    theme_ggstatsplot()+
    theme(plot.title = element_text(size=12,hjust = 0.5))
  
  p2
})
wrap_plots(pca_list, byrow = T, nrow = 2 )




# ============================================ DESeq2-Analysis ===================================================
av <-AggregateExpression(seurat_obj,
                         group.by = c("orig.ident","group"),
                         assays = "RNA",
                         slot = "counts",
                         return.seurat = FALSE)  # 返回总的计数 
av=as.data.frame(av[[1]])
head(av)[1:3,1:3]       # 可以看到是整数矩阵

## ----------------------------------------------- DESeq2 ---------------------------------------------------
library(tibble)
library(DESeq2)
### Get counts matrix
counts_res <- av
head(counts_res)

### generate sample level metadata
colData <- data.frame(samples = colnames(counts_res))
colData <- colData %>%
  mutate(condition = ifelse(grepl('Normal', samples), 'Normal', 'Tumor')) %>%
  column_to_rownames(var = 'samples')

### !!Create DESeq2 object!! 
dds <- DESeqDataSetFromMatrix(countData = counts_res,
                              colData = colData,
                              design = ~ condition) # condition 表示差异分析将基于colData的condition 变量进行

### filter Counts >= 10 
keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

### !!run DESeq2!!
dds <- DESeq(dds)
resultsNames(dds)                  # Check the coefficients for the comparison

### Generate results object
res <- results(dds, name = "condition_Normal_vs_Tumor")
res
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG_deseq2 = na.omit(DEG)

#添加上下调信息
DEG_deseq2 <- DEG_deseq2 %>%
  mutate(Type = if_else(padj > 0.05, "stable",
                        if_else(abs(log2FoldChange) < 1, "stable",
                                if_else(log2FoldChange >= 1, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>% 
  rownames_to_column("Symbol")

top_genes <- DEG_deseq2 %>%
  filter(Type %in% c("up", "down")) %>% # 只选显著基因
  group_by(Type) %>%
  slice_head(n = 5) %>% # 每个分组选Top10
  ungroup()


ggplot(DEG_deseq2, aes(log2FoldChange,-log10(padj))) +
  geom_point(size = 3.5, alpha = 0.8,
             aes(color = Type),show.legend = T)  +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) +
  xlim(-10, 10) +
  labs(x = "Log2(fold change)", y = "-log10(padj)") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'black',lwd=0.8) + 
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'black',lwd=0.8)+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_label_repel(
    data = top_genes,                        # 仅标注top基因
    aes(label = Symbol),                     # 基因名（Symbol列）
    size = 3,                                # 标签字体大小
    fill = "white",                          # 标签背景色
    alpha = 0.8,                             # 背景透明度
    label.padding = unit(0.2, "lines"),      # 标签内边距
    max.overlaps = Inf,                      # 允许所有标签显示（即使重叠，也会尽量分开）
    box.padding = unit(0.3, "lines"),        # 标签与点的距离
    point.padding = unit(0.4, "lines"),      # 点与标签框的距离
    color = "black"                          # 标签字体颜色
  )

