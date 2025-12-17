# ============================== Prepare Environment ===============================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("workspace")
getwd()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(scRNAtoolVis)
library(CytoTRACE2)
library(monocle)
library(ggpubr)
library(argparse)
library(tidyverse)
library(igraph)   # Version-2.0.3  !!!
set.seed(1234)
# 查看工作路径下的文件
list.files()
outdir <- getwd()
outputDir <- file.path(outdir,"monocle2")
if (!dir.exists(paste0(outputDir))){
  dir.create(paste0(outputDir))
}

# ============================================================ Load scData for TraceAnalysis ==================================================
seurat_obj <- qread("seurat_obj.qs")
ncol(seurat_obj)
Idents(seurat_obj)

clusterCornerAxes(object = seurat_obj,
                  reduction = 'umap',
                  clusterCol = "subcelltype",  #replace your explore subtype
                  noSplit = T,
                  cellLabel = T,
                  cellLabelSize = 5)
DimPlot(seurat_obj, redcution="umap", goup.by="subcelltype",pt.size=2 )
Idents(seurat_obj)<-"subcelltype"
levels(Idents(seurat_obj))

# ============================================================= Color =============================================================================
MYCOLOR<-c("#FF6700", "#9370DB", "#F8D568","#00AD43", "#89CFF0", "#BA160C",
           "#FF91AF", "#A6A6A6", "#006DB0","#C154C1", "#D99A6C", "#96C8A2",
           "#FBEC5D", "#77B5FE", "#E29CD2","#007F5C", "#ACBF60", "#7B68EE",
           "#00FFCD", "#D3AF37", "#50C878","#1974D2", "#FF0000", "#9F00FF",
           "#91A3B0", "#8B4513", "#4166F5","#C19A6B", "#6EAEA1", "#39FF14",
           "#0247FE", "#AF6E4D", "#FF66CC","#9400D3", "#A57164", "#00FFFF",
           "#DF00FF", "#FFAA1D", "#20B2AA","#CB6D51", "#5218FA", "#BC8F8F","#000000")
MYCOLOR2=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF",
            "#F39B7FFF","#8491B4FF","#91D1C2FF","#DC0000FF",
            "#d1cc78","#b2387e","#E64B354C","#4DBBD54C",
            "#00A0874C","#3C54884C","#F39B7F4C","#8491B44C",
            "#91D1C24C","#DC00004C","#fcec7c","#ff92e0",
            "#E64B357F","#4DBBD57F","#00A0877F","#3C54887F",
            "#F39B7F7F","#8491B47F","#91D1C27F","#DC00007F",
  "#eae8c1","#ef4db1","#E64B35B2","#4DBBD5B2",
  "#00A087B2","#3C5488B2","#F39B7FB2","#8491B4B2",
  "#91D1C2B2","#DC0000B2","#edd600","#c64ea4" )
#其他配色
cols <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3', 'darkmagenta', 'black', 'grey','lightpink', 'red', 'orange', 'gold', 'darkgreen')

#state 配色
colsp <-c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF','#FD7446FF','#D5E4A2FF','#197EC0FF','#F05C3BFF','#46732EFF',
          '#71D0F5FF','#370335FF','#075149FF','#C80813FF','#91331FFF','#1A9993FF','#FD8CC1FF','#FF6700','#9370DB',
          '#F8D568','#00AD43','#89CFF0','#BA160C','#FF91AF','#A6A6A6','#006DB0','#C154C1','#D99A6C','#96C8A2','#FBEC5D')


# ============================================================= Prepare InputData =============================================================================
## 提取表达数据
count_matrix <- GetAssayData(smc, assay = "RNA", layer = "counts")
count_count<- as(as.matrix(count_matrix),'spareMatrix')
## 提取表型数据
pd <- new('AnnotatedDataFrame', data = seurat_obj@meta.data)
## 提取基因信息
fData <- data.frame(gene_short_name = row.names(seurat_obj), row.names = row.names(seurat_obj))
fd <- new('AnnotatedDataFrame', data = fData)

## ------------------------------------------构建CDS文件--------------------------------------------------
### 将FPKM/TPM数据转换为UMI数据（read count）rpc_matrix <- relative2abs(cd)
cds <- newCellDataSet(count_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size()       # expressionFamily: 数据为TPM/FPKM时设置为tobit(Lower = 0.1)，数据为count时设置为negbinomial.size())
                     )

## ------------------------------------------估计size factor和离散度--------------------------------------
### size factors有助于消除细胞间mRNA捕获的差异；dispersions用于后续的差异表达分析
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# ============================================================= 轨迹构建 =============================================================================
## ------------------------------Selecting genes with high dispersion across cells---------------
## ---------挑选细胞间高度离散的基因 高度离散基因的筛选标准，可根据数据情况设置mean_expression的值------------

### Method-1(使用seurat选择的高变基因）
express_genes <- VariableFeatures(seurat_obj)
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

### Method-2(使用clusters差异表达基因）
deg.cluster <- FindAllMarkers(seurat_obj)
express_genes <- subset(deg.cluster, p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

### Method-3（使用monocle选择的高变基因）
disp_table <- dispersionTable(cds)
ordering_genes_temp <- subset(disp_table, mean_expression >= 0.1  & dispersion_empirical >= 1 * dispersion_fit) 
ordering_genes<-ordering_genes_temp$gene_id
cds <- setOrderingFilter(gbm_cds, ordering_genes)
write.table(ordering_genes, file = paste0(outputDir,"/","dispersion.ordering_genes.xls"), row.names = T, quote = F,sep='\t',col.names=T)
pdf(file = paste(outputDir,"ordering_genes.pdf",sep='/'), width = 9, height = 5)
print(plot_ordering_genes(cds))
dev.off()



## ----------------------------------------------------逆时间轴轨迹构建和在逆时间内排列细胞-----------------------------------------------
### 用DDRTree构建trajectory
tryCatch({
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
}, error = function(e) {
  print("Error in reduceDimension(). Try to apply auto_param_selection = F")
  cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree',auto_param_selection = F)
})

### order cells along the trajectory
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = 1)  # 自定义root_state

cds$barcode <- colnames(cds)
rep <- as.data.frame(pData(cds))
write.table(pData(cds), "cell_Pseudotime.txt", row.names = T, quote = F,sep=',')

### 对state_levels画图的大小进行定义
state_levels <- levels(pData(cds)$State)
if (length(state_levels) <= 2){
  widths_state = 8
  heights_state = 5
  nrows_state = 1
}else if(1 < length(state_levels) & length(state_levels) <= 8){
  widths_state = length(state_levels)*1.5
  heights_state = length(state_levels)*1.5
  nrows_state = 2
}else{
  widths_state = length(state_levels)*1.5
  heights_state = length(state_levels)*1.5
  nrows_state = 3
}


## --------------------------------------------------------------Visulazition----------------------------------------------------------------
### 对"celltype" 画图的大小进行定义!!!需要修改！！！
celltype_levels <-levels(cds$celltype)
if (length(celltype_levels) <= 2){
  widths_celltype = 8
  heights_celltype = 5
  nrows_celltype = 1
}else if(1 < length(celltype_levels) & length(celltype_levels) <= 8){
  widths_celltype = length(celltype_levels)*1.5
  heights_celltype = length(celltype_levels)*1.5
  nrows_celltype = 2
}else{
  widths = length(celltype_levels)*1.5
  heights = length(celltype_levels)*1.5
  nrows = 3
}

### 对"sample" 画图的大小进行定义
sample_levels <-levels(cds$orig.ident)
if (length(celltype_levels) <= 2){
  widths_orig = 8
  heights_orig = 5
  nrows_orig = 1
}else if(1 < length(celltype_levels) & length(celltype_levels) <= 8){
  widths_orig = length(celltype_levels)*1.5
  heights_orig = length(celltype_levels)*1.5
  nrows_orig = 2
}else{
  widths_orig = length(celltype_levels)*1.5
  heights_orig = length(celltype_levels)*1.5
  nrows_orig = 3
}

### Pseudotime_Score
p1 <-plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")
p1
ggsave(p1,file = paste0(outputDir,"/","1_Pseudotime_trajectory.pdf"), width = 8, height = 6, dpi=300, limitsize = FALSE)
p1 <-plot_cell_trajectory(cds, color_by = "Pseudotime")  + theme(plot.title = element_text(hjust = 0.5),legend.position = "top")+facet_wrap(~orig.ident)
p1
ggsave(p1,file = paste0(outputDir,"/","2_Pseudotime_trajectory_splited_by_orig.pdf"), width = widths_orig, height = heights_orig,limitsize = FALSE)

### state_Score
p2 <-plot_cell_trajectory(cds, color_by = "State") +scale_color_manual(values=colsp)+
  ggtitle("State") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p2
p2 <-plot_cell_trajectory(cds, color_by = "State", cell_size = 0.5,cell_link_size = 0.5) +
  facet_wrap(~State, nrow = nrows_state, scales = "free")+
  scale_color_manual(values=colsp)+ggtitle("State") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p2
p2 <-plot_cell_trajectory(cds, color_by = "State", cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~orig.ident, nrow = nrows_orig,scales = "free")+scale_color_manual(values=colsp)+ggtitle("State") + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p2
p2 <- ggplot(rep, aes(Pseudotime, fill = State)) +geom_density() +facet_wrap(~State,nrow =nrows_state, scales = "free_y") +
  theme_bw() +RotatedAxis() +
  theme(
    strip.text = element_blank(),
    strip.background = element_rect(color = "white", fill = "white"),
    panel.grid = element_blank()
  ) +scale_fill_manual(values = colsp)
p2

### celltype
p3 <-plot_cell_trajectory(gbm_cds, color_by = 'celltype')+scale_color_manual(values=MYCOLOR)+
  ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p3
ggsave(p3 ,file =paste0(outputDir,"/","7_Celltype_trajectory.pdf"), width = 8, height = 6,limitsize = FALSE)

# celltype splited by orig.ident
p3 <-plot_cell_trajectory(gbm_cds, color_by = 'celltype')+scale_color_manual(values=MYCOLOR)+facet_wrap(~orig.ident)+
  ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p3
ggsave(p3 ,file =paste0(outputDir,"/","8_Celltype_trajectory_splited_by_orig.pdf"), width = widths_orig, height = heights_orig,limitsize = FALSE)

p3 <-plot_cell_trajectory(gbm_cds, color_by = 'celltype', cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~celltype, nrow = nrows_celltype, scales = "free")+scale_color_manual(values=MYCOLOR)+ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p3
ggsave(p3 ,file =paste0(outputDir,"/","9_Celltype_trajectory_splited_by_celltype.pdf"), width = widths_celltype, height = heights_celltype,limitsize = FALSE)

p3 <- ggplot(rep, aes(Pseudotime, fill = celltype)) +geom_density() +facet_grid(~celltype, scales = "free_y") +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()
) +scale_fill_manual(values = MYCOLOR)
p3
ggsave(p3 ,file =paste0(outputDir,"/","10_Density_trajectory_splited_by_celltype.pdf"), height = heights_celltype,  width = widths_celltype,limitsize = FALSE)

p3 <- ggplot(rep, aes(Pseudotime, fill = celltype)) +geom_density() +facet_wrap(~orig.ident) +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()) +scale_fill_manual(values = MYCOLOR)
p3
ggsave(p3 ,file =paste0(outputDir,"/","11_Density_trajectory_splited_by_orig.pdf"), height = heights_orig,  width = widths_orig,limitsize = FALSE)

### Sample
p4 <-plot_cell_trajectory(gbm_cds, color_by = 'orig.ident') +scale_color_manual(values=MYCOLOR)+ggtitle('celltype') + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p4
ggsave(p4 ,file =paste0(outputDir,"/","12_Sample_trajectory.pdf"), width = 8, height = 6,limitsize = FALSE)

p4 <-plot_cell_trajectory(gbm_cds, color_by = 'orig.ident', cell_size = 0.5,cell_link_size = 0.5) +facet_wrap(~celltype, nrow = nrows_celltype, scales = "free")+scale_color_manual(values=MYCOLOR) + theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
p4
ggsave(p4 ,file =paste0(outputDir,"/","12_Sample_trajectory_splited_by_orig.pdf"), width = widths_celltype, height = heights_celltype,limitsize = FALSE)

p4 <- ggplot(rep, aes(Pseudotime, fill = orig.ident)) +geom_density() +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()) +scale_fill_manual(values = MYCOLOR2)
p4
ggsave(p4 ,file =paste0(outputDir,"/","13_Density_trajectory_splited_by_orig.pdf"), height = 8, width = 6,limitsize = FALSE)

p4 <- ggplot(rep, aes(Pseudotime, fill = orig.ident)) +geom_density() +facet_wrap(~celltype) +theme_bw() +RotatedAxis() +theme(
  strip.text = element_blank(),
  strip.background = element_rect(color = "white", fill = "white"),
  panel.grid = element_blank()) +scale_fill_manual(values = MYCOLOR2)
p4
ggsave(p4 ,file =paste0(outputDir,"/","14_Density_trajectory_splited_by_orig.pdf"), height = heights_celltype, width = widths_celltype,limitsize = FALSE)

### 树形图
p2 <- plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                   color_by = "celltype") +
  scale_color_manual(values = color) +
theme(legend.title = element_blank())
p2

### 沿时间轴的细胞密度图
df <- pData(cds) # 取出cds对象中cds@phenoData@data的内容
View(df)
ggplot(df, aes(Pseudotime, colour = celltype, fill = celltype))+
geom_density(bw=0.5, size=1,alpha = 0.5) +theme_classic2()

###！！！提取感兴趣的细胞！！！
pdata <- Biobase::pData(cds)
s.cells <- subset(pdata,State == "7") %>% rownames()
save(s.cells,file="Celltype_state7.rad")
write.csv(pData(cds),"pseudotime.csv")
save(cds,file="cds.rda")


# ============================================================= 指定基因的可视化 =============================================================================
### 选择前4个top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
## 可视化
p1 <- plot_genes_in_pseudotime(cds_subset,color_by="State")
p2 <- plot_genes_in_pseudotime(cds_subset,color_by="celltype")
p3 <- plot_genes_in_pseudotime(cds_subset,color_by="Pseudotime")

## 指定基因
s.genes <- c("CD4","CD8",...)
p1 <- plot_genes_jitter(cds[s.genes,],grouping = "State", color_by ="State")
p2 <- plot_genes_violin(cds[s.genes,],grouping = "State", color_by ="State")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by ="State")

## 逆时序展示单个基因表达量
colnames(pData(cds))
pData(cds)$CCL5 = log2(exprs(cds)['CCL5',]+1)
p1 <- plot_cell_trajectory(cds, color_by = "CCL5") = scale_color_gsea()
p1


# ============================================================= 逆时差异基因 =============================================================================
## disp_table <- dispersionTable(cds)
## ordering_genes <- subset(disp_table,mean_expression >= 0.1 &dispersion_empirical >= 1 * dispersion_fit)
cds_expressed_genes <- rownames(fData(cds))
diff_test_res <- differentialGeneTest(cds[cds_expressed_genes,],cores = 1，
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")

sig_gene_names <- diff_test_res %>% 
                  dplyr::filter(qval < 0.01) %>% 
                  arrange(qval) %>% row.names()  
diff_test_res1 <- diff_test_res %>% rownames_to_column("gene")

write.table(diff_test_res1, file = paste0(outputDir,'/','diff_test_ordering_genes.xls'), row.names = F,col.names = T,quote = F,sep="\t")
#write.table(ordering_genes,file=paste0(outputDir,'/',"ordering_genes.csv"),row.names = F,col.names = T,quote = F,sep="\t")
write.table(sig_gene_names, file = paste0(outputDir,'/','ordering_genes_sig_gene_names.all.xls'), row.names = F,col.names = F,quote = F)

##-----------------------------------------筛选top50 展示热图------------------------------------
tryCatch({
  p <-plot_pseudotime_heatmap(cds[sig_gene_names[1:50],],
                              num_clusters = 3,
                              cores = 4,
                              show_rownames = T,
                              return_heatmap= T)
  ggsave(p,file=paste0(outputDir,'/',"pseudotime_dependent_gene_heatmap_top50.pdf"),width = 9, height = 5,limitsize = FALSE)
  # 拆分拟时序分析轴上每个cluster的基因（拆分top50个基因）
  clusters <- cutree(p$tree_row, k = 3)
  clustering <- as.data.frame(clusters) %>% rownames_to_column(var = "gene")
  colnames(clustering) <- c("gene","cluster")
  write.table(clustering,file=paste0(outputDir,'/',"pseudotime_dependent_gene_heatmap_top50_gene_by_cluster.xls"),sep="\t",quote = F,col.names = T,row.names = F)
}, error = function(e) {
  print("Fewer than 50 significant genes")
})

## ----------------------------------------手动提取基因单独分析-------------------------------------
p$tree_row
clusters <- cutree(p$tree_row, k = 2)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.csv(clustering, "Time_clustering_all.csv",row.names=F)

### 逆时差异基因热图绘制
Time_genes <- top_n(diff_test_res, n =100, desc(qval) )%>% 
                    pull(gene_short_name) %>% 
                    as.character()
p = plot_pseudotime_heatmap(cds[Time_genes,], 
                            num_clusters = 4, 
                            show_rownames = T,
                            return_heatmap =T)

### 显著差异基因按热图结果排序
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- diff_test_res[hp.genes, c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"Time_diff_sig.csv",row.names = F)

### 手动选择基因绘制热图
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_short_name %in% c("A","B",...)))
diff_test_res <- differentialGeneTest(cds[cds_expressed_genes,],cores = 1，
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval <0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                              num_clusters = 3,
                              cores = 4,
                              show_rownames = T,
                              return_heatmap= T)

# ============================================================= 差异表达基因做GO分析 =============================================================================
# TYPE的可选类型：SYMBOL, ENSEMBL, ENTREZID
# OrgDb的可选类型：'org.Dr.eg.db','org.Mm.eg.db','org.Hs.eg.db'
library(org.Hs.eg.db) # 确保加载了org.Hs.eg.db包
library(enrichplot)   # 确保加载了enrichplot包，用于可视化
library(clusterProfiler)

dir.create('GO')
# 初始化或清空错误信息文件
#这行代码的作用是打开（或创建）一个名为 "error_messages.txt" 的文件，以追加模式进行写入。
#这意味着任何写入这个文件的数据都会被添加到文件的末尾，而不会覆盖原有内容。
#如果 "error_messages.txt" 文件不存在，这行代码将会创建这个文件。

error_file <- file("GO/error_messages.txt", "a")
writeLines(c("Error log for GO enrichment analysis:", 
             "------------------------------------------------",
             'starting time at:', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             "------------------------------------------------"), 
           error_file, useBytes = FALSE)

for (i in unique(clustering$cluster)) {
  # 尝试转换基因ID
  gene <- tryCatch({
    bitr(geneID = clustering$gene[clustering$cluster == i],
         fromType = 'SYMBOL',
         toType = 'ENTREZID',
         OrgDb = 'org.Hs.eg.db')}, 
    error = function(e) {
      message_text <- paste("Error in gene ID conversion for :", i, ": ", e$message)
      writeLines(message_text, error_file, useBytes = FALSE)
      return(NULL)})
  
  # 如果基因ID转换失败，跳过当前循环的剩余部分
  if (is.null(gene)) {
    next  }
  
  # 进行GO富集分析
  ego <- tryCatch({
    enrichGO(gene = gene$ENTREZID,
             OrgDb = 'org.Hs.eg.db',
             ont = "BP",
             pvalueCutoff = 0.05,
             qvalueCutoff = 1,
             readable = TRUE)
  }, error = function(e) {
    message_text <- paste("Error in GO enrichment analysis for :", i, ": ", e$message)
    writeLines(message_text, error_file, useBytes = FALSE)
    return(NULL)})
  
  # 检查是否富集出任何通路
  # 如果没有富集出任何通路，跳过画图代码
  
  if (is.null(ego) || nrow(summary(ego)) == 0) {
    message_text <- paste("No enriched GO terms found for :", i)
    writeLines(message_text, error_file, useBytes = FALSE)
    next}
  
  # 保存富集分析结果
  file_ego <- paste0("GO/",'Cluster_',i, "_Upregulated_GO_terms.csv", sep = "")
  write.csv(summary(ego), file_ego)
  
  # 绘制图形并设置标题
  file_name <- paste0("GO/Barplot_GO_terms_enriched_in_Cluster_", i, ".pdf")
  p <- barplot(ego, title = paste0("GO terms enriched in ", i))
  ggsave(file_name, plot = p, width = 8, height = 6)
  
  # 由于dotplot函数的使用可能有误，这里假设你使用的是 enrichplot 包的 dotplot 函数
  # 如果不是，请根据实际使用的包和函数进行调整
  file_name <- paste0("GO/Dotplot_GO_terms_enriched_in_Cluster_", i, ".pdf")
  p<- dotplot(ego, x = 0, title = paste0("GO terms enriched in cluster ", i, sep = ""))
  ggsave(file_name, plot = p,width = 3.83, height = 5.00)
  }

close(error_file)



# ============================================================= 单细胞轨迹的分支分析 =============================================================================
## ------------------------------------------BEAM进行统计分析----------------------------------------------------
### dpFeature
expressed_genes <- row.names(fData(cds))
diff <- differentialGeneTest(cds[expressed_genes,],
                             fullModelFormulaStr = "~celltype",  # 理论上可以为p_data的任意列名
                             cores=1)
head(diff)

#### 差异表达基因作为轨迹构建的基因，差异基因的选择标准是qval<0.01, decreasing = F表示按数值增加排序
deg <- subset(diff, qval< 0.01)
deg <- deg[order(deg$qval, decreasing = F),]
head(deg)
write.table(deg, file="train.monocle.DEG.xls",col.names=T, row.names=F,sep="\t",quote=F)

ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)   # 将基因列表嵌入cds对象
## 以上操作储存在cds@featureData@data[["use_for_ordering"]]
## 通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)
# ordergene <- row.names(deg)[order(deg$qval)][1:400] #选择用于排序的基因数目一般在2000 
BEAM_res <- BEAM(cds[ordergene,], branch_point = 1, cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name","pval","qval")]
head(BEAM_res)
write.csv(BEAM_res,"BEAM_res.csv", row.names = F)
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
### 该热图显示的是同一时间点两个谱系的变化，热图的列是伪时间的点，行是基因。这张图最上面的条条日澳，灰色的代表分叉前，红色代表左边cellfate
### 热图基因按照等级聚类

### 选前100个基因可视化
BBEAM_genes <- top_n(BEAM_res, n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,], 
                                 branch_point = 1,
                                 num_clusters = 3,
                                 show_rownames = T,
                                 return_heatmap = T)
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes, c("gene_short_name","pval","qval")]
write.csv(BEAM_sig, "BEAM_sig.csv", row.names = F)
head(BEAM_sig)
genes <- row.names(susbet(fData(cds),
                          gene_short_name %in% c("A","B",...)))
plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)




















