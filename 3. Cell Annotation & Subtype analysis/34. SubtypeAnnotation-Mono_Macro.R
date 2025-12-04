# =========================== Prepare Environment ==============================
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("worksapce")  #replace your workspace
getwd()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)
library(SCP)
set.seed(1234)
# 查看工作路径下的文件
list.files()

# =================================== Subcelltype Annotation =====================================
## Loading Major-subtype and Subset Monocyte
seurat_obj <- qread("seurat_obj_annotation.qs")
DimPlot(seurat_obj,reduction = "umap",group.by = "seurat_clusters",label = T,pt.size = 0.25)+NoLegend()
seurat_obj <- subset(seurat_obj, subset=celltype_major=="Mono/Macro")
print(seurat_obj)


## 用细胞总 UMI 计数的中位数作为缩放因子消除细胞间测序差异
seurat_obj <- NormalizeData(seurat_obj, normalization.method ="LogNormalize", 
                            scale.factor = median(seurat_obj@meta.data$nCount_RNA))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000) 

## 计算细胞周期评分
cc.genes.updated.2019 <- cc.genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes                        
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes) 
## 回归掉不感兴趣的变量
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo1",
                                                        "percent_ribo2","percent_mt","percent_RBC"))
## 使用HVG去跑PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
ElbowPlot(seurat_obj,ndims = 50) 

## Perform Harmony batch correction
seurat_obj <- seurat_obj %>% RunHarmony(
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony",    
  plot_convergence = TRUE,
  max.iter = 30,
  verbose = FALSE
)
## Save Harmony convergence plot
ElbowPlot(seurat_obj, reduction = "harmony", ndims = 50) + 
  theme_minimal() + 
  ggtitle("Elbow Plot for Harmony-Corrected Dimensions")
ggsave(file.path(out_dir, paste0(obj_name, "_Harmony_convergence.pdf")), width = 8, height = 6, dpi = 300)

seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30, 
                                reduction.name = "umap", verbose = T)
seurat_obj <- RunTSNE(seurat_obj, reduction = "harmony", dims = 1:30, 
                                reduction.name = "tsne", verbose = T)



seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "harmony")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2)
table(seurat_obj@meta.data$seurat_clusters)
table(seurat_obj@meta.data$orig.ident)

plot1 <- DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T,pt.size = 0.25)
plot2 <- DimPlot(seurat_obj,reduction = "umap",group.by = "seurat_clusters",label = T,pt.size = 0.25)+NoLegend()
plot1|plot2

# ==== 加载我们要用到的marker ====
known_markers = list(
  "cMo"   = c("CD14","CCR2"),
  "ncMo"  = c("NR4A1","CX3CR1"),
  "intMo" = c("HLA-DR","CD86"),
  "MDSC"  = c("S100A8","S100A9","ARG1","STAT3"),
  "M1 Macro" =
  "M2 Macro" =
  "TAM"      =
  "cDC1"     =
  "cDC2"     =
  "pDC"      =
  "migDC"    =
  )

tcell_markers <- list(
  # 基本谱系
  "CD4 lineage"   = c("CD4"),
  "CD8 lineage"   = c("CD8A","CD8B"),
  # 初始/naive & 中央记忆
  "Naive/Tcm"     = c("CCR7","SELL","TCF7","LEF1","IL7R","LTB"),
  # 效应型 (Teff) / 效应记忆 (Tem)
  "Effector/Tem"  = c("GZMA","GZMB","GZMK","PRF1","GNLY","NKG7","IFNG","CX3CR1"),
  # 组织驻留 (Trm)
  "Trm"           = c("ITGAE","ZNF683","HOPX","CD69"),
  # 耗竭 (Tex)
  "Exhaustion"    = c("PDCD1","HAVCR2","LAG3","TOX"),
  # 调节性T细胞 (Treg)
  "Treg"          = c("FOXP3","IL2RA","CTLA4","TNFRSF4","IKZF2"),
  # 辅助 Tfh
  "Tfh"           = c("CXCR5","BCL6","ICOS","IL21"),
  # 增殖/活化
  "Proliferation" = c("MKI67","TOP2A","STMN1","HMGB2"),
  "Activation"    = c("TIGIT","TNFRSF9")
)

# ==== 经典气泡图，看不同cluster的marker gene表达情况 ====
DotPlot(object = seurat_obj,
        features = known_markers,
        scale=T,
        group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()

DotPlot(object = seurat_obj,
        features = tcell_markers,
        scale=T,
        group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()

# ==== 细胞注释 ====
meta_supp = data.frame(seurat_cluster = 0:(length(unique(seurat_obj$seurat_clusters)) - 1), Celltype = NA)
meta_supp[meta_supp$seurat_cluster %in% c(0), 'Celltype'] = 'Naive/Central Memory T'
meta_supp[meta_supp$seurat_cluster %in% c(1), 'Celltype'] = 'CD8+ Effector T'
meta_supp[meta_supp$seurat_cluster %in% c(2), 'Celltype'] = 'Cytotoxic NK'
meta_supp[meta_supp$seurat_cluster %in% c(3), 'Celltype'] = 'CD8+ Effector-Memory T'
meta_supp[meta_supp$seurat_cluster %in% c(4), 'Celltype'] = 'CD4+ Regulatory T'
meta_supp[meta_supp$seurat_cluster %in% c(5), 'Celltype'] = 'Proliferative T'
meta_supp[meta_supp$seurat_cluster %in% c(6), 'Celltype'] = 'CD8+ Residue-Memory T'

for (i in 1:nrow(meta_supp)) {
  seurat_obj@meta.data[which(seurat_obj$seurat_clusters == meta_supp$seurat_cluster[i]), 'subtype'] = meta_supp$Celltype[i]
}
Idents(seurat_obj) <- 'subtype'

# 看看注释情况
CellDimPlot(
  seurat_obj,
  group.by = "celltype_major",
  theme_use = "theme_blank",
  xlab = "UMAP1",
  ylab = "UMAP2",
  label = TRUE,           
  label_insitu = TRUE,
  show_stat = F,      
  legend.position = "none" 
)
DotPlot(object = seurat_obj,
              features = known_markers,
              group.by = "Celltype_major")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
DotPlot(object = seurat_obj,
        features = tcell_markers,
        scale=T,
        group.by = "Celltype_major")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
table(seurat_obj@meta.data$Celltype_major)

# ==== 不要忘记保存注释后的对象 ====
qsave(seurat_obj, file = 'T_NK_post_annotation.qs')

# ==== 画个饼图，看看不同的患者中不同细胞的占比 ====
pie_data=data.frame(
  row.names = rownames(seurat_obj@meta.data),
  sample=seurat_obj@meta.data$orig.ident,
  Celltype=seurat_obj@meta.data$celltype_major
)
# 统计每个患者中不同细胞类型的数量
plot_data <- pie_data %>%
  group_by(sample, Celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count),
         label = ifelse(freq > 0.05, paste0(round(freq * 100, 1), "%"), ""))  # 只标 >5%

# 画饼图 + 数字
ggplot(plot_data, aes(x = "", y = freq, fill = Celltype)) +
  geom_col(width = 1, color = "white") +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), size = 4) +  # 堆叠居中
  coord_polar(theta = "y") +
  facet_wrap(~ sample) +
  theme_void() +
  theme(
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "right"
  )
