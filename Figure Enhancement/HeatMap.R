# ========================ComplexHeatMap============================
markers<-FindAllMarkers(seurat_obj, only.pos=TRUE,
                                min.pct=0.25,
                                logfc.threshold=0.25)
#get top 5 genes
top_5 <-markers%>%
  dplyr::group_by(cluster)%>%
  dplyr::top_n(n=5,wt=avg_log2FC)
































































# ======================geom_sc_heatmap.R==================================
source("geom_sc_heatmap.R")
#读取 Seurat object
pbmc=readRDS("pbmc.rds")
markers=FindAllMarkers(pbmc)
strip.col=c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148","#B09C85")

## ==========普通热图==============
ggscplot(object = pbmc,features = markers$gene,
         mapping = aes(y = gene_name,x = seurat_clusters,
                       fill = mean_exp)) +
    geom_sc_heatmap() +
    scale_fill_gradient(low = "grey90",high = "red") +
    theme_bw()


ggscplot(object = pbmc,features = markers$gene,
         featuresAnno = markers$cluster,
         mapping = aes(x = gene_name,y = seurat_clusters,
                       fill = mean_exp,exp = mean_exp,
                       celltype = featureAnno)) +
    geom_sc_heatmap(branch.label.hjust=1,branch.label.rot=-90) +
    scale_fill_gradient(low = "grey90",high = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,hjust = 1),
          axis.text = element_text(colour = "black"),
          plot.margin = margin(t = 0.3,unit = "npc"))
ggscplot(object = pbmc,features = markers$gene,
         featuresAnno = markers$cluster,
         mapping = aes(y = gene_name,x = seurat_clusters,
                       fill = mean_exp)) +
    geom_sc_heatmap() +
    facet_new(facet_row = "featureAnno",scales = "free_y",
              strip.col = strip.col,
              switch = "y")
ggscplot(object = pbmc,features = markers$gene,
         featuresAnno = markers$cluster,
         mapping = aes(x = gene_name,y = seurat_clusters,
                       fill = mean_exp)) +
    geom_sc_heatmap() +
    facet_new(facet_col = "featureAnno",scales = "free_x",
              strip.col = strip.col,
              x.angle = 90,x.label.hjust = 1)

ggscplot(object = pbmc,features = markers$gene,
         featuresAnno = markers$cluster,
         mapping = aes(x = gene_name,y = seurat_clusters,
                       fill = mean_exp,exp = mean_exp,
                       celltype = featureAnno)) +
    geom_sc_heatmap(branch.label.hjust=1,branch.label.rot=-90,add.tree = T,add.tree.y = T,tree.y.side = "right",new.ylabel.x = 0) +
    scale_fill_gradient(low = "grey90",high = "red") +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90,hjust = 1),
          axis.text = element_text(colour = "black"),
          plot.margin = margin(t = 0.3,r = 0.2,l = 0.1,unit = "npc"),
          legend.position = c(1.15,0.5)) +
    ylab("")

ggscplot(object = pbmc,features = markers$gene,
         featuresAnno = markers$cluster,
         mapping = aes(y = gene_name,x = seurat_clusters,
                       fill = mean_exp,exp = mean_exp)) +
    geom_sc_heatmap(add.tree = T,add.tree.y = T,tree.y.side = "right",new.ylabel.x = 0) +
    facet_new(facet_row = "featureAnno",scales = "free_y",
              strip.col = circlize::rand_color(9),
              switch = "y",
              r = 0.3,no.ylabel = T,
              legend.position = c(1.4,0.5)) +
    theme(strip.switch.pad.grid = unit(0.2,"npc"))








