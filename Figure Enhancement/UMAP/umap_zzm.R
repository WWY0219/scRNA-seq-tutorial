# ================================================ 张泽民云雾状UMAP图 ====================================================

## ===============================================Single-UMAP============================================================
# 加载必需包
library(ggplot2)
library(ggh4x)

## ----------------------UMAP-Data Subset--------------------------
Idents(seurat_obj) <- "celltype_major"
umap = seurat_obj@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = seurat_obj@meta.data$celltype) # 注释后的label信息 ，改为cell_type
head(umap)
min_umap1 <- min(umap$umap_1)
min_umap2 <- min(umap$umap_2)
x_end <- min_umap1 + 3
y_end <- min_umap2 + 3

## --------------------- Color Vector ------------------------------
cell_colors <- c("#919ac2","#ffac98","#70a4c8","#a5a9af","#63917d","#dbd1b4","#6e729a","#9ba4bd","#c5ae5f","#b9b8d6")
cell_major_colors <- c("Tumor"= "#919ac2",
                       "Fibroblast" = "#ffac98",
                       "Monocyte" = "#70a4c8",
                       "T/NK"="#a5a9af",
                       "SMC"="#cea5c7",
                       "Endothelial" = "#AECDE1",
                       "Mastcell" = "#ece399",
                       "Epithelial" = "#de8b36",
                       "VSMC" = "#6e729a",
                      "Proliferating_cell"="#b9b8d6")

## --------------------- Draw Picture ------------------------------
p <- ggplot(umap, aes(x = umap_1, y = umap_2, color = cell_type)) +  
  # 1. 绘制细胞点
  geom_point(size = 0.6, shape = 16, stroke = 0) +  
  # 2. 颜色映射
  scale_color_manual(values = cell_major_colors) +
  # 3. 主题美化（合并theme，避免重复）
  theme(
    panel.grid.major = element_blank(),    # 主网格线
    panel.grid.minor = element_blank(),    # 次网格线
    panel.border = element_blank(),        # 边框
    axis.title = element_blank(),          # 轴标题
    axis.text = element_blank(),           # 轴文本
    axis.ticks = element_blank(),          # 轴刻度
    panel.background = element_rect(fill = 'white'),  # 面板背景
    plot.background = element_rect(fill = "white"),   # 绘图区背景
    legend.title = element_blank(),        # 图例标题
    legend.key = element_rect(fill = 'white'), # 图例背景
    legend.text = element_text(size = 20), # 图例文本大小
    legend.key.size = unit(1, 'cm')        # 图例符号大小
  ) +  
  # 4. 调整图例中点的大小
  guides(color = guide_legend(override.aes = list(size = 5)),
        title.theme = element_blank()) + 
  # 5. 绘制X轴箭头（用annotate替代geom_segment，消除长度警告）
  annotate(
    "segment",
    x = min_umap1, y = min_umap2,          # 起点
    xend = x_end, yend = min_umap2,        # 终点
    colour = "black", size = 1,
    arrow = arrow(length = unit(0.3, "cm"))
  ) + 
  # 6. 绘制Y轴箭头（用annotate替代geom_segment）
  annotate(
    "segment",
    x = min_umap1, y = min_umap2,          # 起点
    xend = min_umap1, yend = y_end,        # 终点
    colour = "black", size = 1,
    arrow = arrow(length = unit(0.3, "cm"))
  ) +
  # 7. 标注X轴文字
  annotate(
    "text", 
    x = min_umap1 + 1.5, y = min_umap2 - 1, 
    label = "UMAP_1", color = "black", size = 3, fontface = "bold"
  ) + 
  # 8. 标注Y轴文字
  annotate(
    "text", 
    x = min_umap1 - 1, y = min_umap2 + 1.5, 
    label = "UMAP_2", color = "black", size = 3, fontface = "bold", angle = 90
  )

# 打印图形
print(p)
ggsave("sft_umap.pdf",plot=p,width=8,height=6,dpi=300)

## ===============================================分面绘制UMAP（按group分面，保留cell_type着色）============================================================

p <- ggplot(umap, aes(x = umap_1, y = umap_2, color = cell_type)) +  
  geom_point(size = 0.05, shape = 16, stroke = 0) +  
  scale_color_manual(values = cell_major_colors) +
  facet_wrap(~group, ncol = 3) +                               # 按分组拆分子图
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"),
    legend.title = element_blank(),
    legend.key = element_rect(fill = 'white'),
    legend.text = element_text(size = 12),                     
    legend.key.size = unit(0.8, 'cm'),
    strip.background = element_rect(fill = "lightgray"),        # 分面标题背景
    strip.text = element_text(size = 14)                        # 分面标题文字
  ) +  
  guides(color = guide_legend(override.aes = list(size = 5))) + 
  # 若保留箭头，需将annotate改为geom_segment（annotate在分面中会只绘制一次）
  geom_segment(
    aes(x = min_umap1, y = min_umap2, xend = x_end, yend = min_umap2),
    colour = "black", size = 1,
    arrow = arrow(length = unit(0.3, "cm")),
    inherit.aes = FALSE                                        # 不继承color映射
  ) + 
  geom_segment(
    aes(x = min_umap1, y = min_umap2, xend = min_umap1, yend = y_end),
    colour = "black", size = 1,
    arrow = arrow(length = unit(0.3, "cm")),
    inherit.aes = FALSE
  ) +
  annotate(
    "text", 
    x = min_umap1 + 1.5, y = min_umap2 - 1, 
    label = "UMAP_1", color = "black", size = 3, fontface = "bold"
  ) + 
  annotate(
    "text", 
    x = min_umap1 - 1, y = min_umap2 + 1.5, 
    label = "UMAP_2", color = "black", size = 3, fontface = "bold", angle = 90
  )

print(p)
ggsave("t_umap_facet_by_group.pdf", plot = p, width = 20, height = 8, dpi = 300)
