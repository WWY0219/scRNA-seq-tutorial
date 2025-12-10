# ==================================================Cell Proportion by Group=====================================================
## 统计每组celltype 的细胞数量
celltype_counts <- table(seurat_obj$group, seurat_obj$celltype) %>% 
  as.data.frame() %>% 
  rename(group = Var1, celltype = Var2, count = Freq)

# 计算每个 group 中各细胞类型的占比（百分比）
celltype_proportions <- celltype_counts %>%
  group_by(group) %>%
  mutate(percentage = count / sum(count) * 100)

group_prop <- ggplot(celltype_proportions, aes(x = group, y = percentage, fill = celltype)) +
  geom_col(position = "fill") +                            # 百分比堆积
  scale_y_continuous(labels = scales::percent_format()) +  # y轴显示百分比
  labs(x = "Samples", y = "celltype percentage", fill = "Celltype") +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 15),       
    axis.text.y  = element_text(size = 15),                                 
    axis.title   = element_text(size = 18, face = "bold"),           
    legend.title = element_text(size = 18, face = "bold"),         
    legend.text  = element_text(size = 15),                        
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5))+
scale_fill_manual(values = cell_colors)
print(group_prop)
ggsave("group_prop.pdf",plot=group_prop,width = 6,height=8,dpi=300)

# ===================================================Cell Counts by Group=======================================================
celltype_counts <- as.data.frame(table(seurat_obj$group, seurat_obj$celltype)) %>%
  rename(
    orig.ident = Var1,    # 样本名
    celltype = Var2,      # 细胞类型
    count = Freq          # 绝对数量
  )

group_count <- ggplot(celltype_counts, aes(x = orig.ident, y = count, fill = celltype)) +
geom_col(position = "stack") +  
  labs(
    x = "Samples", 
    y = "Number of Cells", 
    fill = "Celltype"
  ) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),  
    axis.text.y = element_text(size = 15),                         
    axis.title = element_text(size = 18, face = "bold"),            
    legend.title = element_text(size = 18, face = "bold"),         
    legend.text = element_text(size = 15),                         
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )+
  scale_fill_manual(values = cell_colors) +
  scale_y_continuous(labels = scales::comma_format())
print(group_count)
ggsave("group_count.pdf",plot=group_count,width = 6,height=8,dpi=300)
