# ==================================================Cell Proportion by Group=====================================================
## 统计每组celltype 的细胞数量
celltype_counts <- table(seurat_obj$group, seurat_obj$celltype) %>% 
  as.data.frame() %>% 
  rename(group = Var1, celltype = Var2, count = Freq)

# 计算每个 orig.ident 中各细胞类型的占比（百分比）
celltype_proportions <- celltype_counts %>%
  group_by(group) %>%
  mutate(percentage = count / sum(count) * 100)

group_prop <- ggplot(celltype_proportions, aes(x = group, y = percentage, fill = celltype)) +
  geom_col(position = "fill") +  # 百分比堆积
  scale_y_continuous(labels = scales::percent_format()) +  # y轴显示百分比
  labs(x = "Samples", y = "celltype percentage", fill = "Celltype") +
  theme(
      panel.background = element_blank(),
    # 2. 移除绘图区边框
    panel.border = element_blank(),
    # 3. 移除网格线（x/y轴都移除）
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),  # X轴基因名旋转45度
    axis.text.y = element_text(size = 15),                         # Y轴字体大小
    axis.title = element_text(size = 18, face = "bold"),            # 轴
    legend.title = element_text(size = 18, face = "bold"),         # 图例标题加粗
    legend.text = element_text(size = 15),                         # 图例文字大小
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5))+
scale_fill_manual(values = cell_colors)
print(group_prop)
ggsave("group_prop.pdf",plot=group_prop,width = 6,height=8,dpi=300)
