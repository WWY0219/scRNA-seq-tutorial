# =====================================================Theme=======================================================
## ----------------------------ggplot2::theme----------------------------------
theme(
  # 基本几何元素
line = element_line()       # 所有线条
rect = element_rect()       # 所有矩形区域
text = element_text()       # 所有文本
title = element_text()      # 所有标题
point = element_point()     # 所有点
polygon = element_polygon() # 所有多边形

# 布局相关
aspect.ratio = NULL        # 图形宽高比
spacing = unit()           # 间距单位
margins = margin()         # 边距
theme(
  # 全局设置
  text = element_text(family = "Arial", size = 12),
  line = element_line(linewidth = 0.5),
  
  # 坐标轴
  axis.title = element_text(face = "bold"),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.line = element_line(color = "black"),
  axis.ticks = element_line(color = "black"),
  
  # 面板
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_line(color = "grey80"),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "black", fill = NA),
  
  # 图例
  legend.position = "top",
  legend.background = element_rect(fill = "grey95"),
  legend.key = element_rect(fill = "white"),
  
  # 绘图区域
  plot.title = element_text(hjust = 0.5, face = "bold"),
  plot.background = element_rect(fill = "white"),
  plot.margin = margin(20, 20, 20, 20)
)
