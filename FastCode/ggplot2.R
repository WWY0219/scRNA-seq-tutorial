geom_boxplot(
  aes(x = group, y = value, fill = group),
  # 离群点设置
  outlier.color = "red",       # 离群点边框颜色
  outlier.fill = "yellow",     # 离群点填充颜色
  outlier.shape = 21,          # 形状（21=带填充的圆形）
  outlier.size = 2,           # 离群点大小
  outlier.stroke = 1,         # 离群点边框粗细
  outlier.alpha = 0.7,        # 离群点透明度
  
  # 箱体设置
  color = "black",            # 箱线图边框颜色
  size = 0.5,                 # 箱线图边框粗细
  alpha = 0.8,                # 箱线图透明度
  linetype = "solid",         # 线条类型
  
  # 异常值显示控制
  show.outliers = TRUE       # 是否显示离群点
)
