# ======================================================geom_boxplot==========================================
geom_boxplot(
  # 几何外观参数
  aes(x = group, y = value, fill = group),
  ## 离群点设置
  outlier.color = "red",       # 离群点边框颜色
  outlier.fill = "yellow",     # 离群点填充颜色
  outlier.shape = 21,          # 形状（21=带填充的圆形）
  outlier.size = 2,            # 离群点大小
  outlier.stroke = 1,          # 离群点边框粗细
  outlier.alpha = 0.7,         # 离群点透明度
  ## 箱体设置
  color = "black",             # 箱线图边框颜色
  size = 0.5,                  # 箱线图边框粗细
  alpha = 0.8,                 # 箱线图透明度
  linetype = "solid",          # 线条类型
  ## 异常值显示控制
  show.outliers = TRUE         # 是否显示离群点

  #箱线图形态参数
  ## 凹口箱线图（用于比较中位数）
  notch = TRUE,              # 是否显示凹口
  notchwidth = 0.5,          # 凹口宽度（0-1）
  ## 宽度设置
  width = 0.7,              # 箱线图宽度（0-1）
  varwidth = FALSE,         # 是否根据样本量调整宽度
  ## 位置调整
  position = position_dodge(width = 0.8)  # 分组并排

  #统计与位置参数
  stat = "boxplot",          # 使用的统计变换
  position = "dodge2",       # 位置调整方式，可选：
  # - "identity": 不调整
  # - "dodge": 分组并排
  # - "dodge2": 改进的分组并排
  # - "jitter": 添加抖动
  # - position_dodge(): 自定义并排
  # - position_dodge2(): 改进的自定义并排
  orientation = "x",         # 方向："x"或"y"
  na.rm = FALSE              # 是否移除NA值

  show.legend = c(fill = TRUE, color = FALSE),  # 控制哪些美学显示图例
  inherit.aes = TRUE         # 是否继承父图层的aes()
)

