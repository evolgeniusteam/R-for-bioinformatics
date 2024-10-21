library(ggplot2)
library(dplyr)
library(scales)

# 创建示例数据框
data <- data.frame(
  category = c("A", "B", "C", "D"),
  count = c(10, 15, 30, 45)
)

# 计算每个类别的百分比和标签位置
data <- data %>%
  mutate(
    percentage = count / sum(count),                             # 计算每个类别所占百分比
    label = percent(percentage),                             # 转换为百分比格式，用于显示标签
    label_position = sum(count) - (cumsum(count) - count /2)
  )

# 创建 Donut Plot
ggplot(data, aes(x = 2, y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +       # 使用 stat = "identity"，设置边框颜色为白色
  coord_polar(theta = "y") +                                      # 转换为极坐标
  xlim(0.5, 2.5) +                                                # 控制内外圈大小，形成甜甜圈
  geom_text(aes( y = label_position, label = label),       # 使用 label_position 确定标签位置
            color = "black", size = 5) +                          # 设置文本颜色和字体大小
  theme_void() +                                                  # 移除背景和网格线
  theme(legend.position = "right") +                              # 调整图例位置
  ggtitle("Donut Plot with Correct Percentage Labels")            # 设置图表标题
