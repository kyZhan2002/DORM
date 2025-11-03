library(ggplot2)
library(reshape2)

df_mix0 <- data.frame(
  Index = seq(0.05,0.5,0.05),
  DORM_condA = c(72.95,72.99,73.50, 75.11, 77.89, 76.55, 77.00, 76.92, 76.82, 77.79),
  #BestSingle = c(16.66, 26.98, 37.12, 47.70, 55.64, 58.32, 59.80, 61.25, 62.53, 63.47),
  SimpleAve = c(79.34, 80.18, 80.66, 80.62, 81.38, 81.28, 81.51, 82.02, 82.03, 81.58),
  # RA = c(19.78, 10.60578116, 15.10092016, 19.89342186, 24.52878898),
  Maximin = c(111.22, 109.71, 108.12, 105.76, 101.95, 101.77,97.03, 96.76, 94.93, 93.37),
  RhoAve = c(74.20, 75.57, 76.21, 78.02, 78.58, 79.43, 79.05, 82.45, 82.35, 84.82)
)

df_mix0[-1] = df_mix0[-1]/100

# 将数据框转换为长格式
mix0_long <- melt(df_mix0, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")
p_mix0 = NULL
p_mix0 <- ggplot(mix0_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), linewidth = 0.5) + 
  geom_point(size = 1) +
  ylim(0.6,1.2) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM_condA" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("DORM_condA" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13))


p_mix0
ggsave("CondA_plot.png", plot = p_mix0, dpi = 300, width = 8, height = 6, path = here('simu', 'simu_10242025'))


