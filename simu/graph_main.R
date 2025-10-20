
# 导入必要的库
library(ggplot2)
library(reshape2)

# 创建数据框,from 24-1-4 AARESULT manually selected
data <- data.frame(
  Index = seq(0.05,0.5,0.05),
  REMIX = c(21.20,27.20,33.35, 37.89, 41.23,44.89, 47.56, 49.02, 50.01, 50.86),
  REMIX_true = c(20.30, 26.70, 32.92, 37.70, 41.43, 45.06, 47.31,48.93, 50.60, 51.67),
  REMIX_P0.05 = c(21.56, 27.93, 34.09, 38.82, 43.19, 46.06, 48.55, 50.41,51.52, 53.08),
  REMIX_P0.1 = c(21.56,30.16,35.59,41.21, 44.62, 48.24, 50.21, 51.41,52.77, 53.64),
  # BestSingle = c(16.66, 26.98, 37.12, 47.70, 55.64, 58.32, 59.80, 61.25, 62.53, 63.47),
  SimpleAve = c(40.02, 42.08, 44.02, 46.10, 48.35, 50.07, 52.13, 54.34, 56.51, 58.29),
  # RA = c(19.78, 10.60578116, 15.10092016, 19.89342186, 24.52878898),
  Maximin = c(54.72, 54.77, 54.91, 54.93, 55.04, 55.07, 55.22, 55.37, 55.50, 55.37),
  RhoAve = c(19.92, 26.70, 33.16, 40.52, 47.42, 53.95, 61.51, 68.77, 75.30, 82.92)
)

data[-1] = data[-1]/50

# 将数据框转换为长格式
data_long <- melt(data, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

# 使用ggplot2创建折线图
p <- ggplot(data_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.25,1.85) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX_true" = "tan", 
                                "REMIX" = "#FF0000", 
                                "REMIX_P0.05" = "#FF00FF", 
                                "REMIX_P0.1" = "#777777", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX_true" = "solid", 
                                   "REMIX" = "solid", 
                                   "REMIX_P0.05" = "solid", 
                                   "REMIX_P0.1" = "solid", 
                                   #"BestSingle" = "dashed", 
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) 
  
# 保存图像
p
ggsave("Main_plot.png", plot = p, dpi = 300, width = 8, height = 6)

#######################################

df_mix0 <- data.frame(
  Index = seq(0.05,0.5,0.05),
  REMIX = c(21.20,27.20,33.35, 37.89, 41.23,44.89, 47.56, 49.02, 50.01, 50.86),
  #BestSingle = c(16.66, 26.98, 37.12, 47.70, 55.64, 58.32, 59.80, 61.25, 62.53, 63.47),
  SimpleAve = c(40.02, 42.08, 44.02, 46.10, 48.35, 50.07, 52.13, 54.34, 56.51, 58.29),
  # RA = c(19.78, 10.60578116, 15.10092016, 19.89342186, 24.52878898),
  Maximin = c(54.72, 54.77, 54.91, 54.93, 55.04, 55.07, 55.22, 55.37, 55.50, 55.37),
  RhoAve = c(19.82, 26.60, 33.06, 40.42, 47.32, 53.95, 61.41, 68.67, 75.20, 82.82)
)

df_mix0[-1] = df_mix0[-1]/50

# 将数据框转换为长格式
mix0_long <- melt(df_mix0, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_mix0 <- ggplot(mix0_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.75) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
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
ggsave("Main_plot_mix0.png", plot = p_mix0, dpi = 300, width = 8, height = 6)

#################################################
## mix1

df_mix1 <- data.frame(
  Index = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
  REMIX = c(21.48710136, 27.4897697, 32.70937432, 37.08811541, 40.71248378, 
            43.54266802, 45.20542468, 46.23549406, 47.36447967, 48.79567508),
  #BestSingle = c(18.63936757, 29.24077058, 35.94571352, 40.02834565, 43.16543054, 
  #       46.29769877, 49.46821983, 52.59578095, 55.73632207, 58.84139175),
  SimpleAve = c(22.74117946, 26.6493208, 30.50039511, 34.38464916, 38.2671985, 
         42.12616003, 46.02840709, 49.90890444, 53.78866405, 57.63340971),
  Maximin = c(34.47404217, 36.06907744, 37.62741921, 39.21226796, 40.79701611, 
         42.36005111, 43.94550913, 45.52581984, 47.13200544, 48.64908587),
  RhoAve = c(21.33602029, 27.99031054, 34.62651952, 41.2395427, 47.87275795, 
         54.48348736, 61.16704483, 67.8066522, 74.40791174, 81.03383266)
)


df_mix1[-1] = df_mix1[-1]/48

# 将数据框转换为长格式
mix1_long <- melt(df_mix1, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_mix1 <- ggplot(mix1_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.75) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  theme_bw()+ theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13))
   

p_mix1
ggsave("Main_plot_mix1.png", plot = p_mix1, dpi = 300, width = 8, height = 6)
## mix2

df_mix2 <- data.frame(
  Index = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
  REMIX = c(21.16466798, 26.34179251, 30.65817614, 34.40882817, 
           37.41439831, 39.77900474, 41.23376723, 41.98555327, 
           42.57390626, 42.93941966),
  # BestSingle = c(15.85685161, 23.46374795, 31.13392738, 38.76102827, 45.04068209, 47.4856588, 49.90115332, 52.27474033, 54.74613582, 57.10854388),
  SimpleAve = c(31.48969021, 33.2841731, 35.03174573, 36.82241361, 
         38.60598492, 40.3850516, 42.15557034, 43.87909277, 
         45.71582718, 47.46628445),
  Maximin = c(41.88019327, 42.02231802, 42.10652467, 42.22873982, 
         42.35588998, 42.47783872, 42.58489009, 42.67836754, 
         42.82902815, 42.91555856),
  RhoAve = c(19.72087579, 24.13174846, 28.39639391, 32.7484552, 
         37.0610322, 41.44059118, 45.75535914, 50.00525975, 
         54.42398381, 58.70577186)
)


df_mix2[-1] = df_mix2[-1]/40

# 将数据框转换为长格式
mix2_long <- melt(df_mix2, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_mix2 <- ggplot(mix2_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.5) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) 
  

p_mix2
ggsave("Main_plot_mix2.png", plot = p_mix2, dpi = 300, width = 8, height = 6)

## mix3
df_mix3 <- data.frame(
  Index = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
  REMIX = c(59.19474442, 59.47648462, 59.42879475, 59.1501148, 
           58.99058709, 58.69947816, 58.22301669, 57.5936125, 
           56.89771865, 56.26263566),
  #BestSingle = c(48.18344518, 52.34781629, 56.51092638, 60.70497767,  64.83988587, 68.95163644, 73.23634304, 77.29836587, 74.43841607, 72.22440807),
  SimpleAve = c(59.20775805, 60.24151914, 61.24544815, 62.28646484, 
         63.31270851, 64.24725215, 65.37065198, 66.3433663, 
         67.38854758, 68.45080664),
  Maximin = c(62.42888968, 61.74333487, 61.03759438, 60.36500051, 
         59.69293531, 58.91832835, 58.3103133, 57.62326497, 
         56.92969695, 56.28737647),
  RhoAve = c(59.59803586, 60.66489209, 61.72162603, 62.90758141, 
         64.08362061, 65.11266704, 66.38380803, 67.49576102, 
         68.64577974, 69.92453767)
)


df_mix3[-1] = df_mix3[-1]/60

# 将数据框转换为长格式
mix3_long <- melt(df_mix3, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_mix3 <- ggplot(mix3_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.7,1.4) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_mix3
ggsave("Main_plot_mix3.png", plot = p_mix3, dpi = 300, width = 8, height = 6)


# L2
df_L2 <- data.frame(
  Index = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
  REMIX = c(7.458512385, 7.698833652, 7.840798662, 7.901392027, 7.880774587, 
          7.795818893, 7.647269421, 7.553678783, 7.731320353, 7.902286826),
  #BestSingle = c(9.177873146, 10.23763637, 11.27108406, 12.29182429, 13.24834858, 13.78827215, 13.86132623, 13.79113505, 13.84748372, 13.86942549),
  SimpleAve = c(7.131033246, 7.120277725, 7.144185423, 7.142200448, 7.153670215, 
         7.147016647, 7.149041038, 7.122826714, 7.156976174, 7.176336815),
  Maximin = c(6.502270644, 6.646841063, 6.816203918, 6.969525436, 7.124157247, 
         7.269028794, 7.422463766, 7.53459566, 7.726081136, 7.895536308),
  RhoAve = c(8.059515737, 8.982816495, 9.876406009, 10.75552399, 11.60063427, 
         12.47848556, 13.31836654, 14.09183514, 15.01169042, 15.9350087)
)


df_L2[-1] = df_L2[-1]/9

# 将数据框转换为长格式
L2_long <- melt(df_L2, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_L2 <- ggplot(L2_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.7,1.8) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_L2
ggsave("Main_plot_L2.png", plot = p_L2, dpi = 300, width = 8, height = 6)


## L10
df_L10 <- data.frame(
  Index = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5),
  REMIX = c(19.60830629, 26.38578721, 32.32611946, 37.33531571, 
          41.57633021, 44.98795116, 47.29386591, 48.84558115, 
          50.11606535, 51.08579419),
  # BestSingle = c(15.63163768, 25.8625062, 36.58864116, 46.97821569, 55.77011521, 58.09853198, 59.2579017, 60.7383798,62.09905749, 63.59316113),
  SimpleAve = c(39.63903168, 41.74901345, 43.672971, 45.7927782, 
         47.8819245, 49.99454149, 51.87060984, 53.97662975, 
         56.10476182, 58.16418298),
  Maximin = c(54.21166159, 54.25918514, 54.0779366, 54.05377184, 
         53.98711305, 53.86184268, 53.61866736, 53.38584162, 
         53.38555791, 53.35861974),
  RhoAve = c(18.98232904, 25.93669119, 33.016942, 39.97411933, 
         46.97556501, 54.17398205, 60.983652, 68.32889652, 
         75.37084591, 82.24804733)
)

df_L10[-1] = df_L10[-1]/52

# 将数据框转换为长格式
L10_long <- melt(df_L10, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_L10 <- ggplot(L10_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.7) +
  labs(x = "Violation level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "RhoAve" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "RhoAve" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_L10
ggsave("Main_plot_L10.png", plot = p_L10, dpi = 300, width = 8, height = 6)


############################
df_high <- data.frame(
  Index = c(50, 100, 150, 200),
  REMIX = c(7.7129, 7.6297, 7.7256, 7.6415),
  SimpleAve = c(11.75,11.76,11.765,11.76),
  TransLasso = c(20.8009, 20.1519, 16.875, 17.02),
  TransGLM = c(25.147, 15.903, 13.108, 10.896), 
  PTL = c(10.113865, 9.364381 , 9.616546,  8.728799 )
)

df_high[-1] = df_high[-1]/8.2521

# 将数据框转换为长格式
high_long <- melt(df_high, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_high <- ggplot(high_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  labs(x = "Target sample size", y = "Standarized MSE") +
  scale_x_continuous(breaks = seq(50, 200, 50)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                "SimpleAve" = 'blue2',
                                "TransLasso" = "orange", 
                                "TransGLM" = "purple", 
                                "PTL" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   "TransLasso" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "TransGLM" = "dashed", 
                                   "PTL" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_high
ggsave("highdimTL.png", plot = p_high, dpi = 300, width = 7, height = 6)

#################

df_highw <- data.frame(
  Index = c(50, 100, 150, 200),
  REMIX = c(14.03, 14.02, 14.03, 14.04),
  SimpleAve = c(19.63,19.76,19.76,19.73),
  TransLasso = c(31.5156, 31.173, 30.01, 30.44),
  TransGLM = c(31.52, 25.96, 22.31, 20.04), 
  PTL = c(16.35652 ,15.91864 ,16.13212, 15.69885 )
)

df_highw[-1] = df_highw[-1]/8.2521

# 将数据框转换为长格式
highw_long <- melt(df_highw, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

p_highw <- ggplot(highw_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  labs(x = "Target sample size", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(50, 200, 50)) +
  scale_color_manual(values = c("REMIX" = "#FF0000", 
                                "SimpleAve" = 'blue2',
                                "TransLasso" = "orange", 
                                "TransGLM" = "purple", 
                                "PTL" = "green3")) +
  scale_linetype_manual(values = c("REMIX" = "solid", 
                                   "TransLasso" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "TransGLM" = "dashed", 
                                   "PTL" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_highw
ggsave("highdimTL_worst.png", plot = p_highw, dpi = 300, width = 7, height = 6)




















data_long <- melt(data, id.vars = "Index", variable.name = "Benchmarks", value.name = "Value")

# 使用ggplot2创建折线图
p <- ggplot(data_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), linewidth = 0.5) + 
  geom_point(size = 1) +
  labs(x = "Level of shift", y = "Worst Case LoBestSingle") +
  scale_x_continuous(breaks = seq(0, 3, 0.5)) +
  scale_color_manual(values = c("true" = "#FF0000", 
                                "BestSingle" = "#FFB90F", 
                                "SimpleAve" = "#8EBA42", 
                                "MI" = "#0000CC",
                                "PA" = "#00DD00")) +
  scale_linetype_manual(values = c("true" = "solid", 
                                   "BestSingle" = "dashed", 
                                   "SimpleAve" = "dashed", 
                                   "MI" = "dashed", 
                                   "PA" = "dashed")) 
# 保存图像
ggSimpleAveve("Main_plot_center.png", plot = p, dpi = 300, width = 8, height = 6)













##
data = read.csv("Acenter.csv")
data = data[,-5]
colnames(data) = c("Index","true","BestSingle","SimpleAve","MI","PA")
data




