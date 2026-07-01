library(ggplot2)
library(reshape2)

# This script is for replicating (most of) the figures in DORM paper
# Processed and summarized simulation results are stored in .rds files

# Read data
data_long <- readRDS("data_long.rds") 

# Figure 1
p <- ggplot(data_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.25,1.85) +
  labs(x = "Decoupling level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM_true" = "tan", 
                                "DORM" = "#FF0000", 
                                "DORM_P0.05" = "#FF00FF", 
                                "DORM_P0.1" = "#777777", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue", 
                                "Maximin" = "purple",
                                "Joint_mixing" = "green3")) +
  scale_linetype_manual(values = c("DORM_true" = "solid", 
                                   "DORM" = "solid", 
                                   "DORM_P0.05" = "solid", 
                                   "DORM_P0.1" = "solid", 
                                   #"BestSingle" = "dashed", 
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "Joint_mixing" = "dashed")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) 

p
ggsave("Main_plot.png", plot = p, dpi = 300, width = 8, height = 6)

#######################################

# Different rho
# Mix 0: 0.5,0.5

mix0_long <- readRDS("mix0_long.rds") 

p_mix0 <- ggplot(mix0_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.75) +
  labs(x = "Decoupling level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "Joint_mixing" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "Joint_mixing" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13))


p_mix0
ggsave("Main_plot_mix0.png", plot = p_mix0, dpi = 300, width = 8, height = 6)

#################################################
## mix1: 0.8,0.2

mix1_long <- readRDS("mix1_long.rds") 

p_mix1 <- ggplot(mix1_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.75) +
  labs(x = "Decoupling level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "Joint_mixing" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "Joint_mixing" = "dashed")) +
  theme_bw()+ theme(axis.title = element_text(size = 15),
                    axis.text = element_text(size = 13),
                    legend.title = element_text(size = 14),
                    legend.text = element_text(size = 13))


p_mix1
ggsave("Main_plot_mix1.png", plot = p_mix1, dpi = 300, width = 8, height = 6)

## mix2: all 1/3

mix2_long <- readRDS("mix2_long.rds") 

p_mix2 <- ggplot(mix2_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.5) +
  labs(x = "Decoupling level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "Joint_mixing" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "Joint_mixing" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 


p_mix2
ggsave("Main_plot_mix2.png", plot = p_mix2, dpi = 300, width = 8, height = 6)

## mix3: all 1/5

mix3_long <- readRDS("mix3_long.rds") 

p_mix3 <- ggplot(mix3_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.7,1.4) +
  labs(x = "Decoupling level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "Joint_mixing" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "Joint_mixing" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_mix3
ggsave("Main_plot_mix3.png", plot = p_mix3, dpi = 300, width = 8, height = 6)


## L = 10

L10_long <- readRDS("L10_long.rds") 

p_L10 <- ggplot(L10_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  ylim(0.3,1.7) +
  labs(x = "Decoupling level s*", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(0, 0.5, 0.05)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                #"BestSingle" = "orange", 
                                "SimpleAve" = "blue2", 
                                "Maximin" = "purple",
                                "Joint_mixing" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
                                   #"BestSingle" = "dashed",
                                   "SimpleAve" = "dashed", 
                                   "Maximin" = "dashed", 
                                   "Joint_mixing" = "dashed")) +
  #geom_hline(yintercept = 0, linetype = "solid", color = "gray", size = 1) +
  theme_bw() + theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 13),
                     legend.title = element_text(size = 14),
                     legend.text = element_text(size = 13)) 

p_L10
ggsave("Main_plot_L10.png", plot = p_L10, dpi = 300, width = 8, height = 6)


############################
# High dimension: on target

high_long <- readRDS("high_long.rds") 

p_high <- ggplot(high_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  labs(x = "Target sample size", y = "Standarized MSE") +
  scale_x_continuous(breaks = seq(50, 200, 50)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                "SimpleAve" = 'blue2',
                                "TransLasso" = "orange", 
                                "TransGLM" = "purple", 
                                "PTL" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
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
# High dimension: worst case

highw_long <- readRDS("highw_long.rds") 

p_highw <- ggplot(highw_long, aes(x = Index, y = Value, color = Benchmarks, group = Benchmarks)) +
  geom_line(aes(linetype = Benchmarks), size = 0.5) + 
  geom_point(size = 1) +
  labs(x = "Target sample size", y = "Worst Case Standarized MSE") +
  scale_x_continuous(breaks = seq(50, 200, 50)) +
  scale_color_manual(values = c("DORM" = "#FF0000", 
                                "SimpleAve" = 'blue2',
                                "TransLasso" = "orange", 
                                "TransGLM" = "purple", 
                                "PTL" = "green3")) +
  scale_linetype_manual(values = c("DORM" = "solid", 
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
