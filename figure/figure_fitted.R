# 加载所需包
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
options(scipen = 200)
####################################figure4##############################
# 读取两个 sheet 数据
df_A <- read_excel("fitting_phase_all_5401.xlsx", sheet = "流感拟合结果") %>%
  mutate(disease = "IAV")
df_B <- read_excel("fitting_phase_all_5401.xlsx", sheet = "新冠拟合结果") %>%
  mutate(disease = "SARS-Cov-2")

df_A$Date=as.Date(df_A$Date)
df_A[, 2:5] <- df_A[, 2:5] / 21

df_B$Date=as.Date(df_B$Date)
df_B[, 2:5] <- df_B[, 2:5] / 21


# 流感
p1<-ggplot(df_A, aes(x = Date)) +
  geom_line(aes(y = Observed), color = "black", size = 0.5) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = "#f4a582", alpha = 0.3) +
  geom_line(aes(y =Fitted), color = "#ca0020", size = 0.5) +
  theme_bw() + theme(  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  scale_x_date(breaks = "3 months")+
  labs(x =NULL, y = "IAV incidence\n(1/1000000)")+
  theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 45),
        axis.text.y = element_text(size = 10, color="black",vjust = 0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text( size=12))



# 新冠
p2<-ggplot(df_B, aes(x = Date)) +
  geom_line(aes(y = Observed), color = "black", size = 0.5) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = "#92c5de", alpha = 0.3) +
  geom_line(aes(y =Fitted), color = "#0571b0", size = 0.5) +
  theme_bw() + theme(  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  labs(x = NULL, y = "SARS-Cov-2 incidence\n(1/1000000)")+
  scale_x_date(breaks = "3 months")+
  theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 45),
        axis.text.y = element_text(size = 10, color="black",vjust = 0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text( size=12))



library(ggpubr)
p<-ggarrange(  p1,p2,   nrow =2 , align = 'v',labels = c("A","B"))
p
ggsave("figure 4.pdf", plot = p, width=16, height=12, units="cm")
ggsave("figure 4.tiff", plot = p, width=16, height=15, units="cm", dpi=300, compression = "lzw")



####################################figure simulation##############################
# 读取两个 sheet 数据
df_A <- read_excel("fitting_phase_all_simulation.xlsx", sheet = "流感拟合结果") %>%
  mutate(disease = "IAV")
df_B <- read_excel("fitting_phase_all_simulation.xlsx", sheet = "新冠拟合结果") %>%
  mutate(disease = "SARS-Cov-2")

df_A$Date=as.Date(df_A$Date)
df_A[, 2:5] <- df_A[, 2:5] / 21

df_B$Date=as.Date(df_B$Date)
df_B[, 2:5] <- df_B[, 2:5] / 21


# 流感
p1<-ggplot(df_A, aes(x = Date)) +
  geom_line(aes(y = Observed), color = "black", size = 0.5) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = "#f4a582", alpha = 0.3) +
  geom_line(aes(y =Fitted), color = "#ca0020", size = 0.5) +
  theme_bw() + theme(  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  scale_x_date(breaks = "3 months")+
  labs(x =NULL, y = "IAV incidence\n(1/1000000)")+
  theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 45),
        axis.text.y = element_text(size = 10, color="black",vjust = 0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text( size=12))



# 新冠
p2<-ggplot(df_B, aes(x = Date)) +
  geom_line(aes(y = Observed), color = "black", size = 0.5) +
  geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), fill = "#92c5de", alpha = 0.3) +
  geom_line(aes(y =Fitted), color = "#0571b0", size = 0.5) +
  theme_bw() + theme(  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  labs(x = NULL, y = "SARS-Cov-2 incidence\n(1/1000000)")+
  scale_x_date(breaks = "3 months")+
  theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 45),
        axis.text.y = element_text(size = 10, color="black",vjust = 0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text( size=12))



library(ggpubr)
p<-ggarrange(  p1,p2,   nrow =2 , align = 'v',labels = c("A","B"))
p

ggsave("figure 4S.tiff", plot = p, width=16, height=15, units="cm", dpi=300, compression = "lzw")





####################################MCMC plot##############################
chain1 <- read_excel("fitting_phase_all_5401.xlsx", sheet = "MCMC链条") %>%
  mutate(iter = 1:n(), run = "Run1")

chain2 <- read_excel("fitting_phase_all_5402.xlsx", sheet = "MCMC链条") %>%
  mutate(iter = 1:n(), run = "Run2")

chain3 <- read_excel("fitting_phase_all_5403.xlsx", sheet = "MCMC链条") %>%
  mutate(iter = 1:n(), run = "Run3")



long_chain <- bind_rows(chain1, chain2, chain3) %>%
  pivot_longer(cols = -c(iter, run),
               names_to = "parameter",
               values_to = "value")

# 指定参数顺序
param_order <- c("IS"  ,   "RS"   ,  "SI" ,    "SR"   ,  
                 "RR", "B2","c2","d2",
                 "sigma1", "sigma2" ,"p1"  ,   "p2"  )
long_chain$parameter <- factor(long_chain$parameter, levels = param_order)
param_labels <- c(

  IS     = "IS",
  RS     = "RS",
  SI     = "SI",
  SR     = "SR",
  RR     = "RR",
  B2="B[2]",
  c2="c[2]",
  d2="d[2]",
  sigma1 = "sigma[1]",
  sigma2 = "sigma[2]",
  p1     = "1/rho[1]",
  p2     = "1/rho[2]"
)

y_limits <- tibble::tibble(
  parameter = factor(c("IS", "RS", "SI", "SR", "RR", "B2", "c2", "d2", "sigma1", "sigma2", "p1", "p2"),
                     levels = param_order),
  ymin = c(750, 16795000, 2600, 200, 800,  1.1,  0.07,  -8,  -0.4,   -1,  0,   30),
  ymax = c(1250, 16805000, 3400, 800, 1200, 1.2,  0.1,  -7.8,  0.4,   -0.7,  15,   50)
)

# 第二步：合并 y 限制到原数据
long_chain_limited <- long_chain %>%
  left_join(y_limits, by = "parameter")

plot<-ggplot(long_chain_limited, aes(x = iter, y = value, color = run)) +
  geom_point(alpha = 0.6,size=0.1) +
  # 利用 geom_blank 设置 y 轴范围
  geom_blank(aes(y = ymin)) +
  geom_blank(aes(y = ymax)) +
  scale_x_continuous(expand = c(0, 0),breaks=c(0,nrow(chain1)))+
  facet_wrap(~parameter, scales = "free_y", ncol = 3,
             labeller = as_labeller(param_labels, label_parsed)) +
  theme_bw() + theme(  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
  scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fb9a99")) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 12)) +
  labs(x = "Iteration", y = NULL, color = "Chain")+
  theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 0),
        axis.text.y = element_text(size = 10, color="black",vjust = 0.5),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text( size=12))

#ggsave("mcmcplot.pdf", plot = plot, width=16, height=12, units="cm")
ggsave("mcmcplotnew.tiff", plot = plot, width=16, height=20, units="cm", dpi=300, compression = "lzw")

