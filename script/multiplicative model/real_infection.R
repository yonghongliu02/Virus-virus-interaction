library(readxl)
library(tsibble)
library(reshape2)
library(ggplot2)
library(dplyr)
library(zoo)
#周数据
df_all<-read.csv("../beijing-all-1128.csv")
#日数据
df_ILI<-read_excel("../beijing.xlsx")
df_ILI$date=as.Date(df_ILI$date) 
df_ILI$week=yearweek(df_ILI$date)
df_all$week=yearweek(df_all$week)

#————————————————————————————————————————————————————————————————————————
#感染人数1:=法定报告病例数/（显性感染率*就诊率*诊断率*采样检测成功率）
#感染人数2:=C*ILI占比*阳性率=一段时间的患病率


#—————————————————————#方法1乘法模型:估算新冠的感染人数———————————————————————####

report_cases=df_ILI$flu_reported

#有症状比例：p_sym=58.3-74.5
A_min <- 0.583 
A_max <- 0.745

#显性中ILI比例：p_ILI=26.0-42.0
B_min <- 0.26
B_max <- 0.42

#就诊率：p_visit=48.2-76.2
C_min <- 0.482
C_max <- 0.762

#检测率：p_sampling=80-90
D_min <- 0.8
D_max <- 0.9

#灵敏度：p_test=95—100
E_min <- 0.95
E_max <- 1

# 蒙特卡罗模拟的次数
n_simulations <- 10000 


flu_cases_simulated <- matrix(NA,nrow=nrow(df_ILI),ncol=n_simulations)


set.seed(123) 
for(i in 1:nrow(df_ILI)){
for (j in 1:n_simulations) {
  A <- runif(1, A_min, A_max)
  B <- runif(1, B_min, B_max)
  C <- runif(1, C_min, C_max)
  D <- runif(1, D_min, D_max)
  E <- runif(1, E_min, E_max)
  
  flu_cases_simulated[i,j] <- report_cases[i] /(A * B * C * D *E)
}
}

row_means <- rowMeans(flu_cases_simulated)
ci_95_rows <- apply(flu_cases_simulated, 1, function(x) {
  t_test <- t.test(x)
  return(c(t_test$conf.int))
})


flu_cases_simulated_final=cbind(data.frame(row_means),t(data.frame(ci_95_rows)))
colnames(flu_cases_simulated_final)=c("flu_estimated","flu_lower","flu_upper")

#plot
df_ILI=cbind(df_ILI,flu_cases_simulated_final)


df_ILI_melt=melt(df_ILI[,c(1,2,3,7)],"date")
#df_ILI_melt=df_ILI_melt[which(df_ILI_melt$date>="2023-1-1"),]

ggplot(df_ILI_melt,aes(x=date,y=value,fill=variable))+
  geom_line()+
  facet_wrap(~variable,scales = "free_y")


#write.csv(df_ILI[,c(1,2,3,7)],"estimated_cases.csv",row.names = F)








#———————————————————————————#感染人数2:=C*ILI占比*阳性率—————————————————————————####
#阳性率差值
df_all$date=as.Date(df_all$week)
# 构造从开始到结束的每日序列
daily_dates <- data.frame(date = seq(min(df_all$date), max(df_all$date), by = "day"))


# 样条插值
interp <- spline(x = as.numeric(df_all$date), y = df_all$IAV, xout = as.numeric(daily_dates$date))

non_negative_values <- pmax(interp$y, 0) 

daily_df <- data.frame(
  date = as.Date(interp$x, origin = "1970-01-01"),
  IAV = non_negative_values
)




#合并数据集
df_ILI=merge(df_ILI,daily_df,by="date",all.x = T)

#插值
df_ILI$IAV[is.na(df_ILI$IAV)]=0

df_ILI$flu_infection=df_ILI$ILI/df_ILI$visits*df_ILI$IAV

#平滑
df_ILI$flu_infection <- rollmean(df_ILI$flu_infection, k = 7, fill = NA, align = "center")


beta=95029/6.260538465#按照乘法模型的病例数放大系数

df_ILI$flu_infection_estimated=df_ILI$flu_infection*beta



#plot
df_ILI_melt=melt(df_ILI[,c(1,2,3,7,10,11,12)],"date")
#df_ILI_melt=df_ILI_melt[which(df_ILI_melt$date>="2023-1-1"),]

ggplot(df_ILI_melt,aes(x=date,y=value,fill=variable))+
  geom_line()+
  facet_wrap(~variable,scales = "free_y")



df_ILI$flu_infection_estimated[is.na(df_ILI$flu_infection_estimated)]=0
write.csv(df_ILI[,c(1,2,3,7,12)],"estimated_cases0731.csv",row.names = F)

