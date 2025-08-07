library(readr)
library(readxl)
library(data.table)
library(dbplyr)
library(plyr)
library(gtsummary)
library(tidyverse)
library(tsibble)
library(reshape2)
library(imputeTS)
#####1.load data######### 
data<-read_excel("country_epi_climate_data.xlsx",sheet="all0229")

unique(data$subregion)
table(data$subregion)
table(data$yearweek)

data$yearweek=yearweek(data$yearweek)
#去除1:固定时间段
data=data[which(data$yearweek>yearweek("2021 W39") & 
                  data$yearweek<yearweek("2024 W06")),]

unique(data$subregion)

data[,4:21] <- apply(data[,4:21], 2, as.numeric)

summary(data[,4:21])
##去除2:香港不合理
data=data[complete.cases(data$MA3_cov_test),]


#write.csv(data,"data-lmy.csv")

#####2.climate######### 
library(readr)
library(purrr)
df_cli <-  list.files(path="../data/6_Climate/清理后数据", pattern="*.csv", full.names=TRUE) |>
  set_names() |>
  map_dfr(read_delim, .id = "file")

pattern <- "(?<=-).*?(?=.csv)"
df_cli$country <- str_extract(df_cli$file, pattern)
final_cli=df_cli[,c("country","region","week","TEMP","RH")]

#改名
# 使用ifelse()进行条件赋值
final_cli$region <- ifelse(is.na(final_cli$region), final_cli$country, final_cli$region)
unique(final_cli$region)
colnames(final_cli)=c("country","subregion","yearweek","TEMP","RH")

unique(final_cli$yearweek)
final_cli$yearweek=gsub("-"," ",final_cli$yearweek)

final_cli$country=gsub("USall","Usall",final_cli$country)
final_cli$subregion=gsub("USall","Usall",final_cli$subregion)

#合并
final_cli=final_cli[,-1]
final_cli=aggregate(cbind(TEMP,RH)~subregion+yearweek,final_cli,mean)
final_cli$yearweek=yearweek(final_cli$yearweek)
data_cli=merge(data,final_cli,by=c("subregion","yearweek"),all.x = T)




#####3.NPI######### 
df_npi<-read.csv("../data/4_NPI/All-NPI-0223.csv")
df_npi=df_npi[,-1]
table(df_npi$subregion)
df_npi$subregion=gsub("USall","Usall",df_npi$subregion)
df_npi$subregion=gsub("Hong Kong","HK",df_npi$subregion)
unique(df_npi$subregion)
df_npi$yearweek=yearweek(df_npi$yearweek)

#插值
wk1=yearweek("2023 W01")
ts=wk1+0:56
subregion=unique(df_npi$subregion)
newts=merge(ts,subregion)
colnames(newts)=c("yearweek","subregion")
final_npi=merge(df_npi,newts,by=c("subregion","yearweek"),all = T)


library(zoo)
library(imputeTS)
final_npi <- final_npi %>%
  group_by(subregion) %>%
  arrange(yearweek) %>%
  mutate(new_npi = na.approx(NPI, method = "constant", rule = 3))


ggplot(final_npi)+
  geom_line(aes(x=yearweek,y=NPI),color="red")+
  geom_line(aes(x=yearweek,y=new_npi),color="blue")+
  facet_wrap(~subregion,ncol=4)

unique(final_npi$subregion)
final_npi=final_npi[,-3]

data_cli$yearweek=yearweek(data_cli$yearweek)
data_cli_npi=merge(data_cli,final_npi,
                   by=c("subregion","yearweek"),all.x = T)




#####5.variant######### 
df_variant<-read.csv("../data/7_variant/变异株4.csv")
# 或者使用apply函数
par(mfrow=c(2,4))  # 设置图形布局为一行三列
apply(df_variant[,4:11], 2, function(column) hist(column,  xlab="Value", col="lightblue", border="black"))


df_variant$BA.2.75._new=df_variant$BA.2.75.+df_variant$CH.1.
df_variant$BA.2.86._new=df_variant$BA.2.86.+df_variant$JN.1.
df_variant$XBB_new=df_variant$XBB+df_variant$EG.5.
df_variant$yearweek=yearweek(as.Date(df_variant$Week.prior.to))

final_variant=df_variant[,c("yearweek","Country","Other",
                            "VOC.Omicron.GRA..B.1.1.529.BA...",
                            "BA.2.75._new",
                            "BA.2.86._new","XBB_new")]


final_variant.melt=melt(final_variant,id=c("yearweek","Country"))

##plot
ggplot(final_variant.melt)+
  geom_line(aes(x=yearweek,y=value,color=variable))+
  facet_wrap(~Country)



##合并
colnames(final_variant)=c("yearweek","country","other","omicron","BA.2.75","BA.2.86","XBB")
final_variant$country=gsub("USA","Usall",final_variant$country)
final_variant$country=gsub("Hong Kong","HK",final_variant$country)
final_variant$country=gsub("United Kingdom","UK",final_variant$country)
x=final_variant[final_variant$country=="Usall",]
x$country="US"

final_variant=rbind(final_variant,x)

data_cli_npi_variant=merge(data_cli_npi,final_variant,
                                 by=c("country","yearweek"),all.x = T)




#####6. influenza season######### 
data_cli_npi_variant$season=ifelse((data_cli_npi_variant$MA3_ABflu>=5), 40, 0)

data_cli_npi_variant$week=as.numeric(substr(data_cli_npi_variant$yearweek,7,8))
data_cli_npi_variant$season2=ifelse(data_cli_npi_variant$week>=40|data_cli_npi_variant$week<=13, 40, 0)

#香港全是流行季
data_cli_npi_variant$season2[data_cli_npi_variant$country=="HK" & (year(data_cli_npi_variant$yearweek)=="2023"|
                              year(data_cli_npi_variant$yearweek)=="2024")]=40



#####6. 导出数据######### 
#write.csv(data_cli_npi_variant,"country_epi_climate_data_0229.csv")


##plot
#不画美国的
data=data_cli_npi_variant
data=data[!grepl("Region", data$subregion), c("country", "yearweek","season2","subregion" ,
                                              "MA3_cov" ,"MA3_Aflu", "MA3_rsv"    )]

data <- data %>%
  arrange(subregion, yearweek)

data_s=data[which(data$season2==40),]

#lag
data$lag_cov_8 <- tsModel::Lag(data$MA3_cov, group = data$subregion,k = 8)
data$lag_iav_5 <- tsModel::Lag(data$MA3_Aflu, group = data$subregion,k = 5)
data$lag_iav_6 <- tsModel::Lag(data$MA3_Aflu, group = data$subregion,k = 6)
data$lag_rsv_8 <- tsModel::Lag(data$MA3_rsv, group = data$subregion,k = 8)


#plot1:iav_lag5 cov
ggplot(data)+
  geom_line(aes(x=yearweek,y=MA3_cov,colour="MA3_cov"))+
  geom_line(aes(x=yearweek,y=lag_iav_5,colour="lag_iav_5"))+
  scale_color_manual(name = "Y series", 
                     values = c("MA3_cov" = "blue", "lag_iav_5" = "red"),
                     labels = c("MA3_cov" = "cov", "lag_iav_5" = "iav_lag5"))+
  geom_line(aes(x=yearweek,y=season2))+
  ylab("Positive rate (%)")+xlab(NULL)+
   theme(legend.position = "bottom")+  
  facet_wrap(~subregion,nrow=5)



#plot2:iav cov
ggplot(data)+
  geom_line(aes(x=yearweek,y=MA3_cov,colour="MA3_cov"))+
  geom_line(aes(x=yearweek,y=MA3_Aflu,colour="MA3_iav"))+
  scale_color_manual(name = "Y series", 
                     values = c("MA3_cov" = "blue", "MA3_iav" = "red"),
                     labels = c("MA3_cov" = "cov", "MA3_iav" = "iav"))+
  geom_line(aes(x=yearweek,y=season2))+
  ylab("Positive rate (%)")+xlab(NULL)+
  theme(legend.position = "bottom")+  
  facet_wrap(~subregion,nrow=5)



#plot3:rsv_lag8 cov
ggplot(data)+
  geom_line(aes(x=yearweek,y=MA3_cov,colour="MA3_cov"))+
  geom_line(aes(x=yearweek,y=lag_rsv_8,colour="lag_rsv_8"))+
  scale_color_manual(name = "Y series", 
                     values = c("MA3_cov" = "blue", "lag_rsv_8" = "red"),
                     labels = c("MA3_cov" = "cov", "lag_rsv_8" = "rsv_lag8"))+
  geom_line(aes(x=yearweek,y=season2))+
  ylab("Positive rate (%)")+xlab(NULL)+
  theme(legend.position = "bottom")+  
  facet_wrap(~subregion,nrow=5)

cor.test(data_s$lag_rsv_8,data_s$MA3_cov)


#plot4:rsv cov
ggplot(data)+
  geom_line(aes(x=yearweek,y=MA3_cov,colour="MA3_cov"))+
  geom_line(aes(x=yearweek,y=MA3_Aflu,colour="MA3_rsv"))+
  scale_color_manual(name = "Y series", 
                     values = c("MA3_cov" = "blue", "MA3_rsv" = "red"),
                     labels = c("MA3_cov" = "cov", "MA3_rsv" = "rsv"))+
  geom_line(aes(x=yearweek,y=season2))+
  ylab("Positive rate (%)")+xlab(NULL)+
  theme(legend.position = "bottom")+  
  facet_wrap(~subregion,nrow=5)

cor.test(data_s$MA3_cov,data_s$MA3_rsv)



#plot2:rsv_lag8 iav
ggplot(data)+
     geom_line(aes(x=yearweek,y=MA3_Aflu,colour="MA3_Aflu"))+
     geom_line(aes(x=yearweek,y=lag_rsv_8,colour="lag_rsv_8"))+
     scale_color_manual(name = "Y series", 
                        values = c("MA3_Aflu" = "blue", "lag_rsv_8" = "red"),
                       labels = c("MA3_Aflu" = "IAV", "lag_rsv_8" = "rsv_lag8"))+
     geom_line(aes(x=yearweek,y=season2))+
     ylab("Positive rate (%)")+xlab(NULL)+
  theme(legend.position = "bottom")+ 
     facet_wrap(~subregion,nrow=5)
cor.test(data_s$MA3_rsv,data_s$lag_rsv_8)


#plot2:iav_lag6 rsv
ggplot(data)+
  geom_line(aes(x=yearweek,y=MA3_rsv,colour="MA3_rsv"))+
  geom_line(aes(x=yearweek,y=lag_iav_6,colour="lag_iav_6"))+
  scale_color_manual(name = "Y series", 
                     values = c("MA3_rsv" = "blue", "lag_iav_6" = "red"),
                     labels = c("MA3_rsv" = "RSV", "lag_iav_6" = "iav_lag6"))+
  geom_line(aes(x=yearweek,y=season2))+
  ylab("Positive rate (%)")+xlab(NULL)+
  theme(legend.position = "bottom")+ 
  facet_wrap(~subregion,nrow=5)

cor.test(data_s$MA3_rsv,data_s$lag_iav_6)
