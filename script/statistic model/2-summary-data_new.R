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
#####1.data old######### 
data<-read.csv("country_epi_climate_data_0229.csv")
unique(data$country)
unique(data$subregion)
unique(data$yearweek)
data$yearweek=yearweek(data$yearweek)
##1去掉香港
data=data[-which(data$country=="HK"|data$subregion=="Usall"),]
#2去掉变量：X,TEMP,RH，season,week,season2
data=data[,!(names(data) %in% c("X","TEMP","RH","season","week","season2"))]
#3去掉2023+2024
data=data[which(data$yearweek<yearweek("2023 W01")),]

#####2.data new#########
data_new<-read.csv("../data/country_epi_climate_data_ 2024-07-02 .csv")
unique(data_new$country)
unique(data_new$subregion)
unique(data_new$yearweek)
data_new$yearweek=yearweek(data_new$yearweek)
#1去掉TEMP，RH
data_new=data_new[,!(names(data_new) %in% c("TEMP","RH"))]

#####3. combine 1#########
colnames(data)
colnames(data_new)
data1=rbind(data,data_new)

#####4. combine climate#########
data<-read.csv("country_epi_climate_data_0229.csv")
data$yearweek=yearweek(data$yearweek)
data=data[which(data$yearweek<yearweek("2024 W01")),]
data=data[-which(data$country=="HK"|data$subregion=="Usall"),]
data=data[,c("subregion","yearweek","TEMP","RH")]


data_new<-read.csv("../data/country_epi_climate_data_ 2024-07-02 .csv")
data_new$yearweek=yearweek(data_new$yearweek)
data_new=data_new[complete.cases(data_new$TEMP),c("subregion","yearweek","TEMP","RH")]

data2=rbind(data,data_new)

#####5. combine all#########
finaldf=merge(data1,data2,by=c("subregion","yearweek"),all = T)
summary(finaldf[,3:ncol(finaldf)])

write.csv(finaldf,"combine-data.csv")

library(zoo)
library(imputeTS)
finaldf <- finaldf %>%
  group_by(subregion) %>%
  arrange(yearweek) %>%
  mutate(other = na.approx(other, method = "constant", rule = 3),
         omicron = na.approx(omicron, method = "constant", rule = 3),
        BA.2.75 = na.approx(BA.2.75, method = "constant", rule = 3),
        BA.2.86 = na.approx(BA.2.86, method = "constant", rule = 3),
         XBB = na.approx(XBB, method = "constant", rule = 3))

#####6. plot#########

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=MA3_cov))+
  facet_wrap(~subregion,)

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=MA3_Aflu))+
  facet_wrap(~subregion,)

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=MA3_Bflu))+
  facet_wrap(~subregion,)

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=MA3_rsv))+
  facet_wrap(~subregion,)

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=TEMP))+
  facet_wrap(~subregion,)

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=RH))+
  facet_wrap(~subregion,)

ggplot(finaldf)+
  geom_line(aes(x=yearweek,y=med365_cov))+
  facet_wrap(~subregion,)


data=finaldf

##add season
data$season=ifelse((data$MA3_ABflu>=5), 40, 0)
data$week=as.numeric(substr(data$yearweek,7,8))
data$season2=ifelse(data$week>=40|data$week<=13, 40, 0)



##export
library(writexl)

#write.csv(finaldf,"country_epi_climate_data_0704.csv")



###标准化####
#1. 4个病毒的阳性率
#2. 气象要素
#3. 免疫水平
#4. NPI

standardized_df <- data %>%
  group_by(subregion) %>%
  mutate(
    MA3_cov = (MA3_cov - min(MA3_cov)) / (max(MA3_cov)-min(MA3_cov))*100,
    MA3_Aflu= (MA3_Aflu - min(MA3_Aflu)) / (max(MA3_Aflu)-min(MA3_Aflu))*100,
    MA3_Bflu = (MA3_Aflu - min(MA3_Bflu)) / (max(MA3_Bflu)-min(MA3_Bflu))*100,
    MA3_rsv = (MA3_rsv - min(MA3_rsv)) / (max(MA3_rsv)-min(MA3_rsv))*100)



ggplot(standardized_df )+
  geom_line(aes(x=yearweek,y=MA3_cov))+
  facet_wrap(~subregion,)


ggplot(data )+
  geom_line(aes(x=yearweek,y=MA3_cov))+
  facet_wrap(~subregion,)


data=standardized_df


