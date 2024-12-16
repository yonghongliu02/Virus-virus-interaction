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

##read date
finaldf<-read.csv("../data/combine-data.csv")
finaldf$yearweek=yearweek(finaldf$yearweek)

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
  geom_line(aes(x=yearweek,y=MA3_ABflu))+
  geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
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
data$week=as.numeric(substr(data$yearweek,7,8))
data$season2=ifelse(data$week>=40|data$week<=13, 40, 0)




###standardized####
standardized_df <- data %>%
  group_by(subregion) %>%
  mutate(
    MA3_cov = (MA3_cov - min(MA3_cov)) / (max(MA3_cov)-min(MA3_cov))*100,
    MA3_Aflu= (MA3_Aflu - min(MA3_Aflu)) / (max(MA3_Aflu)-min(MA3_Aflu))*100,
    MA3_Bflu = (MA3_Bflu - min(MA3_Bflu)) / (max(MA3_Bflu)-min(MA3_Bflu))*100,
    MA3_rsv = (MA3_rsv - min(MA3_rsv)) / (max(MA3_rsv)-min(MA3_rsv))*100)



ggplot(standardized_df )+
  geom_line(aes(x=yearweek,y=MA3_cov))+
  facet_wrap(~subregion,)


ggplot(data )+
  geom_line(aes(x=yearweek,y=MA3_cov))+
  facet_wrap(~subregion,)


data=standardized_df


