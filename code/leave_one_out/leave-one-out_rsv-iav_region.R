rm(list = ls())
library(INLA)
library(readxl)
library(tsModel)
library(ggplot2)
library(reshape2)
library(dplyr)
library(zoo)
library(ggcorrplot)
library(corrplot)
library(GGally)
library(dlnm)
library(data.table)
library(PerformanceAnalytics)
library(imputeTS)
library(dbplyr)
library(tsibble)
library(RColorBrewer)
library(stats)
library(geofacet)
library(ModelMetrics)

my_rbind <- function(df1, df2){    # Create own merging function
  rbind(df1, df2)
}

#####################1. load data###########################################
source("1-summary-data.R")
data=data[!grepl("Region", data$subregion), ]
allcountry=data

region=unique(allcountry$subregion)


for (j in 5) {
  
  
data=allcountry[allcountry$subregion==region[j],]
data$region=as.integer(factor(data$subregion))



data$MA3_cov_log=log(data$MA3_cov+0.001)
data$MA3_Aflu_log=log(data$MA3_Aflu+0.001)
data$MA3_Bflu_log=log(data$MA3_Bflu+0.001)
data$MA3_ABflu_log=log(data$MA3_ABflu+0.001)
data$MA3_rsv_log=log(data$MA3_rsv+0.001)






########################3.Create lagged variables###########################
data <- data %>%
  arrange(subregion, yearweek)
# set maximum lag
nlag = 8

#virus
lag_cov <- tsModel::Lag(data$MA3_cov, group = data$subregion,k = 0:nlag)
lag_iav <- tsModel::Lag(data$MA3_Aflu,  group = data$subregion,k = 0:nlag)
lag_ibv <- tsModel::Lag(data$MA3_Bflu, group = data$subregion, k = 0:nlag)
lag_iabv <- tsModel::Lag(data$MA3_ABflu, group = data$subregion, k = 0:nlag)
lag_rsv <- tsModel::Lag(data$MA3_rsv, group = data$subregion, k = 0:nlag)

# Remove lag missing value
lag_cov <- lag_cov[complete.cases(lag_cov),]
lag_iav <- lag_iav[complete.cases(lag_iav),]
lag_ibv <- lag_ibv[complete.cases(lag_ibv),]
lag_iabv <- lag_iabv[complete.cases(lag_iabv),]
lag_rsv <- lag_rsv[complete.cases(lag_rsv),]


#data <- data[data$week > yearweek("2022 W08"),]


data <- data %>%
  group_by(subregion) %>%
  filter(row_number() > 8)
#########################4.define cross-basis matrix###########################
# set lag knots
lagknot = equalknots(0:nlag, 2)

# COVID-19
var <- lag_cov
basis_cov <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(data$MA3_cov, 2)),
                        arglag = list(fun = "ns", knots = lagknot))

# Influenza
var <- lag_iav
basis_iav<- crossbasis(var,
                       argvar = list(fun = "ns", knots = equalknots(data$MA3_Aflu, 2)),
                       arglag = list(fun = "ns", knots = lagknot))
var <- lag_ibv
basis_ibv<- crossbasis(var,
                       argvar = list(fun = "ns", knots = equalknots(data$MA3_Bflu, 2)),
                       arglag = list(fun = "ns", knots = lagknot))
var <- lag_iabv
basis_iabv<- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(data$MA3_ABflu, 2)),
                        arglag = list(fun = "ns", knots = lagknot))

# Rsv
var <- lag_rsv
basis_rsv <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(data$MA3_rsv, 2)),
                        arglag = list(fun = "ns", knots = lagknot))




# assign unique column names to cross-basis matrix for inla() model
colnames(basis_cov) = paste0("basis_cov.", colnames(basis_cov))
colnames(basis_iav) = paste0("basis_iav.", colnames(basis_iav))
colnames(basis_ibv) = paste0("basis_ibv.", colnames(basis_ibv))
colnames(basis_iabv) = paste0("basis_iabv.", colnames(basis_iabv))
colnames(basis_rsv) = paste0("basis_rsv.", colnames(basis_rsv))


########################5.create dataframe for model testing###########################

# set data for models
Y_COV  <- data$MA3_cov_log # response variable
Y_IAV  <- data$MA3_Aflu_log # response variable
Y_IBV  <- data$MA3_Bflu_log# response variable
Y_IABV  <- data$MA3_ABflu_log# response variable
Y_RSV  <- data$MA3_rsv_log # response variable
T1 <- as.numeric(substring(data$yearweek,7,8)) # for random effect to account for annual cycle (seasonality)
T2 <-as.numeric(substring(data$yearweek,1,4))-2020 # for random effect to account for inter-annual variability
S<-as.numeric(gsub("Region ","",data$subregion)) # for random effect to account for spatial

temp<-data$TEMP
rh<-data$RH
npi<-data$new_npi



med365_cov<-data$med365_cov
med365_Aflu<-data$med365_Aflu
med365_Bflu<-data$med365_Bflu
med365_ABflu<-data$med365_ABflu
med365_rsv<-data$med365_rsv
week=data$yearweek

df <- data.frame(Y_COV,Y_IAV,Y_IBV,Y_IABV,Y_RSV,T1,T2,S,week,
                 temp,rh,npi,
                 med365_rsv,med365_ABflu,med365_Bflu,med365_Aflu,med365_cov)

nyear=length(unique(year(data$yearweek)))
nregion=length(unique(data$region))

# include formula and set defaults for data, family (to allow other prob dist models e.g. Poisson) and config (to allow for sampling)
mymodel <- function(formula, data = df, family = "gaussian", config = FALSE)
  
{
  model <- inla(formula = formula, data = data, family = family, 
                control.inla = list(strategy = 'adaptive'), 
                control.compute = list(dic = TRUE, config = config, 
                                       cpo = TRUE, return.marginals = FALSE),
                control.fixed = list(correlation.matrix = TRUE, 
                                     prec.intercept = 1, prec = 1),
                control.predictor = list(link = 1, compute = TRUE), 
                verbose = FALSE)
  model <- inla.rerun(model)
  return(model)
}

# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))


}

######################## cross validation###########################


#best model
baseformula <- Y_RSV~ 1 + 
  f(T1,  model = "rw1", cyclic = TRUE, constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior) 
formula0.11<- update.formula(baseformula, ~. + rh +npi+med365_rsv)
formula1.1 <- update.formula(formula0.11, ~. + basis_iav)

df$Y_RSV_1=df$Y_RSV



for(i in 1:nrow(df)){
  
  # Step 1: rerun the selected model (fitted with config = TRUE for sampling) 
  # Step 2: produce cross-validated posterior predictive samples leaving out one year and one week at a time
  
  # define number of samples
  s <- 1000
  # replace data in testing period with NA for out of sample prediction
  casestopred <- df$Y_RSV_1 # response variable
  idx.pred <- i
  casestopred[idx.pred] <- NA # replace cases in year and month of interest to NA
  mpred <- length(idx.pred)
  
  # set response variable and year indicator
  df$Y_RSV<- casestopred
  
  mymodel <-function(formula, data = df, family = "gaussian", config = FALSE){
    
    
    model <- inla(formula = formula, data = df, family = family, 
                  control.inla = list(strategy = 'adaptive'), 
                  control.compute = list(dic = TRUE, config = TRUE, 
                                         cpo = TRUE, return.marginals = FALSE),
                  control.fixed = list(correlation.matrix = TRUE, 
                                       prec.intercept = 1, prec = 1),
                  control.predictor = list(link = 1, compute = TRUE), 
                  verbose = FALSE)
    model <- inla.rerun(model)
    return(model)
  }
  
  
  ##最优模型
  model <- mymodel( formula1.1 , df)
  
  
  #生成1000个后验分布的样本
  xx <- inla.posterior.sample(s, model) 
  
  #从后验分布提取信息
  xx.s <- inla.posterior.sample.eval(function(...) c(theta[1], Predictor[idx.pred]), xx)
  
  #存储分布的随机数
  y.pred <- matrix(NA, mpred, s)
  
  #
  for(s.idx in 1:s) {
    xx.sample <- xx.s[, s.idx]##提取第N次抽样的后验参数信息
    y.pred[, s.idx] <- rnorm(mpred, mean = (xx.sample[-1]), sd = xx.sample[1])
  }
  
  
  preds <- list(year = 2020 + df[i,"T2"], week = df[i,"T1"], idx.pred = idx.pred, 
                mean = apply(y.pred, 1, mean), median = apply(y.pred, 1, median),
                lci = apply(y.pred, 1, quantile, probs = c(0.025)),
                uci = apply(y.pred, 1, quantile, probs = c(0.975)))
  save(preds, file = paste0("../output_new/preds_rsv_iav_",region[j],"/preds_",2020 + df[i,"T2"], "_",  df[i,"T1"], ".RData"))
  
}



df$pred.mean <- NA
df$pred.median <- NA
df$pred.lci <- NA
df$pred.uci <- NA


for(i in 1:nrow(df)){    
  
  load(paste0("../output_new/preds_rsv_iav_",region[j],"/preds_",2020 +  df[i,"T2"],"_", df[i,"T1"],".RData"))
  
  df$pred.mean[preds$idx.pred] <- preds$mean
  df$pred.median[preds$idx.pred] <- preds$median
  df$pred.lci[preds$idx.pred] <- preds$lci
  df$pred.uci[preds$idx.pred] <- preds$uci
  
}




# plot observed v fitted 
ggplot(df) + 
  geom_ribbon(aes(x = week, ymin = pred.lci, ymax = pred.uci), 
              fill = "#D37295", alpha = 0.5) + 
  geom_line(aes(x =week, y = pred.mean, col = "#D37295")) +
  geom_line(aes(x = week, y = Y_RSV_1, col = "#499894")) +
  scale_colour_identity(name = "",
                        breaks = c("#499894", "#D37295"),
                        labels = c("Observed", "Posterior predicted"),
                        guide = "legend") +
  theme_bw() + 
  ylab("RSV positivity rate (standardized, log)")+
  xlab(NULL)+
  theme(axis.text.x = element_text(size = 10, color="black",vjust = 0.5, angle = 90))



#######
write.csv(df,paste0("../output_new/preds_rsv_iav_",region[j],"/rsv-data.csv"))



###combine 5 results
df1<-read.csv("../output_new/preds_rsv_iav_Denmark/rsv-data.csv")
df2<-read.csv("../output_new/preds_rsv_iav_Ireland/rsv-data.csv")
df3<-read.csv("../output_new/preds_rsv_iav_Portugal/rsv-data.csv")
df4<-read.csv("../output_new/preds_rsv_iav_Slovenia/rsv-data.csv")
df5<-read.csv("../output_new/preds_rsv_iav_England/rsv-data.csv")

df1$subregion="Denmark"
df2$subregion="Ireland"
df3$subregion="Portugal"
df4$subregion="Slovenia"
df5$subregion="England"

finaldf=rbind(df1,df2,df3,df4,df5)

finaldf$yearweek=yearweek(finaldf$week)
# plot observed v fitted 
p<-ggplot(finaldf) + 
  geom_ribbon(aes(x = yearweek, ymin = pred.lci, ymax = pred.uci), 
              fill = "#D37295", alpha = 0.5) + 
  geom_line(aes(x =yearweek, y = pred.mean, col = "#D37295",group=subregion)) +
  geom_line(aes(x = yearweek, y = Y_RSV_1, col = "#499894",group=subregion)) +
  scale_colour_identity(name = "",
                        breaks = c("#499894", "#D37295"),
                        labels = c("Observed", "Posterior predicted"),
                        guide = "legend") +
  theme_bw() + 
  ylab("RSV positivity rate (standardized, log)")+
  xlab(NULL)+
  theme(axis.text.x = element_text(size = 10, color="black",vjust = 0.5, angle = 90))+
  facet_wrap( ~ subregion, scales = "free") +
  theme(legend.position = "bottom") 


ggsave("../figure/rsv_iav_5countries.tiff", p, width=20, height=15, units="cm", dpi=300, compression = "lzw")


