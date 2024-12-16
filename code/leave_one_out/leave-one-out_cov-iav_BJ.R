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
data=read.csv("../data/beijing-all-1128.csv")
data$yearweek=yearweek(data$week)

data$RSV=0

p<-ggplot(data)+
  geom_line(aes(x=yearweek,y=COV,color="COV"))+
  geom_line(aes(x=yearweek,y=IAV,color="IAV"))+
  geom_line(aes(x=yearweek,y=IBV,color="IBV"))+
  geom_line(aes(x=yearweek,y=RSV,color="RSV"))+
  theme_minimal() +
  scale_color_manual(breaks=c("COV","IAV","IBV","RSV"),
                     labels = c('SARS-CoV-2', 'IAV',"IBV","RSV"), 
                     values = c('#1f78b4', '#e31a1c', '#eeb479','#33a02c'))+
  ylab("Positivity rate (%)")+xlab(NULL)+
  scale_y_continuous(limits=c(0,70),expand = c(0, 0))+
  scale_x_yearweek(date_breaks="10 weeks",expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color="black",hjust = 1, angle = 60),
        axis.text.y = element_text(size = 9, color="black",vjust = 0.1),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text( size=10))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1.5),'lines'))+
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(col = guide_legend(ncol = 2))+
  theme(legend.title = element_blank())+  
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)
  
p

ggsave("../figure/BJ2.tiff", p, width=18, height=4.3, units="cm", dpi=300, compression = "lzw")


#############

data <- data %>%
  mutate(
    COV = (COV - min(COV)) / (max(COV)-min(COV))*100,
    IAV= (IAV - min(IAV)) / (max(IAV)-min(IAV))*100)
   



hist(data$COV)
hist(data$IAV)

data$COV_log=log(data$COV+0.001)
data$IAV_log=log(data$IAV+0.001)


########################3.Create lagged variables###########################

data <- data %>%
  arrange( yearweek)
# set maximum lag
nlag = 8

#virus
lag_cov <- tsModel::Lag(data$COV, group = data$subregion,k = 0:nlag)
lag_iav <- tsModel::Lag(data$IAV,  group = data$subregion,k = 0:nlag)

# Remove lag missing value
lag_cov <- lag_cov[complete.cases(lag_cov),]
lag_iav <- lag_iav[complete.cases(lag_iav),]


data <- data %>%
  filter(row_number() > 8)
#########################4.define cross-basis matrix###########################
# set lag knots
lagknot = equalknots(0:nlag, 2)

# COVID-19
var <- lag_cov
basis_cov <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(data$COV, 2)),
                        arglag = list(fun = "ns", knots = lagknot))

# Influenza
var <- lag_iav
basis_iav<- crossbasis(var,
                       argvar = list(fun = "ns", knots = equalknots(data$IAV, 2)),
                       arglag = list(fun = "ns", knots = lagknot))



# assign unique column names to cross-basis matrix for inla() model
colnames(basis_cov) = paste0("basis_cov.", colnames(basis_cov))
colnames(basis_iav) = paste0("basis_iav.", colnames(basis_iav))


########################5.create dataframe for model testing###########################

# set data for models
Y_COV  <- data$COV_log # response variable
T1 <- as.numeric(substring(data$yearweek,7,8)) # for random effect to account for annual cycle (seasonality)
T2 <-as.numeric(substring(data$yearweek,1,4))-2022 # for random effect to account for inter-annual variability
week=data$yearweek

rh<-data$RH
npi<-data$NPI

BA.2.86<-data$BA.2.86
BA.2.75<-data$BA.2.75
omicron<-data$Omicron

med365_cov<-data$med365_COV


df <- data.frame(Y_COV,T1,T2,week,
                 rh,npi,
                 BA.2.75,BA.2.86,omicron,
                 med365_cov)

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


precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

######################## 6. cross validation###########################


#best model
baseformula <- Y_COV ~ 1 + 
  f(T1,  model = "rw1", cyclic = TRUE, constr = TRUE,
    scale.model = TRUE,  hyper = precision.prior) 
formula0.11<- update.formula(baseformula, ~. + rh +
                               npi+omicron+BA.2.75+BA.2.86+med365_cov)
formula1.1 <- update.formula(formula0.11, ~. + basis_iav)

df$Y_COV_1=df$Y_COV



for(i in 1:nrow(df)){

  # Step 1: rerun the selected model (fitted with config = TRUE for sampling) 
  # Step 2: produce cross-validated posterior predictive samples leaving out one year and one week at a time
  
  # define number of samples
  s <- 1000
  # replace data in testing period with NA for out of sample prediction
  casestopred <- df$Y_COV_1 # response variable
  idx.pred <- i
  casestopred[idx.pred] <- NA # replace cases in year and month of interest to NA
  mpred <- length(idx.pred)
  
  # set response variable and year indicator
  df$Y_COV <- casestopred
  
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
  
  
  preds <- list(year = 2022 + df[i,"T2"], week = df[i,"T1"], idx.pred = idx.pred, 
                mean = apply(y.pred, 1, mean), median = apply(y.pred, 1, median),
                lci = apply(y.pred, 1, quantile, probs = c(0.025)),
                uci = apply(y.pred, 1, quantile, probs = c(0.975)))
  save(preds, file = paste0("../output_new/preds_cov_iav_BJ/preds_",2022 + df[i,"T2"], "_",  df[i,"T1"], ".RData"))
  
}



df$pred.mean <- NA
df$pred.median <- NA
df$pred.lci <- NA
df$pred.uci <- NA


for(i in 1:nrow(df)){    
  
  load(paste0("../output_new/preds_cov_iav_BJ/preds_",2022 +  df[i,"T2"],"_", df[i,"T1"],".RData"))
  
  df$pred.mean[preds$idx.pred] <- preds$mean
  df$pred.median[preds$idx.pred] <- preds$median
  df$pred.lci[preds$idx.pred] <- preds$lci
  df$pred.uci[preds$idx.pred] <- preds$uci
  
}




# plot observed v fitted 
p<-ggplot(df) + 
  geom_ribbon(aes(x = week, ymin = pred.lci, ymax = pred.uci), 
              fill = "#D37295", alpha = 0.5) + 
  geom_line(aes(x =week, y = pred.mean, col = "#D37295")) +
  geom_line(aes(x = week, y = Y_COV_1, col = "#499894")) +
  scale_colour_identity(name = "",
                        breaks = c("#499894", "#D37295"),
                        labels = c("Observed", "Posterior predicted"),
                        guide = "legend") +
  theme_bw() + 
  ylab("SARS-CoV-2 positivity rate\n (standardized, log)")+
  theme(axis.title.y = element_text( size=9))+
  xlab(NULL)+
  theme(axis.text.x = element_text(size = 10, color="black",vjust = 0.5, angle = 90))+
  theme(legend.position = "bottom") 



ggsave("../figure/cov_iav_BJ.tiff", p, width=10, height=8, units="cm", dpi=300, compression = "lzw")



