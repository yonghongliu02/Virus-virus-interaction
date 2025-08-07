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
########################1. load map###########################
library(sf)
map <- read_sf("gadm41_USA_shp/gadm41_USA_1.shp") 
region<-read_excel("../data/US/region_us.xlsx")
colnames(region)=c("region","NAME_1")
#add region to map
map.us=merge(map,region,id="NAME_1")


map.us.region <- map.us %>%
  group_by(region) %>%
  summarise(geometry = st_union(geometry))


#plot(map.us.region)

#library(spdep)
#nb <- poly2nb(map.us.region)
#head(nb)
#nb2INLA("map.adj", nb)

#####################0. load data###########################################
g <- inla.read.graph(filename = "map.adj")
source("2-summary-data_new.R")
data=data[grepl("Region", data$subregion), ]

##add popu
popu=read_excel("meta.xlsx")
popu=popu[grepl("Region", popu$country), c("country","GDP(US$ millions)","popdensity(2021)", "order than 65 year")]
colnames(popu)=c("subregion","GDP","popu_density","above65")
popu$GDP=popu$GDP/10^7

data=merge(data,popu,by="subregion",all.x=T)
data$region=as.integer(factor(data$subregion))

##log

data$MA3_cov_log=log(data$MA3_cov+0.001)
data$MA3_Aflu_log=log(data$MA3_Aflu+0.001)
data$MA3_Bflu_log=log(data$MA3_Bflu+0.001)
data$MA3_ABflu_log=log(data$MA3_ABflu+0.001)
data$MA3_rsv_log=log(data$MA3_rsv+0.001)

#########################3.data distribution###########################

#method = c("pearson", "kendall", "spearman"),
chart.Correlation(data[,20:35], histogram = TRUE, pch = 19,
                  method = "pearson")
qqplot(qnorm(ppoints(length(data$MA3_cov_log))),data$MA3_cov_log)



########################4.Create lagged variables###########################
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


data <- data %>%
  group_by(subregion) %>%
  filter(row_number() > 8)


##season
lag_cov <- lag_cov[-which(data$season2==0),]
lag_iav <- lag_iav[-which(data$season2==0),]
lag_ibv <- lag_ibv[-which(data$season2==0),]
lag_iabv <- lag_iabv[-which(data$season2==0),]
lag_rsv <- lag_rsv[-which(data$season2==0),]

data=data[-which(data$season2==0),]
#########################5.define cross-basis matrix###########################
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

########################6.create dataframe for model testing###########################

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
popu_density<-data$popu_density
above65<-data$above65
npi<-data$new_npi

BA.2.86<-data$BA.2.86
BA.2.75<-data$BA.2.75
omicron<-data$omicron

med365_cov<-data$med365_cov
med365_Aflu<-data$med365_Aflu
med365_Bflu<-data$med365_Bflu
med365_ABflu<-data$med365_ABflu
med365_rsv<-data$med365_rsv

df <- data.frame(Y_COV,Y_IAV,Y_IBV,Y_RSV,T1,T2,S,
                 temp,rh,popu_density,above65,npi,
                 BA.2.75,BA.2.86,omicron,
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




#########################7.best model for COVID-19###########################
# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

baseformula <- Y_RSV ~ 1 + 
  f(T1, replicate = S, model = "rw1", cyclic = TRUE, constr = TRUE,
  scale.model = TRUE,  hyper = precision.prior) +
  f(S, model = "bym2", replicate = T2, graph = g,
    scale.model = TRUE, hyper = precision.prior) 

# define formulas by updating the baseline formula with different combinations of climate,virus,flow,NPI cross-basis functions
formula0.1 <- update.formula(baseformula, ~. + temp)
formula0.2 <- update.formula(baseformula, ~. + rh)
formula0.3 <- update.formula(baseformula, ~. + popu_density)
formula0.4 <- update.formula(baseformula, ~. + above65)
formula0.5 <- update.formula(baseformula, ~. + npi)

formula0.9 <- update.formula(baseformula, ~. + med365_rsv)
formula0.10 <- update.formula(baseformula, ~. + temp +
                                popu_density + above65+
                                npi+med365_rsv)
formula0.11<- update.formula(baseformula, ~. + rh +
                               popu_density + above65+
                               npi+med365_rsv)



# create a list of formulas
formulas <- list(baseformula, formula0.1, formula0.2, formula0.3, formula0.4, 
                 formula0.5,  formula0.9,
                 formula0.10, formula0.11)
# create model label string
lab <- c("basemodel", "model0.1", "model0.2","model0.3", "model0.4",
         "model0.5", "model0.9",
         "model0.10", "model0.11")


# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas), 
                 function(i) {
                   model <- mymodel(formulas[[i]], df)
                   save(model, file = paste0("output_new2/select_model_us_rsv/rsv_", lab[i],".RData"))})

# create table to store DIC and select best model 
table0 <- data.table(Model  = c("base", "temp", "rh","pop_density","above65",
                                "npi","med365_rsv",
                         
                                "tem +
                                popu_density + above65+
                                npi+med365_rsv",
                                "rh +
                                popu_density + above65+
                                npi+med365_rsv"), 
                     DIC = NA,
                     logscore = NA)

for(i in 1:length(formulas)){
  
  load(paste0("output_new2/select_model_us_rsv/rsv_",lab[i],".RData"))
  table0$DIC[i] <- round(model$dic$dic, 0)
  
  table0$logscore[i] <- round(-mean(log(model$cpo$cpo), na.rm = T), 3)
}

# view table
table0

# define position of best fitting model
which.min(table0$DIC)



#########################8. best.flit+virus###########################



formula1.1 <- update.formula(formula0.11, ~. + basis_cov)
formula1.2 <- update.formula(formula0.11, ~. + basis_ibv)
formula1.3 <- update.formula(formula0.11, ~. + basis_iabv)
formula1.4 <- update.formula(formula0.11, ~. + basis_iav)


# create a list of formulas
formulas.best <- list(formula1.1, formula1.2, formula1.3, formula1.4)
# create model label string
lab.best <- c( "model1.1", "model1.2","model1.3","model1.4")

# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas.best), 
                 function(i) {
                   model <- mymodel(formulas.best[[i]], df)
                   save(model, file = paste0("output_new2/select_model_us_rsv/rsv_", lab.best[i],".RData"))})

# create table to store DIC and select best model 
table0.best <- data.table(Model  = c("cov","ibv","iabv","iav"), 
                          DIC = NA,
                          logscore = NA)

for(i in 1:length(formulas.best)){
  
  load(paste0("output_new2/select_model_us_rsv/rsv_",lab.best[i],".RData"))
  table0.best$DIC[i] <- round(model$dic$dic, 0)
  
  table0.best$logscore[i] <- round(-mean(log(model$cpo$cpo), na.rm = T), 3)
}

# view table
table0.best


######################## 9. compare MAE###########################
#basemodel
load("output_new2/select_model_us_rsv/rsv_basemodel.RData")
basemodel <- model
#best model
load("output_new2/select_model_us_rsv/rsv_model1.1.RData")


# add baseline fitted model result summaries (2.5%, 50%, 97.5% percentiles) to data
data$base.fit <- basemodel$summary.fitted.values$`0.5quant`
data$base.fit.lci <- basemodel$summary.fitted.values$`0.025quant`
data$base.fit.uci <- basemodel$summary.fitted.values$`0.975quant`

# add selected fitted model result summaries (2.5%, 50%, 97.5% percentiles) to data
data$fit <- model$summary.fitted.values$`0.5quant`
data$fit.lci<-model$summary.fitted.values$`0.025quant`
data$fit.uci<-model$summary.fitted.values$`0.975quant`

# compute mean absolute error and compare base model to final model
MAE <- as.data.frame(matrix(NA, nrow = nregion, ncol = 2))
names(MAE) <-c("base", "new")

# calculate the MAE for observed and mean fit cases
for (i in 1:10) {
  
  # cases
  MAE$base[i] <- mae(data$base.fit[df$S == i], 
                               data$MA3_rsv_log[df$S == i], 
                               na.rm = TRUE)
  MAE$new[i] <- mae(data$fit[df$S == i], 
                              data$MA3_rsv_log[df$S == i], 
                              na.rm = TRUE)
  
}

# calculate difference between MAE from the baseline model and MAE from the selected model
MAE$diff <- MAE$base - MAE$new
mn <-min(MAE$diff)
mx <-max(MAE$diff)

MAE$value <- 1
# assign areas where the difference is zero or more as '2', i.e. MAE from new model is smaller (fits better than the base model)
MAE$value[MAE$diff >= 0 ] <- 2
MAE
# plot map to show areas where the new model provided 'added value' over the basemodel (e.g. MAE is smaller for new model) (Appendix Fig S11)
 ggplot(map.us.region) + 
  geom_sf(aes(fill = factor(MAE$value)), lwd = 0) +
  scale_fill_manual(values = c("#17becf","#dc5fbd"), breaks = 1:2, 
                    labels = c("No added value", "Added value")) +
  labs(fill = "") +
  theme_void()



######################## 11. plot model output###########################

table0.best.out=table0.best[1:4,]
table0.best.out$load=c("output_new2/select_model_us_rsv/rsv_model1.1.RData",
                       "output_new2/select_model_us_rsv/rsv_model1.2.RData",
                       "output_new2/select_model_us_rsv/rsv_model1.3.RData",
                       "output_new2/select_model_us_rsv/rsv_model1.4.RData")
table0.best.out$variable=c("basis_cov","basis_ibv","basis_iabv","basis_iav")
table0.best.out$cen=c(round(mean(data$MA3_cov),2),
                      round(mean(data$MA3_Bflu),2),
                      round(mean(data$MA3_ABflu),2),
                      round(mean(data$MA3_Aflu),2)
                      )
table0.best.out$print=c(1,1,1,1)


dev.off()


table2d.list<-NULL
coef.list<-NULL
vcov.list<-NULL



##for plot

for (i in 1:nrow(table0.best.out)) {

load(as.character(table0.best.out[i,"load"]))
  
  
##plot1  
# extract full coef and vcov and create indicators for each term
coef <- model$summary.fixed$mean
vcov <- model$misc$lincomb.derived.covariance.matrix
# find position of the terms associated with influenza crossbasis
indt <- grep(table0.best.out[i,"variable"], model$names.fixed)
# extract predictions from the influenza DLNM centred on overall mean 
predt <- crosspred(
                   eval(parse(text = table0.best.out[i,"variable"])),
                   coef = coef[indt], vcov=vcov[indt,indt],
                   bylag = 0.25, 
                   cen=0)
                  # cen =  table0.best.out[i,"cen"])


coef.list[[table0.best.out$Model[i]]]<-as.data.frame(coef[indt])
vcov.list[[table0.best.out$Model[i]]]<-as.data.frame(vcov[indt,indt])


# contour and scenario plots for each virus
y <- predt$predvar
x <- seq(0, nlag, 0.25)
z <- t(predt$matfit)
library(RColorBrewer)
pal <- rev(brewer.pal(11, "PRGn"))
levels <- pretty(z, 20)
col1 <- colorRampPalette(pal[1:6])
col2 <- colorRampPalette(pal[6:11])
cols <- c(col1(sum(levels < 0)), col2(sum(levels > 0)))

filled.contour(x,(y),z,
               xlab = "Lag", ylab = table0.best.out[i,"Model"], main = "RSV",
               col = cols,levels = levels,
               plot.axes = { axis(1, at = 0:nlag, c(0:nlag)) 
                 axis(2)})




#plot2：累计效应
vars <- predt$predvar
allrr <- as.vector(predt$allfit)
allrr.lci <- as.vector(predt$alllow)
allrr.uci <- as.vector(predt$allhigh)
cumulative.df=as.data.frame(cbind(vars,allrr,allrr.lci,allrr.uci))

cumulative.df.select <- cumulative.df[(cumulative.df$allrr < 0 & cumulative.df$allrr.lci < 0 & cumulative.df$allrr.uci < 0) |
                                        (cumulative.df$allrr > 0 &cumulative.df$allrr.lci > 0 & cumulative.df$allrr.uci > 0), ]


p<-ggplot( cumulative.df) + 
  geom_ribbon(aes(x = vars, ymin = (allrr.lci), ymax = (allrr.uci)), 
              fill = "#D37295", alpha = 0.5) + 
  geom_line(aes(x =vars, y = (allrr), col = "#D37295"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("cumulative risk")

print(p)


table2d.list[[table0.best.out$Model[i]]]<-cumulative.df




}




########export
library(writexl)
write_xlsx(table2d.list, "output_new2/select_model_us_rsv/table2d.list.xlsx")
write_xlsx(coef.list, "output_new2/select_model_us_rsv/coef.list.xlsx")
write_xlsx(vcov.list, "output_new2/select_model_us_rsv/vcov.list.xlsx")

