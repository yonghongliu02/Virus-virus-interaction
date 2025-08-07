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
source("2-summary-data_new.R")
data=data[!grepl("Region", data$subregion),]
allcountry=data

region=unique(allcountry$subregion)

table2d.list.Denmark<-NULL
coef.list.Denmark<-NULL
vcov.list.Denmark<-NULL

table2d.list.Ireland<-NULL
coef.list.Ireland<-NULL
vcov.list.Ireland<-NULL

table2d.list.Portugal<-NULL
coef.list.Portugal<-NULL
vcov.list.Portugal<-NULL

table2d.list.Slovenia<-NULL
coef.list.Slovenia<-NULL
vcov.list.Slovenia<-NULL

table2d.list.England<-NULL
coef.list.England<-NULL
vcov.list.England<-NULL



for (j in 1:5) {
  

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


##season
lag_cov <- lag_cov[(data$season2!=0),]
lag_iav <- lag_iav[(data$season2!=0),]
lag_ibv <- lag_ibv[(data$season2!=0),]
lag_iabv <- lag_iabv[(data$season2!=0),]
lag_rsv <- lag_rsv[(data$season2!=0),]

data=data[(data$season2!=0),]
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

BA.2.86<-data$BA.2.86
BA.2.75<-data$BA.2.75
omicron<-data$omicron

med365_cov<-data$med365_cov
med365_Aflu<-data$med365_Aflu
med365_Bflu<-data$med365_Bflu
med365_ABflu<-data$med365_ABflu
med365_rsv<-data$med365_rsv

df <- data.frame(Y_COV,Y_IAV,Y_IBV,Y_IABV,Y_RSV,T1,T2,S,
                 temp,rh,npi,
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




#########################6.best model for COVID-19###########################
# define priors
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

baseformula <- Y_IAV ~ 1 + 
  f(T1,  model = "rw1", cyclic = TRUE, constr = TRUE,
  scale.model = TRUE,  hyper = precision.prior) 

# define formulas by updating the baseline formula with different combinations of climate,virus,flow,NPI cross-basis functions
formula0.1 <- update.formula(baseformula, ~. + temp)
formula0.2 <- update.formula(baseformula, ~. + rh)
formula0.5 <- update.formula(baseformula, ~. + npi)

formula0.9 <- update.formula(baseformula, ~. + med365_Aflu)
formula0.10 <- update.formula(baseformula, ~. + temp +
                                npi+med365_Aflu)
formula0.11<- update.formula(baseformula, ~. + rh +
                               npi+med365_Aflu)



# create a list of formulas
formulas <- list(baseformula, formula0.1, formula0.2, 
                 formula0.5,  formula0.9,
                 formula0.10, formula0.11)
# create model label string
lab <- c("basemodel", "model0.1", "model0.2",
         "model0.5", "model0.9",
         "model0.10", "model0.11")


# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas), 
                 function(i) {
                   model <- mymodel(formulas[[i]], df)
                   save(model, file = paste0("output_new2/select_model_region_iav/iav_",region[j], lab[i],".RData"))})

# create table to store DIC and select best model 
table0 <- data.table(Model  = c("base", "temp", "rh",
                                "npi","med365_Aflu",
                                
                                "tem +
                                npi+med365_Aflu",
                                "rh +
                                npi+med365_Aflu"), 
                     DIC = NA,
                     logscore = NA)


for(i in 1:length(formulas)){
  
  load(paste0("output_new2/select_model_region_iav/iav_",region[j],lab[i],".RData"))
  table0$DIC[i] <- round(model$dic$dic, 0)
  
  table0$logscore[i] <- round(-mean(log(model$cpo$cpo), na.rm = T), 3)
}

# view table
table0

# define position of best fitting model
which.min(table0$DIC)


table0

#########################7. best.flit+virus###########################


formula1.1 <- update.formula(formula0.11, ~. + basis_cov)
formula1.2 <- update.formula(formula0.11, ~. + basis_ibv)
formula1.3 <- update.formula(formula0.11, ~. + basis_iabv)
formula1.4 <- update.formula(formula0.11, ~. + basis_rsv)
# create a list of formulas
formulas.best <- list(formula1.1, formula1.2, formula1.3, 
                      formula1.4)
# create model label string
lab.best <- c( "model1.1", "model1.2","model1.3",
               "model1.4")

# create a function to run a model for each formula in the list and save the model output to file
# WARNING: this may take a long time to run
models <- lapply(1:length(formulas.best), 
                 function(i) {
                   model <- mymodel(formulas.best[[i]], df)
                   save(model, file = paste0("output_new2/select_model_region_iav/iav_",region[j], lab.best[i],".RData"))})

# create table to store DIC and select best model 
table0.best <- data.table(Model  = c("cov","ibv","iabv","rsv" ), 
                          DIC = NA,
                          logscore = NA)

for(i in 1:length(formulas.best)){
  
  load(paste0("output_new2/select_model_region_iav/iav_",region[j],lab.best[i],".RData"))
  table0.best$DIC[i] <- round(model$dic$dic, 0)
  
  table0.best$logscore[i] <- round(-mean(log(model$cpo$cpo), na.rm = T), 3)
}

# view table
table0.best


table0.best



######################## 9. plot model output###########################
   
table0.best.out=table0.best[1:4,]
table0.best.out$load=c(paste0("output_new2/select_model_region_iav/iav_",region[j],"model1.1",".RData"),
                       paste0("output_new2/select_model_region_iav/iav_",region[j],"model1.2",".RData"),
                       paste0("output_new2/select_model_region_iav/iav_",region[j],"model1.3",".RData"),
                       paste0("output_new2/select_model_region_iav/iav_",region[j],"model1.4",".RData"))
table0.best.out$variable=c("basis_cov","basis_ibv","basis_iabv","basis_rsv")
table0.best.out$cen=c(round(mean(data$MA3_cov),2),
                      round(mean(data$MA3_Bflu),2),
                      round(mean(data$MA3_ABflu),2),
                      round(mean(data$MA3_rsv),2)
)
table0.best.out$print=c(1,1,1,1)


#dev.off()

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
                   #cen =  table0.best.out[i,"cen"])

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


all.df=as.data.frame(as.table(z))
all.df$lag=as.numeric(gsub("lag","",all.df$Var1))
all.df$Var2=as.numeric(as.character(all.df$Var2))

all.plot<-ggplot(all.df, aes(x=lag,y=(Var2),  z=Freq)) +
  geom_contour_filled(breaks =levels) + 
  theme_bw()+
  theme(plot.title = element_text(size = 10,hjust = 0.5)) +
  scale_fill_manual(values = cols,drop=FALSE) +
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0))+
  ylab( table0.best.out[i,"Model"])+
  labs(title = paste0("IAV-",region[j]) )+
  guides(fill = guide_colorsteps(direction = "vertical",
                                 barwidth = unit(0.7, "cm"),
                                 barheight = unit(5, "cm")))



#plot2：累计效应
vars <- predt$predvar
allrr <- as.vector(predt$allfit)
allrr.lci <- as.vector(predt$alllow)
allrr.uci <- as.vector(predt$allhigh)
cumulative.df=as.data.frame(cbind(vars,allrr,allrr.lci,allrr.uci))

cumulative.df.select <- cumulative.df[(cumulative.df$allrr < 0 & cumulative.df$allrr.lci < 0 & cumulative.df$allrr.uci < 0) |
                                        (cumulative.df$allrr > 0 &cumulative.df$allrr.lci > 0 & cumulative.df$allrr.uci > 0), ]


cumu.plot<-ggplot( cumulative.df) + 
  geom_ribbon(aes(x =(vars), ymin = (allrr.lci), ymax = (allrr.uci)), 
              fill = "#D37295", alpha = 0.5) + 
  geom_line(aes(x =(vars), y = (allrr), col = "#D37295"))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  ylab("cumulative risk")+
  xlab(table0.best.out[i,"Model"])+
  labs(title=paste0(table0.best.out[i,"Model"],"-",region[j]))+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="none")




# 三次条件判断输出各地区的结果
if (j==1) {
  table2d.list.Denmark[[table0.best.out$Model[i]]]<-cumulative.df
  coef.list.Denmark[[table0.best.out$Model[i]]]<-as.data.frame(coef[indt])
  vcov.list.Denmark[[table0.best.out$Model[i]]]<-as.data.frame(vcov[indt,indt])
} else if (j == 2) {
  table2d.list.Ireland[[table0.best.out$Model[i]]]<-cumulative.df
  coef.list.Ireland[[table0.best.out$Model[i]]]<-as.data.frame(coef[indt])
  vcov.list.Ireland[[table0.best.out$Model[i]]]<-as.data.frame(vcov[indt,indt])
} else if (j == 3) {
  table2d.list.Portugal[[table0.best.out$Model[i]]]<-cumulative.df
  coef.list.Portugal[[table0.best.out$Model[i]]]<-as.data.frame(coef[indt])
  vcov.list.Portugal[[table0.best.out$Model[i]]]<-as.data.frame(vcov[indt,indt])
} else if (j == 4) {
  table2d.list.Slovenia[[table0.best.out$Model[i]]]<-cumulative.df
  coef.list.Slovenia[[table0.best.out$Model[i]]]<-as.data.frame(coef[indt])
  vcov.list.Slovenia[[table0.best.out$Model[i]]]<-as.data.frame(vcov[indt,indt])
} else if (j == 5) {
  table2d.list.England[[table0.best.out$Model[i]]]<-cumulative.df
  coef.list.England[[table0.best.out$Model[i]]]<-as.data.frame(coef[indt])
  vcov.list.England[[table0.best.out$Model[i]]]<-as.data.frame(vcov[indt,indt])
} 

}




}




###############10. save result ###########################


########export
library(writexl)
write_xlsx(table2d.list.Denmark, "output_new2/select_model_region_iav/table2d.list.Denmark.xlsx")
write_xlsx(coef.list.Denmark, "output_new2/select_model_region_iav/coef.list.Denmark.xlsx")
write_xlsx(vcov.list.Denmark, "output_new2/select_model_region_iav/vcov.list.Denmark.xlsx")

write_xlsx(table2d.list.Ireland, "output_new2/select_model_region_iav/table2d.list.Ireland.xlsx")
write_xlsx(coef.list.Ireland, "output_new2/select_model_region_iav/coef.list.Ireland.xlsx")
write_xlsx(vcov.list.Ireland, "output_new2/select_model_region_iav/vcov.list.Ireland.xlsx")

write_xlsx(table2d.list.Portugal, "output_new2/select_model_region_iav/table2d.list.Portugal.xlsx")
write_xlsx(coef.list.Portugal, "output_new2/select_model_region_iav/coef.list.Portugal.xlsx")
write_xlsx(vcov.list.Portugal, "output_new2/select_model_region_iav/vcov.list.Portugal.xlsx")

write_xlsx(table2d.list.Slovenia, "output_new2/select_model_region_iav/table2d.list.Slovenia.xlsx")
write_xlsx(coef.list.Slovenia, "output_new2/select_model_region_iav/coef.list.Slovenia.xlsx")
write_xlsx(vcov.list.Slovenia, "output_new2/select_model_region_iav/vcov.list.Slovenia.xlsx")

write_xlsx(table2d.list.England, "output_new2/select_model_region_iav/table2d.list.England.xlsx")
write_xlsx(coef.list.England, "output_new2/select_model_region_iav/coef.list.England.xlsx")
write_xlsx(vcov.list.England, "output_new2/select_model_region_iav/vcov.list.England.xlsx")



