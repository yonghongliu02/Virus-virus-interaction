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
library(metR)
my_rbind <- function(df1, df2){    # Create own merging function
  rbind(df1, df2)
}




source("1-summary-data.R")

####1. meta variable########
meta<-read_excel("../data/meta.xlsx")
meta=meta[!grepl("Region", meta$country), ]
meta=meta[!is.na(meta$`order than 65 year`),]
colnames(meta)=c("country","year","GDP","UHC","pop_density",
                 "under5","above65","SDI","popu","area")


meta.new=meta[,c("country","GDP","UHC","pop_density",
                 "above65","SDI")]

meta.new$GDP=meta.new$GDP/10^5
chart.Correlation(meta.new[,2:6], histogram = TRUE, pch = 19,
                  method = "spearman")


meta.new$country=gsub( "United States","US",meta.new$country)
meta.new$country=gsub( "Hong Kong SAR, China","HK",meta.new$country)

meta.new <- as.data.frame(meta.new)
rownames(meta.new)=meta.new$country
meta.new=meta.new[-3,-1]


####2. basic lag variable########

#data=df_us
data=data
# 美国的平均
us_data <- data %>% filter(country == "US")

us_data1 <- us_data%>%
  group_by(yearweek)%>%
  dplyr::summarise(MA3_cov = mean(MA3_cov),
                   MA3_Aflu = mean(MA3_Aflu),
                   MA3_rsv= mean(MA3_rsv))
# 非美国的平均
other_data <- data %>% filter(country != "US")

other_data1 <- other_data%>%
  group_by(yearweek)%>%
  dplyr::summarise(MA3_cov = mean(MA3_cov),
                   MA3_Aflu = mean(MA3_Aflu),
                   MA3_rsv= mean(MA3_rsv))

data=rbind(us_data1,other_data1)
data <- data%>%
  group_by(yearweek)%>%
  dplyr::summarise(MA3_cov = mean(MA3_cov),
                   MA3_Aflu = mean(MA3_Aflu),
                   MA3_rsv= mean(MA3_rsv))

data$subregion="all"


#标准化
data <- data %>%
  mutate(
    MA3_cov = (MA3_cov - min(MA3_cov)) / (max(MA3_cov)-min(MA3_cov))*100,
    MA3_Aflu= (MA3_Aflu - min(MA3_Aflu)) / (max(MA3_Aflu)-min(MA3_Aflu))*100,
    MA3_rsv = (MA3_rsv - min(MA3_rsv)) / (max(MA3_rsv)-min(MA3_rsv))*100)

data <- data %>%
  arrange(subregion, yearweek)


##season2
data$week=as.numeric(substr(data$yearweek,7,8))
data$season2=ifelse(data$week>=40|data$week<=13, 40, 0)

# set maximum lag
nlag = 8

#virus
lag_cov <- tsModel::Lag(data$MA3_cov, group = data$subregion,k = 0:nlag)
lag_iav <- tsModel::Lag(data$MA3_Aflu,  group = data$subregion,k = 0:nlag)
lag_rsv <- tsModel::Lag(data$MA3_rsv, group = data$subregion, k = 0:nlag)

# Remove lag missing value
lag_cov <- lag_cov[complete.cases(lag_cov),]
lag_iav <- lag_iav[complete.cases(lag_iav),]
lag_rsv <- lag_rsv[complete.cases(lag_rsv),]



data <- data %>%
  group_by(subregion) %>%
  filter(row_number() > 8)



##season
lag_cov <- lag_cov[(data$season2!=0),]
lag_iav <- lag_iav[(data$season2!=0),]
lag_rsv <- lag_rsv[(data$season2!=0),]

data=data[(data$season2!=0),]



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

# Rsv
var <- lag_rsv
basis_rsv <- crossbasis(var,
                        argvar = list(fun = "ns", knots = equalknots(data$MA3_rsv, 2)),
                        arglag = list(fun = "ns", knots = lagknot))


# assign unique column names to cross-basis matrix for inla() model
colnames(basis_cov) = paste0("basis_cov.", colnames(basis_cov))
colnames(basis_iav) = paste0("basis_iav.", colnames(basis_iav))
colnames(basis_rsv) = paste0("basis_rsv.", colnames(basis_rsv))



####3. model coef vcov########
library(readxl)
y=c("cov","iav","rsv")

vcov.list<-NULL
cumulative.df.list<-NULL
mvall.result.list<-NULL

allplot.list<-NULL
cumuplot.list<-NULL
bestplotmin.list<-NULL
bestplotmax.list<-NULL

for (j in 1:3) {
  
  x=y[y!=y[j]]
  
  for (i in 1:2) {
    
    #1.data
    df.us<-read_excel(paste0("../output_new2/select_model_us_",y[j],"/coef.list.xlsx"),sheet=x[i])
    df.slovenia<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.Slovenia.xlsx"),sheet=x[i])
   # df.hk<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.HK.xlsx"),sheet=x[i])
    df.denmark<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.Denmark.xlsx"),sheet=x[i])
    df.ireland<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.Ireland.xlsx"),sheet=x[i])
    df.portugal<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.Portugal.xlsx"),sheet=x[i])
    df.england<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.England.xlsx"),sheet=x[i])
    # df.usall<-read_excel(paste0("../output_new2/select_model_region_",y[j],"/coef.list.Usall.xlsx"),sheet=x[i])
    
    
    df.coef=t(cbind(df.us,df.slovenia,
                    #df.hk,
                    df.denmark,
                    df.ireland,
                    df.portugal,
                    df.england))
    
    
    
    rownames(df.coef)=c("US","Slovenia",
                       # "HK",
                        "Denmark","Ireland","Portugal","England")
    
    
    vcov.list[[1]]=read_excel(paste0("../output_new2/select_model_us_",y[j],"/vcov.list.xlsx"),sheet=x[i])
    vcov.list[[2]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.Slovenia.xlsx"),sheet=x[i])
   # vcov.list[[3]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.HK.xlsx"),sheet=x[i])
    vcov.list[[3]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.Denmark.xlsx"),sheet=x[i])
    vcov.list[[4]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.Ireland.xlsx"),sheet=x[i])
    vcov.list[[5]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.Portugal.xlsx"),sheet=x[i])
    vcov.list[[6]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.England.xlsx"),sheet=x[i])
    #vcov.list[[7]]=read_excel(paste0("../output_new2/select_model_region_",y[j],"/vcov.list.Usall.xlsx"),sheet=x[i])
    
    
    
    
    
    #2.meta
    library(mvmeta)
    method <- "reml" #限制性极大似然估计
    
    #不加协变量
    mvall <- mvmeta(df.coef~1,vcov.list,
                    # control=list(showiter=F,igls.iter=10000)，
                    method=method)
    #加协变量
    mvall.all <- mvmeta(df.coef~1+GDP+UHC+pop_density+above65+SDI ,vcov.list,
                        data=meta.new,
                        # control=list(showiter=F,igls.iter=10000)，
                        method=method)
    
    
    #记录异质性检验结果
    mvall.list<-NULL
    mvall.list <- lapply(meta.new, function(x) mvmeta(df.coef~x,vcov.list,method= method))
    mvall.list[["0-base"]]<-mvall
    mvall.list[["0-all"]]<-mvall.all
    mvall.result<-as.data.frame(cbind(variable=NA,Q=NA,df=NA,p=NA,I2=NA))
    
    
    for (z in 1: length(mvall.list)) {
      
      model.qt=qtest(mvall.list[[z]])
      
      mvall.result[z,"variable"]=names(mvall.list)[z]
      mvall.result[z,"Q"]=model.qt$Q[1]
      mvall.result[z,"df"]=model.qt$df[1]
      mvall.result[z,"p"]=model.qt$pvalue[1]
      mvall.result[z,"I2"]=(model.qt$Q[1]-model.qt$df[1])/model.qt$Q[1]*100
      
    }
    
    mvall.result.list[[paste0(y[j],"-",x[i])]]=mvall.result
    
    
    
    ##预测
    predlat <- predict(mvall, 
                       data.frame(GDP=apply(meta.new, 2, median)[1],
                                  UHC=apply(meta.new, 2, median)[2],
                                  pop_density=apply(meta.new, 2, median)[3],
                                  above65=apply(meta.new, 2, median)[4],
                                  SDI=apply(meta.new, 2, median)[5]),
                       vcov=T)
    
    
    
    
    predt <- crosspred(
      eval(parse(text =  paste0("basis_",x[i]))),
      coef=predlat$fit,
      vcov=predlat$vcov,
      bylag = 0.25, 
      model.link = "log",
      cen=0)
    
    
    ##plot1:3d
    zz <- (t(predt$matRRfit))
    library(RColorBrewer)
    pal <-rev(brewer.pal(11, "RdBu"))
    levels <-seq(0,2,0.1)
   # levels<-pretty(zz,20)
    col1 <- colorRampPalette(pal[1:6])
    col2 <- colorRampPalette(pal[6:11])
    cols <-c(col1(sum(levels <1)), col2(sum(levels >1)))
    
    all.df=as.data.frame(as.table(zz))
    all.df$lag=as.numeric(gsub("lag","",all.df$Var1))
    all.df$Var2=as.numeric(as.character(all.df$Var2))
    
    all.df$x=x[i]
    all.df$y=y[j]
    
    #字符替换
    all.df$x <- gsub("cov", "SARS-CoV-2", all.df$x)
    all.df$x <- gsub("iav", "IAV", all.df$x)
    all.df$x <- gsub("rsv", "RSV", all.df$x)
    all.df$y <- gsub("cov", "SARS-CoV-2", all.df$y)
    all.df$y <- gsub("iav", "IAV", all.df$y)
    all.df$y <- gsub("rsv", "RSV", all.df$y)
    
    all.plot<-ggplot(all.df, aes(x=lag,y=(Var2),  z=Freq)) +
      geom_contour_fill() + 
      theme_bw()+
      scale_fill_gradientn(limits = c(0.2, 1.8), colours = cols)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+

      ylab(paste0("Percentile of ",unique(all.df$x)," (%)"))+
      
      theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5))+
      theme(axis.text.y = element_text(size = 9, color="black",vjust = 0.5))+
     labs(title = paste0("RR of ",unique(all.df$y))) + 
      labs(fill =NULL) +
      theme(plot.title = element_text(size = 10))  +# 调整主标题大小
      theme(plot.title = element_text(hjust = 0.5)) +  # 设置图例标签名称
      theme(legend.title = element_text(size =8),
            legend.key.width = unit(0.3, "cm"),  # Adjusted width
            legend.key.height = unit(0.6, "cm") ) # Adjusted height)
  
    
    all.plot
    
    
    allplot.list[[paste0(y[j],"-",x[i])]]<- all.plot
    
    #plot2：累计效应
    vars <- predt$predvar
    allrr <- as.vector(predt$allRRfit)
    allrr.lci <- as.vector(predt$allRRlow)
    allrr.uci <- as.vector(predt$allRRhigh)
    cumulative.df=as.data.frame(cbind(vars,allrr,allrr.lci,allrr.uci))
  
    
    cumulative.df$x=x[i]
    cumulative.df$y=y[j]
    
    #字符替换
    cumulative.df$x <- gsub("cov", "SARS-CoV-2", cumulative.df$x)
    cumulative.df$x <- gsub("iav", "IAV", cumulative.df$x)
    cumulative.df$x <- gsub("rsv", "RSV", cumulative.df$x)
    cumulative.df$y <- gsub("cov", "SARS-CoV-2", cumulative.df$y)
    cumulative.df$y <- gsub("iav", "IAV", cumulative.df$y)
    cumulative.df$y <- gsub("rsv", "RSV", cumulative.df$y)
    
    cumu.plot<-ggplot( cumulative.df) + 
      geom_ribbon(aes(x =(vars), ymin = (allrr.lci), ymax = (allrr.uci)), 
                  fill = "#66c2a5", alpha = 0.2) + 
      geom_line(aes(x =(vars), y =(allrr)), color = "#66c2a5",size=0.9)+
      geom_hline(yintercept = 1, linetype = "dashed")+
      ylab(paste0( "RR of ",unique(all.df$y)))+
      xlab(paste0("Percentile of ",unique(all.df$x), " (%)"))+
       theme_classic()+
       theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5))+
       theme(axis.text.y = element_text(size = 9, color="black",vjust = 0.5))+
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(legend.position="none")
    
   cumuplot.list[[paste0(y[j],"-",x[i])]]<-  cumu.plot
    
  
   
   
   
   
   
   
   
  
  #plot3：specify x min##
  rr <- (predt$matRRfit)
  rr.lci <- (predt$matRRlow)
  rr.uci <- (predt$matRRhigh)
  
  mx<-which(round(vars, 2) ==  round(cumulative.df$vars[which(cumulative.df$allrr==min(cumulative.df$allrr),)],2))
 print( cumulative.df[which.min(cumulative.df$allrr),])
  specify.df=as.data.frame(t(rbind(
    lagbylag =seq(0, nlag, 0.25),
    low=rr.lci[ mx,],
    high=rr.uci[ mx,],
    best=rr[ mx,])))
  
  

  specify.df$y=y[j]
  
  #字符替换
  specify.df$y <- gsub("cov", "SARS-CoV-2", specify.df$y)
  specify.df$y <- gsub("iav", "IAV", specify.df$y)
  specify.df$y <- gsub("rsv", "RSV", specify.df$y)
  
  best.plot<- ggplot( specify.df) + 
    geom_ribbon(aes(x = lagbylag, ymin = (low), ymax = (high)), 
                fill = "#fc8d62", alpha = 0.2) + 
    geom_line(aes(x =lagbylag, y = best,  colour = "#fc8d62"),size=0.9)+
    geom_hline(yintercept = 1, linetype = "dashed")+
    ylab(paste0("RR of ",unique(all.df$y)))+
    xlab("Lag")+
    theme_classic()+
    theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5))+
    theme(axis.text.y = element_text(size = 9, color="black",vjust = 0.5))+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position = c(0.3,0.9)) +
    scale_colour_discrete(labels=paste0("Percentile = ",vars[mx],"%"))+
    theme(legend.background = element_blank(),legend.key = element_blank())+ 
    theme(legend.title=element_blank())
  
   bestplotmin.list[[paste0(y[j],"-",x[i])]]<-  best.plot
   
   
   #plot3：specify x max##
   rr <- (predt$matRRfit)
   rr.lci <- (predt$matRRlow)
   rr.uci <- (predt$matRRhigh)
   
   mx<-which(round(vars, 2) ==  round(cumulative.df$vars[which(cumulative.df$allrr==max(cumulative.df$allrr),)],2))
   print(cumulative.df[which.max(cumulative.df$allrr),])
   specify.df=as.data.frame(t(rbind(
     lagbylag =seq(0, nlag, 0.25),
     low=rr.lci[ mx,],
     high=rr.uci[ mx,],
     best=rr[ mx,])))
   
   best.plot<- ggplot( specify.df) + 
     geom_ribbon(aes(x = lagbylag, ymin = (low), ymax =(high)), 
                 fill = "#fc8d62", alpha = 0.2) + 
     geom_line(aes(x =lagbylag, y = best,colour = "#fc8d62"),size=0.9)+
     geom_hline(yintercept = 1, linetype = "dashed")+
     ylab(paste0("RR of ",unique(all.df$y)))+
     xlab("Lag")+
     theme_classic()+
     theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5))+
     theme(axis.text.y = element_text(size = 9, color="black",vjust = 0.5))+
     theme(plot.title = element_text(hjust = 0.5)) + 
     theme(legend.position = c(0.3,0.9)) +
     scale_colour_discrete(labels=paste0("Percentile = ",vars[mx],"%"))+
     theme(legend.background = element_blank(),legend.key = element_blank())+ 
        theme(legend.title=element_blank())
   
   bestplotmax.list[[paste0(y[j],"-",x[i])]]<-  best.plot
  
  }
}




library(ggpubr)



arranged_plot<-ggarrange(
  allplot.list[[2]], cumuplot.list[[2]],
  allplot.list[[3]], cumuplot.list[[3]],
  allplot.list[[5]], cumuplot.list[[5]],
  
  
  ncol = 2,  nrow = 3, align = 'hv',hjust=-2,
  labels = c("A", "B", "C",
             "D", "E", "F"))

arranged_plot
ggsave("../figure/figS5.tiff", plot = arranged_plot, width=16, height=18, units="cm", dpi=300, compression = "lzw")


