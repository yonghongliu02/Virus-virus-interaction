####figure1#####
library(data.table)
library(readxl)
library(tsibble)
library(ggplot2)

data<-read.csv("../data/combine-data.csv")
data$yearweek=yearweek(data$yearweek)
unique(data$yearweek)
data=data[which(data$yearweek>yearweek("2021 W39")&
                  data$yearweek<yearweek("2024 W21")),]


##1. us all positive####
df_us=data[data$country=="US",
           c("yearweek","subregion","MA3_cov","MA3_Aflu","MA3_Bflu","MA3_rsv")]



us.melt1=melt(df_us,id=c("yearweek","subregion"))



maxt=max(data$yearweek)
mint=min(data$yearweek)

p_us<-ggplot(us.melt1,aes(x=yearweek,y=value,color=variable,group = interaction(variable,subregion))) +
   geom_line(size=0.2,alpha=0.7)+
 
  scale_color_manual(breaks=c("MA3_cov","MA3_Aflu","MA3_Bflu","MA3_rsv"),
    labels = c('COV2', 'Aflu', 'Bflu', 'RSV'), 
                     values = c('#1f78b4', '#e31a1c', '#eeb479', '#33a02c'))+
  theme_bw()+
  xlab(NULL) +  
  ylab("Positive rate (%)") + 
  scale_y_continuous(limits=c(0,80),expand = c(0, 0))+
  scale_x_yearweek(date_breaks="20 weeks",limits=as.Date(c(mint,maxt)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color="black",hjust = 1, angle =30),
        axis.text.y = element_text(size = 9, color="black",vjust = 0.1),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text( size=10))+
  theme(plot.margin=unit(c(0.5,1,0.5,2.5),'lines'))+
  theme(legend.position = c(0.4,0.8))+
  theme(legend.background = element_blank(),legend.key = element_blank())+
  theme(legend.title = element_blank())+
  labs(title="US")+
  theme(legend.text = element_text(size = 8))+
  guides(color = guide_legend(ncol = 4))+
  theme(plot.title = element_text(hjust = 0.5)) 
  
p_us



##2. US HHS####
geo<-read_excel("../data/region_us.xlsx",sheet="map")
colnames(us.melt1)=c( "yearweek","name" ,"variable", "value")
us.melt1.merge=merge(us.melt1,geo,by="name")
us.melt1$name=factor(us.melt1$name,levels = c("Region 1", "Region 2" , "Region 3" ,
                                              "Region 4" , "Region 5",  "Region 6"  ,"Region 7" ,
                                              "Region 8" ,"Region 9",  "Region 10" ))

library(geofacet)
p1<-ggplot(us.melt1,aes(x=yearweek,y=value,color=variable,group = interaction(variable,name))) +
  geom_line(size=0.4)+
  scale_color_manual(breaks=c("MA3_cov","MA3_Aflu","MA3_Bflu","MA3_rsv"),
                     labels = c('Cov-2', 'IAV', 'IBV', 'RSV'), 
                     values = c('#1f78b4', '#e31a1c', '#eeb479', '#33a02c'))+
  theme_minimal() +
  scale_y_continuous(limits=c(0,65),expand = c(0, 0))+
  scale_x_yearweek(date_breaks="20 weeks",limits=as.Date(c(mint,maxt)),expand = c(0, 0))+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color="black",hjust = 1, angle =60),
        axis.text.y = element_text(size = 9, color="black",vjust = 0.1),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text( size=10))+
  # organise by state name in grid file
  facet_wrap(~name,ncol = 2)+
  theme(legend.background = element_blank(),legend.key = element_blank(),
       strip.background =element_rect(color="white",fill="white") ,
       strip.text.x = element_text(size = 12),
       panel.grid=element_blank())+
  theme(legend.title = element_blank())+  
  xlab(NULL) +  
  ylab("Positivity rate (%)") +  
  theme(legend.text = element_text(size = 7))+
  theme(legend.key.height = unit(0.8, "lines"))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(legend.position="none")

p1



##3. other####

p_list=NULL
df_nonus=data[!grepl("Region", data$subregion), 
              c("yearweek","subregion","MA3_cov","MA3_Aflu","MA3_Bflu","MA3_rsv")]



df.melt=melt(df_nonus,id=c("yearweek","subregion"))

df.melt$subregion=factor(df.melt$subregion,levels = c("Denmark" ,"Ireland","Portugal" ,"Slovenia", "England"),
                         labels = c("Denmark" ,"Ireland","Portugal" ,"Slovenia", "England" ))

p2<-ggplot(df.melt,aes(x=yearweek,y=value,color=variable)) +
  geom_line(size=0.6)+
  theme_minimal() +
  scale_color_manual(breaks=c("MA3_cov","MA3_Aflu","MA3_Bflu","MA3_rsv"),
                     labels = c('SARS-CoV-2', 'IAV', 'IBV', 'RSV'), 
                     values = c('#1f78b4', '#e31a1c', '#eeb479', '#33a02c'))+
  xlab(NULL) +  
  ylab("Positivity rate (%)") + 
  scale_y_continuous(limits=c(0,65),expand = c(0, 0))+
  scale_x_yearweek(date_breaks="20 weeks",limits=as.Date(c(mint,maxt)),expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, color="black",hjust = 1, angle = 60),
        axis.text.y = element_text(size = 9, color="black",vjust = 0.1),
        axis.title.x=element_text(size=10),
        axis.title.y=element_text( size=10))+
  theme(plot.margin=unit(c(0.5,0.5,0.5,1.5),'lines'))+
  theme(legend.background = element_blank(),legend.key = element_blank(),
        strip.background =element_rect(color="white",fill="white") ,
        strip.text.x = element_text(size = 12),
        panel.grid=element_blank())+
  theme(legend.title = element_blank())+  
  theme(legend.position = c(0.4,0.98))+
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~subregion,ncol = 1)+
 
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  guides(color = guide_legend(nrow = 1))
p2

library(ggpubr)
arranged_plot<-ggarrange(  p1,p2,   ncol=2 , align = 'v',
            labels = c("A","B"))
arranged_plot
ggsave("../figure/fig1.tiff", plot = arranged_plot, width=24, height=18, 
       units="cm", dpi=300, compression = "lzw")


