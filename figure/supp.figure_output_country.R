####load data####
library(readxl)
library(ggplot2)
library(gridExtra) 
library(RColorBrewer)
y=c("cov","iav","rsv")
df.list<-NULL

for (j in 1:3) {
  
  x=y[y!=y[j]]

for (i in 1:2) {
  


df.us<-read_excel(paste0("../output_new/select_model_us_",y[j],"/table2d.list.xlsx"),sheet=x[i])

df.denmark<-read_excel(paste0("../output_new/select_model_region_",y[j],"/table2d.list.Denmark.xlsx"),sheet=x[i])
df.ireland<-read_excel(paste0("../output_new/select_model_region_",y[j],"/table2d.list.Ireland.xlsx"),sheet=x[i])
df.portugal<-read_excel(paste0("../output_new/select_model_region_",y[j],"/table2d.list.Portugal.xlsx"),sheet=x[i])
df.slovenia<-read_excel(paste0("../output_new/select_model_region_",y[j],"/table2d.list.Slovenia.xlsx"),sheet=x[i])
df.england<-read_excel(paste0("../output_new/select_model_region_",y[j],"/table2d.list.England.xlsx"),sheet=x[i])


df.us$region="US"
df.denmark$region="Denmark"
df.ireland$region="Ireland"
df.portugal$region="Portugal"
df.slovenia$region="Slovenia"
df.england$region="England"



df=rbind(df.us,
           df.denmark,
         df.ireland,
         df.portugal,
         df.slovenia,
         df.england)

df$x=x[i]
df$y=y[j]

df.list[[paste0(y[j],"-",x[i])]]  <-df
}
  

}



my_rbind <- function(df1, df2){    # Create own merging function
  rbind(df1, df2)
}
alldf=Reduce(my_rbind, df.list)  


#order
alldf$x=factor(alldf$x,level=c("cov","iav","rsv"))
alldf$y=factor(alldf$y,level=c("cov","iav","rsv"))



#plot
ggplot(alldf,aes(x=vars,y=allrr))+
  geom_ribbon(aes(ymin=allrr.lci, ymax=allrr.uci,fill=region), alpha = 0.2) +
  geom_line(size=0.8,aes(color=region)) +
   theme(strip.text.x = element_text(size = 18, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
  theme(strip.text.y = element_text(size = 18, colour = "black")) + # 设置分面的字字体大小、颜色、背景、边框，
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(y~x,scales="free",drop=F) +
  xlab("Explanatory (%)")+
  ylab("Risk")



################################################################################

alldf$x <- gsub("cov", "SARS-CoV-2", alldf$x)
alldf$x <- gsub("iav", "IAV", alldf$x)
alldf$x <- gsub("rsv", "RSV", alldf$x)
alldf$y <- gsub("cov", "SARS-CoV-2", alldf$y)
alldf$y <- gsub("iav", "IAV", alldf$y)
alldf$y <- gsub("rsv", "RSV", alldf$y)
alldf$region <- gsub("US", "USA", alldf$region)

virus=unique(alldf[,c("x","y")])
colors <- brewer.pal(6, "Accent")
plot.list<-NULL

for (i in 1:6) {
  
  alldf_new=alldf[which(alldf$x==as.vector(virus[i,1])&
                  alldf$y==as.vector(virus[i,2])),]
  

  p<-ggplot(alldf_new)+
    geom_ribbon(aes(x=vars,ymin=(allrr.lci), ymax=(allrr.uci),fill=region), alpha = 0.1) +
    geom_line(aes(x=vars,y=allrr,color=region),size=0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    scale_color_manual(values = colors,breaks =unique(alldf$region) ) +  # 手动指定颜色
      theme_classic()+
    theme(axis.text.x = element_text(size = 9, color="black",vjust = 0.5, angle = 0))+
    theme(axis.text.y = element_text(size = 9, color="black",vjust = 0.5, angle = 0))+
    xlab(paste0(unique(alldf_new$x)," percentile (%)"))+
    ylab(paste0("RR of ",unique(alldf_new$y)))+
    theme(legend.position = "none") 
  
  plot.list[[paste0(virus[i,1],virus[i,2])]]=p
  
}


library(ggpubr)


arranged_plot <- ggarrange(plot.list[[1]], plot.list[[3]], 
                           plot.list[[2]],plot.list[[5]],
                           plot.list[[4]],plot.list[[6]], 
                           
                           nrow=3,ncol=2, align = 'hv',
                           labels = c("A", "B","C", "D","E", "F"))
arranged_plot 


# plot1 with legend 
plot1_legend <- ggplot(alldf_new)+
  geom_ribbon(aes(x=vars,ymin=(allrr.lci), ymax=(allrr.uci),fill=region), alpha = 0.1) +
  geom_line(aes(x=vars,y=(allrr),color=region),size=0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  guides(fill = guide_legend(ncol = 6))+
   theme(legend.position = "bottom",
         legend.title=element_blank()) 

# function to extract legend from plot 
get_only_legend <- function(plot) { 
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  legend <- plot_table$grobs[[legend_plot]] 
  return(legend) 
} 

# extract legend from plot1 using above function 
legend <- get_only_legend(plot1_legend)    

# final combined plot with shared legend 
p<-grid.arrange(arranged_plot , legend, nrow = 2, heights = c(10, 1))

ggsave("../figure/figS2.tiff", plot = p, width=18, height=20, units="cm", dpi=300, compression = "lzw")
