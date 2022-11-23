suppressMessages({
library(corrgram)
library(ggplot2)
library(RColorBrewer)
library(ggunchained) 
library(ggsci)
library(stringr)
})
args <- commandArgs(T)
args <- c("/home/zhl/projects/batch_effect_correction/new_results/protein_after_batch")

melt_data <- function(df_cor,sample,method.name){
  mydata<- reshape2::melt(df_cor)
  mydata$group <- sample
  mydata$method <- method.name
  mydata$value[mydata$value==1] <- NA
  mydata <- na.omit(mydata)
  mydata <- mydata[,c(4,3,5)]
  return(mydata)
}

path=args[1]
files <- dir(path, pattern = ".csv")
# files <- files[-c(6:7)]
Data <- data.frame()
if (length(files)!=0){
  for (i in 1:length(files)){
    #method type
    method.name <-strsplit(files[i],'protein_after_batch_')[[1]][2]
    method.name <- strsplit(method.name,'.csv')[[1]][1]
    #input data
    df <-read.csv(paste(path,files[i],sep='/'), header=T,check.names = F,quote = "",row.names = 1)
    df <-  df[,-1]
    df[df==0] <- NA
    df <- log10(df)
    df <- data.frame(df)
    df_A <-  df[,colnames(df)[grepl("_A_",colnames(df))]]
    df_B <-  df[,colnames(df)[grepl("_B_",colnames(df))]]
    df_A_cor <- round(cor(df_A, use="pair"),3)
    df_B_cor <- round(cor(df_B, use="pair"),3)
    #
    mydata_A <- melt_data(df_cor=df_A_cor,sample="A",method.name=method.name) 
    mydata_B <- melt_data(df_cor=df_B_cor,sample="B",method.name=method.name)
    mydata <- rbind(mydata_A,mydata_B)
    Data <- rbind(Data,mydata)
    }
  }

method <- c("raw","Combat","Limma","Median-Norm","QC-Mnorm","RUVg","RUVIII","RUVs")
#method <- c("raw","Combat","Median-Norm","QC-Mnorm","RUVg","RUVIII","RUVs")
Data$method <- factor(Data$method,levels = method)

p1 <- ggplot(Data, aes(x = method, y = value,fill= group))+
  #geom_split_violin(trim = T,color = "white",alpha=0.85)+
  #geom_point(stat = "summary",fun=median,position = position_dodge(width = 0.9))+
  #geom_boxplot(width=0.15,outlier.shape = NA,coef=0.45,color="black",lwd=0.1,fatten=8)+
  geom_boxplot(outlier.shape = NA,alpha=0.85)+
  scale_fill_manual(values = brewer.pal(3,'Set1'))+theme_bw(base_rect_size = 1)+
  theme(#legend.position = c(0.038,0.948),
        legend.position = c(0.964,0.943),#c(0.964,0.948)
        legend.title=element_blank(),
        panel.grid=element_blank(),
        legend.background = element_rect(color = "black", linetype = "solid", size = 0.25),
        axis.text.x=element_text(color='black',size = 15,face = "bold",hjust=1,angle = 60),
        axis.text.y=element_text(color='black',size = 15,face = "bold"),
        axis.title = element_text(size = 16,color = "black"))+labs(y="Pearson Correlation Coefficient",x="")+
  scale_y_continuous(breaks = seq(0.7,1,0.05),limits = c(0.6,1))
ggsave(paste(path,"Pearson_boxplot.pdf",sep='/'),p1,width = 8.5,height = 8)
ggsave(paste(path,"Pearson_boxplot.png",sep='/'),p1,width = 8.5,height = 8)


# subData <-  subset(Data,method!="RUVs")
p2 <- ggplot(Data, aes(x = method, y = value,fill= group))+
  geom_split_violin(trim = T,color = "grey50",alpha=0.8,size=0.1)+
  #geom_point(stat = "summary",fun=median,position = position_dodge(width = 0.9))+
  geom_boxplot(width=0.15,outlier.shape = NA,coef=0.45,color="black",lwd=0.1,fatten=8)+
  #geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = brewer.pal(3,'Set1'))+theme_bw(base_rect_size = 1)+
  theme(legend.position = c(0.964,0.943),
    legend.title=element_blank(),
    legend.background = element_rect(color = "black", linetype = "solid", size = 0.25),
    panel.grid.major = element_blank(),
    panel.grid=element_blank(),
    axis.text.x=element_text(color='black',size = 15,face = "bold",hjust=1,angle = 60),
    axis.text.y=element_text(color='black',size = 15,face = "bold"),
    axis.title = element_text(size = 16,color = "black"))+labs(y="Pearson Correlation Coefficient",x="")+
  scale_y_continuous(breaks = seq(0.7,1,0.05),limits = c(0.6,1))
ggsave(paste(path,"Pearson_violin.png",sep='/'),p2,width = 8.5,height = 8)
ggsave(paste(path,"Pearson_violin.pdf",sep='/'),p2,width = 8.5,height = 8)

