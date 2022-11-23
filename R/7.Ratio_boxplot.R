library(ggplot2)
library(ggsci)
library(stringr)
args=commandArgs(T)
# args <- c("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/uniprot-Ecoli_K12_201709/Ratio/Peptide_Ratio.csv",
#           "/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/uniprot-Ecoli_K12_201709/Ratio/Peptide_Ratio.png")
df <- read.csv(args[1],row.names = 1,check.names = F)

if (grepl("Ecoli",args[1])) {
  ylabs = "E.coli: log2(A/B)"
  yratio_line = c(2)
  #limits = c(0.2,4)
  #breaks=seq(0.2,4,0.15)
  limits = c(-4,4)
  breaks=seq(-5,5,1)
} else if (grepl("Saccharomyces",args[1])) {
  ylabs = "Yeast: log2(A/B)"
  yratio_line = c(-1)
  # limits = c(0.2,1.5)
  # breaks=seq(0.2,1.5,0.15)
  limits = c(-4,4)
  breaks=seq(-4,4,1)
}else{
  ylabs = "Human: log2(A/B)"
  yratio_line = c(0)
  # limits = c(0.2,1.5)
  # breaks=seq(0.2,1.5,0.15)
  limits = c(-4,4)
  breaks=seq(-4,4,1)
}



df_t <- reshape2::melt(df,variable.name = "Method")
df_t <- na.omit(df_t)
df_t$value <- log2(df_t$value+0.00001)
method <- c("raw","Combat","Limma","QC-Mnorm","Median-Norm","QC-NORM","RUVs","RUVg","RUVIII")
df_t$Method <- factor(df_t$Method,levels = method)
p <- ggplot(df_t,aes(y=value,x=Method,fill=Method))+#geom_violin(adjust = 1,trim = FALSE,alpha=0.8)+
  geom_boxplot(position=position_dodge(2),alpha=0.8,outlier.size = 0.5)+
  scale_fill_lancet()+theme_test()+
  geom_hline(yintercept = yratio_line,lty=2,col=c("red"))+
  theme(axis.text=element_text(color='black'),
        axis.text.x = element_text(color='black',size = 15,face = "bold",hjust=1,angle = 60),
        axis.text.y=element_text(color='black',size = 15,face = "bold"),
        axis.title=element_text(size=17,color='black'),
        legend.position = "none",
        panel.grid=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(colour=NA),
        panel.background=element_rect(fill="white",colour="black"))+
  labs(y=ylabs,x="")+scale_y_continuous(limits = limits,breaks=breaks)

ggsave(args[2],p,width = 8.5,height = 7)
ggsave(args[3],p,width = 8.5,height = 7)
