library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(stringr)
options(warn = -1)
args=commandArgs(T)
args = c("/home/zhl/projects/batch_effect_correction/new_results/QC/RSD_/RSD_method.csv",
         "/home/zhl/projects/batch_effect_correction/new_results/QC/RSD_/")
data <- read.csv(args[1],header=T,comment.char = "",check.names=F)
out1 <- paste(args[2],"RSD_method.png",sep = "/")
out2<- paste(args[2],"RSD_method.pdf",sep = "/")
colnames(data) <- c("Sample","RSD")
data <- data[data$Sample!="Limma",]
if (grepl("protein", args[1])) {
#protein
Method <- c("raw","Combat","Limma","Combat_seq", "EigenMS","QC-KNN","QC-XGB","QC-CTB","QC-BT","QC-DT","QC-Mnorm","Median-Norm",
  "RUVIII","RUVg","RUVs","RUVr","WaveICA","WaveICA2")
} else {
#peptides
Method <- c("raw","Combat","Limma","Combat_seq", "EigenMS","QC-KNN","QC-XGB","QC-CTB","QC-BT","QC-DT","QC-Mnorm","Median-Norm",
            "RUVIII","RUVIII_C","RUVg","RUVs","RUVr","WaveICA","WaveICA2")}
Method <- rev(c("raw","RUVIII","QC-Mnorm","RUVs","RUVg","Median-Norm","Combat"))
data$Sample <- factor(data$Sample,as.vector(Method))
data <- data[data$Sample!="Limma",]
w=length(unique(data$Sample))*1.125+2.5
p=ggplot(data,aes(y=RSD,x=Sample))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(aes(fill=Sample),outlier.colour = NA)+ylim(0,1)+
  #scale_fill_lancet()+
  scale_fill_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2')))+
  scale_color_manual(values = c(brewer.pal(9,'Set1'),brewer.pal(8,'Set2')))+
  theme_bw()+guides(fill=guide_legend(title="Method type"))+
  theme(axis.title=element_text(size=20,color='black'),
        axis.text.x = element_text(size=15,color='black',face = "bold",hjust=1,angle = 60),
        axis.text.y = element_text(size=15,color='black',face = "bold"),
        #legend.position="none",
        panel.background=element_rect(fill="white",colour="black"),
        #panel.grid.major=element_line(color="gray"),
        #panel.grid.minor=element_line(color="gray"),
        panel.grid=element_blank())+labs(y = "RSD",x = "")
p
ggsave(out1,p,height=7.5,width=9.5)
ggsave(out2,p,height=7.5,width=9.5)
