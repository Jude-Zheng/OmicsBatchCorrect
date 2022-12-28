library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(stringr)
options(warn = -1)
args=commandArgs(T)
data <- read.csv(args[1],header=T,comment.char = "",check.names=F)
filename <- str_split_fixed(args[1],"rsd_protein_after_batch_",2)[,2]
filename <- str_split_fixed(filename,".csv",2)[,1]
out1 <- paste0(filename,"_RSD.png")
out2 <- paste0(filename,"_RSD.pdf")

colnames(data) <- c("Sample","RSD")
labs <- c("Lab1","Lab2","Lab3_1","Lab3_2","Lab4","Lab5","Lab6","Lab7","Lab8","Lab9","Lab10","Lab11")
data$Sample <- factor(data$Sample,levels = as.vector(labs))
w=length(unique(data$Sample))*1.125+2.25
p=ggplot(data,aes(y=RSD,x=Sample))+
  stat_boxplot(geom="errorbar",width=0.3)+
  geom_boxplot(aes(fill=Sample),outlier.colour = NA)+ylim(0,1)+
  #scale_fill_lancet()+
  scale_fill_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  scale_color_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  theme_bw()+
  guides(fill=guide_legend(title="Lab type"))+
  theme(axis.text.x = element_text(size=15,color='black',face = "bold",hjust=1,angle = 60),
        axis.text.y = element_text(size=15,color='black',face = "bold"),
        axis.title= element_text(size=22,color='black'),
        legend.title=element_text(size=22,color='black'),
        legend.text=element_text(size=15,color='black'),
        panel.background=element_rect(fill="white",colour="black"),
        #panel.grid.major=element_line(color="gray"),
        #panel.grid.minor=element_line(color="gray"),
        panel.grid=element_blank())+xlab("")
#ggsave(paste(args[2],"RSD.pdf", sep='/'),p,height=6,width=w)
ggsave(paste(args[2],out1, sep='/'),p,height=6,width=9)
ggsave(paste(args[2],out2, sep='/'),p,height=6,width=9)