library(ggpubr)
args <- c("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/Correction/peptides_after_MVI_Combat.csv",
          "/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/group.xls",
          "/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/")

data <- read.csv(args[1],row.names = 1)
data <- as.data.frame(t(data))
group <- readxl::read_xls("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/group.xls")
colnames(group)[1:3]=c("ID","Type","injection")
ds=merge(group,data,by.x = "ID",by.y = "row.names")
rownames(ds) <- ds$ID
ds <- ds[order(ds$order),]
#QC data
ds_rsd <- subset(ds, Rep == "QC")
ds_rsd <- ds_rsd[,-c(1:7)]
file4=paste(args[3],"RSD_QC_after_raw.png",sep="/")
file5=paste(args[3],"barplot_QC_after_raw.png",sep="/")
#============================================
#QC RSD
#============================================
qc_list <- list()
for(i in seq(2, nrow(ds_rsd), by=1)){
    #print(i)
    #x=i-1
    a = 1:i
    x=c(i/2)
    qc_list[[i]]=a
  }

qc_list = qc_list[-1]

RSD_QC_results <- list()
RSD_QC_func <- function (x) { 
  for (i in (1:length(qc_list))){
    RSD_QC_results[[i]] <- apply(x[qc_list[[i]],], 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100)}
  
  for (i in (1:length(RSD_QC_results))){
    RSD_QC_results[[i]] <- data.frame(colnames(x),RSD_QC_results[[i]])}
  
  for (i in (1:length(RSD_QC_results))){
    colnames(RSD_QC_results[[i]]) <- c("name", "rsd")}
  
  n <- c(paste(c(2:c(length(qc_list)+1)), "number of QC"))
  names(RSD_QC_results) <- n
  RSD_QC_results
}

rsd <- RSD_QC_func(ds_rsd) 

cutoff1 <- 30 # adjust to your data
RSD_QC_func_cutoff1 <- c()
for (i in (1:length(rsd))){
  RSD_QC_func_cutoff1[i] <- round(length(which(rsd[[i]]$rsd <= cutoff1))/length(rsd[[i]]$rsd)*100,0)
  print(RSD_QC_func_cutoff1[i])}
n <- c(paste(c(2:(length(qc_list)+1)), "number of QC")) # names
RSD_QC_func_cutoff1 <- data.frame("QC" = n,RSD_QC_func_cutoff1)
colnames(RSD_QC_func_cutoff1)[2] <- c(paste(cutoff1, "% RSD"))

ggbarplot(RSD_QC_func_cutoff1, 'QC', '30 % RSD',position=position_dodge(0.7), label = T)+theme_test()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),legend.title=element_blank())


ggsave(file4,width=10,height=8.5)
#========================================
# QC boxplot
#========================================``
group_QC <- group[group$Type=="QC",]
group_QC$order <- factor(group_QC$order,levels = as.vector(unique(group_QC$order)))
ID = unique(as.vector(group_QC$ID))
data_qc <- log2(subset(data,select=ID))
#BOX
box.qc=reshape2::melt(data_qc,measure.vars=colnames(data_qc),variable.name = "Sample")
box.qc <- na.omit(box.qc)
box.qc <- merge(box.qc,group,by.x = "Sample",by.y = "ID")

ggplot(box.qc,aes(y=value,x=Sample,fill=batch))+#geom_violin(trim = FALSE,alpha=0.8)+
  geom_boxplot(position=position_dodge(2),outlier.colour = NA,alpha=0.8,lwd=0.1)+ #outlier.size=0.5,
  scale_fill_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  scale_color_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  #scale_fill_d3()+#scale_color_npg()
  #geom_vline(xintercept = batch_nu,lty=2,col="grey") + 
  geom_hline(yintercept = median(box.qc$value),lty=2,col="black",size=0.3) + 
  theme_test()+
  theme(axis.text=element_text(color='black'),
        axis.text.x = element_text(hjust=1,angle = 60),
        legend.key.size = unit(12, "pt"),
        #axis.title=element_text(size=15,color='black'),
        legend.position = "top")+
  labs(y="Log2 (Intensity)",x="order")#+scale_y_continuous(limits = c(15,38))
#file4=paste(args[3],"barplot_QC.pdf",sep="/")

#ggsave(file4,width=10,height=8.5)
ggsave(file5,width=10,height=8.5)