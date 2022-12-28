suppressMessages({
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(ggsci)
library(data.table)
library(variancePartition)
library(factoextra)
library(rafalib)#as.fumeric
})
#===============raw data============
#setwd("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/Correction/maxLFQ")
# args <- c("/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/protein_after_batch/protein_after_batch_raw.csv",
#           "/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/group.xls",
#           "/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/QC/Protein")
# args <- c("/home/zhl/projects/batch_effect_correction/new_results/protein_after_batch/protein_after_batch_Limma.csv",
#           "/home/zhl/projects/batch_effect_correction/new_results/group_new.xls",
#           "/home/zhl/projects/batch_effect_correction/new_results/QC")
#===================================
args=commandArgs(T)
data <- read.csv(args[1])
group <- readxl::read_xls(args[2])
colnames(group)[1:3]=c("ID","Type","injection")
#==============================================================

if (grepl("peptides_raw_matrix", args[1])) {
  filename <- "peptides_raw"
  pca_name <- paste0(args[3],"/","PCA_",filename,".png")
  box_name <- paste0(args[3],"/","boxplot_",filename,".png")
  point_name <- paste0(args[3],"/","pointplot_",filename,".png")
  variancePartition <- paste0(args[3],"/","variancePartition_",filename,".png")
  hcluster <- paste0(args[3],"/","HCluster_",filename,".png")
  
  pca_name_ <- paste0(args[3],"/","PCA_",filename,".pdf")
  box_name_ <- paste0(args[3],"/","boxplot_",filename,".pdf")
  point_name_ <- paste0(args[3],"/","pointplot_",filename,".pdf")
  variancePartition_ <- paste0(args[3],"/","variancePartition_",filename,".pdf")
  hcluster_ <- paste0(args[3],"/","HCluster_",filename,".pdf")
} else {
filename <- str_split_fixed(args[1],"protein_after_batch_",2)[,2]
filename <- str_split_fixed(filename,".csv",2)[,1]
pca_name <- paste0(args[3],"/","PCA/","PCA_",filename,".png")
box_name <- paste0(args[3],"/","boxplot/","boxplot_",filename,".png")
point_name <- paste0(args[3],"/","pointplot/","pointplot_",filename,".png")
variancePartition <- paste0(args[3],"/","variancePartition/",filename,".png")
hcluster <- paste0(args[3],"/","Hcluster/","HCluster_",filename,".png")

pca_name_ <- paste0(args[3],"/","PCA/","PCA_",filename,".pdf")
box_name_ <- paste0(args[3],"/","boxplot/","boxplot_",filename,".pdf")
point_name_ <- paste0(args[3],"/","pointplot/","pointplot_",filename,".pdf")
variancePartition_ <- paste0(args[3],"/","variancePartition/",filename,".pdf")
hcluster_ <- paste0(args[3],"/","Hcluster/","HCluster_",filename,".pdf")
}
#===============================================================
group <- merge(group,data.frame(ID=colnames(data)))
group$Batch <- paste0(group$Batch,"b")
group <- group[order(group$order),]
subdata <- subset(data,select = as.vector(group$ID))
subdata[subdata<1] <- NA
subdata <- na.omit(subdata)
subdata_plot=log2(subdata+1)

#==========boxplot============
group$Type <- factor(group$Type,levels = as.vector(unique(group$Type)))
batch_nu <- cumsum(table(group$Type))[1:(length(cumsum(table(group$Type)))-1)]

box.data=reshape2::melt(subdata_plot,measure.vars=colnames(subdata_plot),variable.name = "Sample")
box.data <- na.omit(box.data)
box.data <- merge(box.data,group,by.x = "Sample",by.y = "ID")
p3 = ggplot(box.data,aes(y=value,x=Sample,fill=Type))+#geom_violin(trim = FALSE,alpha=0.8)+
  geom_boxplot(position=position_dodge(2),outlier.colour = NA,alpha=0.8,lwd=0.15)+ #outlier.size=0.5,
  scale_fill_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  scale_color_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  #scale_fill_d3()+#scale_color_npg()
  geom_vline(xintercept = batch_nu,lty=2,col="grey") +
  theme_test()+
  theme(axis.text=element_text(color='black'),
        axis.text.x = element_text(hjust=0.8,angle = 60,size = 4.5),
        axis.text.y = element_text(size = 16,face = "bold"),
        axis.title = element_text(size=18),
        legend.key.size = unit(12, "pt"),
        #axis.title=element_text(size=15,color='black'),
        legend.position = "top")+xlab("Samples injection")+
  ylab("Log2 (Intensity)")+scale_x_discrete(labels = group$injection)+scale_y_continuous(limits = c(11,23),breaks = seq(11,23,2))

ggsave(box_name,p3,height = 8.5,width = 10.5)
ggsave(box_name_,p3,height = 8.5,width = 10.5)
#==========point PLOT=============
point_data <- as.data.frame(apply(subdata_plot, 2, function(x) median(x,na.rm = TRUE)))
point_data$ID <- row.names(point_data)
colnames(point_data)[1] <- c("intensity")
point_data <- merge(point_data,group)
point_data$Type <- factor(group$Type,levels = as.vector(unique(group$Type)))
p4 = ggplot(point_data,aes(x=ID,y=intensity,colour=Type,group=Type))+geom_point(size = 2.5)+#geom_line(size=1)+
  theme_test()+geom_vline(xintercept = batch_nu,lty=2,col="grey") +
  #scale_color_D3()+
  scale_color_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  geom_hline(yintercept=c(16.9,17.2),linetype="dashed",color="grey50")+
  theme(axis.text=element_text(color='black'),
        axis.text.x = element_text(hjust=0.8,angle = 60,size = 4.5),
        axis.text.y = element_text(size = 16,face = "bold"),
        axis.title = element_text(size=18),
        #legend.text=element_text(size=6),
        legend.key.size = unit(13, "pt"),
        legend.position = "top")+
  scale_x_discrete(labels = group$injection)+
  #scale_y_continuous(limits = c(min(point_data$intensity)-1.5,max(point_data$intensity)+1.5))+
  scale_y_continuous(limits = c(16.5,17.5),breaks = seq(16.5,17.5,0.2))+
  labs(y="Log2 (Median Intensity)",x="Samples injection")
ggsave(point_name,p4,height = 7,width = 10.5)
ggsave(point_name_,p4,height = 7,width = 10.5)
#======================PCA===============================
subdata.pca=t(na.omit(subdata_plot))
data.pca = prcomp(subdata.pca,scale. = TRUE)
pca=as.data.frame(data.pca$x)
pca$ID=rownames(pca)
pca=merge(pca,group,by="ID")
percent1=sprintf("%.1f%%",data.pca$sdev[1]^2/sum(data.pca$sdev^2)*100)
percent2=sprintf("%.1f%%",data.pca$sdev[2]^2/sum(data.pca$sdev^2)*100)
xlabel=paste("PC1 (",percent1,")")
ylabel=paste("PC2 (",percent2,")")
p2=ggscatter(pca,x="PC1",y="PC2",color="Type",ellipse=F,ellipse.level = 0.95,shape ="Sample", #label="ID",
          ellipse.type='confidence',ellipse.alpha=0.3,
          show.legend.text=F,size=2.5,font.label=c(10,"plain","black"),repel=T)+
  scale_fill_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  scale_color_manual(values = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2')))+
  theme_bw()+xlab(xlabel)+ylab(ylabel)+guides(color=guide_legend(title="Lab type"),fill="none")+
  theme(axis.text=element_text(size=18,color='black'),
        axis.title=element_text(size=18,color='black'),
        #legend.title=element_text(size=18,color='black'),
        #legend.text=element_text(size=15,color='black')
        legend.key.size = unit(12, "pt"))
ggsave(pca_name,p2,height = 7.5,width = 8.5)
ggsave(pca_name_,p2,height = 7.5,width = 8.5)
##########################################################
####### HIERARCHICAL CLUSTER ANALYSIS
##########################################################
# type <- as.character(group$Type) # 4 for batch; 2 for lab
# # number of groups
# k <- length(unique(type)) # groups in HC
# mat <- as.data.frame(t(subdata))
# 
# Cols = function(vec, ord){
#   cols = c(pal_lancet(palette = c("lanonc"), alpha = 1)(9),rev(brewer.pal(5,'Set2')))[1:length(unique(vec))]
#   return(cols[as.fumeric(vec)[ord]])}
# # Cols = function(vec, ord){
# #   cols = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2'))
# #   return(cols[as.fumeric(vec)[ord]])}
# #mat_1 <- data.frame(scale(mat, center = T, scale = T))
# rownames(mat) = paste0(group$Type,".",group$injection," ",group$Sample)
# res.dist <- dist(mat, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
# res.hc <- hclust(d = res.dist, method = "ward.D2") #{ward (ward.D)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}
# 
# hca <- fviz_dend(res.hc, #k = k, # Cut in k groups
#                  k=2,k_colors = c("grey50", "black"),
#                  cex = 0.62, # label size 0.3/0.7
#                  #k_colors = unique(Cols(type,res.hc$order)),
#                  color_labels_by_k = TRUE, # 数据标签也根据颜色设定
#                  label_cols = Cols(type,res.hc$order),
#                  rect = T, # Add rectangle around groups
#                  rect_fill = T,
#                  #rect_border = unique(Cols(type,res.hc$order)),
#                  horiz = F,
#                  type = "circular",
#                  lwd = 0.6, # lines size 0.3/0.7
#                  show_labels = T,
#                  main = "",
#                  ylab = "") 
# 
# hp <- hca + theme(axis.text.y = element_blank())+scale_shape_manual(values=rep(0:length(unique(type)))) # number of groups
# ggsave(hcluster,hp,height = 8,width = 8)
# ggsave(hcluster_,hp,height = 8,width = 8)

#==================variancePartition=====================
# varpart <- function(femData,info){
#   femData[is.na(femData)] <- 0
#   form <- ~ (1|Batch)
#   varPart <- fitExtractVarPartModel( femData, form, info )
#   vp <- sortCols( varPart )
#   plotPercentBars( vp[1:30,] )
#   plotVarPart(vp)
# }
# p1=varpart(subdata_plot,group)
# p1_ = p1+theme(axis.text=element_text(size=19,color='black'),
#     axis.text.x=element_text(size=20,color='black'),
#     axis.title=element_text(size=21,color='black'))
# ggsave(variancePartition,p1_,height = 8,width = 7)
# ggsave(variancePartition_,p1_,height = 8,width = 7)
