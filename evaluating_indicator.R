args <- c("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/Correction/peptides_after_MVI_raw.csv",
          "/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/group.xls",
          "/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC")
args <- c("/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/maxLFQ/Protein_raw_maxLFQ.csv",
          "/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/group.xls",                     
          "/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/maxLFQ/")
data <- read.csv(args[1],row.names = 1)
data=as.data.frame(t(data))
group <- readxl::read_xls(args[2])
colnames(group)[1:3]=c("ID","Type","injection")
ds=merge(group,data,by.x = "ID",by.y = "row.names")
rownames(ds) <- ds$ID
ds <- ds[order(ds$order),]

###############################################################################
########################################## PVCA
###############################################################################
library(proBatch)
# generate batch data
batch <- as.numeric(ds$batch)

# generate group data
class <- as.character(ds$Rep)
qc_id <- which(class == "QC")
class[-qc_id] <- "Subject"

# generate lab
lab <- as.character(ds$Type)

# Sample data
s_d <- data.frame(FullRunName = rownames(ds), batch, lab)
rownames(s_d) <- NULL

# Feature data
f_d <- as.matrix(t(ds[,-c(1:7)]))
colnames(f_d) <- rownames(ds)

# Perform
pvca_df <- calculate_PVCA(f_d, s_d, 
                          factors_for_PVCA = c('batch','lab'), # adjust to your data
                          pca_threshold = .6, variance_threshold = .01, fill_the_missing = 0)

p = plot_PVCA(f_d, s_d,
          technical_factors = "batch",
          biological_factors = 'lab')
ggsave(paste(args[3],"PVCA.png",sep = "/"),p,width=8.5,height = 7.5)
##########################################################
####### HIERARCHICAL CLUSTER ANALYSIS
##########################################################
library(factoextra)
library(RColorBrewer)
library(rafalib)
library(ggsci)
type <- ds[,2] # 4 for batch; 2 for lab
# number of groups
k <- length(unique(type)) # groups in HC
mat <- ds[,-c(1:7)] 

Cols = function(vec, ord){
  cols = c(pal_lancet(palette = c("lanonc"), alpha = 1)(9),rev(brewer.pal(5,'Set2')))[1:length(unique(vec))]
  return(cols[as.fumeric(vec)[ord]])}

# Cols = function(vec, ord){
#   cols = c(brewer.pal(7,'Set1'),brewer.pal(5,'Set2'))
#   return(cols[as.fumeric(vec)[ord]])}


mat_1 <- mat
#mat_1 <- data.frame(scale(mat, center = T, scale = T))
rownames(mat_1) = paste(ds$Type,ds$injection,ds$Sample)
res.dist <- dist(mat_1, method = "manhattan") #{euclidean}, {maximum}, {manhattan}, {canberra}, {binary}, {minkowski}
res.hc <- hclust(d = res.dist, method = "ward.D2") #{ward (ward.D)}, {single}, {complete}, {average}, {mcquitty},{median}, {centroid}

hca <- fviz_dend(res.hc, #k = k, # Cut in k groups
                 cex = 0.45, # label size 0.3/0.7
                 #k_colors = unique(Cols(type,res.hc$order)),
                 color_labels_by_k = TRUE, # 数据标签也根据颜色设定
                 label_cols = Cols(type,res.hc$order),
                 rect = T, # Add rectangle around groups
                 rect_fill = T,
                 #rect_border = unique(Cols(type,res.hc$order)),
                 horiz = F,
                 type = "circular",
                 lwd = 0.6, # lines size 0.3/0.7
                 show_labels = T,
                 main = "",
                 ylab = "") 

p <- hca +scale_shape_manual(values=rep(0:length(unique(type)))) # number of groups
ggsave("Cluster_Dendrogram.png",p,height = 8.5,width = 8.5)
#============================================
#QC RSD
#============================================
# after processing from "Generate data for evaluation"
ds_rsd <- subset(ds, Rep == "QC")
ds_rsd <- ds_rsd[,-c(1:7)]
ds_rsd <- ds_rsd[order(ds_rsd$batch, decreasing = F),]

# Select batch by pattern
tn <- ds_rsd$batch # batch variables
tnu <- unique(tn)
vb <- list()
for (i in (1:length(tnu))){
  vb[[i]] <- which(tn == tnu[i])}

# RSD calculation by batches
RSD_by_batch_results <- list()
RSD_by_batch <- function (x) { 
    for (i in (1:length(vb))){
      RSD_by_batch_results[[i]] <- apply(x[vb[[i]],], 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100)}
    # all batches
      RSD_by_batch_results[[(length(vb)+1)]] <- apply(x, 2, function(y) sd(y, na.rm = T)/mean(y, na.rm = T)*100)
      
    for (i in (1:length(RSD_by_batch_results))){
      RSD_by_batch_results[[i]] <- data.frame(colnames(x),RSD_by_batch_results[[i]])}
      
    for (i in (1:length(RSD_by_batch_results))){
      colnames(RSD_by_batch_results[[i]]) <- c("name", "rsd")}
      
    n <- c(paste(c(1:length(vb)), "batch"), "all batches")
    names(RSD_by_batch_results) <- n
    RSD_by_batch_results
    }


ds_rsd1 <- ds_rsd[,-1] # sample in row only metabolites variables
rsd <- RSD_by_batch(ds_rsd1) 

# count in % by cut-off
# cutoff 1
cutoff1 <- 20 # adjust to your data
RSD_by_batch_cutoff1 <- c()
for (i in (1:length(rsd))){
  RSD_by_batch_cutoff1[i] <- round(length(which(rsd[[i]]$rsd <= cutoff1))/length(rsd[[i]]$rsd)*100,0)
  print(RSD_by_batch_cutoff1[i])}
n <- c(paste(c(1:length(vb)), "batch"), "all batches") # names
RSD_by_batch_cutoff1 <- data.frame("batch" = n,RSD_by_batch_cutoff1)
colnames(RSD_by_batch_cutoff1)[2] <- c(paste(cutoff1, "% RSD"))

# cutoff 2 
cutoff2 <- 30 # adjust to your data
RSD_by_batch_cutoff2 <- c()
for (i in (1:length(rsd))){
  RSD_by_batch_cutoff2[i] <- round(length(which(rsd[[i]]$rsd <= cutoff2))/length(rsd[[i]]$rsd)*100,0)}
RSD_by_batch_cutoff2 <- data.frame("batch" = n,RSD_by_batch_cutoff2)
colnames(RSD_by_batch_cutoff2)[2] <- c(paste(cutoff2, "% RSD"))

# result
result_RSD <- rbind(t(RSD_by_batch_cutoff1), t(RSD_by_batch_cutoff2))
result_RSD

