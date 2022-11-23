##############################################################################################################################################################
# Perform correction
##############################################################################################################################################################
# Content:
# 1. Combat (Combat)
# 2. RUVs
# 3. RUVg
# 4. RUVIII (RUVIII)
# 5. QC-norm
###############################################################################
#data input
###############################################################################
library(ggplot2)
library(data.table)
args <- commandArgs(T)
# args <- c("/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/maxLFQ/ds.RData","/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/maxLFQ/Protein_raw_maxLFQ_impute0.4.csv","/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/protein_after_batch/")
load(args[1])
data_name <- as.data.frame((fread(input = args[2], header=T,check.names = F,quote = "")))[,1:2]
##############################
#raw
##############################
raw <-merge(data_name,round(t(ds[,-c(1:6)]),3),by.x = "Protein",by.y = "row.names")
fwrite(raw, paste0(args[3],"protein_after_batch_raw.csv"))
print("raw finish!!")

###############################################################################
#1. Combat (Combat)
###############################################################################
library(sva)
# generate batch data
dat <- as.matrix(t(log2(ds[,-c(1:6)])))
batch <- as.numeric(ds$batch)
model <- model.matrix(~as.factor(ds$Type))
######################### Combat
ds_combat = ComBat(dat=dat, batch=batch, par.prior=TRUE, prior.plots=FALSE) #mod = model
ds_combat = round(as.data.frame(2^ds_combat),3)

ds_combat <-merge(data_name,ds_combat,by.x = "Protein",by.y = "row.names")
# save
fwrite(ds_combat, paste0(args[3],"protein_after_batch_Combat.csv"))
print("Combat finish!!")
###############################################################################
# RUVSeq`s based methods*
###############################################################################
library(RUVSeq)
library(FactoMineR)
library(fpc)
# generate data for RUVSeq
dat <- ds[,-c(1:6)]
dat <- as.matrix(log2(dat))
medians <- colMedians(dat)
dat_stand <- apply(dat, 1, function(x) x-medians)

idxQC <- which(ds$Rep == "QC")
replicates.ind <- matrix(-1, nrow(dat) - length(idxQC) + 1, length(idxQC))
replicates.ind[1,] <- idxQC
replicates.ind[-1,1] <- (1:nrow(dat))[-idxQC]
# generate data for RUVg
QC_data <- as.data.frame(t(dat[idxQC,]))
QC_data$averge <- apply(QC_data, 1, median)
QC_data <- QC_data[order(QC_data$averge,decreasing = TRUE),]
spikes <- rownames(QC_data[1:300,])
######################
#2. RUVs correction*
######################
NumOfComp <- 30 # Max number of components

selectNcomp <- sapply(1:NumOfComp, function(NumOfComp) {
  sprintf("Doing Batch correction using %d number of components\n",
          NumOfComp)
  doRUV <- RUVs(x = dat_stand,
                  cIdx = 1:ncol(dat),
                  k = NumOfComp,
                  scIdx = replicates.ind,
                  round = FALSE,
                  isLog = TRUE)
  doRUV=t(doRUV$normalizedCounts)
  
  pca.ds <- PCA(doRUV[idxQC,], scale.unit = T, graph = F)
  df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components
  
  pca.class <- ds$batch
  pca.class <- pca.class[idxQC]
  scores <- cbind(pca.class, df.pca.x)
  means <- lapply(1:length(unique(scores[,1])), function(y) colMeans(scores[scores$pca.class == unique(scores[,1])[y],-1] , na.rm = T))
  covs <- lapply(1:length(unique(scores[,1])), function(y) cov(scores[scores$pca.class == unique(scores[,1])[y],-1] , ))
  
  batch.dist_b <- matrix(0, length(unique(scores[,1])), length(unique(scores[,1])))
  for (i in 2:length(unique(scores[,1]))) {
    for (j in 1:(i-1)) {
      batch.dist_b[j, i] <- bhattacharyya.dist(means[[j]],
                                               means[[i]],
                                               covs[[j]],
                                               covs[[i]])
    }}
  
  m_b_d <- round(mean(batch.dist_b[col(batch.dist_b) > row(batch.dist_b)]), 2)
  
  return(m_b_d)
})

best_ncomp <- which.min(selectNcomp)


RUVs <- RUVs(x = dat_stand,
             cIdx = spikes,
             k = best_ncomp,
             scIdx = replicates.ind,
             round = F,
             isLog = T)
RUVs <- as.data.frame(RUVs$normalizedCounts)
ruvs <- apply(RUVs, 2, function(x) x+medians)
ruvs <- round(2^ruvs,3)
ruvs <-merge(data_name,ruvs,by.x = "Protein",by.y = "row.names")
# save
fwrite(ruvs, paste0(args[3],"protein_after_batch_RUVs.csv"))
print("RUVs finish!!")
############################
# 3. RUVg correction *
###########################
NumOfComp <- 30 # Max number of components

selectNcomp <- sapply(1:NumOfComp, function(NumOfComp) {
  sprintf("Doing Batch correction using %d number of components\n",
          NumOfComp)
  doRUV <- t(RUVg(x = dat_stand,
                  cIdx = spikes,
                  k = NumOfComp,
                  round = FALSE,
                  isLog = TRUE)$normalizedCounts)
  
  pca.ds <- PCA(doRUV[idxQC,], scale.unit = T, graph = F)
  df.pca.x <- as.data.frame(pca.ds$ind$coord[,1:2]) # or pca.ds$X[,1:n], where n - number of components
  
  pca.class <- ds$batch
  pca.class <- pca.class[idxQC]
  scores <- cbind(pca.class, df.pca.x)
  means <- lapply(1:length(unique(scores[,1])), function(y) colMeans(scores[scores$pca.class == unique(scores[,1])[y],-1] , na.rm = T))
  covs <- lapply(1:length(unique(scores[,1])), function(y) cov(scores[scores$pca.class == unique(scores[,1])[y],-1] , ))
  
  batch.dist_b <- matrix(0, length(unique(scores[,1])), length(unique(scores[,1])))
  for (i in 2:length(unique(scores[,1]))) {
    for (j in 1:(i-1)) {
      batch.dist_b[j, i] <- bhattacharyya.dist(means[[j]],
                                               means[[i]],
                                               covs[[j]],
                                               covs[[i]])
    }}
  
  m_b_d <- round(mean(batch.dist_b[col(batch.dist_b) > row(batch.dist_b)]), 2)
  
  return(m_b_d)
})

best_ncomp <- which.min(selectNcomp)

RUVg <- RUVg(x = dat_stand,
             cIdx = spikes,
             k = best_ncomp,
             round = F,
             isLog = T)
RUVg <- as.data.frame(RUVg$normalizedCounts)
ruvg <- apply(RUVg, 2, function(x) x+medians)
ruvg <- round(2^ruvg,3)
ruvg <-merge(data_name,ruvg,by.x = "Protein",by.y = "row.names")
# save
fwrite(ruvg, paste0(args[3],"protein_after_batch_RUVg.csv"))
print("RUVg finish!!")
###############################################################################
#5. RUVIII
###############################################################################
##############RUVIII 需要QC top pipetides 或者 spiked-in peptides
library(ruv)
#Y,M 
sampleinfo = ds[,c(1:6)]
Y = as.matrix(log2(ds[,-c(1:6)]))
M = replicate.matrix(sampleinfo[,c("Type")])

idxQC <- which(ds$Rep == "QC")
QC_data <- as.data.frame(t(Y[idxQC,]))
QC_data$averge <- apply(QC_data, 1, median)
QC_data <- QC_data[order(QC_data$averge,decreasing = TRUE),]
QC_data$house <- FALSE
QC_data$house[1:300] <- TRUE #top300 peptide
peptideinfo <- data.frame(house=QC_data$house)
rownames(peptideinfo) <- rownames(QC_data)

medians <- colMedians(Y)
Y_stand <- apply(Y, 1, function(x) x-medians)
#RUVIII
ruviii <- RUVIII(t(Y_stand), M, ctl=peptideinfo$house, k=20)
ruviii <- apply(ruviii, 1, function(x) x + medians)
gg_gender_region = list(aes(color=sampleinfo$Type,shape = sampleinfo$Sample),
                        labs(color="lab"))

ruv_svdplot(t(ruviii)) + gg_gender_region

ruviii = round(as.data.frame(2^ruviii),3)
ruviii <-merge(data_name,ruviii,by.x = "Protein",by.y = "row.names")
# save
fwrite(ruviii, paste0(args[3],"protein_after_batch_RUVIII.csv"))
print("RUVIII finish!!")
#############RUVIII_C
# library(RUVIIIC)
# #Y,M
# sampleinfo = ds[,c(1:6)]
# #ds[,-c(1:6)] <- as.numeric(as.character(ds[,-c(1:6)]))
# #ds[,-c(1:6)] <- apply(ds[,-c(1:6)], 2, as.numeric)
# Y = data.matrix(log2(ds[,-c(1:6)]))
# M = replicate.matrix(sampleinfo[,c("Type")])
# #M <- model.matrix(~ Type - 1, data = sampleinfo)
# 
# #QC
# QCinfo <- sampleinfo[sampleinfo$Rep=="QC",]
# QC_data <- as.data.frame(t(subset(ds,Rep=="QC")[,-c(1:6)]))
# QC_data$averge <- apply(QC_data, 1, function(x) mean(x,na.rm=T))
# QC_data <- QC_data[order(QC_data$averge,decreasing = TRUE),]
# QC_data$house <- FALSE
# QC_data$house[1:300] <- TRUE #top300peptide
# 
# QC_peptide <- QC_data[QC_data$house==TRUE,]
# potentialControlsAlwaysFound <- names(which(apply(Y[,rownames(QC_peptide)], 2, function(x) sum(is.na(x))) == 0))
# actualControls <- head(potentialControlsAlwaysFound, 300)
# 
# try(RUVIIIC::omp_set_num_threads(1L), silent=TRUE)
# 
# ruviii_c <- RUVIII_C(k = 16, Y = Y, M = M, toCorrect = colnames(Y),
#                      controls = actualControls)
# 
# gg_gender_region = list(aes(color=sampleinfo$Type,shape = sampleinfo$Sample),
#                         labs(color="lab"))
# ruv_svdplot(na.omit(ruviii_c)) + gg_gender_region
# # ruviii_c <- RUVIII_C_Varying(k = 11, Y = Y, M = M, toCorrect =colnames(Y), 
# #                              potentialControls = actualControls)
# ruviii_c <- as.data.frame(t(ruviii_c))
# ruviii_c <-merge(data_name,ruviii_c,by.x = "Protein",by.y = "row.names")
# # save
# fwrite(ruviii_c, paste0(args[3],"protein_after_batch_RUVIII_C.txt"),sep = "\t")
# print("RUVIII_C finish!!")
###############################################################################
########################################## Median Between Batches (Median-norm)
###############################################################################

# generate batch data
batch <- ds$batch
batch <- as.numeric(batch)

# generate data
ds_batch_dat <- cbind(batch = ds$batch, ds[,-c(1:6)])
ds_batch_dat$batch <- batch

# perform
Median.Norm <- function(data){
library(dplyr)
ds_batch_dat <- data
b_b_c_subsets <- lapply(1:length(unique(ds_batch_dat[,1])), function(y) filter(ds_batch_dat[,-1], ds_batch_dat$batch == unique(ds_batch_dat[,1])[y])) # list of subsets by batches
b_b_c_factor <- lapply(1:length(unique(ds_batch_dat[,1])), function(y) sapply(1:ncol(b_b_c_subsets[[1]]), function(z) median(b_b_c_subsets[[y]][,z], na.rm = T)/median(ds_batch_dat[,(z+1)], na.rm = T))) # calculate factor for every feature
b_b_c_result_l <- lapply(1:length(unique(ds_batch_dat[,1])), function(y) sapply(1:ncol(b_b_c_subsets[[1]]), function(z) b_b_c_subsets[[y]][,z]/b_b_c_factor[[y]][z])) # results by batch
b_b_c_result <- data.frame(do.call("rbind",b_b_c_result_l)) # results in data frame
rownames(b_b_c_result) <- rownames(ds_batch_dat)
colnames(b_b_c_result) <- colnames(ds_batch_dat[,-c(1)])
return(b_b_c_result)
}

ds_median_norm <- Median.Norm(data = ds_batch_dat)
ds_median_norm <- round(as.data.frame(t(ds_median_norm)),3)
ds_median_norm <-merge(data_name,ds_median_norm,by.x = "Protein",by.y = "row.names")
# save
fwrite(ds_median_norm, paste0(args[3],"protein_after_batch_Median-Norm.csv"))
print("Median-Norm finish!!")
###############################################################################
########################################## Median Between QC'Batches (QC-Mnorm)
###############################################################################

# generate qc data
ds_qc <- subset(ds,Rep=="QC")
qc_batch <- as.numeric(ds_qc$batch)
ds_qc_dat <- cbind(batch = ds_qc$batch, ds_qc[,-c(1:6)])
ds_qc_dat$batch <- qc_batch
#generate data
batch <- ds$batch
batch <- as.numeric(batch)
ds_batch_dat <- cbind(batch = ds$batch, ds[,-c(1:6)])
ds_batch_dat$batch <- batch

# perform
QC.Mnorm <- function(qc_data,data){
  library(dplyr)
  ds_qc_dat <- qc_data
  ds_batch_dat <- data
  ds_qc_dat_subsets <- lapply(1:length(unique(ds_qc_dat[,1])), function(y) filter(ds_qc_dat[,-1], ds_qc_dat$batch == unique(ds_qc_dat[,1])[y])) # list of subsets by batches
  ds_qc_dat_factor <- lapply(1:length(unique(ds_qc_dat[,1])), function(y) sapply(1:ncol(ds_qc_dat_subsets[[1]]), function(z) median(ds_qc_dat_subsets[[y]][,z], na.rm = T)/median(ds_qc_dat[,(z+1)], na.rm = T))) # calculate factor for every feature
  
  ds_batch_dat_subsets <- lapply(1:length(unique(ds_batch_dat[,1])), function(y) filter(ds_batch_dat[,-1], ds_batch_dat$batch == unique(ds_batch_dat[,1])[y]))
  
  ds_batch_dat_result <- lapply(1:length(unique(ds_batch_dat[,1])), function(y) sapply(1:ncol(ds_batch_dat_subsets[[1]]), function(z) ds_batch_dat_subsets[[y]][,z]/ds_qc_dat_factor[[y]][z])) # results by batch
  ds_batch_dat_result <- data.frame(do.call("rbind",ds_batch_dat_result)) # results in data frame
  rownames(ds_batch_dat_result) <- rownames(ds_batch_dat)
  colnames(ds_batch_dat_result) <- colnames(ds_batch_dat[,-c(1)])
  return(ds_batch_dat_result)
}

ds_qc_norm <- QC.Mnorm(qc_data = ds_qc_dat,data=ds_batch_dat)
ds_qc_norm <- round(as.data.frame(t(ds_qc_norm)),3)
ds_qc_norm <- merge(data_name,ds_qc_norm,by.x = "Protein",by.y = "row.names")
# save
fwrite(ds_qc_norm, paste0(args[3],"protein_after_batch_QC-Mnorm.csv"))
print("QC-Mnorm finish!!")

###############################################################################
#limma 
###############################################################################
library(limma)
# generate batch data
dat <- as.matrix(t(log2(ds[,-c(1:6)])))
batch <- as.numeric(ds$batch)
model <- model.matrix(~as.factor(ds$Type))
######################### Combat
ds_limma = removeBatchEffect(dat, batch=batch) #mod = model
ds_limma = round(as.data.frame(2^ds_limma),3)

ds_limma <-merge(data_name,ds_limma,by.x = "Protein",by.y = "row.names")
# save
fwrite(ds_limma, paste0(args[3],"protein_after_batch_Limma.csv"))
print("Limma finish!!")
################################################################
######################QC Ratio_G
################################################################
# #data
# batch <- factor(ds$batch)
# batches = levels(batch)
# nbatches = length(batches)
# ds_dat <- ds[,-c(1:6)]
# #qc data
# ds_qc <- subset(ds,Rep=="QC")
# qc_batch <- as.factor(ds_qc$batch)
# qc_batches = levels(qc_batch)
# qc_nbatches = length(qc_batches)
# qc_ds_dat <- ds_qc[,-c(1:6)]
# ######
# means = as.list(rep(0,nbatches))
# 
# for (i in 1:nbatches) {
#   means[[i]] <- apply(qc_ds_dat[qc_batch==qc_batches[i],], 2, function(x) exp(mean(log(x))))
#   ds_dat[batch==batches[i],] = scale(ds_dat[batch==batches[i],],center=rep(0,ncol(ds_dat)),scale=means[[i]])
# }
# ds_norm <- round(as.data.frame(t(ds_dat)),3)
# Protein <- rownames(ds_norm)
# ds_norm <-cbind(Protein,ds_norm)
