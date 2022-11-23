#setwd("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC")
args <- c("/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/maxLFQ/Protein_raw_maxLFQ_impute0.4_T.csv",
          "/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/group.xls",
          "/home/boot/Script/test/batch_pipeline/NC_batch_pipeline2/Peptide/maxLFQ/")

library(data.table)
library(dplyr)
library(stringr)
args = commandArgs(T)
##############################################################################################################################################################
# Metadata generating
##############################################################################################################################################################
dsr <-as.data.frame((fread(input = args[1], header=T,check.names = F,quote = ""))) # 第一列为样本名,列名为肽段,行为样本
rownames(dsr) <- dsr[,1]
dsr <- dsr[,-1]
########## METADATA GENERATING
# dataset with all metadata
group <- readxl::read_xls(args[2]) #xls
colnames(group)[1:3]=c("ID","Type","injection")
subid <- data.frame(ID=rownames(dsr))
all_meta_data <- merge(group,subid)
ds <- data.frame(cbind(all_meta_data, dsr),check.names = F)[,-1] # 1-2 => drop ID
save(ds,file = paste0(args[3],"ds.RData"))
rm(list = ls()[ls()!="ds"])
