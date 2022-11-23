#setwd("/home/boot/Script/test/batch_pipeline/batch_effect_correction/NC/Correction/long_format")
library(iq)
library(dplyr)
library(data.table)
args = commandArgs(T)
proteinGroup <- as.data.frame(fread(args[1],check.names = F))
# proteinGroup <- as.data.frame(fread("peptides_after_MVI_EigenMS.csv",check.names = F,quote = ""))
proteinGroup$PEP.Quantity <- as.numeric(proteinGroup$PEP.Quantity)
primary_id = "PG.ProteinAccessions"
sample_id  <- "R.FileName" 
intensity_col = "PEP.Quantity"
secondary_id <- c("EG.PrecursorId","PG.FASTAName")
annotation_columns <- "PG.FASTAName"

#==========================
#MaxLFQ
#==========================
norm_data <- iq::preprocess(proteinGroup,
                            primary_id = primary_id,
                            sample_id  = sample_id,
                            intensity_col = intensity_col,
                            secondary_id = secondary_id,
                            pdf_out = NULL,
                            median_normalization = FALSE)
#============fast============
#fast_MaxLFQ = iq::create_protein_list + iq::create_protein_table
result_faster <- iq::fast_MaxLFQ(norm_data)
#annotation
extra_names <- iq::extract_annotation(rownames(result_faster$estimate),
                                      proteinGroup,
                                      primary_id = "PG.ProteinAccessions",
                                      annotation_columns = annotation_columns)
#save
write.csv(cbind(Protein = rownames(result_faster$estimate),
                PG.FASTAName=extra_names[, annotation_columns],
                2^result_faster$estimate),quote = FALSE,
                args[2], row.names = FALSE)

# #==================
# #top3
# #==================
# 
# result_top3 <- iq::create_protein_table(protein_list,method = "topN")
# 
# extra_names <- iq::extract_annotation(rownames(result_top3$estimate), 
#                                       proteinGroup, 
#                                       annotation_columns = annotation_columns)
# 
# write.table(cbind(Protein = rownames(result_top3$estimate),
#                   extra_names[, annotation_columns],
#                   2^result_top3$estimate),quote = FALSE,
#             "Protein-top3.txt", sep = "\t", row.names = FALSE)
# # result_top3 <- iq::create_protein_table(protein_list,method = "topN")
# # 
# # extra_names <- iq::extract_annotation(rownames(result_top3$estimate), 
# #                                       proteinGroup, 
# #                                       annotation_columns = annotation_columns)
# # 
# # write.table(cbind(Protein = rownames(result_top3$estimate),
# #                   extra_names[, annotation_columns],
# #                   2^result_top3$estimate),quote = FALSE,
# #             "Protein-top3.txt", sep = "\t", row.names = FALSE)
