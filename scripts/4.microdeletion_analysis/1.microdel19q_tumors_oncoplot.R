##### Script purpose: Oncoplot of PMN-negative samples with microdeletions in 19q
##### Author: Ege Ulgen
##### Date: Dec 2023

library(maftools)
library(TCGAbiolinks)

annotated_df <- read.csv("output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv")
microdel_19q_pts <- unique(unlist(strsplit(annotated_df$altered_ids, ";")))

metadata_df <- readRDS("data/selected_data/meta.RDS")

# oncoplot ----------------------------------------------------------------
som_df <- readRDS("data/selected_data/maf.RDS")
som_df$Tumor_Sample_Barcode <- som_df$patient_barcode
table(som_df$Variant_Classification)
high_conseq <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                 "Translation_Start_Site","Nonsense_Mutation", 
                 "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", 
                 "Missense_Mutation")
som_df <- som_df[som_df$Variant_Classification %in% high_conseq, ]
som_df <- som_df[som_df$Tumor_Sample_Barcode %in% microdel_19q_pts, ]

scna_df <- readRDS("data/selected_data/CN_gene_level.RDS")
scna_df <- scna_df[abs(scna_df$Segment_Mean) >= .3, ]
scna_df <- scna_df[scna_df$patient_barcode %in% microdel_19q_pts, ]
scna_df$alt_type <- ifelse(scna_df$Segment_Mean >= .3, "amp", "del")
cn_table <- scna_df[, c("symbol", "patient_barcode", "alt_type")]

clin_df <- metadata_df[metadata_df$patient %in% microdel_19q_pts, ]
clin_df$Tumor_Sample_Barcode <- clin_df$patient

maf_df_del19q <- read.maf(maf = som_df,
                          cnTable = cn_table,
                          clinicalData = clin_df,
                          verbose = TRUE, isTCGA = TRUE)

pdf("output/23.microdel19q_oncoplot.pdf", width = 10, height = 8)
oncoplot(maf_df_del19q, top = 25)
dev.off()
