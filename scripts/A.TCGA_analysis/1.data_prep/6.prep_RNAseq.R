##### Script purpose: Prepare RNAseq data
##### Author: Ege Ulgen 
##### Date: Dec 2023

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)

metadata_df <- readRDS("data/selected_data/meta.RDS")
raw_counts <- readRDS("data/selected_data/expr_raw_counts.RDS")
raw_exp_mat <- assays(raw_counts)[[1]]  # unstranded

# keep only genes with 10 or more reads in at least 25% of all samples
raw_exp_mat <- raw_exp_mat[rowSums(raw_exp_mat >= 10) >= ncol(raw_exp_mat) *.25, ]
colnames(raw_exp_mat) <- substr(colnames(raw_exp_mat), 1, 12)

metadata <- metadata_df[match(colnames(raw_exp_mat), metadata_df$patient), ]
rownames(metadata) <- metadata$patient

dds <- DESeqDataSetFromMatrix(countData = raw_exp_mat,
                              colData = metadata, 
                              rowData = rowData(raw_counts)[rownames(raw_exp_mat), ],
                              design = ~ 1)
dds

# Normalize and transform the data in the `DESeqDataSet` object using the `vst()`
# function from the `DESEq2` R package
dds_norm <- vst(dds)

saveRDS(dds_norm, "data/selected_data/processed_expr.RDS")
