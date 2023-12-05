##### Script purpose: Identify gene whose expressions are correlated with 
## MYC expr. and whose alterations are associated with increased MYC expression
## in non-PMN-hit samples
##### Author: Ege Ulgen
##### Date: Dec 2023

library(DESeq2)
library(ggplot2)
library(maftools)
library(pathfindR)

PMN_genes <- read.csv("data/PMN_gene_list.csv")
excluded_genes <- c(PMN_genes$Gene.name, "IDH1", "IDH2", "TP53", "ATRX", "CDKN2A", "CDKN2B")

metadata_df <- readRDS("data/selected_data/meta.RDS")
metadata_df$PMN_hit <- as.character(metadata_df$PMN_hit)
metadata_df$PMN_hit <- sub("\\n.+", "", metadata_df$PMN_hit)
metadata_df$PMN_hit <- sub("PMN-", "", metadata_df$PMN_hit)
metadata_df$PMN_hit <- factor(metadata_df$PMN_hit, levels = c("no hit", "hit"))

# expression correlation with MYC -----------------------------------------
dds_norm <- readRDS("data/selected_data/processed_expr.RDS")
dds_norm$PMN_hit <- metadata_df$PMN_hit[match(dds_norm$patient, metadata_df$patient)]

# Retrieve the normalized data from the `DESeqDataSet` and Transpose this data
normalized_counts <- assay(dds_norm)
gene_syms <- rowData(dds_norm)[rownames(normalized_counts), "gene_name"]
normalized_counts <- t(limma::avereps(normalized_counts, ID = gene_syms))
norm_no_hit <- normalized_counts[metadata_df$patient[metadata_df$PMN_hit == "no hit"], ]

# Determine genes whose expression is correlated with MYC expression
all_other_genes <- setdiff(colnames(norm_no_hit), excluded_genes)
cor_res <- c()
for (gene in all_other_genes) {
    res <- cor.test(norm_no_hit[, "MYC"], norm_no_hit[, gene])
    cor_res <- rbind(cor_res,
                     data.frame(gene = gene,
                                cor = res$estimate,
                                p = res$p.value))
}
cor_res$adj_p <- p.adjust(cor_res$p, "fdr")
sum(cor_res$adj_p < 0.2 & abs(cor_res$cor) > .2)

# MYC expr. - som. alteration assoc ---------------------------------------
donor_no_PMN <- metadata_df[metadata_df$PMN_hit == "no hit", ]

myc_expr <- assay(dds_norm[rowData(dds_norm)$gene_name == "MYC", ])
myc_expr <- myc_expr[, match(donor_no_PMN$patient, colnames(myc_expr))]

scna_df <- readRDS("data/selected_data/CN_gene_level.RDS")
scna_df <- scna_df[abs(scna_df$Segment_Mean) >= .3, ]
scna_df <- scna_df[scna_df$patient_barcode %in% donor_no_PMN$patient, ]
scna_df$alt_type <- ifelse(scna_df$Segment_Mean >= .3, "amp", "del")

som_df <- readRDS("data/selected_data/maf.RDS")
som_df$Tumor_Sample_Barcode <- som_df$patient_barcode
table(som_df$Variant_Classification)
high_conseq <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                 "Translation_Start_Site","Nonsense_Mutation", 
                 "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", 
                 "Missense_Mutation")
som_df <- som_df[som_df$Variant_Classification %in% high_conseq, ]
som_df <- som_df[som_df$Tumor_Sample_Barcode %in% donor_no_PMN$patient, ]

cn_table <- scna_df[, c("symbol", "patient_barcode", "alt_type")]

clin_df <- donor_no_PMN
clin_df$Tumor_Sample_Barcode <- clin_df$patient

maf_df_no_pmn <- read.maf(maf = som_df,
                          cnTable = cn_table,
                          clinicalData = clin_df,
                          verbose = TRUE, isTCGA = TRUE)

gene_summary <- getGeneSummary(maf_df_no_pmn)

# keep only altered genes altered in >0.1 and < 0.9 of samples
min_samples <- round(nrow(clin_df) * 0.1)
all_genes <- gene_summary$Hugo_Symbol[gene_summary$AlteredSamples > min_samples & gene_summary$AlteredSamples < nrow(clin_df) - min_samples]
all_genes <- setdiff(all_genes, excluded_genes)

# determine altered genes associated with increased MYC expression
assoc_df <- c()
for (gene in all_genes) {
    
    tmp <- maf_df_no_pmn@data
    tmp <- tmp[tmp$Hugo_Symbol == gene, ]
    alt_seen <- paste(unique(tmp$Variant_Classification), collapse = ";")
    
    res <- wilcox.test(myc_expr[names(myc_expr) %in% tmp$Tumor_Sample_Barcode], 
                       myc_expr[!names(myc_expr) %in% tmp$Tumor_Sample_Barcode],
                       alternative = "greater")
    
    assoc_df <- rbind(assoc_df, data.frame(gene = gene,
                                           alt_seen = alt_seen,
                                           median_WT = median(myc_expr[!names(myc_expr) %in% tmp$Tumor_Sample_Barcode]),
                                           median_alt = median(myc_expr[names(myc_expr) %in% tmp$Tumor_Sample_Barcode]),
                                           n_WT = sum(!names(myc_expr) %in% tmp$Tumor_Sample_Barcode),
                                           n_alt = sum(names(myc_expr) %in% tmp$Tumor_Sample_Barcode),
                                           p = res$p.value))
}

assoc_df$median_diff <- assoc_df$median_alt - assoc_df$median_WT
assoc_df$adj_p <- p.adjust(assoc_df$p, "fdr")
sum(assoc_df$adj_p < 0.2)

# enrichment --------------------------------------------------------------
filtered_genes <- intersect(assoc_df$gene[assoc_df$adj_p < 0.2], cor_res$gene[cor_res$adj_p < 0.2 & abs(cor_res$cor) > .2])
input_df <- assoc_df[assoc_df$gene %in% filtered_genes, c("gene", "adj_p")]

res_pf_df <- run_pathfindR(input_df, output_dir = "output/PMN_neg_analysis", list_active_snw_genes = TRUE, p_val_threshold = 0.2)
saveRDS(res_pf_df, "output/PMN_neg_analysis/enr_result.RDS")

g <- enrichment_chart(res_pf_df)
ggsave("output/PMN_neg_analysis/enr_chart.pdf", g, width = 6, height = 9)

write.csv(filtered_genes, "output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv", row.names = FALSE)
