##### Script purpose: Annotate altered genes associated with increased MYC expression in PMN-negative samples
##### Author: Ege Ulgen
##### Date: Dec 2023

library(biomaRt)
library(ggpubr)

# Annotate localization ---------------------------------------------------
assoc_genes <- read.csv("output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv")
assoc_genes <- assoc_genes[, 1]

# annotate localization of genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
annotated_df <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "band"),
                      filters = "hgnc_symbol",
                      values = assoc_genes,
                      mart = ensembl)
annotated_df <- annotated_df[!grepl("CHR", annotated_df$chromosome_name), ]
all(assoc_genes %in% annotated_df$hgnc_symbol)

annotated_df$chromosome_name <- paste0("chr", annotated_df$chromosome_name)
annotated_df$cytoband <- paste0(annotated_df$chromosome_name, annotated_df$band)
annotated_df$arm <- sub("\\d+(\\.\\d+)*$", "", annotated_df$cytoband)

annotated_df <- annotated_df[, c("hgnc_symbol", "chromosome_name", "arm", "cytoband")]

table(annotated_df$arm)
table(annotated_df$cytoband)

write.csv(annotated_df, "output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv", row.names = FALSE)

# Annotate observed event frequencies -------------------------------------
metadata_df <- readRDS("data/selected_data/meta.RDS")
dds_norm <- readRDS("data/selected_data/processed_expr.RDS")

metadata_df$PMN_hit <- as.character(metadata_df$PMN_hit)
metadata_df$PMN_hit <- sub("\\n.+", "", metadata_df$PMN_hit)
metadata_df$PMN_hit <- sub("PMN-", "", metadata_df$PMN_hit)
metadata_df$PMN_hit <- factor(metadata_df$PMN_hit, levels = c("no hit", "hit"))


donor_no_PMN <- metadata_df[metadata_df$PMN_hit == "no hit", ]

scna_df <- readRDS("data/selected_data/CN_gene_level.RDS")
scna_df <- scna_df[abs(scna_df$Segment_Mean) >= .3, ]
scna_df <- scna_df[scna_df$patient_barcode %in% donor_no_PMN$patient, ]
scna_df$alt_type <- ifelse(scna_df$Segment_Mean >= .3, "amp", "del")

som_df <- readRDS("data/selected_data/maf.RDS")
som_df$Tumor_Sample_Barcode <- som_df$patient_barcode
high_conseq <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                 "Translation_Start_Site","Nonsense_Mutation", 
                 "Nonstop_Mutation", "In_Frame_Del","In_Frame_Ins", 
                 "Missense_Mutation")
som_df <- som_df[som_df$Variant_Classification %in% high_conseq, ]
som_df <- som_df[som_df$Tumor_Sample_Barcode %in% donor_no_PMN$patient, ]

annotated_df$n_amp <- annotated_df$n_del <- annotated_df$n_mut <- NA
annotated_df$altered_ids <- NA
for (gene in annotated_df$hgnc_symbol) {
    
    mut_ids <- unique(som_df$Tumor_Sample_Barcode[som_df$Hugo_Symbol == gene])
    amp_ids <- unique(scna_df$patient_barcode[scna_df$symbol == gene & scna_df$alt_type == "amp"])
    del_ids <- unique(scna_df$patient_barcode[scna_df$symbol == gene & scna_df$alt_type == "del"])
    
    annotated_df$n_mut[annotated_df$hgnc_symbol == gene] <- length(mut_ids)
    annotated_df$n_amp[annotated_df$hgnc_symbol == gene] <- length(amp_ids)
    annotated_df$n_del[annotated_df$hgnc_symbol == gene] <- length(del_ids)
    
    if (length(del_ids) > length(mut_ids) & length(del_ids) > length(amp_ids)) {
        sel_ids <- del_ids
    } else if (length(amp_ids) > length(mut_ids)) {
        sel_ids <- amp_ids
    } else {
        sel_ids <- mut_ids 
    }
    
    annotated_df$altered_ids[annotated_df$hgnc_symbol == gene] <- paste(sel_ids, collapse = ";")
    
}

# deletion is the most frequent alt. for all genes
table(apply(annotated_df, 1, function(x) which.max(x[5:7])))

write.csv(annotated_df, "output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv", row.names = FALSE)


altered_patients <- unique(unlist(strsplit(annotated_df$altered_ids, ";")))
table(metadata_df$mol_grade[match(altered_patients, metadata_df$patient)])
