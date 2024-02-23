##### Script purpose: Obtain and filter for primary LGG samples with complete data (mutation + SCNA + RNAseq) and classify acc. to IDH/codel
##### Author: Ege Ulgen
##### Date: 10 Dec 2023

library(TCGAbiolinks)
library(SummarizedExperiment)
import::from("scripts/utils.R", determine_chr_arm_level_CN_ratios, produce_genome_SCNA_plot)

dir.create("output")
cytoband_df <- read.delim("data/hg38_cytoBand.txt", header = FALSE)

# clinical ----------------------------------------------------------------
query.clin <- GDCquery(project = "TCGA-LGG", data.category = "Clinical", data.type = "Clinical Supplement", data.format = "BCR Biotab")
GDCdownload(query.clin, directory = "data")

# mutation ----------------------------------------------------------------
query.mut <- GDCquery(project = "TCGA-LGG", data.category = "Simple Nucleotide Variation", access = "open",
    data.type = "Masked Somatic Mutation", workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query.mut, directory = "data")

# SCNA --------------------------------------------------------------------
query.SCNA.seg <- GDCquery(project = "TCGA-LGG", data.category = "Copy Number Variation", data.type = "Masked Copy Number Segment",
    access = "open")
GDCdownload(query.SCNA.seg, directory = "data")

# RNAseq ------------------------------------------------------------------
query.exp <- GDCquery(project = "TCGA-LGG", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts", sample.type = "Primary Tumor")
GDCdownload(query.exp, directory = "data")

# prepare raw data --------------------------------------------------------
# clinical data
clinical.BCRtab.all <- GDCprepare(query.clin, directory = "data")
patients_df <- clinical.BCRtab.all$clinical_patient_lgg[-c(1:2), ]
saveRDS(patients_df, "data/all_LGG_patients_df.RDS")

# somatic mutation
maf <- GDCprepare(query.mut, directory = "data")
maf$patient_barcode <- substr(maf$Tumor_Sample_Barcode, 1, 12)
table(substr(maf$Tumor_Sample_Barcode, 14, 16))
prim_maf <- maf[substr(maf$Tumor_Sample_Barcode, 14, 15) == "01", ]

# SCNA
scna_segs <- GDCprepare(query.SCNA.seg, directory = "data")
scna_segs$patient_barcode <- substr(scna_segs$Sample, 1, 12)
prim_scna_segs <- scna_segs[substr(scna_segs$Sample, 14, 15) == "01", ]

# RNAseq
raw.counts <- GDCprepare(query = query.exp, directory = "data")
prim_raw.counts <- subset(raw.counts, select = raw.counts$shortLetterCode == "TP")

# determine IDHmut --------------------------------------------------------
### obtain IDH-mutant patient IDs
idh_mut_df <- prim_maf[prim_maf$Hugo_Symbol %in% c("IDH1", "IDH2"), ]
table(idh_mut_df$callers)   # all supported by >1 caller
table(idh_mut_df$Hugo_Symbol, idh_mut_df$HGVSp_Short)  # all have valid HGVSp consequence
table(idh_mut_df$dbSNP_RS, idh_mut_df$Variant_Classification, exclude = FALSE)  # all have non-syn. mut.s

prim_idh_mut_patients <- unique(substr(idh_mut_df$Tumor_Sample_Barcode, 1, 12))

codel_check_df <- prim_scna_segs[prim_scna_segs$Chromosome %in% c("1", "19"), ]

# determine codel. --------------------------------------------------------
# process subsetted segments df
codel_check_df <- as.data.frame(codel_check_df)
codel_check_df$chr <- paste0("chr", codel_check_df$Chromosome)
codel_check_df$Ratio <- 2^codel_check_df$Segment_Mean ## because it's log2(R)

chr_arm_df <- determine_chr_arm_level_CN_ratios(codel_check_df, cytoband_df)
chr_arm_df <- chr_arm_df[chr_arm_df$arm %in% c("chr1p", "chr19q"), ]

codel_df <- data.frame(Case_ID = unique(prim_scna_segs$Sample), codel_stat = "non-codel")
for (case in codel_df$Case_ID) {
    tmp <- chr_arm_df[chr_arm_df$Case_ID == case, ]
    if (all(tmp$avg_ratio < 1)) {
        codel_df$codel_stat[codel_df$Case_ID == case] <- "codel"
    }
}
table(codel_df$codel_stat)

prim_codel_patients <- unique(substr(codel_df$Case_ID[codel_df$codel_stat == "codel"], 1, 12))

# Subset for primary samples with complete data ---------------------------
all_prim_patients <- intersect(unique(prim_maf$patient_barcode), intersect(unique(prim_scna_segs$patient_barcode), unique(prim_raw.counts$patient)))

prim_maf_final <- prim_maf[prim_maf$patient_barcode %in% all_prim_patients, ]
prim_scna_segs_final <- prim_scna_segs[prim_scna_segs$patient_barcode %in% all_prim_patients, ]
prim_raw.counts_final <- subset(prim_raw.counts, select = prim_raw.counts$patient %in% all_prim_patients)

meta_df <- colData(prim_raw.counts_final)

meta_df$IDH_mut_WXS <- ifelse(meta_df$patient %in% prim_idh_mut_patients, "mut", "WT")
meta_df$Codel_GDC <- ifelse(meta_df$patient %in% prim_codel_patients, "codel", "non-codel")
meta_df$manual_subtype <- ifelse(meta_df$IDH_mut_WXS == "WT", "IDH-WT", paste0(meta_df$IDH_mut_WXS, "--", meta_df$Codel_GDC))

# Final class -------------------------------------------------------------
# we manually check IDH-mut and 1p/19q codel. status to confirm
table(meta_df$manual_subtype, meta_df$paper_IDH.codel.subtype)

meta_df$final_subtype <- meta_df$manual_subtype

# IDH-WT in paper -- has INS in IDH2s
meta_df$patient[meta_df$paper_IDH.codel.subtype == "IDHwt" & meta_df$manual_subtype == "mut--non-codel"]
tmp <- prim_maf_final[prim_maf_final$patient_barcode == "TCGA-P5-A72U", ]
tmp[tmp$Hugo_Symbol %in% c("IDH1", "IDH2"), ]

# mutant in paper -- will assume IDHmut (published as such)
meta_df[meta_df$paper_IDH.codel.subtype == "IDHmut-non-codel" & meta_df$manual_subtype == "IDH-WT", ]
produce_genome_SCNA_plot("TCGA-WY-A85A", prim_scna_segs_final, cytoband_df) # visual inspection reveals non-codel
meta_df$final_subtype[meta_df$paper_IDH.codel.subtype == "IDHmut-non-codel" & meta_df$manual_subtype == "IDH-WT"] <- "mut--non-codel"

# non-codel in paper / current subtype - codel
codel_mismatch_df <- meta_df[meta_df$paper_IDH.codel.subtype == "IDHmut-non-codel" & meta_df$manual_subtype == "mut--codel", ]

pdf("output/mistaken_codels_for_manual_inspection.pdf", width = 25, height = 6)
for (patient_id in codel_mismatch_df$patient) {
    produce_genome_SCNA_plot(patient_id, prim_scna_segs_final, cytoband_df)
}
dev.off()

true_codels <- c("TCGA-HT-8113", "TCGA-S9-A6U8", "TCGA-FG-7637", "TCGA-VM-A8CA",
                 "TCGA-E1-A7YY", "TCGA-DU-A7TG", "TCGA-CS-5394", "TCGA-HT-7602",
                 "TCGA-FG-8189", "TCGA-S9-A6WI", "TCGA-S9-A6WQ", "TCGA-WY-A85E")

non_codels <- setdiff(codel_mismatch_df$patient, true_codels)
meta_df$final_subtype[meta_df$patient %in% non_codels] <- "mut--non-codel"

table(meta_df$final_subtype, meta_df$paper_IDH.codel.subtype, exclude = FALSE)

# save data ---------------------------------------------------------------
dir.create("data/all_primary_complete")
saveRDS(meta_df, "data/all_primary_complete/meta_data.RDS")
saveRDS(prim_maf_final, "data/all_primary_complete/maf.RDS")
saveRDS(prim_scna_segs_final, "data/all_primary_complete/CN_segments.RDS")
saveRDS(prim_raw.counts, "data/all_primary_complete/expr_raw_counts.RDS")
