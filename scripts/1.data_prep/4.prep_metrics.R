##### Script purpose: Calculate genomic instability metrics
##### Author: Ege Ulgen 
##### Date: Dec 2023

source("scripts/utils.R")

EXOME_LEN <- 44
cytoband_df <- read.delim("data/hg38_cytoBand.txt", header = FALSE)

metadata_df <- readRDS("data/selected_data/meta.RDS")

# SNV/indel burden --------------------------------------------------------
som_df <- readRDS("data/selected_data/maf.RDS")

metadata_df$indel_burden <- metadata_df$snv_burden <- NA
for (donor_id in metadata_df$patient) {
    donor_som_df <- som_df[som_df$patient_barcode == donor_id, ]
    # Filter for depth greater than 20X in tumor samples & greater than 10X in normal samples
    donor_som_df <- donor_som_df[donor_som_df$t_depth >= 20 & donor_som_df$n_depth >= 10, ]

    # subset for tumor_f > 0.05
    donor_som_df$VAF <- donor_som_df$t_alt_count / donor_som_df$t_depth
    donor_som_df <- donor_som_df[donor_som_df$VAF > 0.05, ]
    
    metadata_df$snv_burden[metadata_df$patient == donor_id] <- sum(donor_som_df$Variant_Type == "SNP") / EXOME_LEN
    metadata_df$indel_burden[metadata_df$patient == donor_id] <- sum(donor_som_df$Variant_Type != "SNP") / EXOME_LEN
}

# CNA burden --------------------------------------------------------------
SCNA_df <- readRDS("data/selected_data/CN_segments.RDS")

hist(SCNA_df$Segment_Mean, breaks = 100)
abline(v = .3, col = "red")
abline(v = -.3, col = "red")
SCNA_df <- SCNA_df[abs(SCNA_df$Segment_Mean) > .3, ]

all_cases <- unique(SCNA_df$patient_barcode)
chr_list <- 1:22

chr_sizes_hg38 <- tapply(cytoband_df$V3, cytoband_df$V1, function(x) max(x) + 1)[paste0("chr", chr_list)]

tcga_wgii_df <- c()
for (case in all_cases) {
    segs_df <- SCNA_df[SCNA_df$patient_barcode == case, ]
    segs_df$width <- segs_df$End - segs_df$Start + 1
    
    # WGII
    ratio_by_chr <- c()
    for (chrm in chr_list) {
        seg_chr_df <- segs_df[segs_df$Chromosome == chrm, ]
        chr_width <- chr_sizes_hg38[chrm]
        ratio_by_chr <- c(ratio_by_chr, sum(seg_chr_df$width) / chr_width)
    }
    wgii <- mean(ratio_by_chr)
    
    tcga_wgii_df <- rbind(tcga_wgii_df,
                          data.frame(case_ID = case,
                                     WGII = wgii))
}

metadata_df$WGII <- tcga_wgii_df$WGII[match(metadata_df$patient, tcga_wgii_df$case_ID)]
metadata_df$WGII[is.na(metadata_df$WGII)] <- 0

# CAER --------------------------------------------------------------------
SCNA_df$chr <- paste0("chr", SCNA_df$Chromosome)
SCNA_df$Ratio <- 2 ^ SCNA_df$Segment_Mean

chr_arm_df <- determine_chr_arm_level_CN_ratios(SCNA_df, cytoband_df,  paste0("chr", 1:22))
chr_arm_df$log_ratio <- log2(chr_arm_df$avg_ratio)

hist(chr_arm_df$log_ratio, breaks = 100)
abline(v = -.3, col = "red")
abline(v = .3, col = "red")

# 1 for amp, 0 for neutral, -1 for del
chr_arm_df$SCNA_event <- ifelse(chr_arm_df$log_ratio < -.3, -1,
                                ifelse(chr_arm_df$log_ratio > .3, 1, 0))
write.csv(chr_arm_df, "output/chr_arm_scna.csv", row.names = FALSE)

binary_mat <- reshape2::acast(chr_arm_df, Case_ID ~ arm, value.var = "SCNA_event")
binary_mat <- ifelse(binary_mat != 0 , 1, 0)
chr_arm_event_ratio <- rowSums(binary_mat) / 44

hist(chr_arm_event_ratio)

metadata_df$chr_arm_event_ratio <- chr_arm_event_ratio[match(metadata_df$patient, substr(names(chr_arm_event_ratio), 1, 12))]
metadata_df$chr_arm_event_ratio[is.na(metadata_df$chr_arm_event_ratio)] <- 0

# CN amplitude ------------------------------------------------------------
max_cn_df <- data.frame()
for (case_ID in all_cases) {
    segs_df <- SCNA_df[SCNA_df$patient_barcode == case_ID, ]
    
    max_cn_df <- rbind(max_cn_df,
                       data.frame(Case_ID = case_ID,
                                  max_CN = round(2 ^ max(segs_df$Segment_Mean) * 2)))
}

metadata_df$maxCN <- max_cn_df$max_CN[match(metadata_df$patient, max_cn_df$Case_ID)]
metadata_df$maxCN[is.na(metadata_df$maxCN)] <- 2


saveRDS(metadata_df, "data/selected_data/meta.RDS")
