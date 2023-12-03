##### Script purpose: Select primary grade 2-3 IDH-mutant Non-codel samples
##### Author: Ege Ulgen 
##### Date: Dec 2023

library(survival)
library(survminer)
library(ggplot2)


meta_df <- readRDS("data/all_primary_complete/meta_data.RDS")
gene_scna <- readRDS("data/all_primary_complete/CN_gene_level.RDS")
expr_raw_counts <- readRDS("data/all_primary_complete/expr_raw_counts.RDS")
maf_df <- readRDS("data/all_primary_complete/maf.RDS")
SCNA_segs <- readRDS("data/all_primary_complete/CN_segments.RDS")

# finalize selection ------------------------------------------------------
patients_df <- readRDS("data/all_LGG_patients_df.RDS")
selected_meta_df <- meta_df[meta_df$final_subtype == "mut--non-codel", ]

final_meta_df <- merge(selected_meta_df, patients_df, by.x = "patient", by.y = "bcr_patient_barcode")
final_meta_df <- final_meta_df[final_meta_df$tumor_grade.y != "[Discrepancy]", ]
table(final_meta_df$tumor_grade.y)

## check for CDKN2A/2B hom. deletion
tmp <- gene_scna[gene_scna$symbol %in% c("CDKN2A", "CDKN2B"), ]
cdkn_hom_loss <- unique(tmp$patient_barcode[tmp$Segment_Mean < -1])
final_meta_df$CDKN2A_2B_hom_loss <- final_meta_df$patient %in% cdkn_hom_loss

final_meta_df$mol_grade <- ifelse(final_meta_df$CDKN2A_2B_hom_loss, 4, as.numeric(sub("G", "", final_meta_df$tumor_grade.y)))

table(final_meta_df$mol_grade)

selected_scna_genes <- gene_scna[gene_scna$patient_barcode %in% final_meta_df$patient, ]
selected_raw.counts <- subset(expr_raw_counts, select = expr_raw_counts$patient %in% final_meta_df$patient)

final_meta_df$death_stat <- final_meta_df$vital_status.y == "Dead"
final_meta_df$OS_time <- as.numeric(ifelse(final_meta_df$death_days_to != "[Not Applicable]", final_meta_df$death_days_to, final_meta_df$last_contact_days_to))

# add TCGA - CDR survival data
# https://xenabrowser.net/datapages/?dataset=survival%2FLGG_survival.txt&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
TCGA_CDR_df <- read.delim("data/survival_LGG_survival.txt")


# using Accurately Defined - https://www.cell.com/action/showFullTableHTML?isHtml=true&tableId=tbl3&pii=S0092-8674%2818%2930229-0
idx <- match(final_meta_df$patient, TCGA_CDR_df$X_PATIENT)
final_meta_df$CDR_OS <- TCGA_CDR_df$OS[idx]
final_meta_df$CDR_OS_time <- TCGA_CDR_df$OS.time[idx]
final_meta_df$CDR_PFI <- TCGA_CDR_df$PFI[idx]
final_meta_df$CDR_PFI_time <- TCGA_CDR_df$PFI.time[idx]
final_meta_df$CDR_DFI <- TCGA_CDR_df$DFI[idx]
final_meta_df$CDR_DFI_time <- TCGA_CDR_df$DFI.time[idx]


dir.create("data/selected_data")
saveRDS(final_meta_df, "data/selected_data/meta.RDS")
saveRDS(maf_df[maf_df$patient_barcode %in% final_meta_df$patient, ], "data/selected_data/maf.RDS")
saveRDS(SCNA_segs[SCNA_segs$patient_barcode %in% final_meta_df$patient, ], "data/selected_data/CN_segments.RDS")
saveRDS(selected_scna_genes, "data/selected_data/CN_gene_level.RDS")
saveRDS(selected_raw.counts, "data/selected_data/expr_raw_counts.RDS")

# KM by grade -------------------------------------------------------------
plot_KM_by_grade <- function(type) {
    mol_grade_cols <- c("2" = "#FBA465", "3" = "#F86E51", "4" = "#D1193E")
    
    tmp <- final_meta_df
    tmp$Time <- final_meta_df[, paste0(type, "_time")]
    tmp$Stat <- final_meta_df[, type]
    
    fit_grade <- survfit(Surv(Time, Stat) ~ mol_grade, data = tmp)
    g <- ggsurvplot(fit_grade, 
                    data = tmp, 
                    size = 1,
                    conf.int = FALSE,
                    pval = TRUE,
                    legend.labs = c("2", "3", "4"),
                    palette = mol_grade_cols,
                    linetype = "strata",
                    xlab = "Time",
                    newpage = TRUE)
    g <- g$plot
    g <- g + theme_minimal()
    g <- g + theme(title = element_text(face = "bold", size = 12),
                   axis.text = element_text(size = 11),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.title.x = element_blank(),
                   legend.title = element_blank())
    g
}

g_grades <- plot_KM_by_grade("CDR_OS")

ggsave("output/KM_by_grade.pdf", g_grades, width = 6, height = 6)
