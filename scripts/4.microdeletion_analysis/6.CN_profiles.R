##### Script purpose: Prepare input for and run SigProfilerMatrixGenerator
##### Author: Ege Ulgen
##### Date: Dec 2023

library(TCGAbiolinks)

metadata_df <- readRDS("data/selected_data/meta.RDS")

query.SCNA.seg_allele_spec <- GDCquery(project = "TCGA-LGG", data.category = "Copy Number Variation", data.type = "Allele-specific Copy Number Segment",
                                       access = "open", barcode = substr(metadata_df$barcode, 1, 16), workflow.type = "ASCAT2")
GDCdownload(query.SCNA.seg_allele_spec, directory = "data")

ASCAT_df <- GDCprepare(query.SCNA.seg_allele_spec, directory = "data")


final_ASCAT_df <- ASCAT_df[, c("Sample", "Chromosome", "Start", "End", "Major_Copy_Number", "Minor_Copy_Number")]
colnames(final_ASCAT_df) <- c("sample", "chromosome", "startpos", "endpos","nMajor", "nMinor")
final_ASCAT_df$sample <- substr(final_ASCAT_df$sample, 1, 12)

write.table(final_ASCAT_df, "data/TCGA_ASCAT.tsv", row.names = FALSE, quote = FALSE, sep = "\t")

# poetry run python scripts/4.microdeletion_analysis/7.CNV_matrix_generation.py