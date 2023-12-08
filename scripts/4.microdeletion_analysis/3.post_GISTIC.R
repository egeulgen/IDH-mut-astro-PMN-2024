##### Script purpose: Analyze chr19 peaks identified by GISTIC in microdel19q cases
##### Author: Ege Ulgen
##### Date: Dec 2023

library(DESeq2)
library(ggpubr)

metadata_df <- readRDS("data/selected_data/meta.RDS")
dds_norm <- readRDS("data/selected_data/processed_expr.RDS")

GISTIC_peaks <- read.delim("output/PMN_neg_analysis/all_GISTIC_res/microdel19q/all_lesions.conf_90.txt")
GISTIC_peaks$X <- NULL
GISTIC_peaks <- GISTIC_peaks[!grepl("CN values", GISTIC_peaks$Unique.Name), ]

peak_info <- GISTIC_peaks[, 1:9]
peak_info$short_name <- sub("e.* ", "",sub("m.* ", "", peak_info$Unique.Name))
GISTIC_peaks <- GISTIC_peaks[, -c(1:9)]
rownames(GISTIC_peaks) <- peak_info$short_name
GISTIC_peaks[grepl("D", rownames(GISTIC_peaks)), ] <- - GISTIC_peaks[grepl("D", rownames(GISTIC_peaks)), ]

rowSums(GISTIC_peaks != 0)
rowSums(GISTIC_peaks != 0)[rowSums(GISTIC_peaks != 0) > 5]

peak_info[peak_info$short_name == "D6", ]

del_peak_genes <- read.delim("output/PMN_neg_analysis/all_GISTIC_res/microdel19q/del_genes.conf_90.txt", skip = 3)
tmp <- del_peak_genes$chr19.54940787.57124931
tmp <- tmp[tmp != ""]
tmp


sel_expr <- assay(dds_norm[rowData(dds_norm)$gene_name %in% tmp, ])

anyDuplicated(rowData(dds_norm)$gene_name[match(rownames(sel_expr), rowData(dds_norm)$gene_id)])
rownames(sel_expr) <- rowData(dds_norm)$gene_name[match(rownames(sel_expr), rowData(dds_norm)$gene_id)]

pdf("output/PMN_neg_analysis/chr19.54940787.57124931_del_genes.pdf", width = 6, height = 6)
for (gene in rownames(sel_expr)) {
    boxplot(sel_expr[gene, colnames(sel_expr) %in% gsub("\\.", "-", colnames(GISTIC_peaks))],
            sel_expr[gene, !colnames(sel_expr) %in% gsub("\\.", "-", colnames(GISTIC_peaks))],
            main = gene, names = c("19q13.43 del", "intact"))
}
dev.off()

