##### Script purpose: Analyze chr19 peaks identified by GISTIC (on GenePattern) in microdel19q cases
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
cat(tmp, file = "output/PMN_neg_analysis/d6_genes.txt", sep = "\n")


sel_expr <- assay(dds_norm[rowData(dds_norm)$gene_name %in% tmp, ])

anyDuplicated(rowData(dds_norm)$gene_name[match(rownames(sel_expr), rowData(dds_norm)$gene_id)])
rownames(sel_expr) <- rowData(dds_norm)$gene_name[match(rownames(sel_expr), rowData(dds_norm)$gene_id)]

sel_expr_long <- reshape2::melt(sel_expr)
colnames(sel_expr_long) <- c("Gene", "Sample", "MYC expression")
sel_expr_long$microdel <- ifelse(sel_expr_long$Sample %in% gsub("\\.", "-", colnames(GISTIC_peaks)), "19q13.43 del", "intact")

g <- ggviolin(sel_expr_long, x = "microdel", y = "MYC expression", color = "microdel", 
              add = "boxplot", facet.by = "Gene", palette = "npg")
g <- g + stat_compare_means(method = "wilcox.test", method.args = list("alternative" = "less"))
g <- g + rremove("legend")
g <- g + rremove("xlab")

ggsave("output/PMN_neg_analysis/chr19.54940787.57124931_del_genes.pdf", g, width = 15, height = 15)
