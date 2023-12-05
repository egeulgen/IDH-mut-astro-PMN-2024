##### Script purpose: Investigate MYC expression in 19q deleted vs. intact samples
##### Author: Ege Ulgen
##### Date: Dec 2023

library(ggpubr)

metadata_df <- readRDS("data/selected_data/meta.RDS")
dds_norm <- readRDS("data/selected_data/processed_expr.RDS")
myc_expr <- assay(dds_norm[rowData(dds_norm)$gene_name == "MYC", ])

metadata_df$MYC_expr <- myc_expr[match(metadata_df$patient, colnames(myc_expr))]

all_arm_scna_df <- read.csv("output/chr_arm_scna.csv")

chr19q_del_samples <- all_arm_scna_df$Case_ID[all_arm_scna_df$arm == "chr19q" & all_arm_scna_df$SCNA_event == -1]
chr19q_del_samples <- substr(chr19q_del_samples, 1, 12)
metadata_df$del_19q <- ifelse(metadata_df$patient %in% chr19q_del_samples, "del_19q", "WT")

g <- ggviolin(metadata_df, x = "del_19q", y = "MYC_expr", color = "del_19q", palette = "lancet", add = "jitter", xlab = "")
g <- g + stat_compare_means()
g <- g + rremove("legend")
g

g2 <- ggviolin(metadata_df, x = "del_19q", y = "MYC_expr", color = "del_19q", palette = "lancet", add = "jitter", xlab = "",
              facet.by = "PMN_hit")
g2 <- g2 + stat_compare_means()
g2 <- g2 + rremove("legend")
g2

ggsave("output/22.MYC_expr_by_19q_del.pdf", ggarrange(g, g2, nrow = 2), width = 5, height = 8)
 