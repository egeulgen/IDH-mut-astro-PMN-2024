##### Script purpose: identify MYC exon 3 hypomethylation and mir34b/c promoter hypermethylation status
##### Author: Ege Ulgen
##### Date: Dec 2023

library(TCGAbiolinks)
library(ggpubr)

source("scripts/utils.R")

# read data ---------------------------------------------------------------
metadata_df <- readRDS("data/selected_data/meta.RDS")

query_met.hg38 <- GDCquery(project = "TCGA-LGG", 
                           data.category = "DNA Methylation", 
                           platform = "Illumina Human Methylation 450",
                           data.type = "Methylation Beta Value",
                           barcode = substr(metadata_df$barcode, 1, 16))
GDCdownload(query_met.hg38, directory = "data")
methyl_data <- GDCprepare(query_met.hg38, directory = "data", summarizedExperiment = FALSE)

methyl_data <- methyl_data[rownames(methyl_data) %in% c("cg00163372", "cg22879515"), ]
colnames(methyl_data) <- substr(colnames(methyl_data), 1, 12)

all(metadata_df$patient %in% colnames(methyl_data))


# MYC exon 3 hypomethylation ----------------------------------------------
metadata_df$MYC_exon3_hypomet <- metadata_df$patient %in% names(which(methyl_data["cg00163372", ] < 0.5))

# mir34b/c TSS hypermethylation -------------------------------------------
methyl_data_mir34bc <- methyl_data["cg22879515", ]

thr_val <- mean(methyl_data_mir34bc) + 2 * sd(methyl_data_mir34bc)

metadata_df$mir34bc_hypermet <- metadata_df$patient %in% names(which(methyl_data_mir34bc > thr_val))

# assoc? ------------------------------------------------------------------
table(metadata_df$MYC_exon3_hypomet, metadata_df$mir34bc_hypermet)


table(metadata_df$PMN_hit, metadata_df$MYC_exon3_hypomet)
table(metadata_df$PMN_hit, metadata_df$mir34bc_hypermet)


g1 <- plot_cat(df = metadata_df,
               var1 = "PMN_hit", var2 = "MYC_exon3_hypomet",
               val = TRUE,
               y_lbl = "MYC exon3 hypometylation",
               y_lims = c(0, .2))


g2 <- plot_cat(df = metadata_df,
               var1 = "PMN_hit", var2 = "mir34bc_hypermet",
               val = TRUE,
               y_lbl = "miR-34b/c TSS hypermetylation",
               y_lims = c(0, .2))

# effect on MYC expr? -----------------------------------------------------
dds_norm <- readRDS("data/selected_data/processed_expr.RDS")

MYC_expr <- assay(dds_norm[rowData(dds_norm)$gene_name == "MYC", ])

metadata_df$MYC_expr <- MYC_expr[match(metadata_df$patient, colnames(MYC_expr))]

g3 <- ggviolin(metadata_df, x = "MYC_exon3_hypomet", y = "MYC_expr", 
               facet.by = "PMN_hit", color = "PMN_hit", palette = "lancet", 
               add = "boxplot",
               xlab = "MYC exon3 hypometylation", ylab = "Normalized MYC expression")
g3 <- g3 + rremove("legend")
g3 <- g3 + stat_compare_means()

g4 <- ggviolin(metadata_df, x = "mir34bc_hypermet", y = "MYC_expr", 
               facet.by = "PMN_hit", color = "PMN_hit", palette = "lancet", 
               add = "boxplot", 
               xlab = "miR-34b/c TSS hypermetylation", ylab = "Normalized MYC expression")
g4 <- g4 + rremove("legend")
g4 <- g4 + stat_compare_means()

g_comb <- ggpubr::ggarrange(g1, g2, g3, g4, labels = LETTERS[1:4])

ggsave("output/8.methylation_changes.pdf", g_comb, width = 8, height = 8)
