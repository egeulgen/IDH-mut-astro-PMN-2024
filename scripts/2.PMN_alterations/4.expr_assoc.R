##### Script purpose: Assess associations of MYC expression with various factors
##### Author: Ege Ulgen
##### Date: Dec 2023

library(DESeq2)
library(ggpubr)

metadata_df <- readRDS("data/selected_data/meta.RDS")
dds_norm <- readRDS("data/selected_data/processed_expr.RDS")

MYC_expr <- assay(dds_norm[rowData(dds_norm)$gene_name == "MYC", ])
metadata_df$MYC_expr <- MYC_expr[match(metadata_df$patient, colnames(MYC_expr))]

grade_levels <- levels(metadata_df$mol_grade)

for_plot_df <- data.frame(donor_id = metadata_df$patient,
                          grade = metadata_df$mol_grade,
                          PMN_hit = metadata_df$PMN_hit,
                          MYC_amp = metadata_df$MYC_amp,
                          MYC_paralog_amp = metadata_df$MYC_paralog_amp,
                          MYC = metadata_df$MYC_expr)

mol_grade_cols <- c("#FBA465", "#F86E51", "#D1193E")
names(mol_grade_cols) <- grade_levels

### expression level by grade ----------------
comps <- list(c(grade_levels[1], grade_levels[2]),
              c(grade_levels[1], grade_levels[3]),
              c(grade_levels[2], grade_levels[3]))
p <- ggviolin(for_plot_df, "grade", "MYC", color = "grade", add = "jitter", palette = mol_grade_cols) 
p <- p + stat_compare_means(comparisons = comps, step.increase = 0.12) + rremove("legend")
p <- p + scale_y_continuous(breaks = 1:25)
p <- p + ylab("MYC Expression")
p <- p + theme(legend.position = "none",
               title = element_text(face = "bold", size = 12),
               axis.text = element_text(size = 11),
               axis.title.y = element_text(face = "bold", size = 14),
               axis.title.x = element_blank())

pdf("output/14.MYC_expr_by_grade.pdf", width = 4, height = 5)
p
dev.off()

### expression level by PMN hit for MYC -------------
p2 <- ggviolin(for_plot_df, "PMN_hit", "MYC", color = "PMN_hit", add = "jitter", palette = "lancet")
p2 <- p2 + stat_compare_means(label.y = 14.35, label.x = 1.25, label.sep = ",\n", size = 3.5) + rremove("legend")
p2 <- p2 + scale_y_continuous(breaks = 1:25)
p2 <- p2 + ylab("MYC Expression")
p2 <- p2 + theme(legend.position = "none",
                 title = element_text(face = "bold", size = 12),
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
p2


for_plot_df$PMN_hit2 <- sub("\\n\\(.*", "", for_plot_df$PMN_hit)
for_plot_df$PMN_hit2 <- factor(for_plot_df$PMN_hit2, levels = c("PMN-hit", "PMN-no hit"))
p3 <- ggviolin(for_plot_df, "PMN_hit2", "MYC", color = "PMN_hit", add = "jitter", palette = "lancet", facet.by = "grade")
p3 <- p3 + stat_compare_means(label.y = 14.35, label.x = 1.25, label.sep = ",\n", size = 3.5) + rremove("legend")
p3 <- p3 + scale_y_continuous(breaks = 1:25)
p3 <- p3 + ylab("MYC Expression")
p3 <- p3 + theme(legend.position = "none",
                 title = element_text(face = "bold", size = 12),
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())

pdf("output/15.MYC_expr_by_hit.pdf", width = 8, height = 8)
ggarrange(ggarrange(ggplot() + theme_void(), p2, ggplot() + theme_void(), nrow = 1), p3, nrow = 2)
dev.off()

### expression level by event for MYC ---------------------------------------
scna_df <- readRDS('data/selected_data/CN_gene_level.RDS')

tmp <- scna_df[scna_df$symbol == "MYC", ]
tmp <- tmp[tmp$patient_barcode %in% dds_norm$patient, ]
tmp <- tmp[order(tmp$Segment_Mean, decreasing = TRUE), ]
tmp <- tmp[!duplicated(tmp$patient_barcode), ]

for_plot_df <- data.frame(donor_id = metadata_df$patient,
                          grade = metadata_df$mol_grade,
                          PMN_hit = metadata_df$PMN_hit,
                          MYC_amp = metadata_df$MYC_amp,
                          MYC_paralog_amp = metadata_df$MYC_paralog_amp,
                          MYC_expr = metadata_df$MYC_expr,
                          MYC_logR = tmp$Segment_Mean[match(dds_norm$patient, tmp$patient_barcode)])

ggscatter(for_plot_df, "MYC_logR", "MYC_expr", color = "grade", palette = mol_grade_cols, cor.coef = TRUE)
ggscatter(for_plot_df, "MYC_logR", "MYC_expr", color = "grade", facet.by = "PMN_hit", palette = mol_grade_cols, cor.coef = TRUE)

### prep
som_df <- readRDS("data/selected_data/maf.RDS")

# proximal MYC network genes from Schaub 2018
# "FBXW7" from Bai 2016
PMN_genes <- read.csv("data/PMN_gene_list.csv")

pmn_som_df <- som_df[som_df$Hugo_Symbol %in% PMN_genes$Gene.name, ]
table(pmn_som_df$Variant_Classification)
pmn_som_df$Tumor_Sample_Barcode <- pmn_som_df$patient_barcode

pmn_scna_df <- scna_df[scna_df$symbol %in% PMN_genes$Gene.name, ]
pmn_scna_df$SCNA_type <- ifelse(pmn_scna_df$Segment_Mean >= 0.3, "amp",
                                ifelse(pmn_scna_df$Segment_Mean <= -0.3, "del", "CN"))
pmn_scna_df <- pmn_scna_df[pmn_scna_df$SCNA_type != "CN", ]

### summarize
PMN_alterations_df <- c()
all_donors <- metadata_df$patient
for (donor in all_donors) {
    donor_som <- pmn_som_df[pmn_som_df$Tumor_Sample_Barcode == donor, ]
    donor_scna <- pmn_scna_df[pmn_scna_df$patient_barcode == donor, ]
    
    for (i in seq_len(nrow(donor_som))) {
        tmp <- data.frame(donor_id = donor,
                          gene = donor_som$Hugo_Symbol[i],
                          alteration = donor_som$Variant_Classification[i])
        PMN_alterations_df <- rbind(PMN_alterations_df, tmp)
    }
    
    for (i in seq_len(nrow(donor_scna))) {
        tmp <- data.frame(donor_id = donor,
                          gene = donor_scna$symbol[i],
                          alteration = donor_scna$SCNA_type[i])
        PMN_alterations_df <- rbind(PMN_alterations_df, tmp)
    }
    
}
table(PMN_alterations_df$alteration)

g_list <- list()
for (gene in PMN_genes$Gene.name) {
    tmp <- PMN_alterations_df[PMN_alterations_df$gene == gene, ]
    
    for_plot_df[, gene] <- tmp$alteration[match(for_plot_df$donor_id, tmp$donor_id)]
    for_plot_df[, gene][is.na(for_plot_df[, gene])] <- "WT"
    
    g_list[[gene]] <- ggviolin(for_plot_df, gene, "MYC_expr", color = gene, add = "jitter", palette = "lancet", title = gene) + 
        stat_compare_means() + rremove("legend") + scale_y_continuous(breaks = 1:25) +
        ylab("MYC Expression") + theme(legend.position = "none",
                                       title = element_text(face = "bold", size = 12),
                                       axis.text = element_text(size = 11),
                                       axis.title.y = element_text(face = "bold", size = 14),
                                       axis.title.x = element_blank())
}

pdf("output/16.MYC_expr_by_PMN_event.pdf", width = 12, height = 12)
ggarrange(plotlist=g_list, labels = LETTERS[1:14], font.label = list(size = 20, face = "bold"))
dev.off()
