##### Script purpose: Assess alterations in the PMN
##### Author: Ege Ulgen
##### Date: Dec 2023

library(ggplot2)
library(maftools)

# prep --------------------------------------------------------------------
metadata_df <- readRDS("data/selected_data/meta.RDS")
som_df <- readRDS("data/selected_data/maf.RDS")
gene_scna <- readRDS("data/selected_data/CN_gene_level.RDS")

# proximal MYC network genes from Schaub et al. 2018
# "FBXW7" from Bai et al. 2016
PMN_genes <- read.csv("data/PMN_gene_list.csv")

pmn_som_df <- som_df[som_df$Hugo_Symbol %in% PMN_genes$Gene.name, ]
table(pmn_som_df$Variant_Classification)
pmn_som_df$Tumor_Sample_Barcode <- pmn_som_df$patient_barcode


pmn_scna_df <- gene_scna[gene_scna$symbol %in% PMN_genes$Gene.name, ]
hist(pmn_scna_df$Segment_Mean, breaks = 100)
abline(v = -.3, col = "red")
abline(v = .3, col = "red")

pmn_scna_df$SCNA_type <- ifelse(pmn_scna_df$Segment_Mean > .3, "amp", 
                                ifelse(pmn_scna_df$Segment_Mean < -.3, "del", "CN"))
pmn_scna_df <- pmn_scna_df[pmn_scna_df$SCNA_type != "CN", ]

# oncoplot ----------------------------------------------------------------
cn_table <- pmn_scna_df[, c("symbol", "patient_barcode", "SCNA_type")]

clin_df <- as.data.frame(metadata_df)
clin_df$Tumor_Sample_Barcode <- clin_df$patient

maf_df_pmn <- read.maf(maf = som_df,
                       cnTable = cn_table,
                       clinicalData = clin_df,
                       verbose = TRUE, isTCGA = TRUE)

mol_grade_cols <- c("2" = "#FBA465", "3" = "#F86E51", "4" = "#D1193E")

pdf("output/1.PMN_oncoplot.pdf", width = 10, height = 5)
oncoplot(maf_df_pmn, genes = PMN_genes$Gene.name, 
         removeNonMutated = FALSE, top = length(PMN_genes$Gene.name),
         drawColBar = FALSE, 
         clinicalFeatures = "mol_grade", sortByAnnotation = TRUE,
         anno_height = .5,
         annotationColor = list(mol_grade = mol_grade_cols))
dev.off()


#exclusive/co-occurance event analysis. 
pdf("output/2.PMN_mutex.pdf", width = 7, height = 7)
somaticInteractions(maf = maf_df_pmn, genes = PMN_genes$Gene.name, pvalue = c(0.05, 0.1))
dev.off()




maf_df_pmn2 <- subsetMaf(maf_df_pmn, clinQuery = "mol_grade == 2")
maf_df_pmn3 <- subsetMaf(maf_df_pmn, clinQuery = "mol_grade == 3")
maf_df_pmn4 <- subsetMaf(maf_df_pmn, clinQuery = "mol_grade == 4")

sorted_PMN_genes <- c("MYC", "FBXW7", "MAX", "MXI1", "MLXIPL", "MLXIP", "MXD3", "MXD4", "MYCN", "MGA", "MYCL", "MLX", "MNT", "MXD1")
pdf("output/3.G2_PMN_oncoplot.pdf", width = 10, height = 5)
oncoplot(maf_df_pmn2, genes = sorted_PMN_genes, 
         removeNonMutated = FALSE, top = length(sorted_PMN_genes),
         keepGeneOrder = TRUE,
         drawColBar = FALSE, drawRowBar = FALSE)
# anno_height = .4,
# clinicalFeatures = "mol_grade", annotationColor = list(mol_grade = mol_grade_cols))
dev.off()

pdf("output/4.G3_PMN_oncoplot.pdf", width = 9, height = 5)
oncoplot(maf_df_pmn3, genes = sorted_PMN_genes, 
         removeNonMutated = FALSE, top = length(sorted_PMN_genes),
         keepGeneOrder = TRUE,
         drawColBar = FALSE, drawRowBar = FALSE)
# anno_height = .4,
# clinicalFeatures = "mol_grade", annotationColor = list(mol_grade = mol_grade_cols))
dev.off()

pdf("output/5.G4_PMN_oncoplot.pdf", width = 5, height = 5)
oncoplot(maf_df_pmn4, genes = sorted_PMN_genes,
         keepGeneOrder = TRUE,
         removeNonMutated = FALSE, top = length(sorted_PMN_genes),
         drawColBar = FALSE, drawRowBar = FALSE)
# anno_height = .4,
# clinicalFeatures = "mol_grade", annotationColor = list(mol_grade = mol_grade_cols))
dev.off()

# for network figure ------------------------------------------------------
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


metadata_df$PMN_hit <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id, "hit", "no hit")
metadata_df$MYC_amp <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MYC"], "amp", "no amp")
metadata_df$MYC_paralog_amp <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene %in% c("MYC", "MYCN", "MYCL") & PMN_alterations_df$alteration == "amp"], "amp", "no amp")

metadata_df$FBXW7 <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "FBXW7"], 1, 0)
metadata_df$MAX <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MAX"], 1, 0)
metadata_df$MXI1 <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MXI1"], 1, 0)
metadata_df$MLXIPL <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MLXIPL"], 1, 0)
metadata_df$MLXIP <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MLXIP"], 1, 0)
metadata_df$MXD3 <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MXD3"], 1, 0)
metadata_df$MXD4 <- ifelse(metadata_df$patient %in% PMN_alterations_df$donor_id[PMN_alterations_df$gene == "MXD4"], 1, 0)

# MYC paralogs %
sum(metadata_df$MYC_paralog_amp == "amp") / nrow(metadata_df) * 100
sum(metadata_df[metadata_df$mol_grade == 2, ]$MYC_paralog_amp == "amp") / sum(metadata_df$mol_grade == 2) * 100
sum(metadata_df[metadata_df$mol_grade == 3, ]$MYC_paralog_amp == "amp") / sum(metadata_df$mol_grade == 3) * 100
sum(metadata_df[metadata_df$mol_grade == 4, ]$MYC_paralog_amp == "amp") / sum(metadata_df$mol_grade == 4) * 100

### Calculate %s of PMN alterations per Grade
all_PMN_perc_df <- c()
for (gene in PMN_genes$Gene.name) {
    
    for (grade in 2:4) {
        
        sub_df <- PMN_alterations_df[PMN_alterations_df$donor_id %in% metadata_df$patient[metadata_df$mol_grade == grade], ]
        N <- sum(metadata_df$mol_grade == grade)
        
        amp_df <- sub_df[sub_df$gene == gene & sub_df$alteration == "amp", ]
        amp_num <- length(unique(amp_df$donor_id))
        
        del_df <- sub_df[sub_df$gene == gene & sub_df$alteration == "del", ]
        del_num <- length(unique(del_df$donor_id))
        
        mut_df <- sub_df[sub_df$gene == gene & !sub_df$alteration %in% c("amp", "del"), ]
        mut_num <- length(unique(mut_df$donor_id))
        
        all_PMN_perc_df <- rbind(all_PMN_perc_df, data.frame(Gene = gene,
                                                             Grade = grade,
                                                             Amp = round(amp_num / N * 100, 2), 
                                                             Del = round(del_num / N * 100, 2), 
                                                             Mut = round(mut_num / N * 100, 2)))
    }
}

range(all_PMN_perc_df$Amp)
range(all_PMN_perc_df$Del)
range(all_PMN_perc_df$Mut)

# Amp heatmap
ga <- ggplot(all_PMN_perc_df, aes(x = Grade, y = Gene))
ga <- ga + geom_tile(aes(fill = Amp))
ga <- ga + geom_text(aes(label = Amp))
ga <- ga + scale_fill_gradient(low = "#f79292", high = "#ff0000", limits = c(0, 31))

# Del heatmap
gd <- ggplot(all_PMN_perc_df, aes(x = Grade, y = Gene))
gd <- gd + geom_tile(aes(fill = Del))
gd <- gd + geom_text(aes(label = Del))
gd <- gd + scale_fill_gradient(low = "#8a97f9", high = "#0013FF", limits = c(0, 40))

# Mut heatmap
gm <- ggplot(all_PMN_perc_df, aes(x = Grade, y = Gene))
gm <- gm + geom_tile(aes(fill = Mut))
gm <- gm + geom_text(aes(label = Mut))
gm <- gm + scale_fill_gradient(low = "#d4f9bb", high = "#64ff00", limits = c(0, 1))

pdf("output/PMN_alt_percs_tile.pdf", width = 10, height = 10)
ggpubr::ggarrange(ga, gd, gm)
dev.off()

# summarize ---------------------------------------------------------------
metadata_df <- as.data.frame(metadata_df)

tmp <- summary(as.factor(metadata_df$PMN_hit))
pmn_status_with_counts <- paste0("PMN-", names(tmp), "\n(n=", tmp, ")")
names(pmn_status_with_counts) <- names(tmp)
metadata_df$PMN_hit <- pmn_status_with_counts[metadata_df$PMN_hit]
metadata_df$PMN_hit <- factor(metadata_df$PMN_hit, levels = pmn_status_with_counts)

tmp <- summary(as.factor(metadata_df$MYC_amp))
myc_status_with_counts <- paste0("MYC-", names(tmp), "\n(n=", tmp, ")")
names(myc_status_with_counts) <- names(tmp)
metadata_df$MYC_amp <- myc_status_with_counts[metadata_df$MYC_amp]
metadata_df$MYC_amp <- factor(metadata_df$MYC_amp, levels = myc_status_with_counts)

metadata_df$paper_MGMT.promoter.status <- factor(metadata_df$paper_MGMT.promoter.status, levels = c("Unmethylated", "Methylated"))
metadata_df$paper_TERT.promoter.status <- factor(metadata_df$paper_TERT.promoter.status, levels = c("WT", "Mutant"))

tmp <- summary(as.factor(metadata_df$MYC_paralog_amp))
myc_para_status_with_counts <- paste0("MYCparalog-", names(tmp), "\n(n=", tmp, ")")
names(myc_para_status_with_counts) <- names(tmp)
metadata_df$MYC_paralog_amp <- myc_para_status_with_counts[metadata_df$MYC_paralog_amp]
metadata_df$MYC_paralog_amp <- factor(metadata_df$MYC_paralog_amp, levels = myc_para_status_with_counts)

tmp <- summary(as.factor(metadata_df$mol_grade))
mol_grade_with_counts <- paste0("Grade ", names(tmp), "\n(n=", tmp, ")")
names(mol_grade_with_counts) <- names(tmp)
metadata_df$mol_grade <- mol_grade_with_counts[as.character(metadata_df$mol_grade)]
metadata_df$mol_grade <- factor(metadata_df$mol_grade, levels = mol_grade_with_counts)

saveRDS(metadata_df, "data/selected_data/meta.RDS")

# Overall PMN mutex -------------------------------

# proximal MYC network genes from Schaub 2018
# "FBXW7" from Bai 2016
PMN_genes <- read.csv("data/PMN_gene_list.csv")
selected_genes <- c("PMN_hit", "ATRX", "TP53", "EGFR", "PTEN", "CDKN2A", "PIK3CA", "BRAF")

sel_som_df <- som_df[som_df$Hugo_Symbol %in% c(PMN_genes$Gene.name, selected_genes), ]
table(sel_som_df$Variant_Classification)
# sel_som_df <- sel_som_df[sel_som_df$Variant_Classification %in% c("Nonsense_Mutation", "Missense_Mutation"), ]
sel_som_df$Hugo_Symbol[sel_som_df$Hugo_Symbol %in% PMN_genes$Gene.name] <- "PMN_hit"
sel_som_df$Tumor_Sample_Barcode <- sel_som_df$patient_barcode

sel_scna_df <- gene_scna[gene_scna$symbol %in% c(PMN_genes$Gene.name, selected_genes), ]
sel_scna_df$SCNA_type <- ifelse(sel_scna_df$Segment_Mean >= .3, "amp",
                                ifelse(sel_scna_df$Segment_Mean <= -.3, "del", "CN"))
sel_scna_df <- sel_scna_df[sel_scna_df$SCNA_type != "CN", ]
sel_scna_df$symbol[sel_scna_df$symbol %in% PMN_genes$Gene.name] <- "PMN_hit"

cn_table <- sel_scna_df[, c("symbol", "patient_barcode", "SCNA_type")]

clin_df <- metadata_df
clin_df$Tumor_Sample_Barcode <- clin_df$patient


maf_df_pmn <- read.maf(maf = sel_som_df,
                       cnTable = cn_table,
                       clinicalData = clin_df,
                       verbose = TRUE, isTCGA = TRUE)


#exclusive/co-occurance event analysis. 
pdf("output/7.PMN_overall_mutex.pdf", width = 7, height = 7)
somaticInteractions(maf = maf_df_pmn, genes = selected_genes, pvalue = c(0.05, 0.1),
                    geneOrder = selected_genes)
dev.off()
