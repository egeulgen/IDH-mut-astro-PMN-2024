##### Script purpose: Assess associations of PMN status with various factors
##### Author: Ege Ulgen
##### Date: Dec 2023

library(ggpubr)
library(survival)
library(survminer)

source("scripts/utils.R")

metadata_df <- readRDS("data/selected_data/meta.RDS")

# clinical ----------------------------------------------------------------
metadata_df$age_at_diagnosis2 <- metadata_df$age_at_diagnosis / 365
### Age at diagnosis
g1 <- ggviolin(metadata_df, "PMN_hit", "age_at_diagnosis2", color = "PMN_hit", add =  "jitter",
               xlab = "PMN Hit", ylab = "Age at Diagnosis", palette = "lancet")
g1 <- g1 + rremove("legend")
g1 <- g1 + stat_compare_means(label.y = 80, label.x = 1.3)
g1 <- g1 + theme_minimal()
g1 <- g1 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())

### Survival
plot_KM_by_PMN_hit <- function(df, type, title = "") {
    tmp <- df
    tmp$Time <- df[, paste0(type, "_time")]
    tmp$Stat <- df[, type]
    
    tmp[, "PMN_hit"] <- as.character(tmp[, "PMN_hit"])
    tmp[, "PMN_hit"] <- sub("\n\\(n=\\d+\\)", "",  tmp[, "PMN_hit"])
    tmp[, "PMN_hit"] <- as.factor(tmp[, "PMN_hit"])
    
    fit_grade <- survfit(Surv(Time, Stat) ~ PMN_hit, data = tmp)
    g <- ggsurvplot(fit_grade, 
                    data = tmp, 
                    size = 1,
                    conf.int = FALSE,
                    pval = TRUE,
                    legend.labs = levels(tmp[, "PMN_hit"]),
                    linetype = "strata",
                    title = paste(title, type),
                    xlab = "Time",
                    newpage = TRUE)
    g <- g$plot
    g <- g + ggsci::scale_color_lancet()
    g <- g + theme(title = element_text(face = "bold", size = 12),
                   axis.text = element_text(size = 11),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.title.x = element_blank(),
                   legend.title = element_blank())
    g
}

all_PMN_KM_list <- list()
all_PMN_KM_list[["overall"]] <- ggarrange(plot_KM_by_PMN_hit(metadata_df, "CDR_OS"),
                                          plot_KM_by_PMN_hit(metadata_df, "CDR_PFI"),
                                          ncol = 2)

for (grade in levels(metadata_df$mol_grade)) {
    grade_df <- metadata_df[metadata_df$mol_grade == grade, ]
    all_PMN_KM_list[[grade]] <- ggarrange(plot_KM_by_PMN_hit(grade_df, "CDR_OS", title = paste("Mol.Grade", grade)),
                                                       plot_KM_by_PMN_hit(grade_df, "CDR_PFI"),
                                                       ncol = 2)
}

ggsave("output/9.PMN_survival.pdf", ggarrange(plotlist = all_PMN_KM_list, nrow = 4, labels = c("A", "B", "C", "D")), width = 12, height = 20)



g_MGMT <- plot_cat(df = metadata_df, var1 = "paper_MGMT.promoter.status", x_lbl = "MGMT promoter status")
g_Telomere <- plot_cat(df = metadata_df, var1 = "paper_Telomere.Maintenance", x_lbl = "Telomere Maintenance")

pdf("output/10.PMN_clinical_assoc.pdf", width = 10, height = 10)
ggarrange(g1, g_MGMT, g_Telomere , labels = LETTERS[1:3], font.label = list(face = "bold", size = 20))
dev.off()

# grade -------------------------------------------------------------------
g_grade_PMN <- plot_cat(df = metadata_df, var1 = "mol_grade", x_lbl = "Grade")
g_grade_MYC <- plot_cat(df = metadata_df, var1 = "mol_grade", x_lbl = "Grade", var2 = "MYC_amp",
                        y_lbl = "MYC amplification", fisher_position = 3.350352)
g_grade_para <- plot_cat(df = metadata_df, var1 = "mol_grade", x_lbl = "Grade", var2 = "MYC_paralog_amp",
                         y_lbl = "MYC paralog amplification", fisher_position = 3.134271)

pdf("output/11.PMN_grade_assoc.pdf", width = 7, height = 7)
ggarrange(g_grade_PMN, g_grade_MYC, g_grade_para, labels = LETTERS[1:3], font.label =  list(face = "bold", size = 20))
dev.off()


# burden - overall --------------------------------------------------------
range(metadata_df$snv_burden)
### SNV burden
g7 <- ggviolin(metadata_df, "PMN_hit", "snv_burden", color = "PMN_hit",
               xlab = "PMN Hit", ylab = "SNV Burden", palette = "lancet", add =  "jitter")
g7 <- g7 + rremove("legend")
g7 <- g7 + stat_compare_means(label.y = 3)
g7 <- g7 + scale_y_continuous(breaks = 0:3)
g7 <- g7 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())

### indel burden
range(metadata_df$indel_burden)
g8 <- ggviolin(metadata_df, "PMN_hit", "indel_burden", color = "PMN_hit",
               xlab = "PMN Hit", ylab = "InDel Burden", palette = "lancet", add =  "jitter")
g8 <- g8 + rremove("legend")
g8 <- g8 + stat_compare_means()
# g8 <- g8 + scale_y_continuous(breaks = seq(0, 0.5, 0.1))
g8 <- g8 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
g8

### WGII
range(metadata_df$WGII)
g9 <- ggviolin(metadata_df, "PMN_hit", "WGII", color = "PMN_hit",
               xlab = "PMN Hit", ylab = "Copy-number Alteration Frequency (wGII)", palette = "lancet", add =  "jitter")
g9 <- g9 + rremove("legend")
g9 <- g9 + stat_compare_means()
# g9 <- g9 + scale_y_continuous(breaks = seq(0, 0.5, 0.1))
g9 <- g9 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
g9

### CAER
g10 <- ggviolin(metadata_df, "PMN_hit", "chr_arm_event_ratio", color = "PMN_hit",
                xlab = "PMN Hit", ylab = "Degree of Aneuploidy (CAER)", palette = "lancet", add =  "jitter")
g10 <- g10 + rremove("legend")
g10 <- g10 + stat_compare_means()
# g10 <- g10 + scale_y_continuous(breaks = seq(0, 2, 0.1))
g10 <- g10 + theme(legend.position = "none",
                   axis.text = element_text(size = 11),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.title.x = element_blank())
g10

### CN amplitude
g11 <- ggviolin(metadata_df, "PMN_hit", "maxCN", color = "PMN_hit",
                xlab = "PMN Hit", ylab = "Copy-number Amplitude", palette = "lancet", add =  "jitter")
g11 <- g11 + rremove("legend")
g11 <- g11 + stat_compare_means()
# g11 <- g11 + scale_y_continuous(breaks = seq(0, 250, 25))
g11 <- g11 + theme(legend.position = "none",
                   axis.text = element_text(size = 11),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.title.x = element_blank())
g11


### CT
g12 <- plot_cat(df = metadata_df, var1 = "PMN_hit", var2 = "CT", val = "CT+", x_lbl = "PMN Hit", y_lbl = "% with Chromothripsis",
                y_lims = c(0, 0.4), incr = 0.05, fisher_position = .35)

pdf("output/12.PMN_metric_assoc_OVERALL.pdf", width = 8, height = 10)
ggarrange(g7, g8, g9, g10, g11, g12, labels = LETTERS[1:6], font.label = list(face = "bold", size = 18))
dev.off()

# burden by grade ---------------------------------------------------------
burden_by_grade <- function(grade) {
    metadata_df_grade <- metadata_df[metadata_df$mol_grade == grade, ]
    metadata_df_grade$PMN_hit <- as.character(metadata_df_grade$PMN_hit)
    
    tmp <- paste0("PMN-no hit\n(n=", as.character(table(metadata_df_grade$PMN_hit))[1])
    tmp2 <- paste0("PMN-hit\n(n=", as.character(table(metadata_df_grade$PMN_hit))[2])
    metadata_df_grade$PMN_hit <- sub("PMN-no hit\\n\\(n=\\d+", tmp, metadata_df_grade$PMN_hit)
    metadata_df_grade$PMN_hit <- sub("PMN-hit\\n\\(n=\\d+", tmp2, metadata_df_grade$PMN_hit)
    
    metadata_df_grade$PMN_hit <- factor(metadata_df_grade$PMN_hit, 
                                     levels = c(paste0(tmp, ")"),
                                                paste0(tmp2, ")")))
    
    ### SNV burden
    g7 <- ggviolin(metadata_df_grade, "PMN_hit", "snv_burden", color = "PMN_hit",
                   xlab = "PMN Hit", ylab = "SNV Burden", palette = "lancet", add =  "jitter")
    g7 <- g7 + rremove("legend")
    g7 <- g7 + stat_compare_means(label.y = 3)
    g7 <- g7 + scale_y_continuous(breaks = 0:3)
    g7 <- g7 + theme(legend.position = "none",
                     axis.text = element_text(size = 11),
                     axis.title.y = element_text(face = "bold", size = 14),
                     axis.title.x = element_blank())

    ### indel burden
    g8 <- ggviolin(metadata_df_grade, "PMN_hit", "indel_burden", color = "PMN_hit",
                   xlab = "PMN Hit", ylab = "InDel Burden", palette = "lancet", add =  "jitter")
    g8 <- g8 + rremove("legend")
    g8 <- g8 + stat_compare_means()
    g8 <- g8 + scale_y_continuous(breaks = seq(0, 0.4, 0.1))
    g8 <- g8 + theme(legend.position = "none",
                     axis.text = element_text(size = 11),
                     axis.title.y = element_text(face = "bold", size = 14),
                     axis.title.x = element_blank())
    
    ### WGII
    g9 <- ggviolin(metadata_df_grade, "PMN_hit", "WGII", color = "PMN_hit",
                   xlab = "PMN Hit", ylab = "Copy-number Alteration Frequency (wGII)", palette = "lancet", add =  "jitter")
    g9 <- g9 + rremove("legend")
    g9 <- g9 + stat_compare_means()
    g9 <- g9 + scale_y_continuous(breaks = seq(0, 0.3, 0.05))
    g9 <- g9 + theme(legend.position = "none",
                     axis.text = element_text(size = 11),
                     axis.title.y = element_text(face = "bold", size = 14),
                     axis.title.x = element_blank())
    
    ### CAER
    g10 <- ggviolin(metadata_df_grade, "PMN_hit", "chr_arm_event_ratio", color = "PMN_hit",
                    xlab = "PMN Hit", ylab = "Degree of Aneuploidy (CAER)", palette = "lancet", add =  "jitter")
    g10 <- g10 + rremove("legend")
    g10 <- g10 + stat_compare_means()
    g10 <- g10 + scale_y_continuous(breaks = seq(0, 0.3, 0.05))
    g10 <- g10 + theme(legend.position = "none",
                       axis.text = element_text(size = 11),
                       axis.title.y = element_text(face = "bold", size = 14),
                       axis.title.x = element_blank())
    
    ### CN amplitude
    g11 <- ggviolin(metadata_df_grade, "PMN_hit", "maxCN", color = "PMN_hit",
                    xlab = "PMN Hit", ylab = "Copy-number Amplitude", palette = "lancet", add =  "jitter")
    g11 <- g11 + rremove("legend")
    g11 <- g11 + stat_compare_means()
    g11 <- g11 + scale_y_continuous(breaks = seq(0, 100, 25))
    g11 <- g11 + theme(legend.position = "none",
                       axis.text = element_text(size = 11),
                       axis.title.y = element_text(face = "bold", size = 14),
                       axis.title.x = element_blank())
    
    ### CT
    g12 <- plot_cat(df = metadata_df_grade, var1 = "PMN_hit", var2 = "CT", val = "CT+", x_lbl = "PMN Hit", y_lbl = "% with Chromothripsis",
                    y_lims = c(0, 0.3), incr = 0.05, fisher_position = 0.3)
    
    g_final <- ggarrange(g7, g8, g9, g10, g11, g12, labels = LETTERS[1:6], font.label = list(face = "bold", size = 18))
    g_final <- annotate_figure(g_final, text_grob(grade, color = "black", face = "bold", size = 14))
    return(g_final)
}

pdf("output/13.PMN_metric_assoc_by_GRADE.pdf", width = 8, height = 10)
for (grade in levels(metadata_df$mol_grade)) {
    plot(burden_by_grade(grade))
}
dev.off()
