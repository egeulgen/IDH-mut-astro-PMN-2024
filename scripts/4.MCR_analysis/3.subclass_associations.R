##### Script purpose: Investigate various associations of PMN classes (WT, 19q microdel, PMN hit)
##### Author: Ege Ulgen
##### Date: Feb 2024

library(ggpubr)
library(survival)
library(survminer)

metadata_df <- readRDS("data/selected_data/meta.RDS")

metadata_df$PMN_hit <- as.character(metadata_df$PMN_hit)
metadata_df$PMN_hit <- sub("\\n.+", "", metadata_df$PMN_hit)
metadata_df$PMN_hit <- sub("PMN-", "", metadata_df$PMN_hit)
metadata_df$PMN_hit <- factor(metadata_df$PMN_hit, levels = c("no hit", "hit"))

annotated_df <- read.csv("output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv")
altered_patients <- unique(unlist(strsplit(annotated_df$altered_ids, ";")))

# MYC expr. ---------------------------------------------------------------
metadata_df$class <- ifelse(metadata_df$PMN_hit == "hit", "PMN-hit",
                          ifelse(metadata_df$patient %in% altered_patients, "microdel19q", "WT"))
metadata_df$class <- factor(metadata_df$class, levels = c("WT", "microdel19q", "PMN-hit"))
table(metadata_df$class)

### overall violin
g <- ggviolin(metadata_df, x = "class", y = "MYC_expr", color = "class", palette = "lancet", add = c("jitter", "boxplot"), xlab = "", ylab = "Normalized MYC expression")
g <- g + stat_compare_means(comparisons = list(c("WT", "microdel19q"), c("WT", "PMN-hit"),
                                               c("microdel19q", "PMN-hit")),
                            method = "wilcox.test",
                            method.args = list("alternative" = "less")
                            )
g <- g + rremove("legend")
g

ggsave("output/23.PMNneg_alter_MYC_expr.pdf", g, width = 6, height = 7)

# OS - KM -----------------------------------------------------------------
plot_KM_by_class <- function(df, type, title = "") {
    tmp <- df
    tmp$Time <- df[, paste0(type, "_time")]
    tmp$Stat <- df[, type]
    
    fit_grade <- survfit(Surv(Time, Stat) ~ class, data = tmp)
    g <- ggsurvplot(fit_grade, 
                    data = tmp, 
                    size = 1,
                    conf.int = FALSE,
                    pval = TRUE,
                    # legend.labs = c("No PMN hit", "PMN hit"),
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

all_class_KM_list <- list()
all_class_KM_list[["overall"]] <- ggarrange(plot_KM_by_class(metadata_df, "CDR_OS"),
                                          plot_KM_by_class(metadata_df, "CDR_PFI"),
                                          ncol = 2)

for (grade in levels(metadata_df$mol_grade)) {
    grade_df <- metadata_df[metadata_df$mol_grade == grade, ]
    all_class_KM_list[[grade]] <- ggarrange(plot_KM_by_class(grade_df, "CDR_OS", title = paste("Mol.Grade", grade)),
                                                       plot_KM_by_class(grade_df, "CDR_PFI"),
                                                       ncol = 2)
}

g_surv <- ggarrange(plotlist = all_class_KM_list, labels = c("A", "B", "C", "D"),  nrow = 4)

ggsave("output/24.alter_class_survival.pdf", g_surv, width = 8, height = 15)

# by grade ----------------------------------------------------------------
tbl <- table(metadata_df$mol_grade, metadata_df$class)

perc_df <- round(tbl / rowSums(tbl), 4)
perc_df <- as.data.frame(perc_df)
perc_df <- perc_df[perc_df$Freq != 0, ]

p <- ggplot(perc_df, aes(y = Freq, x = Var1, fill = Var2))
p <- p + geom_bar(stat = "identity", position = "stack", aes(color = Var2))
p <- p + geom_text(aes(label = paste0(Freq * 100, "%")),
                   position = position_stack(vjust = 0.5), size = 4, color = "white")
p <- p + geom_text(label = paste0("Fisher's Exact,\np = ", round(fisher.test(tbl)$p.value, 3)),
                   size = 3, position = position_stack(vjust = 1),
                   data = perc_df[perc_df$Var1 == perc_df$Var1[2], ][1, ])
p <- p + ggsci::scale_fill_lancet()
p <- p + ggsci::scale_color_lancet()
p <- p + xlab("Molecular Grade")
p <- p + scale_y_continuous(breaks = seq(0, 1, 0.25),
                            labels = scales::percent_format(accuracy = 1),
                            expand = c(0, 0))
p <- p + theme_minimal()
p <- p + theme(legend.title = element_blank(),
               axis.text = element_text(size = 11),
               axis.title.x = element_text(face = "bold", size = 14),
               axis.title.y = element_blank())
p

ggsave("output/25.all_classes_by_grade.pdf", p, width = 6, height = 6)


non_pmn_metadata_df <- metadata_df[metadata_df$class != "PMN-hit", ]
non_pmn_metadata_df$class <- droplevels(non_pmn_metadata_df$class)
# only 1 grade 4, drop
non_pmn_metadata_df <- non_pmn_metadata_df[non_pmn_metadata_df$mol_grade != 4, ]
tbl <- table(non_pmn_metadata_df$mol_grade, non_pmn_metadata_df$class)

perc_df <- round(tbl / rowSums(tbl), 4)
perc_df <- as.data.frame(perc_df)
perc_df <- perc_df[perc_df$Freq != 0, ]

p2 <- ggplot(perc_df, aes(y = Freq, x = Var1, fill = Var2))
p2 <- p2 + geom_bar(stat = "identity", position = "stack", aes(color = Var2))
p2 <- p2 + geom_text(aes(label = paste0(Freq * 100, "%")),
                     position = position_stack(vjust = 0.5), size = 4, color = "white")
p2 <- p2 + geom_text(label = paste0("Fisher's Exact,\np = ", round(fisher.test(tbl)$p.value, 3)),
                     size = 3, position = position_stack(vjust = 1),
                     data = perc_df[perc_df$Var1 == perc_df$Var1[2], ][1, ])
p2 <- p2 + ggsci::scale_fill_lancet()
p2 <- p2 + ggsci::scale_color_lancet()
p2 <- p2 + xlab("Molecular Grade")
p2 <- p2 + scale_y_continuous(breaks = seq(0, 1, 0.25),
                              labels = scales::percent_format(accuracy = 1),
                              expand = c(0, 0))
p2 <- p2 + theme_minimal()
p2 <- p2 + theme(legend.title = element_blank(),
                 axis.text = element_text(size = 11),
                 axis.title.x = element_text(face = "bold", size = 14),
                 axis.title.y = element_blank())
p2

ggsave("output/26.PMNneg_classes_by_grade.pdf", p2, width = 6, height = 6)

# burden assoc - class ----------------------------------------------------
summary(metadata_df$class)

comps_list <- list(c("WT", "microdel19q"), c("WT", "PMN-hit"),
                   c("microdel19q", "PMN-hit"))

### SNV burden
g6 <- ggviolin(metadata_df, "class", "snv_burden", color = "class",
               xlab = "", ylab = "SNV Burden", palette = "lancet", add =  "jitter")
g6 <- g6 + rremove("legend")
g6 <- g6 + stat_compare_means(comparisons = comps_list,
                              method = "wilcox.test",
                              method.args = list("alternative" = "less"))
g6 <- g6 + scale_y_continuous(breaks = 0:15)
g6 <- g6 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
g6

### indel burden
g7 <- ggviolin(metadata_df, "class", "indel_burden", color = "class",
               xlab = "", ylab = "InDel Burden", palette = "lancet", add =  "jitter")
g7 <- g7 + rremove("legend")
g7 <- g7 + stat_compare_means(comparisons = comps_list,
                              method = "wilcox.test",
                              method.args = list("alternative" = "less"))
g7 <- g7 + scale_y_continuous(breaks = seq(0, 2, 0.1))
g7 <- g7 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
g7

### WGII
g8 <- ggviolin(metadata_df, "class", "WGII", color = "class",
               xlab = "", ylab = "Copy-number Alteration Frequency (wGII)", palette = "lancet", add =  "jitter")
g8 <- g8 + rremove("legend")
g8 <- g8 + stat_compare_means(comparisons = comps_list,
                              method = "wilcox.test",
                              method.args = list("alternative" = "less"))
g8 <- g8 + scale_y_continuous(breaks = seq(0, 2, 0.1))
g8 <- g8 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
g8

### CAER
g9 <- ggviolin(metadata_df, "class", "chr_arm_event_ratio", color = "class",
               xlab = "", ylab = "Degree of Aneuploidy (CAER)", palette = "lancet", add =  "jitter")
g9 <- g9 + rremove("legend")
g9 <- g9 + stat_compare_means(comparisons = comps_list,
                              method = "wilcox.test",
                              method.args = list("alternative" = "less"))
g9 <- g9 + scale_y_continuous(breaks = seq(0, 2, 0.1))
g9 <- g9 + theme(legend.position = "none",
                 axis.text = element_text(size = 11),
                 axis.title.y = element_text(face = "bold", size = 14),
                 axis.title.x = element_blank())
g9

### CN amplitude
g10 <- ggviolin(metadata_df, "class", "maxCN", color = "class",
                xlab = "", ylab = "Copy-number Amplitude", palette = "lancet", add =  "jitter")
g10 <- g10 + rremove("legend")
g10 <- g10 + stat_compare_means(comparisons = comps_list,
                                method = "wilcox.test",
                                method.args = list("alternative" = "less"))
g10 <- g10 + scale_y_continuous(breaks = seq(0, 250, 25))
g10 <- g10 + theme(legend.position = "none",
                   axis.text = element_text(size = 11),
                   axis.title.y = element_text(face = "bold", size = 14),
                   axis.title.x = element_blank())
g10


### CT
tbl <- table(metadata_df$class, metadata_df$CT)

perc_df <- round(tbl / rowSums(tbl), 4)
perc_df <- as.data.frame(perc_df)

perc_df <- perc_df[perc_df$Var2 == "CT+", ]

g11 <- ggplot(perc_df, aes(y = Freq, x = Var1, fill = Var1))
g11 <- g11 + geom_bar(stat = "identity", position = "stack", aes(color = Var1))
g11 <- g11 + geom_text(aes(label = paste0(Freq * 100, "%")),
                   position = position_stack(vjust = 0.5), size = 4, color = "white")
g11 <- g11 + geom_text(label = paste0("Cochran-Armitage Test for Trend,\np = ", round(DescTools::CochranArmitageTest(t(tbl), alternative = "one.sided")$p.value, 3)),
                   size = 3, position = position_stack(vjust = 1),
                   data = perc_df[perc_df$Var1 == perc_df$Var1[2], ])
g11 <- g11 + ggsci::scale_fill_lancet()
g11 <- g11 + ggsci::scale_color_lancet()
g11 <- g11 + ylab("% with Chromothripsis") + xlab("")
g11 <- g11 + scale_y_continuous(breaks = seq(0, 1, 0.05), limits = c(0, .5),
                            labels = scales::percent_format(accuracy = 1),
                            expand = c(0, 0))
g11 <- g11 + theme_minimal()
g11 <- g11 + theme(legend.position = "none",
               axis.text = element_text(size = 11),
               axis.title.y = element_text(face = "bold", size = 14),
               axis.title.x = element_blank())
g11

g_comb2 <- ggarrange(g6, g7, g8, g9, g10, g11, labels = LETTERS[3:8], font.label = list(size = 20, face = "bold"))

ggsave("output/27.alter_class_burden.pdf", g_comb2, width = 16, height = 10)
