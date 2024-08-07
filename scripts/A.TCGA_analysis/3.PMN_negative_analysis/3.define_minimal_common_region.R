##### Script purpose: Further analysis of altered genes associated with 
### increased MYC expression in PMN-negative samples, specifically to 
### define a minimal common region based on overlapping copy-number segments
##### Author: Ege Ulgen
##### Date: Feb 2023

library(ggplot2)
library(GenomicRanges)
library(Gviz)
library(ggpubr)


gene_level_scna_per_case_df <- readRDS("data/selected_data/CN_gene_level.RDS")
assoc_genes_df <- read.csv("output/PMN_neg_analysis/PMN_neg_inc_MYC_genes.csv")
altered_cases <- unique(unlist(sapply(assoc_genes_df$altered_ids, function(x) unlist(strsplit(x, ";")))))
gene_level_scna_per_case_df <- gene_level_scna_per_case_df[gene_level_scna_per_case_df$patient_barcode %in% altered_cases, ]

sum(assoc_genes_df$n_mut) # not many mut.s overall
sum(assoc_genes_df$n_amp) # not many amp.s overall
sum(assoc_genes_df$n_del) # deletion is the most common alteration

# investigate number of alterations per each alteration type --------------
aggregate(n_mut ~ cytoband, data = assoc_genes_df, FUN = sum) # not many mut.s overall
aggregate(n_amp ~ cytoband, data = assoc_genes_df, FUN = sum) # not many amp.s overall
aggregate(n_del ~ cytoband, data = assoc_genes_df, FUN = sum) # genes in chr19q13.42 and chr19q13.43 have very high numbers of del.s

# visualize SCNA segments associated with deleted genes -------------------
assoc_gene_per_case_scna_df <- gene_level_scna_per_case_df[gene_level_scna_per_case_df$symbol %in% assoc_genes_df$hgnc_symbol, ]

assoc_segments_per_case_scna_df <- assoc_gene_per_case_scna_df[, c("patient_barcode", "Segment_Mean", "segment_start", "segment_end")]
assoc_del_segments <- assoc_segments_per_case_scna_df[assoc_segments_per_case_scna_df$Segment_Mean < -.3, ]
assoc_del_segments <- unique(assoc_del_segments)
assoc_del_segments <- assoc_del_segments[order(assoc_del_segments$segment_start, decreasing = FALSE), ]

g <- ggplot(assoc_del_segments, aes(
    x = segment_start, xend = segment_end, y = patient_barcode, yend = patient_barcode, color = patient_barcode)
) 
g <- g + geom_segment(linewidth = 1)
g <- g + scale_x_continuous(
    "Position", 
    limits = c(min(assoc_del_segments$segment_start), max(assoc_del_segments$segment_end)),
    labels = scales::comma
)
g <- g + theme_minimal()
g <- g + scale_color_manual(values = colorRampPalette(c("#00264D","#5C93D1"))(21))
g <- g + theme(
    axis.text.y = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
)
g

# find minimal common region ----------------------------------------------
assoc_del_segments_gr <- GRanges(
    seqnames = "chr19", 
    ranges = IRanges(start=assoc_del_segments$segment_start, end=assoc_del_segments$segment_end),
    patient_barcode=assoc_del_segments$patient_barcode
)

MCR_gr <- intersect(
    assoc_del_segments_gr[assoc_del_segments_gr$patient_barcode == altered_cases[1]], 
    assoc_del_segments_gr[assoc_del_segments_gr$patient_barcode == altered_cases[2]]
)
for (case in altered_cases[-c(1,2)]) {
    MCR_gr <- intersect(MCR_gr, assoc_del_segments_gr[assoc_del_segments_gr$patient_barcode == case])
}

MCR_df <- data.frame(MCR_gr)

g <- g + geom_vline(xintercept = c(MCR_df$start,  MCR_df$end), linetype="dashed", color = "red", linewidth = 1.1)
g

# Gviz visualization ------------------------------------------------------
del_genes_df <- unique(assoc_gene_per_case_scna_df[, c("symbol", "transcript_start", "transcript_end")])
del_genes_gr <- GRanges(
    seqnames = "chr19",
    ranges = IRanges(start=del_genes_df$transcript_start, end=del_genes_df$transcript_end),
    symbol=del_genes_df$symbol
)

mcr_track <- AnnotationTrack(MCR_gr, name = "MCR", fill = "red")
atrack <- AnnotationTrack(
    assoc_del_segments_gr, 
    name = "copy-number loss segments",
    showFeatureId = TRUE, 
    id = assoc_del_segments_gr$patient_barcode
)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg38", chromosome = "chr19")
grTrack <- GeneRegionTrack(
    del_genes_gr, chromosome="chr19", name = "Assoc. Genes", 
    geneSymbol = TRUE, transcriptAnnotation = "symbol"
)

pdf("output/17.GViz_all_overlapping_segments.pdf", width = 10, height = 12)
plotTracks(list(itrack, gtrack, mcr_track, atrack, grTrack))
dev.off()

pdf("output/18.GViz_MCR_region.pdf", width = 8, height = 4)
plotTracks(list(itrack, gtrack, grTrack),
           from = MCR_df$start, to = MCR_df$end, 
           main = paste0(MCR_df$seqnames, ":", MCR_df$start, "-", MCR_df$end,
                         "\n(n=15 genes)"),
           cex.main = 1.2
)
dev.off()

# MCR genes - MYC expression  ---------------------------------------------
MCR_genes_df <- del_genes_df[del_genes_df$transcript_start >= MCR_df$start & del_genes_df$transcript_end <= MCR_df$end, ]
MCR_genes_vec <- MCR_genes_df$symbol

dds_norm <- readRDS("data/selected_data/processed_expr.RDS")

sel_expr <- SummarizedExperiment::assay(dds_norm[SummarizedExperiment::rowData(dds_norm)$gene_name %in% MCR_genes_vec, ])

anyDuplicated(SummarizedExperiment::rowData(dds_norm)$gene_name[match(rownames(sel_expr),SummarizedExperiment:: rowData(dds_norm)$gene_id)])
rownames(sel_expr) <- rowData(dds_norm)$gene_name[match(rownames(sel_expr), rowData(dds_norm)$gene_id)]

sel_expr_long <- reshape2::melt(sel_expr)
colnames(sel_expr_long) <- c("Gene", "Sample", "Expression")
sel_expr_long$microdel <- ifelse(sel_expr_long$Sample %in% altered_cases, "MCR del", "intact")

g <- ggviolin(sel_expr_long, x = "microdel", y = "Expression", color = "microdel", 
              add = "boxplot", facet.by = "Gene", palette = "npg")
g <- g + stat_compare_means(method = "wilcox.test", method.args = list("alternative" = "less"))
g <- g + rremove("legend")
g <- g + rremove("xlab")

ggsave("output/19.MCR_del_genes_expression.pdf", g, width = 10, height = 10)

write.csv(MCR_genes_df, "output/MCR_genes_df.csv")
writeLines(MCR_genes_vec, "output/MCR_genes.txt")
write.csv(MCR_df, "output/MCR_df.csv", row.names = FALSE)
