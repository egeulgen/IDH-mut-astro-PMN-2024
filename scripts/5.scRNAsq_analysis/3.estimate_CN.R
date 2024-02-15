##### Script purpose: Estimate copy-number using copykat
##### Author: Ege Ulgen
##### Date: Feb 2024

library(Seurat)
library(copykat)
library(ggplot2)
library(DescTools)
library(pheatmap)

output_dir <- "output/Venteicher_scRNAseq"
astro_obj <- readRDS(file.path(output_dir, "astro_seurat_obj.RDS"))

basic_cluster_type_vec <- readRDS(file.path(output_dir, "basic_cluster_type_vec.RDS"))
normal_clusters <- names(basic_cluster_type_vec)[basic_cluster_type_vec == "Normal Cell"]
cancer_clusters <- names(basic_cluster_type_vec)[basic_cluster_type_vec == "Cancer Cell"]

# CNV estimation ----------------------------------------------------------
exp.rawdata <- as.matrix(astro_obj[["RNA"]]$counts)
normal_cell_names <- names(Idents(astro_obj))[Idents(astro_obj) %in% normal_clusters]

copykat_res <- copykat(
    rawmat = exp.rawdata,
    id.type = "S", 
    ngene.chr = 5, 
    win.size = 25, 
    KS.cut = 0.1, 
    sam.name = "astro", 
    distance = "euclidean", 
    norm.cell.names = normal_cell_names, 
    output.seg = FALSE, 
    plot.genes = FALSE, 
    genome = "hg20",
    n.cores = 10
)

saveRDS(copykat_res, file.path(output_dir, "copykat_res.RDS"))


# determine copykat segments overlapping the MCR --------------------------
MCR_df <- read.csv("output/MCR_df.csv")
copykat_CNA <- data.frame(copykat_res$CNAmat)
copykat_CNA_chr19 <- copykat_CNA[copykat_CNA$chrom == 19, ]

selected_copykat_intervals <- copykat_CNA_chr19[, "chrompos", drop = FALSE]
selected_copykat_intervals$end <- selected_copykat_intervals$chrompos
for (i in seq(2, nrow(selected_copykat_intervals))) {
    selected_copykat_intervals$end[i - 1] <- selected_copykat_intervals$chrompos[i]
}
selected_copykat_intervals$id <- rownames(selected_copykat_intervals)

selected_copykat_intervals <- selected_copykat_intervals[selected_copykat_intervals$chrompos >= MCR_df$start - 1e6 & selected_copykat_intervals$end <= MCR_df$end + 1.e6, ]

g <- ggplot(selected_copykat_intervals, aes(
    x = chrompos, xend = end, y = id, yend = id, color = id)
) 
g <- g + geom_segment(linewidth = 1)
g <- g + scale_x_continuous(
    "Position", 
    limits = c(min(selected_copykat_intervals$chrompos), max(selected_copykat_intervals$end)),
    labels = scales::comma
)
g <- g + theme_minimal()
g <- g + scale_color_manual(values = colorRampPalette(c("#00264D","#5C93D1"))(220))
g <- g + theme(
    axis.text.y = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
)
g <- g + geom_segment(data = MCR_df, aes(x = start, xend = end, y = min(selected_copykat_intervals$id), yend = min(selected_copykat_intervals$id)), color = "red")
g <- g + geom_vline(xintercept = c(MCR_df$start, MCR_df$end), color = "red", linetype = "dashed")
g

ggsave(file.path(output_dir, "copykat_intervals_overlapping_MCR.pdf"), g, width = 10, height = 10)

selected_copykat_intervals$any_overlap_with_MCR <- FALSE
for (i in seq_len(nrow(selected_copykat_intervals))) {
    selected_copykat_intervals$any_overlap_with_MCR[i] <- c(selected_copykat_intervals$chrompos[i], selected_copykat_intervals$end[i]) %overlaps% c(MCR_df$start, MCR_df$end)
}
MCR_overlapping_copykat_intervals <- selected_copykat_intervals[selected_copykat_intervals$any_overlap_with_MCR, ]

saveRDS(MCR_overlapping_copykat_intervals, file.path(output_dir, "MCR_overlapping_copykat_intervals.RDS"))

# assess MCR deletion status in cells/samples -----------------------------
copy_kat_CNA_MCR_segments <- copykat_CNA_chr19[copykat_CNA_chr19$chrompos %in% MCR_overlapping_copykat_intervals$chrompos, ]

copy_kat_MCR_segments_mat <- as.matrix(copy_kat_CNA_MCR_segments[, -c(1:3)]) 

cell_sample_ids <- sapply(colnames(copy_kat_MCR_segments_mat), function(x) unlist(strsplit(x, "\\."))[1])
copy_kat_MCR_segments_mat <- copy_kat_MCR_segments_mat[, order(cell_sample_ids)]

annotation_col <- data.frame(Sample = as.factor(cell_sample_ids))
sample_ids <- unique(cell_sample_ids)
sample_cols <- RColorBrewer::brewer.pal(length(sample_ids), "Set1")
names(sample_cols) <- sample_ids
annotation_colors <- list(Sample = sample_cols)

pheatmap(
    copy_kat_MCR_segments_mat, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    filename = file.path(output_dir, "5.copykat_MCR_segments_per_sample_heatmap.pdf"),
    width = 8,
    height = 3
)

median_expr_mat <- matrix(, nrow = nrow(copy_kat_CNA_MCR_segments), ncol = 0)
rownames(median_expr_mat) <- rownames(copy_kat_CNA_MCR_segments)
for (sample_id in sample_ids) {
    sample_cells <- names(cell_sample_ids)[cell_sample_ids == sample_id]
    sample_mat <- copy_kat_MCR_segments_mat[, sample_cells]
    
    pheatmap(
        sample_mat, 
        cluster_rows = FALSE, 
        show_colnames = FALSE,
        filename = file.path(output_dir, paste0(sample_id, "_copykat_MCR_segments_heatmap.pdf")),
        width = 8,
        height = 3
    )
    
    median_expr_mat <- cbind(
        median_expr_mat, 
        data.frame(apply(sample_mat, 1, median))
    )
    colnames(median_expr_mat)[ncol(median_expr_mat)] <- sample_id
}

pheatmap(
    median_expr_mat, 
    cluster_rows = FALSE, 
    filename = file.path(output_dir, "copykat_MCR_segments_median_expr_by_sample_heatmap.pdf"),
    width = 8,
    height = 3
)
