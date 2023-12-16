##### Script purpose: Analysis of COSMIC CN signatures by PMN class
##### Author: Ege Ulgen
##### Date: Dec 2023

library(ComplexHeatmap)

metadata_df <- readRDS("data/selected_data/meta.RDS")

CN_sig_counts_df <- read.delim("output/SigProfilerAssignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt")
CN_sig_counts_mat <- as.matrix(CN_sig_counts_df[, -1])
rownames(CN_sig_counts_mat) <- CN_sig_counts_df$Samples

CN_sig_mat <- t(apply(CN_sig_counts_mat, 1, function(x) x / sum(x)))

# order by PMN class
sample_class <- metadata_df$class[match(rownames(CN_sig_mat), metadata_df$patient)]
names(sample_class) <- rownames(CN_sig_mat)

# create annotated heatmap ------------------------------------------------
clu_method <- "average"

hm <- Heatmap(CN_sig_mat,
              heatmap_legend_param = list(title = "rel.freq.", title_position = "lefttop-rot"),
              split = sample_class,
              clustering_method_columns = clu_method,
              clustering_method_rows = clu_method,
              show_row_names = FALSE,
              col = colorRampPalette(c("#FCFBE8", "#6f0000"))(100))
draw(hm, annotation_legend_side = "bottom")

cairo_pdf("output/27.CN_signatures_clustering_heatmap.pdf")
draw(hm, annotation_legend_side = "bottom")
dev.off()
