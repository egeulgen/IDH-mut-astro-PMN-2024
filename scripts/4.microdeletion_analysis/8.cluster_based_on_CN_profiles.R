##### Script purpose: Cluster samples based on CNV profles
##### Author: Ege Ulgen
##### Date: Dec 2023

library(ComplexHeatmap)

metadata_df <- readRDS("data/selected_data/meta.RDS")
# metadata_df <- metadata_df[metadata_df$class %in% c("WT", "microdel19q"), ]
# metadata_df <- metadata_df[metadata_df$microdel_19q13.43, ]

CNV_mat_df <- read.delim("output/PMN-latest.CNV48.matrix.tsv")
CNV_mat <- as.matrix(CNV_mat_df[, -1])
rownames(CNV_mat) <- CNV_mat_df$MutationType
colnames(CNV_mat) <- gsub("\\.", '-', colnames(CNV_mat))
selected_samples <- metadata_df$patient[metadata_df$patient %in% colnames(CNV_mat)]
CNV_mat <- CNV_mat[, selected_samples]

# prepare similarity matrix -----------------------------------------------------
RMSE_mat <- c()
for (i in 1:ncol(CNV_mat)) {
    cur_row <- c()
    for (j in 1:ncol(CNV_mat)) {
        cur_row <- c(cur_row, lsa::cosine(CNV_mat[, i], CNV_mat[, j]))
    }
    RMSE_mat <- rbind(RMSE_mat, cur_row)
}
colnames(RMSE_mat) <- rownames(RMSE_mat) <- colnames(CNV_mat)


# create annotated heatmap ------------------------------------------------
clu_method <- "centroid"
sample_class <- metadata_df$class[match(colnames(RMSE_mat), metadata_df$patient)]
names(sample_class) <- colnames(RMSE_mat)
column_ha <- HeatmapAnnotation(PMN_class = sample_class,
                               col = list(PMN_class = c(WT = "#ffe74c",
                                                       microdel19q = "#38618c",
                                                       `PMN-hit` = "#ff5964")),
                               annotation_legend_param = list(PMN_class = list(nrow = 1, title = "")),
                               show_annotation_name = FALSE)
hm <- Heatmap(RMSE_mat,
              heatmap_legend_param = list(title = "Cosine", title_position = "lefttop-rot"),
              clustering_method_columns = clu_method,
              clustering_method_rows = clu_method,
              show_column_dend = FALSE,
              col = colorRampPalette(c("#FCFBE8", "#6f0000"))(100),
              bottom_annotation = column_ha)
draw(hm, annotation_legend_side = "bottom")

cairo_pdf("output/25.CNMat_clustering_heatmap.pdf")
draw(hm, annotation_legend_side = "bottom")
dev.off()
