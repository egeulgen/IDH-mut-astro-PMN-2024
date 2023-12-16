##### Script purpose: Cluster samples based on CNV profles
##### Author: Ege Ulgen
##### Date: Dec 2023

source("scripts/utils.R")
library(ComplexHeatmap)

metadata_df <- readRDS("data/selected_data/meta.RDS")

CNV_mat_df <- read.delim("output/PMN_project.CNV48.matrix.tsv")
CNV_mat <- as.matrix(CNV_mat_df[, -1])
rownames(CNV_mat) <- CNV_mat_df$MutationType
colnames(CNV_mat) <- gsub("\\.", '-', colnames(CNV_mat))
selected_samples <- metadata_df$patient[metadata_df$patient %in% colnames(CNV_mat)]
CNV_mat <- CNV_mat[, selected_samples]

CNV_mat <- apply(CNV_mat, 2, function(x) x / sum(x))

# create annotated heatmap ------------------------------------------------
clu_method <- "average"
# order by PMN class
sample_class <- metadata_df$class[match(colnames(CNV_mat), metadata_df$patient)]
names(sample_class) <- colnames(CNV_mat)

hm <- Heatmap(t(CNV_mat),
              heatmap_legend_param = list(title = "Rel.Freq.", title_position = "lefttop-rot"),
              split = sample_class,
              clustering_method_columns = clu_method,
              clustering_method_rows = clu_method,
              show_row_names = FALSE,
              col = colorRampPalette(c("#FCFBE8", "#6f0000"))(100))
draw(hm, annotation_legend_side = "bottom")

cairo_pdf("output/25.CNMat_clustering_heatmap.pdf")
draw(hm, annotation_legend_side = "bottom")
dev.off()

# filter features ---------------------------------------------------------
hist(apply(CNV_mat, 1, function(x) sum(x > 0.05)))

CNV_mat_filtered <- CNV_mat[rowSums(CNV_mat > 0.05) > 10, ]

row_ha <- rowAnnotation(PMN_class = sample_class,
                        col = list(PMN_class = c(WT = "#ffe74c",
                                                 microdel19q = "#38618c",
                                                 `PMN-hit` = "#ff5964")),
                        annotation_legend_param = list(PMN_class = list(nrow = 1, title = "")),
                        show_annotation_name = FALSE)

hm2 <- Heatmap(t(CNV_mat_filtered),
              heatmap_legend_param = list(title = "Rel.Freq.", title_position = "lefttop-rot"),
              clustering_method_columns = clu_method,
              clustering_method_rows = clu_method,
              show_row_names = FALSE,
              col = colorRampPalette(c("#FCFBE8", "#6f0000"))(100),
              right_annotation = row_ha)
draw(hm2, annotation_legend_side = "bottom")

cairo_pdf("output/26.CNMat_clustering_heatmap_filtered.pdf")
draw(hm, annotation_legend_side = "bottom")
dev.off()
