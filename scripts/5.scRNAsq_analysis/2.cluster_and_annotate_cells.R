##### Script purpose: Cluster cells and annotate clusters
##### Author: Ege Ulgen
##### Date: Feb 2024

library(Seurat)
library(clustree)
library(ggplot2)

astro_obj <- readRDS("output/Venteicher_scRNAseq/astro_seurat_obj.RDS")
output_dir <- "output/Venteicher_scRNAseq"

# utility -----------------------------------------------------------------
#' Parse stdout log from SCSA annotation to get annotation table
#'
#' @param SCSA_stdout stdout log from SCSA annotation
#'
#' @return data frame of cell annotations
parse_SCSA_annotation_from_stdout <- function(SCSA_stdout) {
    result_flag <- FALSE
    columns <- c("Cluster", "Type", "Celltype", "Score", "Times")
    annot_df <- data.frame() 
    for (out_line in SCSA_stdout) {
        if (result_flag) {
            current_row <- unlist(strsplit(gsub("\\[|\\]|\\'", "", out_line), ", "))
            annot_df <- rbind(annot_df, current_row)
        }
        if (out_line == "#Cluster Type Celltype Score Times") {
            result_flag <- TRUE
        }
    }
    colnames(annot_df) <- columns
    return(annot_df)
}

# cluster cells -----------------------------------------------------------
# cluster cells at different resolutions
astro_obj <- FindClusters(
    astro_obj, 
    reduction = "pca", 
    dims = 1:30,
    resolution = seq(0.1, 2, 0.1), 
    algorithm = 1
)

# The stability index from the SC3 package (Kiselev et al. 2017) measures the 
# stability of clusters across resolutions. The stability index is automatically 
# calculated when a clustering tree is built. Note that each level of clustering 
# corresponds to a different resolution.
clustering_prefix <- "RNA_snn_res."
g_clustree <- clustree(astro_obj, prefix = clustering_prefix, node_colour = "sc3_stability")
ggsave(file.path(output_dir, "1.clustering_tree.pdf"), g_clustree, width = 8, height = 12)

# ad-hoc decision to set final cluster resolution to 0.6 upon visual inspection
final_cluster_resolution <- 0.6
Idents(astro_obj) <- paste0(clustering_prefix, final_cluster_resolution)

(g <- DimPlot(astro_obj, label = TRUE))
ggsave(
    file.path(output_dir, "2.clustered_umap_plot.pdf"), g, width = 5, height = 5
)

# cell type annotation of clusters ----------------------------------------
## cell type annotation using SCSA
# GitHub repo: https://github.com/bioinfo-ibms-pumc/SCSA
# Article for the method: https://www.frontiersin.org/articles/10.3389/fgene.2020.00490/full

#### installation of SCSA 
#### (this is done through the command line)
#### As the required python environment that contains the dependencies has been 
#### set up using `renv`, there is no need to install those
#### in your terminal, run the following to clone the SCSA repository before use:
## $ mkdir tools && cd tools/
## $ git clone https://github.com/bioinfo-ibms-pumc/SCSA.git

all_markers <- FindAllMarkers(astro_obj)
# need to change the logFC column name to the one expected by SCSA
colnames(all_markers)[colnames(all_markers) == "avg_log2FC"] <- "avg_logFC"
all_markers_file <- file.path(output_dir, "all_markers.csv")
write.csv(all_markers, all_markers_file, row.names = FALSE)

# thresholds for marker filtering 
foldchange_threshold <- 1
pvalue_threshold <- 0.01

python_script_path <- "tools/SCSA/SCSA.py"
common_arguments <- c(
    paste0("--input ", all_markers_file), 
    "--db tools/SCSA/whole_v2.db", 
    "--source seurat", 
    "--tissue All", 
    "--Gensymbol", 
    "--species Human",
    paste0("--foldchange ", foldchange_threshold),
    paste0("--pvalue ", pvalue_threshold)
)

cancer_annot_command <- paste(c("python", python_script_path, c(common_arguments, "--celltype cancer")), collapse = " ")
normal_annot_command <- paste(c("python", python_script_path, c(common_arguments, "--celltype normal")), collapse = " ")

SCSA_cancer_annotation <- system(cancer_annot_command, intern = TRUE)
SCSA_normal_annotation <- system(normal_annot_command, intern = TRUE)

(cluster_annot_cancer_types <- parse_SCSA_annotation_from_stdout(SCSA_cancer_annotation))
cluster_annot_normal_types <- parse_SCSA_annotation_from_stdout(SCSA_normal_annotation)

basic_cluster_type_vec <- ifelse(grepl("cancer", cluster_annot_cancer_types$Celltype, ignore.case = TRUE), "Cancer Cell", "Normal Cell")
names(basic_cluster_type_vec) <- cluster_annot_cancer_types$Cluster


saveRDS(cluster_annot_normal_types, file.path(output_dir, "cluster_annot_normal_types.RDS"))
saveRDS(cluster_annot_cancer_types, file.path(output_dir, "cluster_annot_cancer_types.RDS"))
saveRDS(basic_cluster_type_vec, file.path(output_dir, "basic_cluster_type_vec.RDS"))
saveRDS(astro_obj, file.path(output_dir, "astro_seurat_obj.RDS"))
