##### Script purpose: Estimate copy-number using copykat
##### Author: Ege Ulgen
##### Date: Feb 2024

library(Seurat)
library(copykat)

astro_obj <- readRDS("output/Venteicher_scRNAseq/astro_seurat_obj.RDS")
output_dir <- "output/Venteicher_scRNAseq"

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
