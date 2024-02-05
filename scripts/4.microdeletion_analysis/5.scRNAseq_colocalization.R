##### Script purpose: Demonstrate colocalization decreased microdel19q genes + increased MYC expression
## using scRNAseq data (Venteicher et al. 2017 - GSE89567) 
##### Author: Ege Ulgen
##### Date: Feb 2023

library(GEOquery)
library(Seurat)
library(clustree)
library(copykat)
library(scGSVA)
library(ggplot2)

source("scripts/utils.R")

output_dir <- "output/Venteicher_scRNAseq"
dir.create(output_dir, showWarnings = FALSE)

# get data from GEO -------------------------------------------------------
# GSE89567: IDH-mutant astrocytoma tumors, 6341 cells from 10 tumors
# Venteicher AS, Tirosh I, Hebert C, Yizhak K et al. Decoupling genetics, lineages, 
# and microenvironment in IDH-mutant gliomas by single-cell RNA-seq. Science 
# 2017 Mar 31;355(6332). PMID: 28360267
data_dir <- "data/Venteicher_scRNAseq_data"
dir.create(data_dir, showWarnings = FALSE)

data_path <- getGEOSuppFiles(
    GEO = "GSE89567",
    baseDir = data_dir
)

# create seurat object ----------------------------------------------------
counts_mat <- read.delim(rownames(data_path), row.names = 1)
# meta data extracted from Table S1 from Venteicher et al. - https://www.science.org/doi/10.1126/science.aai8478#supplementary-materials 
samples_df <- read.csv(file.path(data_dir, "metadata.csv"))

# fix dim. names
rownames(counts_mat) <- gsub("'", "", rownames(counts_mat))
colnames(counts_mat) <- gsub("_", "\\.", colnames(counts_mat))
colnames(counts_mat) <- sub("X57", "MGH57", colnames(counts_mat))
colnames(counts_mat) <- toupper(colnames(counts_mat))
colnames(counts_mat) <- sub("POS", "", colnames(counts_mat))
colnames(counts_mat) <- sub("NEG", "", colnames(counts_mat))

# only keep LGG
counts_mat <- counts_mat[, sub("\\..*", "", colnames(counts_mat)) %in% samples_df$Designation]

meta_df <- samples_df[match(sub("\\..*", "", colnames(counts_mat)), samples_df$Designation), ]
rownames(meta_df) <- colnames(counts_mat)

astro_obj <- CreateSeuratObject(counts = counts_mat,
                                min.cells = 3, min.features = 200,
                                meta.data = meta_df,
                                names.field = 1, names.delim = "\\.")

# process seurat object ---------------------------------------------------
astro_obj[["percent.mt"]] <- PercentageFeatureSet(astro_obj, pattern = "^MT-")
astro_obj <- subset(astro_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

astro_obj <- NormalizeData(astro_obj, verbose = FALSE)
astro_obj <- FindVariableFeatures(astro_obj, selection.method = "vst", nfeatures = 2000)
astro_obj <- ScaleData(astro_obj, verbose = FALSE)
astro_obj <- RunPCA(astro_obj, npcs = 30, verbose = FALSE)
astro_obj <- RunUMAP(astro_obj, reduction = "pca", dims = 1:30)
astro_obj <- FindNeighbors(astro_obj, reduction = "pca", dims = 1:30)


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

copykat_CNA <- data.frame(copykat_res$CNAmat)
copykat_CNA_chr19 <- copykat_CNA[copykat_CNA$chrom == 19, ]
selected_copykat_df <- copykat_CNA_chr19[copykat_CNA_chr19$chrompos >= 54950535 & copykat_CNA_chr19$chrompos <= 56859293, ]


## using keep.scale = “feature” (default; by row/feature scaling) for the blended feature plots: 
# The plots for each individual feature are scaled to the maximum expression of 
# the feature across the conditions provided to split.by
overall_list <- list()
cancer_cells_list <- list()
for (i in seq_len(nrow(selected_copykat_df))) {
    region <- paste0(selected_copykat_df$chrom[i], ":", selected_copykat_df$chrompos[i])
    
    astro_obj[[region]] <- as.numeric(selected_copykat_df[i, colnames(astro_obj)])
    cancer_cells_obj <- subset(astro_obj, idents = cancer_clusters)
    
    overall_list[[region]] <- FeaturePlot(astro_obj, features = c(region, "MYC"), blend = TRUE, label = TRUE)
    cancer_cells_list[[region]] <- FeaturePlot(cancer_cells_obj, features = c(region, "MYC"), blend = TRUE, label = TRUE)
}
    
g1 <- ggpubr::ggarrange(plotlist = overall_list, ncol = 1)
g2 <- ggpubr::ggarrange(plotlist = cancer_cells_list, ncol = 1)

ggsave(file.path(output_dir, "3.all_cells_D6_CNA_value_vs_MYC_expr_colocalization.pdf"), g1, width = 12, height = 25)
ggsave(file.path(output_dir, "4.cancer_cells_D6_CNA_value_vs_MYC_expr_colocalization.pdf"), g2, width = 12, height = 25)

# coexpr.  MYC + microdel genes -------------------------------------------
del_peak_genes <- read.delim("output/PMN_neg_analysis/all_GISTIC_res/microdel19q/del_genes.conf_90.txt", skip = 3)
gset <- del_peak_genes$chr19.54940787.57124931
gset <- gset[gset != ""]

pdf(file.path(output_dir,"MYC_and_del19q13.43_genes.pdf"), width = 12, height = 4)
for (gene in gset) {
    res <- tryCatch({
        plot(FeaturePlot(astro_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()


pdf(file.path(output_dir, "MYC_and_del19q13.43_genes_cancer_clus.pdf"), width = 12, height = 4)
clu6_obj <- subset(astro_obj, idents = cancer_clusters)
for (gene in gset) {
    res <- tryCatch({
        plot(FeaturePlot(clu6_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()

# score gene set ----------------------------------------------------------
gset_custom <- new("Annot",
                   species = "human", 
                   anntype = "custom", 
                   keytype = "SYMBOL", 
                   annot = data.frame(GeneID = gset,
                                      PATH = "microdel_19q13.43",
                                      Annot = "microdel_19q13.43"))

set.seed(123)
res <- scgsva(as.matrix(counts_mat), gset_custom, cores = 10)

astro_obj <- AddMetaData(astro_obj, res@gsva)

pdf(file.path(output_dir, "MYC_and_PMNneg_gene_set.pdf"), width = 12, height = 4)
FeaturePlot(astro_obj, features = c("microdel_19q13.43", "MYC"), blend = TRUE, label = TRUE)
dev.off()

pdf(file.path(output_dir, "by_clu_MYC_and_PMNneg_gene_set.pdf"), width = 12, height = 4)
for (clu_idx in 0:8) {
    clu_obj <- subset(astro_obj, subset = seurat_clusters == clu_idx)
    res <- tryCatch({
        plot(FeaturePlot(clu_obj, features = c("microdel_19q13.43", "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()
