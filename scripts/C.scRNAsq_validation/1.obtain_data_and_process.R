##### Script purpose: Retrieve scRNAseq data (Venteicher et al. 2017 - GSE89567) 
##### and process it into a Seurat object
##### Author: Ege Ulgen
##### Date: Feb 2024

library(GEOquery)
library(Seurat)

output_dir <- "output/Venteicher_scRNAseq"
dir.create(output_dir, showWarnings = FALSE)

# get data from GEO -------------------------------------------------------
# GSE89567: IDH-mutant astrocytoma tumors, 6341 cells from 10 tumors
# Venteicher AS, Tirosh I, Hebert C, Yizhak K et al. Decoupling genetics, lineages, 
# and microenvironment in IDH-mutant gliomas by single-cell RNA-seq. Science 
# 2017 Mar 31;355(6332). PMID: 28360267
data_dir <- "data/Venteicher_scRNAseq_data"
dir.create(data_dir, showWarnings = FALSE)

data_manifest <- getGEOSuppFiles(
    GEO = "GSE89567",
    baseDir = data_dir
)

data_path <- rownames(data_manifest)

# create seurat object ----------------------------------------------------
# data_path <- "data/Venteicher_scRNAseq_data/GSE89567/GSE89567_IDH_A_processed_data.txt.gz"
counts_mat <- read.delim(data_path, row.names = 1)
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

saveRDS(astro_obj, file.path(output_dir, "astro_seurat_obj.RDS"))
