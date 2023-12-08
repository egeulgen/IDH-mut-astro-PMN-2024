##### Script purpose: Demonstrate colocalization decreased microdel19q genes + increased MYC expression
## using scRNAseq data (Venteicher et al. 2017 - GSE89567) 
##### Author: Ege Ulgen
##### Date: Dec 2023

library(Seurat)
library(scGSVA)

# create seurat object ----------------------------------------------------
counts_mat <- read.delim("data/Venteicher_scRNAseq/GSE89567_IDH_A_processed_data.txt.gz", row.names = 1)
samples_df <- read.csv("data/Venteicher_scRNAseq/metadata.csv")

# fix dimnames
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

astro_obj <- Createastro_object(counts = counts_mat,
                                min.cells = 3, min.genes = 200,
                                meta.data = meta_df,
                                names.field = 1, names.delim = "\\.")

saveRDS(astro_obj, "data/Venteicher_scRNAseq/seurat_obj.RDS")

# process seurat object ---------------------------------------------------
astro_obj <- readRDS("data/Venteicher_scRNAseq/seurat_obj.RDS")

astro_obj[["percent.mt"]] <- PercentageFeatureSet(astro_obj, pattern = "^MT-")
astro_obj <- subset(astro_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

astro_obj <- NormalizeData(astro_obj, verbose = FALSE)
astro_obj <- FindVariableFeatures(astro_obj, selection.method = "vst", nfeatures = 2000)
astro_obj <- ScaleData(astro_obj, verbose = FALSE)
astro_obj <- RunPCA(astro_obj, npcs = 30, verbose = FALSE)
astro_obj <- RunUMAP(astro_obj, reduction = "pca", dims = 1:30)
astro_obj <- FindNeighbors(astro_obj, reduction = "pca", dims = 1:30)
astro_obj <- FindClusters(astro_obj, reduction = "pca", resolution = 0.5)

pdf("output/Venteicher_scRNAseq/dimPlot.pdf", width = 5, height = 5)
DimPlot(astro_obj, label = TRUE)
dev.off()

# cell type ann ------------------------------------------------------------
## cell type annot. - https://www.frontiersin.org/articles/10.3389/fgene.2020.00490/full
# https://github.com/bioinfo-ibms-pumc/SCSA
all_markers <- FindAllMarkers(astro_obj)

colnames(all_markers)[colnames(all_markers) == "avg_log2FC"] <- "avg_logFC"
write.csv(all_markers, "data/Venteicher_scRNAseq/SCSA/all_markers.csv", row.names = FALSE)

# run SCSA
# python3 SCSA.py -d whole.db -s seurat -i all_markers.csv -k All --Gensymbol -g Human -f 1 -p 0.01 --celltype cancer
# #Cluster Type Celltype Score Times
# ['0', '?', 'Cancer stem cell|Regulatory T (Treg) cell', '3.0441142732384363|1.7489064755721913', 1.7405815095072428]
# ['1', '?', 'Exhausted CD8+ T cell|Regulatory T (Treg) cell', '0.8552031961088614|0.7172017101585542', 1.192416560077358]
# ['2', 'Good', 'Cancer stem cell', 4.304530680137537, 5.016102919495587]
# ['3', 'Good', 'Cancer stem cell', 1.4999954157198911, 3.0210686120140826]
# ['4', 'N', '-', '-', '-']
# ['5', 'Good', 'Cancer stem cell', 3.2768228560074557, 7.200487391599209]
# ['6', '?', 'Cancer cell|Cancer stem cell', '0.9252797318468855|0.8543580042003678', 1.0830117202599354] ***
# ['7', '?', 'Macrophage|Dendritic cell', '2.8347821263318003|2.210781179864086', 1.2822535998366316]
# ['8', '?', 'Cancer stem cell|Cancer cell', '2.3935892720884127|1.2156978955877964', 1.9689013864181277] *

# python3 SCSA.py -d whole.db -s seurat -i all_markers.csv -k All --Gensymbol -g Human -f 1 -p 0.01
# #Cluster Type Celltype Score Times
# ['0', '?', 'Microglial cell|Macrophage', '6.63345937089474|6.214304615141349', 1.0674499854307282]
# ['1', 'Good', 'Astrocyte', 5.458093829764412, 2.7443717305370625]
# ['2', '?', 'Microglial cell|Macrophage', '6.745305574315266|6.417575197272179', 1.0510676333301696]
# ['3', '?', 'Oligodendrocyte|Neural stem cell', '4.694209337985199|3.275531014455353', 1.433113995034402]
# ['4', 'Good', 'Oligodendrocyte', 2.041021903528525, 5.225224381424466]
# ['5', '?', 'Macrophage|Monocyte', '7.122096666849755|5.557320603821522', 1.281570234035482]
# ['6', 'Good', 'Astrocyte', 8.840872544795136, 7.433811011558405] ***
# ['7', '?', 'Microglial cell|Macrophage', '7.086374961494544|4.6507538532848525', 1.5237045831805094]
# ['8', '?', 'Oligodendrocyte|Astrocyte', '6.540427581726498|4.151135109474173', 1.5755756941756056] *




## CNV estimation? - https://www.nature.com/articles/s41467-019-13779-x
# https://statbiomed.github.io/SingleCell-Workshop-2021/CNV-analysis.html







# coexpr.  MYC + microdel genes -------------------------------------------
del_peak_genes <- read.delim("output/PMN_neg_analysis/del19q_samples_GISTIC/del_genes.conf_90.txt", skip = 3)
gset <- del_peak_genes$chr19.54940787.57124931
gset <- gset[gset != ""]

pdf("output/Venteicher_scRNAseq/MYC_and_del19q13.43_genes.pdf", width = 12, height = 4)
for (gene in gset) {
    res <- tryCatch({
        plot(FeaturePlot(astro_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()


pdf("output/Venteicher_scRNAseq/MYC_and_del19q13.43_genes_clu6.pdf", width = 12, height = 4)
clu6_obj <- subset(astro_obj, subset = seurat_clusters == 6)
for (gene in gset) {
    res <- tryCatch({
        plot(FeaturePlot(clu6_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()

# score gene set ----------------------------------------------------------
gset_custom <- new("Annot",
                   species = "mouse", 
                   anntype = "custom", 
                   keytype = "SYMBOL", 
                   annot = data.frame(GeneID = gset,
                                      PATH = "microdel_19q13.43",
                                      Annot = "microdel_19q13.43"))

set.seed(123)
res <- scgsva(astro_obj, gset_custom, cores = 10)

astro_obj <- AddMetaData(astro_obj, res@gsva)

pdf("output/Venteicher_scRNAseq/MYC_and_PMNneg_gene_set.pdf", width = 12, height = 4)
FeaturePlot(astro_obj, features = c("microdel_19q13.43", "MYC"), blend = TRUE, label = TRUE)
dev.off()

pdf("output/Venteicher_scRNAseq/by_clu_MYC_and_PMNneg_gene_set.pdf", width = 12, height = 4)
for (clu_idx in 0:8) {
    clu_obj <- subset(astro_obj, subset = seurat_clusters == clu_idx)
    res <- tryCatch({
        plot(FeaturePlot(clu_obj, features = c("microdel_19q13.43", "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()
