##### Script purpose: Demonstrate colocalization decreased microdel19q genes + increased MYC expression
## using scRNAseq data (Venteicher et al. 2017 - GSE89567) 
##### Author: Ege Ulgen
##### Date: Dec 2023

library(Seurat)
library(ggplot2)
library(copykat)
library(scGSVA)

dir.create("output/Venteicher_scRNAseq", showWarnings = FALSE)

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

astro_obj <- CreateSeuratObject(counts = counts_mat,
                                min.cells = 3, min.features = 200,
                                meta.data = meta_df,
                                names.field = 1, names.delim = "\\.")

saveRDS(astro_obj, "data/Venteicher_scRNAseq/seurat_obj.RDS")

# process seurat object ---------------------------------------------------
astro_obj[["percent.mt"]] <- PercentageFeatureSet(astro_obj, pattern = "^MT-")
astro_obj <- subset(astro_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

astro_obj <- NormalizeData(astro_obj, verbose = FALSE)
astro_obj <- FindVariableFeatures(astro_obj, selection.method = "vst", nfeatures = 2000)
astro_obj <- ScaleData(astro_obj, verbose = FALSE)
astro_obj <- RunPCA(astro_obj, npcs = 30, verbose = FALSE)
astro_obj <- RunUMAP(astro_obj, reduction = "pca", dims = 1:30)
astro_obj <- FindNeighbors(astro_obj, reduction = "pca", dims = 1:30)
astro_obj <- FindClusters(astro_obj, reduction = "pca", resolution = 0.5)

g <- DimPlot(astro_obj, label = TRUE)
ggsave("output/Venteicher_scRNAseq/dimPlot.pdf", g, width = 5, height = 5)

# cell type ann ------------------------------------------------------------
## cell type annot. - https://www.frontiersin.org/articles/10.3389/fgene.2020.00490/full
# https://github.com/bioinfo-ibms-pumc/SCSA
all_markers <- FindAllMarkers(astro_obj)
colnames(all_markers)[colnames(all_markers) == "avg_log2FC"] <- "avg_logFC"
write.csv(all_markers, "scripts/tools/SCSA/all_markers.csv", row.names = FALSE)

# run SCSA
# cd scripts/tools/SCSA/
# poetry run python SCSA.py -d whole_v2.db -s seurat -i all_markers.csv -k All --Gensymbol -g Human -f 1 -p 0.01 --celltype cancer
# #Cluster Type Celltype Score Times
# ['0', 'Good', 'Macrophage', 10.360214672412098, 2.1633211942646833]
# ['1', '?', 'Fibroblast|Endothelial cell', '8.032428474469345|5.405594317546298', 1.4859473357807245]
# ['2', 'Good', 'Macrophage', 9.670698887642029, 2.031690502450226]
# ['3', '?', 'Cancer stem cell|Fibroblast', '6.388510769882183|3.2755558467200254', 1.950359288265322]
# ['4', '?', 'Cancer stem cell|Cancer cell', '6.66028175074541|4.243321611549725', 1.5695915512548138]
# ['5', 'Good', 'Macrophage', 9.39813692025507, 2.159848668014317]
# ['6', '?', 'Cancer stem cell|Mesenchymal cell', '5.967879826481657|4.072912594859697', 1.465260961900715]
# ['7', 'Good', 'Macrophage', 8.18091535385208, 2.2778441249043335]
# ['8', 'Good', 'Oligodendrocyte', 7.137451774667723, 2.01133816178245]

# poetry run python SCSA.py -d whole_v2.db -s seurat -i all_markers.csv -k All --Gensymbol -g Human -f 1 -p 0.01
# #Cluster Type Celltype Score Times
# ['0', '?', 'Microglial cell|Macrophage', '14.788337801684538|10.055857271520885', 1.4706193020028708]
# ['1', '?', 'Endothelial cell|Neuron', '10.884506391810927|8.507236881011398', 1.2794408506604207]
# ['2', '?', 'Microglial cell|Monocyte', '12.50486863098912|11.297462144054583', 1.1068741343444066]
# ['3', '?', 'Astrocyte|Oligodendrocyte', '7.171877977081825|6.644286468735449', 1.0794052921753068]
# ['4', '?', 'Neuron|Astrocyte', '11.74834475429808|7.183710734902439', 1.635414507605119]
# ['5', '?', 'Monocyte|Macrophage', '13.168146304774416|11.636039162420976', 1.131669129071122]
# ['6', 'Good', 'Astrocyte', 17.349134067736852, 8.227337010561945]
# ['7', '?', 'Microglial cell|Macrophage', '13.132553406797685|9.373210345126887', 1.4010731567146866]
# ['8', 'Good', 'Oligodendrocyte', 14.515628975024633, 6.069562588017648]

cluster_type_vec <- c("0" = "Normal cell", "1" = "Normal cell", "2" = "Normal cell", "3" = "Cancer cell", "4" = "Cancer cell",
                      "5"= "Normal cell", "6" = "Cancer cell", "7" = "Normal cell", "8" = "Cancer cell")

# CNV estimation ----------------------------------------------------------
exp.rawdata <- as.matrix(astro_obj[["RNA"]]$counts)

copykat_res <- copykat(
    rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, 
    sam.name="astro", distance="euclidean", norm.cell.names="", 
    output.seg="FLASE", plot.genes="FALSE", genome="hg20",n.cores=10
)

copykat_pred <- data.frame(copykat_res$prediction)
copykat_pred <- copykat_pred[copykat_pred$copykat.pred != "not.defined", ]  ##remove undefined cells

copykat_CNA<- data.frame(copykat_res$CNAmat)
copykat_CNA_chr19 <- copykat_CNA[copykat_CNA$chrom == 19, ]

selected_copykat_df <- copykat_CNA_chr19[copykat_CNA_chr19$chrompos >= 54950535 & copykat_CNA_chr19$chrompos <= 56859293, ]
# heatmap(ifelse(as.matrix(selected_copykat_df[, -c(1:3)]) < 0, -1, 1))
# selected_copykat_df <- selected_copykat_df[c("10997", "10998", "10999"), ]

D6_ave_copykat <- apply(selected_copykat_df[, -c(1,3)], 1, min)
astro_obj$D6_CNA_min_value <- D6_ave_copykat[colnames(astro_obj)]

(g1 <- FeaturePlot(astro_obj, features = c("D6_CNA_value", "MYC"), blend = TRUE, label = TRUE))

cancer_cells_obj <- subset(astro_obj, seurat_clusters %in% c(3, 4, 6))

(g2 <- FeaturePlot(cancer_cells_obj, features = c("D6_CNA_value", "MYC"), blend = TRUE, label = TRUE))


ggsave("output/Venteicher_scRNAseq/1.all_cells_D6_CNA_value_vs_MYC_expr_colocalization.pdf", g1, width = 12, height = 4)
ggsave("output/Venteicher_scRNAseq/1.cancer_cells_D6_CNA_value_vs_MYC_expr_colocalization.pdf", g2, width = 12, height = 4)


# coexpr.  MYC + microdel genes -------------------------------------------
del_peak_genes <- read.delim("output/PMN_neg_analysis/all_GISTIC_res/microdel19q/del_genes.conf_90.txt", skip = 3)
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
                   species = "human", 
                   anntype = "custom", 
                   keytype = "SYMBOL", 
                   annot = data.frame(GeneID = gset,
                                      PATH = "microdel_19q13.43",
                                      Annot = "microdel_19q13.43"))

set.seed(123)
res <- scgsva(as.ma, gset_custom, cores = 10)

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
