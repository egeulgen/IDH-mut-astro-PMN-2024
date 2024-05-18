##### Script purpose: Visualize CN estimate of MCR window(s)/genes vs. MYC expression in the same cells
##### Author: Ege Ulgen
##### Date: Feb 2024

library(Seurat)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(scGSVA)

output_dir <- "output/Venteicher_scRNAseq"
astro_obj <- readRDS(file.path(output_dir, "astro_seurat_obj.RDS"))

basic_cluster_type_vec <- readRDS(file.path(output_dir, "basic_cluster_type_vec.RDS"))
normal_clusters <- names(basic_cluster_type_vec)[basic_cluster_type_vec == "Normal Cell"]
cancer_clusters <- names(basic_cluster_type_vec)[basic_cluster_type_vec == "Cancer Cell"]

copykat_res <- readRDS(file.path(output_dir, "copykat_res.RDS"))

copykat_CNA <- data.frame(copykat_res$CNAmat)
copykat_CNA_chr19 <- copykat_CNA[copykat_CNA$chrom == 19, ]
MCR_overlapping_copykat_intervals <- readRDS(file.path(output_dir, "MCR_overlapping_copykat_intervals.RDS"))


# workaround for Seurat:::BlendExpression issue ---------------------------

fixed_BlendExpression <- function(data) {
    if (ncol(x = data) != 2) {
        stop("'BlendExpression' only blends two features")
    }
    features <- colnames(x = data)
    data <- as.data.frame(x = apply(
        X = data,
        MARGIN = 2,
        FUN = function(x) {
            if (min(x) == max(x)) {
                return(rep(0, length(x)))
            }
            return(round(x = 9 * (x - min(x)) / (max(x) - min(x))))
        }
    ))
    data[, 3] <- data[, 1] + data[, 2] * 10
    colnames(x = data) <- c(features, paste(features, collapse = '_'))
    for (i in 1:ncol(x = data)) {
        data[, i] <- factor(x = data[, i])
    }
    return(data)
}

assignInNamespace("BlendExpression", fixed_BlendExpression, ns="Seurat", pos="package:Seurat")
 

# split cells into any MCRdel vs no del -----------------------------------
selected_copykat_df <- copykat_CNA_chr19[copykat_CNA_chr19$chrompos %in% MCR_overlapping_copykat_intervals$chrompos, ]
selected_copykat_df$end <- MCR_overlapping_copykat_intervals$end[match(MCR_overlapping_copykat_intervals$chrompos, selected_copykat_df$chrompos)]

selected_copykat_mat <- selected_copykat_df[, -c(1:3)] 

selected_copykat_del_mask <- selected_copykat_mat < 0
cells_with_any_mcr_del <- names(which(colSums(selected_copykat_del_mask) > 5))

astro_obj[["any_MCR_deletion"]] <- ifelse(colnames(astro_obj) %in% cells_with_any_mcr_del, "with MCR deletion", "without MCR deletion")

# copykat intervals within MCR vs. MYC ------------------------------------

## using keep.scale = “feature” (default; by row/feature scaling) for the blended feature plots
## so that plots are comparable: 
# The plots for each individual feature are scaled to the maximum expression of 
# the feature across the conditions provided to split.by
cancer_plot_list <- all_plot_list <- list()
for (i in seq_len(nrow(selected_copykat_df))) {
    region <- paste0("copykat_region_", rownames(selected_copykat_df)[i])
    
    astro_obj[[region]] <- as.numeric(selected_copykat_df[i, colnames(astro_obj)])
    cancer_cells_obj <- subset(astro_obj, idents = cancer_clusters)
    
    plot_all <- FeaturePlot(astro_obj, features = c(region, "MYC"), blend = TRUE, split.by = "any_MCR_deletion")
    plot_cancer <- FeaturePlot(cancer_cells_obj, features = c(region, "MYC"), blend = TRUE, split.by = "any_MCR_deletion")
    all_plot_list[[region]] <- plot_all
    cancer_plot_list[[region]] <- plot_cancer
}

g1 <- cowplot::plot_grid(plotlist = all_plot_list, ncol = 3)
g2 <- cowplot::plot_grid(plotlist = cancer_plot_list, ncol = 3)

ggsave(file.path(output_dir, "4.MCR_CNA_value_vs_MYC_expr_colocalization.pdf"), g1, width = 36, height = 18)
ggsave(file.path(output_dir, "5.cancer_cells_MCR_CNA_value_vs_MYC_expr_colocalization.pdf"), g2, width = 36, height = 18)

# gene CN estimates vs. MYC -----------------------------------------------
MCR_genes_df <- read.csv("output/MCR_genes_df.csv", row.names = 1)
MCR_genes_gr <- GRanges(
    seqnames = "chr19",
    ranges = IRanges(start=MCR_genes_df$transcript_start, end=MCR_genes_df$transcript_end),
    symbol=MCR_genes_df$symbol
)

MCR_overlapping_copykat_intervals_gr <- makeGRangesFromDataFrame(
    selected_copykat_df, 
    seqnames.field = "chrom", 
    start.field = "chrompos", 
    end.field = "end", 
    keep.extra.columns = TRUE
)
seqlevelsStyle(MCR_overlapping_copykat_intervals_gr) <- "UCSC"

hits <- findOverlaps(MCR_overlapping_copykat_intervals_gr, MCR_genes_gr, type = "any", select = "all")

ranges <- MCR_genes_gr[subjectHits(hits)]
mcols(ranges) <- c(mcols(ranges), mcols(MCR_overlapping_copykat_intervals_gr[queryHits(hits)]))

genes_copykat_df <- data.frame(chr = as.vector(seqnames(ranges)), as.data.frame(mcols(ranges)))
genes_copykat_df <- genes_copykat_df %>%
    mutate(segment_start = as.integer(start(ranges(MCR_overlapping_copykat_intervals_gr[queryHits(hits)])))) %>%
    mutate(segment_end = as.integer(end(ranges(MCR_overlapping_copykat_intervals_gr[queryHits(hits)])))) %>%
    mutate(transcript_start = start(ranges)) %>%
    mutate(transcript_end = end(ranges)) %>%
    mutate(chrom = as.character(seqnames(ranges))) %>%
    rowwise() %>%
    mutate(
        transcript_overlap_percent = round(
            as.numeric((min(transcript_end, segment_end) - max(segment_start, transcript_start))/(transcript_end - transcript_start)) * 100, 
            digits = 2)
    )
genes_copykat_df <- as.data.frame(genes_copykat_df)

genes_copykat_df <- genes_copykat_df[genes_copykat_df$transcript_overlap_percent >= 25, ]

overall_list_genes <- list()
cancer_cells_list_genes <- list()
for (i in seq_len(nrow(genes_copykat_df))) {
    gene <- paste(genes_copykat_df$symbol[i], "_CN_estimate")
    
    astro_obj[[gene]] <- as.numeric(genes_copykat_df[i, colnames(astro_obj)])
    cancer_cells_obj <- subset(astro_obj, idents = cancer_clusters)
    
    overall_list_genes[[gene]] <- FeaturePlot(astro_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE)
    cancer_cells_list_genes[[gene]] <- FeaturePlot(cancer_cells_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE)
}

g3 <- ggpubr::ggarrange(plotlist = overall_list_genes, ncol = 1)
g4 <- ggpubr::ggarrange(plotlist = cancer_cells_list_genes, ncol = 1)

ggsave(file.path(output_dir, "6.all_cells_MCR_gene_CNA_value_vs_MYC_expr_colocalization.pdf"), g3, width = 12, height = 40)
ggsave(file.path(output_dir, "7.cancer_cells_MCR_gene_CNA_value_vs_MYC_expr_colocalization.pdf"), g4, width = 12, height = 40)

# coexpr.  MYC + microdel genes -------------------------------------------
gset <- readLines("output/MCR_genes.txt")

pdf(file.path(output_dir,"MYC_and_MCR_genes_expr.pdf"), width = 12, height = 4)
for (gene in gset) {
    res <- tryCatch({
        plot(FeaturePlot(astro_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()


pdf(file.path(output_dir, "MYC_and_MCR_genes_expr_cancer_cells.pdf"), width = 12, height = 4)
cancer_cells_obj <- subset(astro_obj, idents = cancer_clusters)
for (gene in gset) {
    res <- tryCatch({
        plot(FeaturePlot(cancer_cells_obj, features = c(gene, "MYC"), blend = TRUE, label = TRUE))
    },
    error = function(e) {return(NA)})
}
dev.off()

# score gene set ----------------------------------------------------------
counts_mat <- astro_obj[["RNA"]]$counts
gset_custom <- new("Annot",
                   species = "human", 
                   anntype = "custom", 
                   keytype = "SYMBOL", 
                   annot = data.frame(GeneID = gset,
                                      PATH = "MCR_del_genes",
                                      Annot = "MCR_del_genes"))

set.seed(123)
res <- scgsva(as.matrix(counts_mat), gset_custom, cores = 10)

astro_obj <- AddMetaData(astro_obj, res@gsva)

pdf(file.path(output_dir, "MYC_and_MCR_gene_set.pdf"), width = 12, height = 4)
FeaturePlot(astro_obj, features = c("MCR_del_genes", "MYC"), blend = TRUE, label = TRUE)
dev.off()

pdf(file.path(output_dir, "cancer_cells_MYC_and_MCR_gene_set.pdf"), width = 15, height = 4)
cancer_cells_obj <- subset(astro_obj, idents = cancer_clusters)
plot(FeaturePlot(cancer_cells_obj, features = c("MCR_del_genes", "MYC"), blend = TRUE, label = TRUE))
dev.off()
