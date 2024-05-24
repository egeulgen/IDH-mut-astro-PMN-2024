##### Script purpose: Perform MC simulation to determine whether observed number of PMN SCNA are significant
##### Author: Ege Ulgen
##### Date: May 2024

library(data.table)

# prep --------------------------------------------------------------------
# proximal MYC network genes from Schaub et al. 2018
# "FBXW7" from Bai et al. 2016
PMN_genes <- read.csv("data/PMN_gene_list.csv")

gene_scna <- readRDS("data/selected_data/CN_gene_level.RDS")

gene_scna$SCNA_type <- ifelse(
    gene_scna$Segment_Mean > .3, "amp", 
    ifelse(gene_scna$Segment_Mean < -.3, "del", "CN")
)
pmn_scna_df <- gene_scna[gene_scna$symbol %in% PMN_genes$Gene.name, ]
pmn_scna_df <- pmn_scna_df[pmn_scna_df$SCNA_type != "CN", ]

# MC - permutation test ---------------------------------------------------
# How likely are we to observe this number of SCNA events if selected genes at random?
# Permutation of genes, selecting genes of the same size as the original set PMN genes
# observed to have SCNA events.
# stratifying by patients before selecting genes choosing a subset of the patients
# of the same size as the original group to ensure no patient-specific biases.

set.seed(123)

# use data.table for faster manipulation
gene_scna <- as.data.table(gene_scna)
pmn_scna_df <- as.data.table(pmn_scna_df)

selected_genes <- unique(pmn_scna_df$symbol)
selected_patients <- unique(pmn_scna_df$patient_barcode)

observed_count <- nrow(pmn_scna_df)

num_permutations <- 10000
permuted_counts <- numeric(num_permutations)

all_patients <- unique(gene_scna$patient_barcode)
selected_patient_size <- length(selected_patients)
selected_gene_size <- length(selected_genes)

for (i in seq_len(num_permutations)) {
    
    if (i %% 100 == 0)
        cat("Iteration", i, "out of", num_permutations, "             \r")
    
    random_patients <- sample(all_patients, selected_patient_size)
    patient_group <- gene_scna[patient_barcode %in% random_patients]
    
    random_genes <- sample(unique(patient_group$symbol), selected_gene_size)
    
    gene_group <- patient_group[symbol %in% random_genes & SCNA_type %in% c("amp", "del")]
    permuted_counts[i] <- nrow(gene_group)
}

hist(permuted_counts)
abline(v = observed_count, col = "red")

p_value <- sum(permuted_counts >= observed_count) / num_permutations

message("Observed SCNA count: ", observed_count)
message("P-value: ", p_value)

saveRDS(permuted_counts, "output/MC_results.RDS")
