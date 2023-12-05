# PMN-latest - "Proximal MYC Network in IDH-mutant 1p/19q non-codeleted gliomas" Analysis

### 1. Data Preparation
1. [Obtain TCGA-LGG data and classify samples](scripts/1.data_prep/1.obtain_data.R)
2. [Prepare gene-level SCNA data](scripts/1.data_prep/2.prep_gene_level_SCNA.R)
3. [Select IDHmut 1p/19q non-codeleted Grade 2-3 samples](scripts/1.data_prep/3.astro_selection.R)
4. [Prepare genomic instability metrics](scripts/1.data_prep/4.prep_metrics.R)
5. [Detect chromotripsis events](scripts/1.data_prep/5.prep_CT.R)
6. [Prepare RNAseq data](scripts/1.data_prep/6.prep_RNAseq.R)

### 2. Analysis of Somatic Proximal MYC Network (PMN) Alterations
1. [Determine PMN alterations](scripts/2.PMN_alterations/1.PMN_alterations.R)
2. [Investigate MYC-associated methylation changes](scripts/2.PMN_alterations/2.MYC_assoc_methylation.R)
3. [Investigate associations of PMN alterations with various factors](scripts/2.PMN_alterations/3.PMN_assoc.R)
4. [Investigate associations of MYC expression level with various factors](scripts/2.PMN_alterations/4.expr_assoc.R)

### 3. Analysis of PMN-negative samples with increased MYC expression
1. [Determine genes associated with increased MYC expression](scripts/3.PMN_negative_analysis/1.PMN_neg_analysis.R)
2. [Annotate genes associated with increased MYC expression](scripts/3.PMN_negative_analysis/2.annotate_assoc_gene.R)
3. [Investigate associations of PMN classes with various factors](scripts/3.PMN_negative_analysis/3.PMN_neg_subclass_associations.R)
4. [Investigate MYC expression in 19q deleted vs. intact samples](scripts/3.PMN_negative_analysis/4.del19q_effect.R)
