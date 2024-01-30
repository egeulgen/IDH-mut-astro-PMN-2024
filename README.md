# PMN-latest - "Proximal MYC Network in IDH-mutant 1p/19q non-codeleted gliomas" Analysis

## Prequisite
This repo/R project uses [renv](https://rstudio.github.io/renv/index.html) to create a reproducible environment. You may install the latest version of `renv` from CRAN via:

```r
install.packages("renv")
```

and then initialize the environment using:

```r
renv::init()
```

## Outline of analyses

This will initialize the environment and install all necessary dependencies to enable you to run the scripts.

Below is an outline of the analysis steps with links to related scripts:

### 1. Data Preparation
1. [Obtain TCGA-LGG data and classify samples](scripts/1.data_prep/1.obtain_data.R)
2. [Prepare gene-level SCNA data](scripts/1.data_prep/2.prep_gene_level_SCNA.R)
3. [Select IDHmut 1p/19q non-codeleted Grade 2-3 samples](scripts/1.data_prep/3.astro_selection.R)
4. [Prepare genomic instability metrics](scripts/1.data_prep/4.prep_metrics.R)
5. [Detect chromotripsis events](scripts/1.data_prep/5.prep_CT.R)
6. [Prepare RNAseq data](scripts/1.data_prep/6.prep_RNAseq.R)

### 2. Analysis of Somatic Proximal MYC Network (PMN) Alterations
1. [Determine PMN alterations](scripts/2.PMN_alterations/1.PMN_alterations.R)
2. [Investigate associations of PMN alterations with various factors](scripts/2.PMN_alterations/2.PMN_assoc.R)
3. [Investigate associations of MYC expression level with various factors](scripts/2.PMN_alterations/3.expr_assoc.R)

### 3. Analysis of PMN-negative tumors with increased MYC expression
1. [Determine genes associated with increased MYC expression](scripts/3.PMN_negative_analysis/1.PMN_neg_analysis.R)
2. [Annotate genes associated with increased MYC expression](scripts/3.PMN_negative_analysis/2.annotate_assoc_gene.R)
3. [Investigate associations of PMN classes with various factors](scripts/3.PMN_negative_analysis/3.PMN_neg_subclass_associations.R)
4. [Investigate MYC expression in 19q deleted vs. intact tumors](scripts/3.PMN_negative_analysis/4.del19q_effect.R)

### 4. Further analysis of chr19q microdeletion
1. [Generate oncoplot for PMN-no hit tumors with the chr19q microdeletion](scripts/4.microdeletion_analysis/1.microdel19q_tumors_oncoplot.R)
2. [Prepare and save seg files for GISTIC analyses](scripts/4.microdeletion_analysis/2.prep_GISTIC_inputs.R)
3. [Analyze chr19 peaks identified by GISTIC (on GenePattern) in microdel19q cases](scripts/4.microdeletion_analysis/3.post_GISTIC.R)
4. [Determine deletion status for the 19q13.43 peak identified by GISTIC in all samples](scripts/4.microdeletion_analysis/4.D6_in_all_samples.R)
5. [Demonstrate colocalization of lower CNA value/decreased microdel19q genes + increased MYC expression in single-cell RNA sequencing experiment](scripts/4.microdeletion_analysis/5.scRNAseq_colocalization.R)
