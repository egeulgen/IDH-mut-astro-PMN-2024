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

This will initialize the environment and install all necessary dependencies to enable you to run the scripts.

## Outline of analyses

Below is an outline of the analysis steps with links to related scripts:

### 1. Data Preparation
1. [Obtain TCGA-LGG data and classify tumours (based on IDH mutation and 1p/19q codeletion status)](scripts/1.data_prep/1.obtain_data.R)
2. [Prepare gene-level SCNA data](scripts/1.data_prep/2.prep_gene_level_SCNA.R)
3. [Select lower grade Astrocytoma (IDH-mutant 1p/19q non-codeleted Grade 2-3) cases](scripts/1.data_prep/3.astro_selection.R)
4. [Prepare genomic instability metrics](scripts/1.data_prep/4.prep_metrics.R)
5. [Detect chromotripsis events](scripts/1.data_prep/5.prep_CT.R)
6. [Prepare RNAseq data](scripts/1.data_prep/6.prep_RNAseq.R)

### 2. Analysis of Somatic Proximal MYC Network (PMN) Alterations
1. [Determine PMN alterations](scripts/2.PMN_alterations/1.PMN_alterations.R)
2. [Investigate associations of PMN alterations with various factors](scripts/2.PMN_alterations/2.PMN_assoc.R)
3. [Investigate associations of MYC expression level with various factors](scripts/2.PMN_alterations/3.expr_assoc.R)

### 3. Analysis of PMN-negative tumors with increased MYC expression
1. [Determine altered genes associated with increased MYC expression](scripts/3.PMN_negative_analysis/1.PMN_neg_analysis.R)
2. [Annotate genes associated with increased MYC expression](scripts/3.PMN_negative_analysis/2.annotate_assoc_genes.R)
3. [Analyze SCNA segments overlappping altered genes associated with increased MYC expression to define a Minimal Common Region (MCR)](scripts/3.PMN_negative_analysis/3.define_minimal_common_region.R)
4. [Investigate MYC expression in 19q deleted vs. intact tumors](scripts/3.PMN_negative_analysis/4.del19q_MYC_comparison.R)

### 4. Further analysis related to chr19q MCR deletion
1. [Generate oncoplot for PMN-no hit tumors with the chr19q microdeletion](scripts/4.MCR_analysis/1.MCR_tumors_oncoplot.R)
2. [Investigate associations of PMN sub-classes with various factors](scripts/4.MCR_analysis/2.MCR_in_all_samples.R)
3. [Determine deletion status for the 19q MCR in all samples](scripts/4.MCR_analysis/3.subclass_associations.R)

### 5. single-cell RNAseq analysyis
1. [Obtain Data, create processed Seurat object](scripts/5.scRNAsq_analysis/1.obtain_data_and_process.R)
2. [Cluster cells and annotate clusters](scripts/5.scRNAsq_analysis/2.cluster_and_annotate_cells.R)
3. [Estimate copy-number in in (220KB)(play with `win.size`) windows](scripts/5.scRNAsq_analysis/3.estimate_CN.R)
4. [Visualize CN estimate of MCR window(s)/genes vs. MYC expression in the same cells](scripts/5.scRNAsq_analysis/4.colocalization_visualization.R)
