# ALK_rearranged_pADCs_multiomics

**This repository contains the custom R code used for the analyses reported in:**

**Szeitz B, et al. Spatially resolved proteomic and transcriptomic profiling of ALK-rearranged lung adenocarcinomas reveals key players in inter- and intratumoral heterogeneity. Manuscript submitted for publication.**

## Study summary

In this study, we performed the spatial molecular profiling of 7 ALK-rearranged lung adenocarcinoma (pADC) samples using GeoMx digital spatial profiling ( = 84 smaller transcriptomic regions, tROIs) and label-free shotgun proteomics ( = 23 larger proteomic regions, pROIs). With this, we seek a better understanding of the inter- and intratumoral heterogeneity in ALK-driven pDACs.

## Scripts

1.Data cleaning and transformation scripts:

- **Create_RData_objects_from_TableS1.Rmd**: Import Table S1 (xlsx file) from the manuscript which contains NanoString and proteomic data as well as sample and gene/protein annotations. Create an RData object containing all cleaned tables.
- **Create_RData_objects_from_TCGA_and_CPTAC.Rmd**: Import TCGA and CPTAC data from various sources (see links within the Rmd file). Create an RData object containing all cleaned tables.

2.Data analysis scripts:

- **Overview_of_annotations_and_expr_tables.Rmd**: Miscellaneous checks related to data quality, and performing some descriptive statistics.
- **Compare_pROIs_and_tROIs.Rmd**: Correlation analyses between pROIs and tROIs.
- **Consensus_clustering.Rmd**: Perform consensus clustering on the pROIs and tROIs.
- **Histology_comparisons.Rmd**: Differential expression analyses (tumor vs normal, high vs low immune infiltration, high vs low mucin and stroma score), followed by pre-ranked GSEA.
- **Histology_comparisons_CPTAC_validation.Rmd**: Compare our differential expression analysis results with trends in the CPTAC data.
- **Intratumoral_variability_assessment.Rmd**: Seek out proteins and genes that commonly show homogeneous/heterogeneous expression in different regions of the same tumor.
- **FN1_investigation.Rmd**: Cox regression analyses for FN1 expression, both in our data and in the CPTAC and TCGA datasets.

**Helper scripts**:

- **Load_packages.R** : Load all the packages that are used in any of the scripts.
- **Utility_functions.R**: Custom functions written for the data analysis.
- **Color_list.R**: Specifying the colors for heatmap annotations.


## Data availability

The raw proteomics data have been deposited to the MassIVE repository (https://doi.org/doi:10.25345/C5319S63V). The normalized gene counts and protein intensities can be found in Table S1, which is provided as supplementary material alongside the publication. 