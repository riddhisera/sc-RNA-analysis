# sc-RNA-analysis
 Functional Genomics final project which explores single cell RNA and trajectory analysis of differentiating pancreatic stem cells 

## Goals of this project:
1. Identifying cell type identity through gene expression
2. Observing gene expression changes through differentiation
3. Hypothesize developmental trajectory using gene expression
4. Model developmental trajectories and confirm hypothesis
   
## Dataset of interest:
[Dataset link](https://singlecell.broadinstitute.org/single_cell/study/SCP1526/functional-metabolic-and-transcriptional-maturation-of-human-pancreatic-islets-derived-from-stem-cells?cluster=Beta%20cells&spatialGroups=--&annotation=Cell%20type--group--cluster&subsample=all#study-summary) - This dataset contains scRNA data of human pluripotent stem cells differentiating into pancreatic cells taken at different time points. I found this data to be interesting as it was suitable to perform trajectory analysis and understand how gene expression changes through cell differentiation.

## Method of interest:
Standard single cell RNA analysis and trajectory analysis using PAGA and Slingshot. I chose these methods as they can perform trajectory analysis without the need of spliced/unspliced RNA data like sc-Velo or VeloCyto and also explore how I can incorporate gene expression inferences into lineage identification.

[PAGA](https://github.com/theislab/paga) - PAGA trajectory analysis is a method for mapping cell development paths in single-cell RNA sequencing, highlighting how cells transition between different states or types during differentiation.

[Slingshot](https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html) - Slingshot is a computational tool used for inferring cellular lineages and trajectories in single-cell RNA sequencing data. It identifies differentiation paths in multi-dimensional data, allowing researchers to trace the progression of cell states and types over time.

## Summary of files:
1. Final_Project_Report_Riddhi_Sera.pdf - The final report (duh.)
2. Project Code.ipynb - Code for performing sc-RNA analysis and PAGA
3. R code for slingshot.R - Code for implementing Slingshot using Seurat
