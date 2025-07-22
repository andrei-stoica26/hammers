# hammers
`hammers` is an evolving collection of simple tools designed to address common 
tasks in scRNA-seq data analysis. The included tools are intended to combine 
wide applicability with ease of use.

## Installation

To install `hammers`, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("andrei-stoica26/hammers")
```
## Functions implemented so far

`addLineages`: Adds [slingshot](https://www.bioconductor.org/packages/release/bioc/html/slingshot.html) lineages to Seurat metadata.

`addCurveweights`: Adds slingshot curve weights to Seurat metadata.

`repAnalysis`: Given two gene expression metadata categorical columns, performs
representation analysis for all the pairs consisting of a label from each 
category. Both overrepresentation and underrepresentation analysis are supported.
