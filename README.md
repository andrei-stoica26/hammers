# hammers

`hammers` is a utilities suite for scRNA-seq data analysis compatible with both
`Seurat` and `SingleCellExperiment`. It provides simple tools to address tasks 
such as retrieving aggregate gene statistics, finding and removing rare genes, 
performing representation analysis, computing the center of mass for the 
expression of a gene of interest in low-dimensional space, and calculating 
silhouette and cluster-normalized silhouette.

## Installation

To install `hammers`, run the following R code:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("hammers")
```
