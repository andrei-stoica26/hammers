# hammers

`hammers` is a utilities suite for scRNA-seq data analysis and package 
development. It offers functions—both getters and setters—that can operate 
on both Seurat and SingleCellExperiment objects, intended to help 
developers build tools compatible with both types of input. Additionally, 
`hammers` provides simple tools to address tasks in scRNA-seq data analysis, 
such as retrieving aggregate gene statistics, finding and removing rare genes, 
performing multiple testing correction, computing the center of mass for the 
expression of a gene of interest in low-dimensional space, and calculating 
silhouette and cluster-normalized silhouette.

## Installation

To install `hammers`, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("hammers")
```
