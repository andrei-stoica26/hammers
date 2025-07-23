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

`bfCorrectDF`: Perform multiple testing correction with Bonferroni.

`byCorrectDF`: Perform multiple testing correction with Benjamini-Yekutieli.

`pvalRiverPlot`: Create a river plot using 
[henna](https://github.com/andrei-stoica26/henna) from a data frame with two
categorical columns and a p-value column by converting p-values to weights.

`repAnalysis`: Given two gene expression metadata categorical columns, performs
representation analysis for all the pairs consisting of a label from each 
category. Both overrepresentation and underrepresentation analysis are supported.
