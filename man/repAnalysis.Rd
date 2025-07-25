% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composition.R
\name{repAnalysis}
\alias{repAnalysis}
\title{Find the differential representation of two Seurat columns}
\usage{
repAnalysis(seuratObj, column1, column2, doOverrep = TRUE, pvalThr = 0.05)
}
\arguments{
\item{seuratObj}{A Seurat object.}

\item{column1}{First column.}

\item{column2}{Second column.}

\item{doOverrep}{Whether to perform overrepresentation analysis. If
\code{FALSE}, underrepresentation analysis will be performed.}

\item{pvalThr}{p-value threshold.}
}
\value{
An overrepresentation or underrepresentation data frame.
}
\description{
This function find the differential representation of two Seurat columns.
}
