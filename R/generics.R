#' Saves plot or list of plots
#'
#' This function saves a plot or list of plots as a pdf. Can also take as input
#' a function that returns a ggplot object together with its arguments.
#'
#' @param plotObject A function, ggplot object, or list of ggplot objects.
#' @param ... Additional arguments.
#'
#' @return No value. This function is called for its side effect.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = c(1, 2), y = c(3, 5))
#' p <- ggplot(df) + geom_point(aes(x, y))
#' devPlot(p)
#'
#' simplePlot <- function(df, title)
#'     return(ggplot(df) + geom_point(aes(x, y)) + ggtitle(title))
#'
#' devPlot(simplePlot, df, 'Plot title')
#'
#' @export
#'
#'
devPlot <- function(plotObject, ...)
    UseMethod(generic='devPlot', object=plotObject)

#' Extract metadata from object as a data frame
#'
#' This function extracts the metadata from a Seurat or
#' SingleCellExperiment object as a data frame.
#'
#' @param scObj A \code{Seurat} or \code{SingleCellExperiment} object.
#'
#' @return A metadata data frame.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' df <- metadataDF(scObj)
#'
#' @export
#'
metadataDF <- function(scObj)
    UseMethod(generic='metadataDF', object=scObj)


#' Return metadata names.
#'
#' This function extracts metadata names from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataDF
#'
#' @return The names of the metadata columns.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' colNames <- metadataNames(scObj)
#'
#' @export
#'
metadataNames <- function(scObj)
    UseMethod(generic='metadataNames', object=scObj)


#' Extract a metadata/coldata column from object.
#'
#' This function extracts a metadata/coldata column from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#' @param col Column name.
#'
#' @return A vector.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' v <- scCol(scObj, 'label')
#'
#' @export
#'
scCol <- function(scObj, col)
    UseMethod(generic='scCol', object=scObj)

#' Extracts a dimensionality reduction matrix from object.
#'
#' This function extracts a dimensionality reduction matrix (PCA or UMAP)
#' from a Seurat or SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#'
#' @return A PCA or UMAP matrix.
#'
#' @noRd
#'
scDimredMat <- function(scObj, dimred = c('pca', 'umap'))
    UseMethod(generic='scDimredMat', object=scObj)

#' Extracts the expression matrix from object.
#'
#' This function extracts an expression matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams scGeneExp
#' @param genes Selected genes. If \code{NULL}, all genes will be retained
#' @param densify Whether to convert to dense matrix.
#'
#' @return An expression matrix.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' mat <- scExpMat(scObj)
#'
#' @export
#'
scExpMat <- function(scObj, dataType = c('data',
                                         'counts',
                                         'logcounts'),
                     genes = NULL,
                     densify = TRUE)
    UseMethod(generic='scExpMat', object=scObj)

#' Extracts the expression of a single gene
#'
#' This function extracts the expression of a single gene from a Seurat,
#' SingleCellExperiment, dgCMatrix or matrix object.
#
#' @param scObj A \code{Seurat}, \code{SingleCellExperiment},
#' \code{dgCMatrix} or \code{matrix} object.
#' @param gene Selected gene.
#' @param dataType Expression data type. Ignored if \code{scObj} is of class
#' \code{dgCMatrix} or \code{matrix}.
#'
#' @return A gene expression vector.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' v <- scGeneExp(scObj, 'AURKA')
#'
#' @export
#'
scGeneExp <- function(scObj, gene, dataType = c('data',
                                                'counts',
                                                'logcounts'))
    UseMethod(generic='scGeneExp', object=scObj)
