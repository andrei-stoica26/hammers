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
#' SingleCellExperiment object as a data frame
#'
#' @param scObj A \code{Seurat} or \code{SingleCellExperiment} object.
#'
#' @return A metadata data frame.
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
#' @export
#'
metadataNames <- function(scObj)
    UseMethod(generic='metadataNames', object=scObj)


#' Extract a metadata column from object.
#'
#' This function extracts a metadata column from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#' @param colStr Column name.
#'
#' @return A column vector.
#'
#' @export
#'
scCol <- function(scObj, colStr)
    UseMethod(generic='scCol', object=scObj)

#' Extracts the PCA matrix from object.
#'
#' This function extracts the PCA matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#'
#' @return A PCA matrix.
#'
#' @export
#'
scPCAMat <- function(scObj)
    UseMethod(generic='scPCAMat', object=scObj)

#' Extracts the expression matrix from object.
#'
#' This function extracts an expression matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#' @param dataType Expression data type.
#' @param genes Selected genes. If \code{NULL}, all genes will be retained
#' @param densify Whether to convert to dense matrix.
#'
#' @return A PCA matrix.
#'
#' @export
#'
scExpMat <- function(scObj, dataType = c('counts',
                                         'data',
                                         'logcounts'),
                     genes = NULL,
                     densify = TRUE)
    UseMethod(generic='scExpMat', object=scObj)
