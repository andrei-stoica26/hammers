#' @importFrom scLang dimPlot metadataDF metadataNames scCol scCol<- scColCounts scColPairCounts scColPairPercs scExpMat scPCAMat scUMAPMat
#'
NULL

#' Calculate the center of mass of columns
#'
#' This function calculates the center of mass based on the columns of
#' a data frame or matrix and a vector of weights.
#'
#' @param obj A data frame or matrix.
#' @param weights A vector of weights.
#'
#' @return A vector containing the center of mass.
#'
#' @examples
#' obj <- matrix(data=c(2, 3, 1, 3, 6, 8), nrow=3, ncol=2)
#' weights <- c(0.8, 6, 16)
#' centerOfMass(obj, weights)
#'
#' @export
#'
centerOfMass <- function(obj, weights){
    #Less fast than matrixStats::colWeightedMeans,
    #but able to handle negative weights.
    totalWeight <- sum(weights)
    return(apply(obj, 2, function(x) sum(x * weights) / totalWeight))
}

#' Calculate the centers of mass of the expression of input genes
#'
#' This function calculates the centers of mass of the expression of input
#' genes.
#'
#' @param scObj A \code{Seurat}, \code{SingleCellExperiment},
#' \code{dgCMatrix} or \code{matrix} object.
#' @param genes A character vector of genes.
#'
#' @return A data frame containing the centers of mass.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' geneCenters(sceObj, c('Gene_0480', 'Gene_0481', 'Gene_0482'))
#'
#' @export
#'
geneCenters <- function(scObj, genes){
    genesExp <- scExpMat(scObj, genes=genes)
    dimMat <- scUMAPMat(scObj)
    centersDF <- data.frame(t(apply(genesExp, 1,
                       function(x) centerOfMass(dimMat, x))))
    return(centersDF)
}

#' Calculate the centers of mass of metadata/coldata columns
#'
#' This function calculates the centers of mass of  selected metadata/coldata
#' columns from a \code{Seurat} or \code{SingleCellExpression}
#' object.
#'
#' @inheritParams geneCenters
#' @param columns Numeric columns.
#'
#' @return A data frame containing the coordinates of centers of mass.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' colCenters(sceObj, c('sizeFactor'))
#'
#' @export
#'
colCenters <- function(scObj, columns){
    cols <- metadataDF(scObj)[, columns, drop=FALSE]
    dimMat <- scUMAPMat(scObj)
    centersDF <- data.frame(t(apply(cols, 2,
                                    function(x) centerOfMass(dimMat, x))))
    return(centersDF)
}
