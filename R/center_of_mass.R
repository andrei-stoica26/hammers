#' Calculate the coordinates of the center of mass
#'
#' This function calculates the coordinates of the center of mass based on a
#' matrix of cell embeddings and a vector of weights.
#'
#' @param dimMat A matrix of cell embeddings.
#' @param weights A vector of weights.
#'
#' @return A vector containing the coordinates of the center of mass.
#'
#' @examples
#' dimMat <- matrix(data=c(2, 3, 1, 3, 6, 8), nrow=3, ncol=2)
#' weights <- c(0.8, 6, 16)
#' centerOfMass(dimMat, weights)
#'
#' @export
#'
centerOfMass <- function(dimMat, weights){
  totalWeight <- sum(weights)
  return(apply(dimMat, 2, function(x) sum(x * weights) / totalWeight))
}

#' Calculate the coordinates of centers of mass of gene expression
#'
#' This function calculates the coordinates of the center of mass of the
#' expression of input genes.
#'
#' @inheritParams scExpMat
#'
#' @return A data frame containing the coordinates of centers of mass.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runUMAP(scObj)
#' geneCenters(scObj, c('Gene_0980', 'Gene_0981', 'Gene_0982'))
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

#' Calculate the coordinates of centers of mass of metadata/coldata columns
#'
#' This function calculates the coordinates of the center of mass of the
#' selected metadata/coldata columns from a Seurat or SingleCellExpression
#' object.
#'
#' @inheritParams metadataDF
#' @param columns Numeric columns.
#'
#' @return A data frame containing the coordinates of centers of mass.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runUMAP(scObj)
#' colCenters(scObj, c('sizeFactor'))
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
