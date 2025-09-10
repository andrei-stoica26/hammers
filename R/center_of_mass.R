#' Calculate the coordinates of the center of mass
#'
#' This function calculates the coordinates of the center of mass based on a
#' matrix of cell embeddings and a vector of weights
#'
#' @param dimMat A matrix of cell embeddings
#' @param weights A vector of weights
#'
#' @return A vector containing the coordinates of the center of mass
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
#' sceObj <- scRNAseq::BaronPancreasData('human')
#' sceObj <- scuttle::logNormCounts(sceObj)
#' sceObj <- scater::runUMAP(sceObj)
#' geneCenters(sceObj, c('AURKA', 'MKI67', 'TOP2A'))
#'
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
#' @export
#'
colCenters <- function(scObj, columns){
    cols <- metadataDF(scObj)[, columns]
    dimMat <- scUMAPMat(scObj)
    centersDF <- data.frame(t(apply(cols, 2,
                                    function(x) centerOfMass(dimMat, x))))
    return(centersDF)
}


