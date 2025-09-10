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
#' @export
#'
geneCenters <- function(scObj, genes){
    genesExp <- scExpMat(scObj, genes=genes)
    dimMat <- scUMAPMat(scObj)
    centersDF <- data.frame(t(apply(genesExp, 1,
                       function(x) centerOfMass(dimMat, x))))
    return(centersDF)
}
