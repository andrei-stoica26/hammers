#' @importFrom abdiv euclidean
#' @importFrom withr with_seed
#'
NULL

#' Extract count information from Seurat column
#'
#' This function extracts count information from Seurat column.
#'
#' @param seuratObj A Seurat object.
#' @param column Column
#'
#' @return A frequency vector with the unique column values as names.
#'
#' @noRd
#'
countsVector <- function(seuratObj, column){
    df <- dplyr::count(seuratObj[[]], {{column}})
    v <- setNames(df[, 2], as.factor(df[, 1]))
    return(v)
}


#' Get nearest neighbors from distance matrix
#'
#' This function gets the nearest neighbors from a distance matrix.
#'
#' @param A distance matrix.
#'
#' @return A named character vector.
#'
#' @examples
#' df <- data.frame(v = c(1, 2, 4, 5, 6),
#' w = c(2, 3, 1, 5, 8),
#' x = c(2, 8, 7, 1, 1),
#' y = c(2, 3, 2, 2, 4),
#' z = c(1, 9, 9, 7, 6))
#' distMat <- as.matrix(stats::dist(df))
#' rownames(distMat) <- c('v', 'w', 'x', 'y', 'z')
#' colnames(distMat) <- c('v', 'w', 'x', 'y', 'z')
#' nearestNeighbors(distMat)
#'
#' @export
#'
nearestNeighbors <- function(distMat)
    return(apply(distMat, 1, function(x) {
        names(x) <- colnames(distMat)
        x <- x[x > 0]
        return(names(x)[which.min(x)])
    }))

#' Prepare dataframe for alluvial plot
#'
#' This function extracts the relevant information from dataframe and adjusts
#' p-values to be used as weights for the alluvia.
#'
#' @param repDF A representation data frame.
#' @param pvalCol Name of p-value column.
#' @param colIndices A vector respresenting the indices of the two categorical
#' columns from the data frame that will be used.
#'
#' @return A data frame with weight scores in lieu of p-values.
#'
#' @noRd
#'
prepAlluvial <- function(repDF, pvalCol = 'pvalAdj', colIndices = c(1, 2)){
    pvals <- sort(repDF[[pvalCol]])
    pvals[-1] <- -log(pvals[-1])
    if (pvals[1])
        pvals[1] <- -log(pvals[1]) else
            if (length(pvals) > 2)
                pvals[1] <- 2 * pvals[2] - pvals[3] else
                    pvals[1] <- 1


    resDF <- repDF[, colIndices]
    resDF$weight <- sqrt(pvals)
    return(resDF)
}

#' Compute proximity between two vectors based on Euclidean distance
#'
#' This functions computes proximity between two vectors based on Euclidean
#' distance and an input maximum distance.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param maxDist Maximum distance.
#'
#' @return A number between 0 and 1.
#'
#' @export
#'
proximity <- function(x, y, maxDist)
    return(1 - euclidean(x, y) / maxDist)

#' Perform min-max normalization when possible; otherwise return a single-value
#' vector.
#'
#' This function min-max-normalizes a vector when possible, and otherwise returns
#' a single-value vector.
#'
#' @param scores Numeric vector.
#' @param safeVal Value to replace all values with when all values in the vector
#' are the same.
#'
#' @return Min-max-normalized scores or a single-value vector.
#'
#' @export
#'
safeMinmax <- function(scores, safeVal = 0){
    if(length(unique(scores)) < 2)
        return(rep(safeVal, length(scores)))
    return(liver::minmax(scores))
}

#' Replaces genes from vector
#'
#' This function removes and adds genes from vector at random.
#'
#' @inheritParams metadataDF
#' @param genes A character vector.
#' @param lossFrac Fraction of genes than be removed. Must be in \code{[0, 1]}.
#' @param noiseFrac Amount of noise (random genes) in the final gene vector.
#' Must be in \code{[0, 1)}
#' @param geneCountThresh Minimum number of cells in which newly added genes
#' must be expressed.
#' @param seed Random seed.
#'
#' @return Genes vector after changes.
#'
#' @export
#'
shuffleGenes <- function(scObj, genes, lossFrac, noiseFrac,
                         geneCountThresh = 10, seed = 1){
    nGenes <- length(genes)
    nRemovedGenes <- round(lossFrac * nGenes)
    nRetainedGenes <- nGenes - nRemovedGenes

    expression <- scExpMat(scObj, 'counts')
    freq <- rowSums(expression != 0)

    suitableGenes <- setdiff(names(freq[freq >= geneCountThresh]), genes)

    if(lossFrac > 0){
        genes <- c(with_seed(seed, sample(genes, nRetainedGenes)))
        message('Removed ', nRemovedGenes, ' genes.')
    } else
        message('No genes were removed.')

    if(noiseFrac > 0){
        nAddedGenes <- round(noiseFrac * nRetainedGenes / (1 - noiseFrac))
        genes <- c(genes, with_seed(seed, sample(suitableGenes, nAddedGenes)))
        message('Added ', nAddedGenes, ' random genes.')
    } else
        message('No random genes were added')

    return(genes)
}
