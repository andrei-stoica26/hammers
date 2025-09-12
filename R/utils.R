#' @importFrom abdiv euclidean
#' @importFrom withr with_seed
#'
NULL

#' Get nearest neighbors from distance matrix
#'
#' This function gets the nearest neighbors from a distance matrix.
#'
#' @param distMat A distance matrix.
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
#' @param df A data frame.
#' @param pvalCol Name of p-value column.
#' @param colIndices A vector respresenting the indices of the two categorical
#' columns from the data frame that will be used.
#' @param weightExp Exponent used in constructing weight from p-values.
#' @param offset Offset used to avoid zeros inside the logarithm function.
#'
#' @return A data frame with weight scores in lieu of p-values.
#'
#' @keywords internal
#'
prepAlluvial <- function(df,
                         pvalCol = 'pvalAdj',
                         colIndices = c(1, 2),
                         weightExp = 1/2,
                         offset = 1e-317){
    pvals <- sort(df[[pvalCol]])
    pvals[-1] <- -log(pvals[-1] + offset)
    if (pvals[1])
        pvals[1] <- -log(pvals[1]) else
            if (length(pvals) > 2)
                pvals[1] <- 2 * pvals[2] - pvals[3] else
                    pvals[1] <- 1
    resDF <- df[, colIndices]
    resDF$weight <- pvals ^ weightExp
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
#' @examples
#' proximity(2, 3, 6)
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
#' @examples
#' safeMinmax(c(0, 3, 2, 1, 4, 5.5, 6.32, 8, 1.1))
#'
#' @export
#'
safeMinmax <- function(scores, safeVal = 0){
    if(length(unique(scores)) < 2)
        return(rep(safeVal, length(scores)))
    return(liver::minmax(scores))
}

#' Message an input if verbose is set to TRUE
#'
#' This function messages an input if \code{verbose} is set to \code{TRUE}.
#'
#' @param msg Message
#' @param verbose Whether the message should be displayed.
#'
#' @return No return value. This function is called for its side effect
#' (messaging the input if \code{verbose} is set to \code{TRUE}).
#'
#' @examples
#' safeMessage('message')
#'
#' @export
#'
safeMessage <- function(msg, verbose = TRUE)
    if(verbose)
        message(msg)

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
#' @param verbose Whether the output should be verbose.
#'
#' @return Genes vector after changes.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' genes <- c('TOP2A', 'BIRC5', 'MKI67', 'RRM2', 'CENPF', 'PTTG2', 'CLSPN')
#' shuffleGenes(scObj, genes, 0.3, 0.9)
#'
#' @export
#'
shuffleGenes <- function(scObj, genes, lossFrac, noiseFrac,
                         geneCountThresh = 10, seed = 1,
                         verbose = TRUE){
    nGenes <- length(genes)
    nRemovedGenes <- round(lossFrac * nGenes)
    nRetainedGenes <- nGenes - nRemovedGenes

    expression <- scExpMat(scObj, 'counts')
    freq <- rowSums(expression != 0)

    suitableGenes <- setdiff(names(freq[freq >= geneCountThresh]), genes)

    if(lossFrac > 0){
        genes <- c(with_seed(seed, sample(genes, nRetainedGenes)))
        safeMessage(paste0('Removed ', nRemovedGenes, ' genes.'), verbose)
    } else
        safeMessage('No genes were removed.', verbose)

    if(noiseFrac > 0){
        nAddedGenes <- round(noiseFrac * nRetainedGenes / (1 - noiseFrac))
        genes <- c(genes, with_seed(seed, sample(suitableGenes, nAddedGenes)))
        safeMessage(paste0('Added ', nAddedGenes, ' random genes.'), verbose)
    } else
        safeMessage('No genes were added.', verbose)

    return(genes)
}

#' Convert a vector to a data frame based on input row and column names
#'
#' This function converts a vector to a data frame based on input row and
#' column names. Optionally, it also calculates the row means.
#'
#' @param v A vector.
#' @param rowNames A character vector.
#' @param colNames A character vector.
#' @param addRowMeans Whether to add the row means to the data frame.
#' @param sortByRowMeans Whether to sort by row means.
#'
#' @return A data frame.

#' @examples
#' v <- c(2, 3, 4, 19, 15, 25, 32, 8)
#' res <- tabulateVector(v, paste0('r', seq(4)), paste0('c', seq(2)))
#'
#' @export
#'
tabulateVector <- function(v,
                           rowNames,
                           colNames,
                           addRowMeans=FALSE,
                           sortByRowMeans=FALSE){
    df <- data.frame(matrix(v, length(rowNames), length(colNames)))
    rownames(df) <- rowNames
    colnames(df) <- colNames
    if(addRowMeans){
        df$avg <- rowMeans(df)
        if(sortByRowMeans)
            df <- df[order(df$avg), ]
    }
    return(df)
}
