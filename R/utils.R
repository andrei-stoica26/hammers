#' Extract count information from Seurat column.
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
