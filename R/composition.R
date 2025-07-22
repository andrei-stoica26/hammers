#' @importFrom dplyr count
#' @import SeuratObject
#' @importFrom sgof BY
#' @importFrom stats phyper setNames
#'
NULL

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

#' Find the differential representation of two Seurat columns
#'
#' This function find the differential representation of two Seurat columns.
#'
#' @param seuratObj A Seurat object.
#' @param column1 First column.
#' @param column2 Second column.
#' @param doOverrep Whether to perform overrepresentation analysis. If
#' \code{FALSE}, underrepresentation analysis will be performed.
#' @param pvalThr p-value threshold.
#'
#' @return An overrepresentation or underrepresentation data frame.
#'
#' @export
#'
#'
repAnalysis <- function(seuratObj, column1, column2, doOverrep=TRUE,
                        pvalThr=0.05){
    nCells <- ncol(seuratObj)
    df <- dplyr::count(seuratObj[[]], {{column1}}, {{column2}})
    colnames(df)[3] <- 'sharedCount'
    v1 <- countsVector(seuratObj, {{column1}})
    v2 <- countsVector(seuratObj, {{column2}})
    newCols <- paste0(colnames(df)[c(1, 2)], 'Count')
    df[[newCols[1]]] <- vapply(df[, 1], function(x) v1[x], integer(1))
    df[[newCols[2]]] <- vapply(df[, 2], function(x) v2[x], integer(1))
    df$expectedCount <- df[, 4] * df[, 5] / nCells
    df$ratio <- df$sharedCount / df$expectedCount
    df$pval <- mapply(function(q, m, n, k)
        phyper(q, m, n, k, lower.tail=1 - doOverrep),
        df[, 3] - doOverrep, df[, 4], nCells - df[, 4], df[, 5])
    df <- df[order(df$pval), ]
    df$pvalAdj <- BY(df$pval, pvalThr)$Adjusted.pvalues
    df <- subset(df, pvalAdj < pvalThr)
    return(df)
}

#' Prepare representation dataframe for alluvial plot
#'
#' This function extracts the relevant information from representation dataframe
#' generated with \code{repAnalysis} and adjust p-values to be used as weights
#' for the alluvia.
#'
#' @param repDF A representation data frame.
#'
#' @return An adjusted representation data frame.
#'
#' @noRd
#'
prepAlluvial <- function(repDF){
    pvals <- unique(repDF$pvalAdj)
    pvals[-1] <- -log(pvals[-1])
    if (pvals[1])
        pvals[1] <- -log(pvals[1]) else
            if (length(pvals) > 2)
                pvals[1] <- 2 * pvals[2] - pvals[3] else
                    pvals[1] <- 1


    resDF <- repDF[, c(1, 2)]
    resDF$weight <- log(pvals)
    return(resDF)
}

