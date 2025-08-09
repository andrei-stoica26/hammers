#' @importFrom dplyr count
#' @import SeuratObject
#' @importFrom sgof BY
#' @importFrom stats phyper setNames
#'
NULL

#' Find the differential representation of two Seurat columns
#'
#' This function find the differential representation of two Seurat columns.
#'
#' @param seuratObj A Seurat object.
#' @param column1 First column.
#' @param column2 Second column.
#' @param doOverrep Whether to perform overrepresentation analysis. If
#' \code{FALSE}, underrepresentation analysis will be performed instead.
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
