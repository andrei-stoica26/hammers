#' @importFrom dplyr count
#' @import SeuratObject
#' @importFrom stats phyper quantile setNames
#'
NULL

#' Find the differential representation of two Seurat columns
#'
#' This function find the differential representation of two Seurat columns.
#'
#' @inheritParams scColPairCounts
#'
#' @param doOverrep Whether to perform overrepresentation analysis. If
#' \code{FALSE}, underrepresentation analysis will be performed instead.
#' @inheritParams mtCorrectDF
#'
#' @return An overrepresentation or underrepresentation data frame.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' repAnalysis(sceObj, 'Cluster', 'Donor')
#'
#' @export
#'
repAnalysis <- function(scObj,
                        col1 = 'seurat_clusters',
                        col2 = 'orig.ident',
                        doOverrep = TRUE,
                        mtMethod = c('by', 'bf', 'bh'),
                        ...){
    mtMethod <- match.arg(mtMethod, c('by', 'bf', 'bh'))
    nCells <- ncol(scObj)
    df <- scColPairCounts(scObj, col1, col2)
    colnames(df)[3] <- 'sharedCount'
    v1 <- scColCounts(scObj, col1)
    v2 <- scColCounts(scObj, col2)
    newCols <- paste0(colnames(df)[c(1, 2)], 'Count')
    df[[newCols[1]]] <- vapply(df[, 1], function(x) v1[x], integer(1))
    df[[newCols[2]]] <- vapply(df[, 2], function(x) v2[x], integer(1))
    df$expectedCount <- df[, 4] * df[, 5] / nCells
    df$ratio <- df$sharedCount / df$expectedCount
    df$pval <- mapply(function(q, m, n, k)
        phyper(q, m, n, k, lower.tail=1 - doOverrep),
        df[, 3] - doOverrep, df[, 4], nCells - df[, 4], df[, 5])
    df <- mtCorrectDF(df, mtMethod, ...)
    return(df)
}
