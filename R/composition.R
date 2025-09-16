#' @importFrom dplyr count
#' @import SeuratObject
#' @importFrom stats phyper setNames
#'
NULL

#' Find the differential representation of two Seurat columns
#'
#' This function find the differential representation of two Seurat columns.
#'
#' @inheritParams scColPairCounts
#' @param doOverrep Whether to perform overrepresentation analysis. If
#' \code{FALSE}, underrepresentation analysis will be performed instead.
#' @param fdrMethod False discovery rate control method. Options are 'BY'
#' (Benjamini-Yekutieli) and 'BH' (Benjamini-Hochberg).
#' @param pvalThr p-value threshold.
#'
#' @return An overrepresentation or underrepresentation data frame.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scCol(scObj, 'Cluster') <- withr::with_seed(1,
#' sample(paste0('Cluster', seq(5)), dim(scObj)[2], replace=TRUE))
#' scCol(scObj, 'Donor') <- rep('Donor1', dim(scObj)[2])
#' for (i in seq(5)){
#' scCol(scObj, 'Donor')[sample(which(scCol(scObj, 'Cluster') ==
#' paste0('Cluster', i)), 15)]<- paste0('Donor', i)
#' scCol(scObj, 'Donor')[sample(which(scCol(scObj, 'Cluster') ==
#' paste0('Cluster', i))
#' , 15)]<- paste0('Donor', i + 1)
#' }
#' repAnalysis(scObj, 'Cluster', 'Donor')
#'
#' @export
#'
repAnalysis <- function(scObj,
                        col1 = 'seurat_clusters',
                        col2 = 'orig.ident',
                        doOverrep = TRUE,
                        fdrMethod = c('BY', 'BH'),
                        pvalThr = 0.05){
    fdrMethod <- match.arg(fdrMethod, c('BY', 'BH'))
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
    fdrControlFun <- eval(as.name(paste0(tolower(fdrMethod), 'CorrectDF')))
    df <- fdrControlFun(df, pvalThr)
    return(df)
}
