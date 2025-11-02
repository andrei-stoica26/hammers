#' @importFrom text2vec sim2
#' @importFrom cluster silhouette
#' @importFrom liver minmax
#' @importFrom stats cor dist weighted.mean
#'
NULL

#' Compute cluster silhouette for single-cell expression object
#'
#' This function computes the silhouette for each cell in the \code{Seurat}
#' or \code{SingleCellExperiment} object.
#'
#' @inheritParams metadataDF
#' @param idClass Identity class. Must be present among the metadata columns of
#' the single-cell expression object.
#' @param silCol The name of the silhouette column to be added.
#'
#' @return The input object (\code{Seurat} or \code{SingleCellExperiment}) with
#' an added metadata silhouette column.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runPCA(scObj)
#' scCol(scObj, 'Cluster') <- withr::with_seed(1,
#' sample(paste0('Cluster', seq(5)), dim(scObj)[2], replace=TRUE))
#' scObj <- computeSilhouette(scObj, 'Cluster')
#'
#' @export
#'
computeSilhouette <- function(scObj, idClass, silCol = 'silhouette'){
    if (!idClass %in% metadataNames(scObj))
        stop(paste0('Column ', idClass, ' not found.'))
    pcaMat <- scPCAMat(scObj)
    message('Computing cosine distance matrix...')
    distMat <- 1 - text2vec::sim2(pcaMat, method='cosine', norm='l2')
    message(paste0('Computing silhouette for identity class: ',
                   idClass, '...'))
    groupVals <- unclass(factor(scCol(scObj, idClass)))
    scObj[[silCol]] <- cluster::silhouette(groupVals,
                                                 dmatrix=distMat)[, 3]
    return(scObj)
}

#' Normalize silhouette by identity class for single-cell expression object
#'
#' This function normalizes the already computed silhouette for each identity
#' class in the single-cell expression object.
#'
#' @inheritParams computeSilhouette
#' @param silCol The name of the silhouette column.
#'
#' @return A data frame with normalized silhouettes for each unique element in
#' the identity class.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runPCA(scObj)
#' scCol(scObj, 'Cluster') <- withr::with_seed(1,
#' sample(paste0('Cluster', seq(5)), dim(scObj)[2], replace=TRUE))
#' scObj <- computeSilhouette(scObj, 'Cluster')
#' df <- normalizeSilhouette(scObj, 'Cluster')
#'
#' @export
#'
normalizeSilhouette <- function(scObj, idClass, silCol='silhouette'){
    df <- metadataDF(scObj)[, c(idClass, silCol)]
    groups <- unique(df[[idClass]])
    res <- data.frame(matrix(0, nrow(df), length(groups)))
    rownames(res) <- rownames(df)
    colnames(res) <- groups
    for (group in groups){
        dfSub <- df[df[, 1] == group, ]
        sortedSil <- dfSub[, silCol]
        auxMin <- 2 * sortedSil[1] - sortedSil[2]
        normSilVals <- liver::minmax(c(dfSub[, silCol],
                                       auxMin))[seq_along(dfSub[, silCol])]
        res[rownames(dfSub), group] <- normSilVals
    }
    return(res)
}

#' Adds normalized silhouette column to a single-cell expression object
#'
#' This function adds a normalized silhouette column to a single-cell
#' expression object.
#'
#' @inheritParams computeSilhouette
#' @param normSilDF Normalized silhouette data frame.
#' @param normSilCol The name of the normalized silhouette column to be added.
#
#' @return The input object (\code{Seurat} or \code{SingleCellExperiment}) with
#' an added metadata normalized silhouette column.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runPCA(scObj)
#' scCol(scObj, 'Cluster') <- withr::with_seed(1,
#' sample(paste0('Cluster', seq(5)), dim(scObj)[2], replace=TRUE))
#' scObj <- computeSilhouette(scObj, 'Cluster')
#' df <- normalizeSilhouette(scObj, 'Cluster')
#' scObj <- addNormSilhouette(scObj, df)
#'
#' @export
#'
addNormSilhouette <- function(scObj, normSilDF, normSilCol='normSilhouette'){
    scCol(scObj, normSilCol) <- apply(normSilDF, 1, max)
    return(scObj)
}


