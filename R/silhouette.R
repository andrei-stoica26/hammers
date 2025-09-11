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
#' @param pcaMat PCA matrix.
#'
#' @return The input object (\code{Seurat} or \code{SingleCellExperiment}) with
#' an added metadata silhouette column.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' scObj <- computeSilhouette(scObj, 'label')
#' head(scCol(scObj, 'silhouette'))
#'
#' @export
#'
computeSilhouette <- function(scObj, idClass, pcaMat = NULL){
  if (!idClass %in% metadataNames(scObj))
    stop(paste0('Column ', idClass, ' not found.'))
  if(is.null(pcaMat))
    pcaMat <- scPCAMat(scObj)
  message('Computing cosine distance matrix...')
  distMat <- 1 - text2vec::sim2(pcaMat, method='cosine', norm='l2')
  message(paste0('Computing silhouette for identity class: ', idClass, '...'))
  groupVals <- unclass(factor(scCol(scObj, idClass)))
  scObj[['silhouette']] <- cluster::silhouette(groupVals, dmatrix=distMat)[, 3]
  return(scObj)
}

#' Normalize silhouette by identity class for single-cell expression object
#'
#' This function normalizes the already computed silhouette for each identity
#' class in the single-cell expression object.
#'
#' @inheritParams computeSilhouette
#'
#' @return A data frame with normalized silhouettes for each unique element in the
#' identity class.
#'
#' @export
#'
normalizeSilhouette <- function(scObj, idClass){
  df <- metadataDF(scObj)[, c(idClass, 'silhouette')]
  groups <- unique(df[[idClass]])
  res <- data.frame(matrix(0, nrow(df), length(groups)))
  rownames(res) <- rownames(df)
  colnames(res) <- groups
  for (group in groups){
    dfSub <- df[df[, 1] == group, ]
    sortedSil <- dfSub$silhouette
    auxMin <- 2 * sortedSil[1] - sortedSil[2]
    normSilVals <- liver::minmax(c(dfSub$silhouette, auxMin))[seq_along(dfSub$silhouette)]
    res[rownames(dfSub), group] <- normSilVals
  }
  return(res)
}
