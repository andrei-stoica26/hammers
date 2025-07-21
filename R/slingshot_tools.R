#' @importFrom slingshot slingPseudotime slingCurveWeights
#' @importFrom tidyr replace_na
#'
NULL

#' Add slingshot information to Seurat object.
#'
#' This function adds slingshot information to Seurat object.
#'
#' @param seuratObj A Seurat object.
#' @param slingshotObj A PseudotimeOrdering object generated using
#' \code{slingshot::slingshot}
#' @param colStr Column name
#' @param fun Function
#'
#' @return A Seurat object with slingshot information added to the metadata.
#'
#' @keywords internal
#'
addSlingshotResults <- function(seuratObj, slingshotObj, colStr, fun){
    nLineages <- ncol(fun(slingshotObj))
    for (i in 1:nLineages)
        seuratObj[[]][[paste0(colStr, i)]] <- tidyr::replace_na(fun(slingshotObj)[, i], NaN)
    return (seuratObj)
}

#' Add slingshot lineages to Seurat object
#'
#' This function adds slingshot lineages to Seurat object.
#'
#' @inheritParams addSlingshotResults
#'
#' @return A Seurat object with slingshot lineages added to the metadata.
#'
addLineages <- function(seuratObj, slingshotObj)
    return(addSlingshotResults(seuratObj, slingshotObj, "Lineage", slingPseudotime))

#' Add slingshot curve weights to Seurat object
#'
#' This function adds slingshot lineages to Seurat object.
#'
#' @inheritParams addSlingshotResults
#'
#' @return A Seurat object with slingshot curve weights added to the metadata.
#'
addCurveweights <- function(seuratObj, slingshotObj)
    return(addSlingshotResults(seuratObj, slingshotObj, "Curveweight", slingCurveWeights))

