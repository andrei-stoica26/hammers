#' Find rare genes in a Seurat or SingleCellExpression object
#'
#' This function finds genes expressed in a low number of cells in a Seurat or
#' SingleCellExpression object.
#'
#' @inheritParams metadataDF
#' @param nCells Minimum number of cells in which a gene must be expressed to
#' be regarded as non-rare.
#'
#' @return A data frame with the rare genes as rownames and a single column
#' representing their frequencies.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' df <- findRareGenes(sceObj)
#'
#' @export
#'
findRareGenes <- function(scObj, nCells = 10)
    return(genePresence(scObj, maxCutoff=nCells - 1))

#' Remove rare genes from a Seurat or SingleCellExpression object
#'
#' This function removes genes expressed in a low number of cells in a Seurat
#' or SingleCellExpression object.
#'
#' @inheritParams metadataDF
#' @param nCells Minimum number of cells in which a gene must be expressed to
#' be retained.
#' @param verbose Logical; whether the output should be verbose.
#'
#' @return A Seurat or SingleCellExpression object with rare genes removed.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' sceObj <- removeRareGenes(sceObj, 30)
#'
#' @export
#'
removeRareGenes <- function(scObj, nCells = 10, verbose = TRUE){
    rareGenes <- findRareGenes(scObj, nCells)[, 1]
    genes <- setdiff(rownames(scObj), rareGenes)
    scObj <- scObj[genes, ]
    if (verbose){
        message('Removed features found in less than ', nCells, ' cells.')
        message(length(rareGenes), ' rare features removed.')
    }
    return(scObj)
}
