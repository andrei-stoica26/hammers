#' Find rare genes in a Seurat or SingleCellExpression object
#'
#' This function finds genes expressed in a low number of cells in a Seurat or
#' SingleCellExpression object.
#'
#' @inheritParams geneCenters
#' @param nCells Minimum number of cells in which a gene must be expressed to
#' be regarded as non-rare.
#'
#' @return A data frame with the rare genes as rownames and a single column
#' representing their frequencies.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' df <- findRareGenes(sceObj)
#'
#' @export
#'
findRareGenes <- function(scObj, nCells = 10)
    return(genePresence(scObj, maxCutoff=nCells - 1))

#' Remove rare genes from a Seurat or SingleCellExpression object
#'
#' This function removes genes expressed in a low number of cells in a
#' \code{Seurat} or \code{SingleCellExpression} object.
#'
#' @inheritParams geneCenters
#' @param nCells Minimum number of cells in which a gene must be expressed to
#' be retained.
#' @param verbose Logical; whether the output should be verbose.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object with
#' rare genes removed.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' sceObj <- removeRareGenes(sceObj, 30)
#'
#' @export
#'
removeRareGenes <- function(scObj, nCells = 10, verbose = TRUE){
    rareGenes <- findRareGenes(scObj, nCells)[, 1]
    genes <- setdiff(rownames(scObj), rareGenes)
    scObj <- scObj[genes, ]
    safeMessage(paste0('Removed features found in less than ', nCells,
                       ' cells.'), verbose)
    safeMessage(paste0(length(rareGenes), ' rare features removed.'),
                verbose)
    return(scObj)
}
