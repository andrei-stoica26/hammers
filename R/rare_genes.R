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
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' df <- findRareGenes(scObj)
#'
#' @export
#'
findRareGenes <- function(scObj, nCells = 10){
    expMat <- scExpMat(scObj, "counts")
    df <- data.frame(Count = rowSums(expMat != 0))
    df <- df[df[, 1] < nCells, ]
    return(df)
}

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
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' scObj <- removeRareGenes(scObj, 30)
#'
#' @export
#'
removeRareGenes <- function(scObj, nCells = 10, verbose = TRUE){
    rareGenes <- rownames(findRareGenes(scObj, nCells))
    genes <- setdiff(rownames(scObj), rareGenes)
    scObj <- scObj[genes, ]
    if (verbose){
        message('Removed features found in less than ', nCells, ' cells.')
        message(length(rareGenes), ' rare features removed.')
    }
    return(scObj)
}
