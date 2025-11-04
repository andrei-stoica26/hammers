#' Check if all genes exist in the single-cell expression object
#'
#' This function checks if all genes exist in the single-cell expression
#' object.
#'
#' @inheritParams scExpMat
#' @param genes A character vector of genes.
#'
#' @return None. This function is called for its side effect.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' checkGenes(sceObj, c('Gene_0480', 'Gene_0481', 'Gene_0482'))
#'
#' @export
#'
checkGenes <- function(scObj, genes){
    extraGenes <- setdiff(genes, rownames(scObj))
    if(length(extraGenes))
        stop('Gene ', extraGenes[1],
             ' not found in the single-cell expression object.')
}
