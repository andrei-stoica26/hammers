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
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
#' checkGenes(scObj, c('Gene_0980', 'Gene_0981', 'Gene_0982'))
#'
#' @export
#'
checkGenes <- function(scObj, genes){
    extraGenes <- setdiff(genes, rownames(scObj))
    if(length(extraGenes))
        stop('Gene ', extraGenes[1],
             ' not found in the single-cell expression object.')
}
