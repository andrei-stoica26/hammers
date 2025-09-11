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
#' scObj <- scRNAseq::BaronPancreasData('human')
#' checkGenes(scObj, c('AURKA', 'TOP2A', 'MKI67'))
#' checkGenes(scObj, c('AURKA', 'TOP2A', 'MKI67', 'DSFDGDG'))
#'
#' @export
#'
checkGenes <- function(scObj, genes){
    extraGenes <- setdiff(genes, rownames(scObj))
    if(length(extraGenes))
        stop('Gene ', extraGenes[1],
             ' not found in the single-cell expression object.')
}
