#' Check if all genes exist in the single-cell expression object
#'
#' This function checks if all genes exist in the single-cell expression
#' object.
#'
#' @inheritParams geneCenters
#'
#' @return None. This function is called for its side effect.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' checkGenes(sceObj, c('Gene_0480', 'Gene_0481', 'Gene_0482'))
#'
#' @export
#'
checkGenes <- function(scObj, genes){
    extraGenes <- setdiff(genes, rownames(scObj))
    if(length(extraGenes))
        stop(paste0(extraGenes, collapse=', '),
             ' gene(s) not found in the single-cell expression object.')
}
