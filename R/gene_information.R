#' Extract gene presence from a Seurat or SingleCellExperiment object
#'
#' This function extracts the number of cells in which a gene from a Seurat or
#' SingleCellExperiment is expressed.
#'
#' @inheritParams metadataDF
#' @param genes Genes for which the number of cells in which the gene is
#' expressed will be computed. If \code{NULL} (as default), all genes will be
#' included in the assessment.
#' @param minCutoff Minimum cutoff for gene counts. Genes with counts below
#' this value will be omitted.
#' @param minCutoff Maximum cutoff for gene counts. Genes with counts above
#' this value will be omitted.
#'
#' @return A data frame with two columns. The first column lists the genes
#' ordered decreasingly by the number of cells in which they appear, the second
#' lists the corresponding numbers of cells.
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=200))
#' df <- genePresence(scObj)
#'
#' @export
#'
genePresence <- function(scObj,
                         genes = NULL,
                         minCutoff = NULL,
                         maxCutoff = NULL){
    expMat <- scExpMat(scObj, 'counts', genes)
    df <- data.frame(Gene = rownames(scObj),
                     nCells = rowSums(expMat != 0))
    df <- df[order(df$nCells, decreasing=TRUE), ]
    if (!is.null(minCutoff))
        df <- df[df$nCells >= minCutoff, ]
    if (!is.null(maxCutoff))
        df <- df[df$nCells <= maxCutoff, ]
    return(df)
}
