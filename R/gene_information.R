#' Extract gene counts from a Seurat or SingleCellExperiment object
#'
#' This function extracts gene counts from a Seurat or
#' SingleCellExperiment object as a data frame.
#'
#' @inheritParams metadataDF
#' @param minCutoff Minimum cutoff for gene counts. Genes with counts below
#' this value will be omitted.
#' @param minCutoff Maximum cutoff for gene counts. Genes with counts above
#' this value will be omitted.
#'
#' @return A data frame of gene counts
#'
#' @examples
#' scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=200))
#' df <- geneCounts(scObj)
#'
#' @export
#'
geneCounts <- function(scObj, minCutoff = NULL, maxCutoff = NULL){
    expMat <- scExpMat(scObj, 'counts')
    df <- data.frame(Gene = rownames(scObj),
                     Count = rowSums(expMat != 0))
    df <- df[order(df$Count, decreasing=TRUE), ]
    if (!is.null(minCutoff))
        df <- df[df$Count >= minCutoff, ]
    if (!is.null(maxCutoff))
        df <- df[df$Count <= maxCutoff, ]
    return(df)
}
