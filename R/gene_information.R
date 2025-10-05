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
#' @param maxCutoff Maximum cutoff for gene counts. Genes with counts above
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
    df <- data.frame(Gene = rownames(expMat),
                     nCells = rowSums(expMat != 0))
    df <- df[order(df$nCells, decreasing=TRUE), ]
    if (!is.null(minCutoff))
        df <- df[df$nCells >= minCutoff, ]
    if (!is.null(maxCutoff))
        df <- df[df$nCells <= maxCutoff, ]
    return(df)
}

#' Generates cell sets for each input gene
#'
#' This function constructs, for each input gene, sets of cells expressing
#' the gene
#'
#' @inheritParams scExpMat
#'
#' @return A named list of character vectors of cell names.
#'
#' @examples
#' mat <- matrix(0, 1000, 500)
#' rownames(mat) <- paste0('G', seq(1000))
#' colnames(mat) <- paste0('C', seq(500))
#' mat[sample(length(mat), 70000)] <- sample(50, 70000, TRUE)
#' mat <- mat[paste0('G', sample(1000, 3)), ]
#' geneCellSets(mat)
#'
#' @export
#'
geneCellSets <- function(scObj, genes = NULL){
    expMat <- scExpMat(scObj, 'counts', genes)
    genes <- rownames(expMat)
    cellSets <- setNames(lapply(genes, function(x){
        geneExp <- expMat[x, ]
        geneExp  <- geneExp[geneExp > 0]
        return(names(geneExp))
    }), genes)
    cellSets <- cellSets[vapply(cellSets, length, numeric(1)) > 0]
    return(cellSets)
}
