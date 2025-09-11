#' @importFrom SeuratObject Embeddings Reductions
#' @importFrom SingleCellExperiment colData reducedDim reducedDims
#' @importFrom SummarizedExperiment assay
#'
NULL

#' @rdname metadataDF
#' @export
#'
metadataDF.default <- function(scObj)
    stop('Unrecognized input type: scObj must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname metadataDF
#' @export
#'
metadataDF.Seurat <- function(scObj)
    return(scObj[[]])

#' @rdname metadataDF
#' @export
#'
metadataDF.SingleCellExperiment <- function(scObj)
    return(data.frame(colData(scObj)))

###############################################################################
#' @rdname metadataNames
#' @export
#'
metadataNames.default <- function(scObj)
    stop('Unrecognized input type: scObj must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname metadataNames
#' @export
#'
metadataNames.Seurat <- function(scObj)
    return(colnames(scObj[[]]))

#' @rdname metadataNames
#' @export
#'
metadataNames.SingleCellExperiment <- function(scObj)
    return(names(colData(scObj)))

###############################################################################
#' @rdname scCol
#' @export
#'
scCol.default <- function(scObj, colStr)
    stop('Unrecognized input type: scObj must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname scCol
#' @export
#'
scCol.Seurat <- function(scObj, colStr)
    return(scObj[[]][[colStr]])

#' @rdname scCol
#' @export
#'
scCol.SingleCellExperiment <- function(scObj, colStr)
    return(scObj[[colStr]])

###############################################################################
#' @rdname scGeneExp
#' @export
#'
scGeneExp.default <- function(scObj, gene, dataType = c('counts',
                                                        'data',
                                                        'logcounts'))
    stop('Unrecognized input type: scObj must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object')

#' @param dataType Expression data type. Must be one of 'counts' and 'data'.
#'
#' @rdname scGeneExp
#' @export
#'
scGeneExp.Seurat <- function(scObj,
                             gene,
                             dataType = c('counts',
                                          'data',
                                          'logcounts')){
    dataType <- match.arg(dataType, c('counts', 'data', 'logcounts'))
    if (dataType == 'logcounts')
        dataType <- 'data'
    return(LayerData(scObj, dataType)[gene, ])
}

#' @param dataType Expression data type. Must be one of 'counts' and 'logcounts'.
#'
#' @rdname scGeneExp
#' @export
#'
scGeneExp.SingleCellExperiment <- function(scObj,
                                           gene,
                                           dataType = c('counts',
                                                        'data',
                                                        'logcounts')){
    dataType <- match.arg(dataType, c('counts', 'data', 'logcounts'))
    if (dataType == 'data')
        dataType <- 'logcounts'
    return(assay(scObj, dataType)[gene, ])
}

#' @rdname scGeneExp
#' @export
#'
scGeneExp.dgCMatrix <- function(scObj,
                                gene,
                                dataType = c('counts',
                                             'data',
                                             'logcounts'))
    return(scObj[gene, ])

#' @rdname scGeneExp
#' @export
#'
scGeneExp.matrix <- function(scObj,
                             gene,
                             dataType = c('counts',
                                          'data',
                                          'logcounts'))
    return(scObj[gene, ])

###############################################################################
#' @rdname scExpMat
#' @export
#'
scExpMat.default <- function(scObj,
                             dataType = c('data',
                                          'counts',
                                          'logcounts'),
                             genes = NULL,
                             densify = TRUE)
    stop('Unrecognized input type: scObj must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object')

#' @rdname scExpMat
#' @export
#'
scExpMat.Seurat <- function(scObj,
                            dataType = c('data',
                                         'counts',
                                         'logcounts'),
                            genes = NULL,
                            densify = TRUE){
    dataType <- match.arg(dataType, c('data', 'counts', 'logcounts'))
    if (dataType == 'logcounts')
        dataType <- 'data'
    mat <- LayerData(scObj, layer=dataType)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.SingleCellExperiment <- function(scObj,
                                          dataType = c('data',
                                                       'counts',
                                                       'logcounts'),
                                          genes = NULL,
                                          densify = TRUE){

    dataType <- match.arg(dataType, c('data', 'counts', 'logcounts'))
    if (dataType == 'data')
        dataType <- 'logcounts'
    mat <- assay(scObj, dataType)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.dgCMatrix <- function(scObj,
                               dataType = c('data',
                                            'counts',
                                            'logcounts'),
                               genes = NULL,
                               densify = TRUE){

    if (!is.null(genes))
        mat <- scObj[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.matrix <- function(scObj,
                               dataType = c('data',
                                            'counts',
                                            'logcounts'),
                               genes = NULL,
                               densify = TRUE){

    if (!is.null(genes))
        scObj <- scObj[genes, ]
    return(scObj)
}

###############################################################################
#' @noRd
#'
scDimredMat.default <- function(scObj, dimred = c('pca', 'umap'))
    stop('Unrecognized input type: scObj must be a Seurat or ',
         'SingleCellExpression object')

#' @noRd
#'
scDimredMat.Seurat <- function(scObj, dimred = c('pca', 'umap'))
{
    dimred <- match.arg(dimred, c('pca', 'umap'))
    reductions <- Reductions(scObj)
    if(!dimred %in% reductions){
        dimred <- toupper(dimred)
        if(!dimred %in% reductions)
            stop(dimred, ' reduction not found in Seurat object.')
    }
    return(as.matrix(Embeddings(scObj, reduction=dimred)))
}

#' @noRd
#'
scDimredMat.SingleCellExperiment <- function(scObj, dimred = c('pca', 'umap')){
    dimred <- match.arg(dimred, c('pca', 'umap'))
    reductions <- names(reducedDims(scObj))
    if(!dimred %in% reductions){
        dimred <- toupper(dimred)
        if(!dimred %in% reductions)
            stop(dimred, ' reduction not found in SingleCellExperiment object.')
    }
    return(reducedDim(scObj, dimred))
}

###############################################################################
#' Extracts the PCA matrix from object.
#'
#' This function extracts the PCA matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#'
#' @return A PCA matrix.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runPCA
#' pcaMat <- scPCAMat(scObj)
#'
#' @export
#'
scPCAMat <- function(scObj)
    return(scDimredMat(scObj, 'pca'))

#' Extracts the UMAP matrix from object.
#'
#' This function extracts the UMAP matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataNames
#'
#' @return A UMAP matrix.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' scObj <- scuttle::logNormCounts(scObj)
#' scObj <- scater::runUMAP
#' umapMat <- scPCAMat(scObj)
#'
#' @export
#'
scUMAPMat <- function(scObj)
    return(scDimredMat(scObj, 'umap'))

