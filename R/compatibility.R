#' @importFrom SeuratObject Embeddings
#' @importFrom SingleCellExperiment colData reducedDim
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
#' @rdname scExpMat
#' @export
#'
scExpMat.default <- function(scObj,
                             dataType = c('counts',
                                          'data',
                                          'logcounts'),
                             genes = NULL,
                             densify = TRUE)
    stop('Unrecognized input type: scObj must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object')

#' @param dataType Expression data type. Must be one of 'counts' and 'data'.
#'
#' @rdname scExpMat
#' @export
#'
scExpMat.Seurat <- function(scObj,
                            dataType = c('counts',
                                         'data',
                                         'logcounts'),
                            genes = NULL,
                            densify = TRUE){
    dataType <- match.arg(dataType, c('counts', 'data', 'logcounts'))
    if (dataType == 'logcounts')
        dataType <- 'data'
    mat <- LayerData(scObj, layer=dataType)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @param dataType Expression data type. Must be one of 'counts'
#' and 'logcounts'.
#'
#' @rdname scExpMat
#' @export
#'
scExpMat.SingleCellExperiment <- function(scObj,
                                          dataType = c('counts',
                                                       'data',
                                                       'logcounts'),
                                          genes = NULL,
                                          densify = TRUE){

    dataType <- match.arg(dataType, c('counts', 'data', 'logcounts'))
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
                               dataType = c('counts',
                                            'data',
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
                               dataType = c('counts',
                                            'data',
                                            'logcounts'),
                               genes = NULL,
                               densify = TRUE){

    if (!is.null(genes))
        scObj <- scObj[genes, ]
    return(scObj)
}

###############################################################################
#' @rdname scPCAMat
#' @export
#'
scPCAMat.default <- function(scObj)
    stop('Unrecognized input type: scObj must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname scPCAMat
#' @export
#'
scPCAMat.Seurat <- function(scObj)
    return(as.matrix(Embeddings(scObj, reduction="pca")))

#' @rdname scPCAMat
#' @export
#'
scPCAMat.SingleCellExperiment <- function(scObj)
    return(reducedDim(scObj, 'PCA'))

