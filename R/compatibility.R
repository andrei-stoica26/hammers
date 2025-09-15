#' @importFrom SeuratObject Embeddings Reductions
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom S4Vectors DataFrame
#'
NULL

#' @rdname metadataDF
#' @export
#'
metadataDF.default <- function(scObj)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname metadataDF
#' @export
#'
`metadataDF<-.default` <- function(scObj, value)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname metadataDF
#' @export
#'
metadataDF.Seurat <- function(scObj)
    return(scObj[[]])

#' @rdname metadataDF
#' @export
#'
`metadataDF<-.Seurat` <- function(scObj, value) {
    if (!is.data.frame(value))
        stop('`value` must be a data.frame.')
    scObj[[]] <- value
    return(scObj)
}

#' @rdname metadataDF
#' @export
#'
metadataDF.SingleCellExperiment <- function(scObj)
    return(as.data.frame(colData(scObj)))

#' @rdname metadataDF
#' @export
#'
`metadataDF<-.SingleCellExperiment` <- function(scObj, value) {
    if (!is.data.frame(value))
        stop('`value` must be a data.frame.')
    colData(scObj) <- DataFrame(value)
    return(scObj)
}

###############################################################################
#' Return metadata names
#'
#' This function extracts metadata names from a Seurat or
#' SingleCellExperiment object. It can also be used to modify metadata names.
#'
#' @inheritParams metadataDF
#'
#' @return The names of the metadata columns.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' colNames <- metadataNames(scObj)
#'
#' @export
#'
metadataNames <- function(scObj)
    return(colnames(metadataDF(scObj)))

#' @rdname metadataNames
#' @export
`metadataNames<-` <- function(scObj, value){
    if (!is.character(value))
        stop('`value` must be a character.')
    colnames(metadataDF(scObj)) <- value
    return(scObj)
}

###############################################################################
#' @rdname scCol
#' @export
#'
scCol.default <- function(scObj, col)
    stop('Unrecognized input type: scObj must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname scCol
#' @export
#'
scCol.Seurat <- function(scObj, col)
    return(scObj[[]][[col]])

#' @rdname scCol
#' @export
#'
scCol.SingleCellExperiment <- function(scObj, col)
    return(scObj[[col]])

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
#' Extract count information from single-cell expression object column
#'
#' This function extracts count information from the column of a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataDF
#' @param col Column as string.
#'
#' @return A frequency vector with the unique column values as names.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' scColCounts(scObj, 'label')
#'
#' @export
#'
scColCounts <- function(scObj, col='orig.ident'){
    df <- dplyr::count(metadataDF(scObj), .data[[col]])
    v <- setNames(df[, 2], as.factor(df[, 1]))
    return(v)
}

#' Extract count information from Seurat column
#'
#' This function extracts count information from Seurat column.
#'
#' @inheritParams metadataDF
#' @param col1 Column as string.
#' @param col2 Column as string.
#'
#' @return A data frame listing the counts of all combinations of pairs from
#' two categorical columns.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' scColPairCounts(scObj, 'donor', 'label')
#'
#' @export
#'
scColPairCounts <- function(scObj, col1='seurat_clusters', col2='orig.ident')
    return(dplyr::count(metadataDF(scObj), .data[[col1]], .data[[col2]]))


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
#' scObj <- scater::runPCA(scObj)
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
#' scObj <- scater::runUMAP(scObj)
#' umapMat <- scUMAPMat(scObj)
#'
#' @export
#'
scUMAPMat <- function(scObj)
    return(scDimredMat(scObj, 'umap'))

