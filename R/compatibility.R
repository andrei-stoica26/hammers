#' @importFrom SeuratObject Embeddings Reductions
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom S4Vectors DataFrame
#'
NULL

#' @rdname dimredNames
#' @export
#'
dimredNames.default <- function(scObj)
    stop('Unrecognized input type: `scObj` must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object.')

#' @rdname dimredNames
#' @export
#'
dimredNames.Seurat <- function(scObj)
    return(Reductions(scObj))


#' @rdname dimredNames
#' @export
#'
dimredNames.SingleCellExperiment <- function(scObj)
    return(names(reducedDims(scObj)))

###############################################################################
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
`metadataDF<-.Seurat` <- function(scObj, value){
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
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' colNames <- metadataNames(sceObj)
#'
#' @export
#'
metadataNames <- function(scObj)
    return(colnames(metadataDF(scObj)))

#' @rdname metadataNames
#' @export
#'
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
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname scCol
#' @export
#'
`scCol<-.default` <- function(scObj, col, value)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname scCol
#' @export
#'
scCol.Seurat <- function(scObj, col)
    return(scObj[[]][[col]])

#' @rdname scCol
#' @export
#'
`scCol<-.Seurat` <- function(scObj, col, value){
    scObj[[]][[col]] <- value
    return(scObj)
}

#' @rdname scCol
#' @export
#'
scCol.SingleCellExperiment <- function(scObj, col)
    return(scObj[[col]])

#' @rdname scCol
#' @export
#'
`scCol<-.SingleCellExperiment` <- function(scObj, col, value){
    scObj[[col]] <- value
    return(scObj)
}

################################################################################
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

#' @param dataType Expression data type. Must be one of 'counts' and
#' 'logcounts'.
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

################################################################################
#' @rdname scExpMat
#' @export
#'
scExpMat.default <- function(scObj,
                             dataType = c('data',
                                          'counts',
                                          'logcounts'),
                             genes = NULL,
                             densify = TRUE)
    stop('Unrecognized input type: `scObj` must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object.')

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
#' @rdname scDimredMat
#' @export
#'
scDimredMat.default <- function(scObj, dimred)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname scDimredMat
#' @export
#'
scDimredMat.Seurat <- function(scObj, dimred)
{
    reductions <- dimredNames(scObj)
    dimred <- dimredName(scObj, dimred)
    return(as.matrix(Embeddings(scObj, reduction=dimred)))
}

#' @rdname scDimredMat
#' @export
#'
scDimredMat.SingleCellExperiment <- function(scObj, dimred){
    reductions <- dimredNames(scObj)
    dimred <- dimredName(scObj, dimred)
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
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' scColCounts(sceObj, 'Mutation_Status')
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
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' scColPairCounts(sceObj, 'Mutation_Status', 'Cell_Cycle')
#'
#' @export
#'
scColPairCounts <- function(scObj, col1='seurat_clusters', col2='orig.ident')
    return(dplyr::count(metadataDF(scObj), .data[[col1]], .data[[col2]]))


#' Extract dimensionality reduction matrix from single-cell expression object
#'
#' This function extracts a dimensionality reduction matrix from a Seurat
#' or SingleCellExpression object.
#'
#' @inheritParams scDimredMat
#'
#' @return A dimensionality reduction matrix.
#'
#' @noRd
#'
dimredName <- function(scObj, dimred){
    lowerDimred <- tolower(dimred)
    upperDimred <- toupper(dimred)
    dimName <- intersect(dimredNames(scObj), c(dimred, lowerDimred,
                                               upperDimred))
    if (!length(dimName))
        stop('No `', lowerDimred, '` or `', upperDimred,
             '` reduction was found in the object.')
    if (length(dimName) > 1){
        warning('Both `', lowerDimred, '` and ', upperDimred, '` reductions were ',
                'found. Returning ', dimred, '.')
        return(dimred)
    }
    return(dimName)
}

#' Extracts the PCA matrix from object.
#'
#' This function extracts the PCA matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataDF
#'
#' @return A PCA matrix.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' pcaMat <- scPCAMat(sceObj)
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
#' @inheritParams metadataDF
#'
#' @return A UMAP matrix.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='hammers')
#' sceObj <- qs::qread(scePath)
#' umapMat <- scUMAPMat(sceObj)
#'
#' @export
#'
scUMAPMat <- function(scObj)
    return(scDimredMat(scObj, 'umap'))

