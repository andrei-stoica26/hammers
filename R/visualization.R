#' @importFrom ggplot2 aes ggplot geom_bar geom_point labs scale_fill_manual theme_classic
#' @importFrom ggrepel geom_text_repel
#' @importFrom henna centerTitle riverPlot
#' @importFrom rlang .data
#' @importFrom Seurat DimPlot
#'
NULL

#' Plot the distribution of cells across two columns
#'
#' This function plots the distribution of cells across two columns.
#'
#' @inheritParams scColPairCounts
#' @param plotTitle Plot title.
#' @param xLab x axis label.
#' @param yLab y axis label.
#' @param legendLab Legend label.
#' @param palette Color palette.
#'
#' @return A ggplot object.
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' distributionPlot(scObj, col1='donor', col2='label')
#'
#' @export
#'
distributionPlot <- function(scObj,
                             plotTitle = 'Distribution plot',
                             col1 = 'seurat_clusters',
                             col2 = 'orig.ident',
                             xLab = 'Column 1',
                             yLab = 'Count',
                             legendLab = 'Column 2',
                             palette = 'Spectral'){
    df <- scColPairCounts(scObj, col1, col2)
    nColors <- length(unique(df[, 2]))
    p <- ggplot(df) +
        geom_bar(position='stack',
                 stat='identity',
                 aes(x=.data[[names(df)[1]]],
                     fill=.data[[names(df)[2]]],
                     y=.data[[names(df)[3]]])) +
        theme_classic() + labs(x=xLab, y=yLab, fill=legendLab) +
        scale_fill_manual(values=hcl.colors(nColors, palette))
    p <- centerTitle(p, plotTitle)
    return(p)
}

#' Plot representation data frame
#'
#' This function plots representation data frame as an alluvial plot.
#'
#' @inheritParams prepAlluvial
#' @param ... Additional parameters passed to \code{henna::riverPlot}
#'
#' @return A ggplot object
#'
#' @examples
#' scObj <- scRNAseq::BaronPancreasData('human')
#' df <- repAnalysis(scObj, 'donor', 'label')
#' pvalRiverPlot(df)
#'
#' @export
#'
pvalRiverPlot <- function(df, weightExp = 1/2, ...){
    resDF <- prepAlluvial(df, weightExp=weightExp)
    p <- riverPlot(resDF, ...)
    return(p)
}

#' Plot Seurat DimPlot with added labeled points
#'
#' This function plots a Seurat DimPlot with added labeled points.
#'
#' @param seuratObj A Seurat object.
#' @param plotTitle Plot title.
#' @param pointsDF A data frame of points with two columns representing the x
#' and y coordinates.
#' @param pointShape Point shape.
#' @param pointSize Point size.
#' @param labelSize Label size.
#' @param maxOverlaps Maximum overlaps.
#' @param ... Additional parameters passed to \code{Seurat::DimPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sceObj <- scRNAseq::BaronPancreasData('human')
#' sceObj <- scuttle::logNormCounts(sceObj)
#' seuratObj <- suppressWarnings(as.Seurat(sceObj))
#' seuratObj <- Seurat::FindVariableFeatures(seuratObj)
#' seuratObj <- Seurat::ScaleData(seuratObj)
#' seuratObj <- Seurat::RunPCA(seuratObj)
#' seuratObj <- suppressWarnings(Seurat::RunUMAP(seuratObj, dims=1:15))
#' pointsDF <- data.frame(x = c(2, 3),
#' y = c(1, 6),
#' row.names = c('P1', 'P2'))
#' pointsDimPlot(seuratObj, pointsDF=pointsDF)
#'
#' @export
#'
pointsDimPlot <- function(seuratObj,
                          plotTitle = 'Dim plot',
                          pointsDF = NULL,
                          pointShape = 4,
                          pointSize = 2,
                          labelSize = 2.5,
                          maxOverlaps = 30,
                          ...){
    p <- DimPlot(seuratObj, ...)
    if(!is.null(pointsDF))
        p <- p + geom_point(aes(.data[[names(pointsDF)[1]]],
                                .data[[names(pointsDF)[2]]]),
                            data=pointsDF,
                            shape=pointShape,
                            size=pointSize) +
            geom_text_repel(aes(.data[[names(pointsDF)[1]]],
                                .data[[names(pointsDF)[2]]]),
                            data=pointsDF,
                            label=rownames(pointsDF),
                            size=labelSize,
                            max.overlaps=maxOverlaps)
    p <- centerTitle(p, plotTitle)
    return(p)
}

#' Plot Seurat DimPlot with added labeled points for genes
#'
#' This function plots a Seurat DimPlot with added labeled points for genes.
#'
#' @param seuratObj A Seurat object.
#' @param genes Genes whose centers of mass will be plotted.
#' @param ... Additional parameters passed to \code{pointsDimPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sceObj <- scRNAseq::BaronPancreasData('human')
#' sceObj <- scuttle::logNormCounts(sceObj)
#' seuratObj <- suppressWarnings(as.Seurat(sceObj))
#' seuratObj <- Seurat::FindVariableFeatures(seuratObj)
#' seuratObj <- Seurat::ScaleData(seuratObj)
#' seuratObj <- Seurat::RunPCA(seuratObj)
#' seuratObj <- suppressWarnings(Seurat::RunUMAP(seuratObj, dims=1:15))
#' genesDimPlot(seuratObj, c('AURKA', 'TOP2A', 'MKI67'))
#'
#' @export
#'
genesDimPlot <- function(seuratObj, genes, ...){
    centersDF <- geneCenters(seuratObj, genes)
    return(pointsDimPlot(seuratObj, pointsDF=centersDF, ...))
}

#' Plot Seurat DimPlot with added labeled points for numeric columns
#'
#' This function plots a Seurat DimPlot with added labeled points for metadata
#' numeric columns.
#'
#' @param seuratObj A Seurat object.
#' @param cols Genes whose centers of mass will be plotted.
#' @param ... Additional parameters passed to \code{pointsDimPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sceObj <- scRNAseq::BaronPancreasData('human')
#' sceObj <- scuttle::logNormCounts(sceObj)
#' seuratObj <- suppressWarnings(as.Seurat(sceObj))
#' seuratObj <- Seurat::FindVariableFeatures(seuratObj)
#' seuratObj <- Seurat::ScaleData(seuratObj)
#' seuratObj <- Seurat::RunPCA(seuratObj)
#' seuratObj <- suppressWarnings(Seurat::RunUMAP(seuratObj, dims=1:15))
#' colsDimPlot(seuratObj, c('nCount_originalexp', 'nFeature_originalexp'))
#'
#' @export
#'
colsDimPlot <- function(seuratObj, cols, ...){
    centersDF <- colCenters(seuratObj, cols)
    return(pointsDimPlot(seuratObj, pointsDF=centersDF, ...))
}

