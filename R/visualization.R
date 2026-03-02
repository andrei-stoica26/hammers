#' @importFrom ggplot2 aes element_text ggplot geom_bar geom_point labs scale_fill_manual theme theme_classic
#' @importFrom ggrepel geom_text_repel
#' @importFrom henna centerTitle hullPlot riverPlot
#' @importFrom rlang .data
#'
NULL

#' Plot the distribution of cells across two columns
#'
#' This function plots the distribution of cells across two columns.
#'
#' @inheritParams repAnalysis
#' @param title Plot title.
#' @param type Whether the plot should display counts ('counts', default) or
#' percentages ('percs').
#' @param xLab x axis label.
#' @param yLab y axis label.
#' @param legendLab Legend label.
#' @param palette Color palette.
#' @param legendPos Legend position.
#' @param legendTitleSize Legend title size.
#' @param legendTextSize Legend text size.
#' @param axisTextSize Axis text size.
#' @param axisTitleSize Axis title size.
#' @param sigDigits Number of significant digits used by percentages displayed
#' on the plot. Ignored if \code{type} is 'counts'.
#'
#' @return A ggplot object.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' p <- distributionPlot(sceObj, col1='Cluster', col2='Donor')
#'
#' @export
#'
distributionPlot <- function(scObj,
                             title = NULL,
                             col1 = 'seurat_clusters',
                             col2 = 'orig.ident',
                             type = c('counts', 'percs'),
                             xLab = col1,
                             yLab = if(type == 'counts') 'Count' else
                                 'Percentage',
                             legendLab = col2,
                             palette = 'Spectral',
                             legendPos = 'right',
                             legendTextSize = 10,
                             legendTitleSize = 10,
                             axisTextSize = 12,
                             axisTitleSize = 12,
                             sigDigits = 2){
    type <- match.arg(type, c('counts', 'percs'))

    if (type == 'counts'){
        df <- scColPairCounts(scObj, col1, col2)
        colIndex <- 3
    } else {
        df <- scColPairPercs(scObj, col1, col2, sigDigits)
        colIndex <- 4
    }

    nColors <- length(unique(df[, 2]))
    p <- ggplot(df) +
        geom_bar(position='stack',
                 stat='identity',
                 aes(x=.data[[names(df)[1]]],
                     fill=.data[[names(df)[2]]],
                     y=.data[[names(df)[colIndex]]])) +
        theme_classic() + labs(x=xLab, y=yLab, fill=legendLab) +
        scale_fill_manual(values=hcl.colors(nColors, palette)) +
        theme(legend.position=legendPos,
              legend.text=element_text(size=legendTextSize),
              legend.title=element_text(size=legendTitleSize),
              axis.text=element_text(size=axisTextSize),
              axis.title=element_text(size=axisTitleSize))
    p <- centerTitle(p, title)
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
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' df <- repAnalysis(sceObj, 'Cluster', 'Donor')
#' pvalRiverPlot(df)
#'
#' @export
#'
pvalRiverPlot <- function(df,
                          pvalCol = 'pvalAdj',
                          colIndices = c(1, 2),
                          weightExp = 1/2,
                          pvalOffset = 1e-317,
                          ...){
    resDF <- prepAlluvial(df, pvalCol, colIndices, weightExp, pvalOffset)
    p <- riverPlot(resDF, ...)
    return(p)
}

#' Create a single-cell dimensionality reduction plot with added labeled points
#'
#' This function creates a single-cell dimensionality reduction plot with
#' added labeled points.
#'
#' @details A wrapper around \code{scLang::dimPlot}.
#'
#' @inheritParams geneCenters
#' @inheritParams distributionPlot
#' @param pointsObj A data frame or matrix of points with two columns
#' representing x and y coordinates.
#' @param alpha Opaqueness level.
#' @param pointShape Point shape.
#' @param pointSize Point size.
#' @param pointColor Point color.
#' @param labelSize Label size. If \code{NULL}, the points will not receive
#' labels.
#' @param maxOverlaps Maximum overlaps.
#' @param ... Additional parameters passed to \code{Seurat::DimPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' pointsObj <- data.frame(x = c(2, 3),
#' y = c(1, 0),
#' row.names = c('P1', 'P2'))
#' pointsDimPlot(sceObj, pointsObj=pointsObj)
#'
#' @export
#'
pointsDimPlot <- function(scObj,
                          title = NULL,
                          pointsObj = NULL,
                          alpha = 1,
                          pointShape = 4,
                          pointSize = 2,
                          pointColor = 'black',
                          labelSize = 2.5,
                          maxOverlaps = 30,
                          ...){
    p <- suppressWarnings(dimPlot(scObj, ...))
    if(!is.null(pointsObj)){
        p <- p + geom_point(aes(.data[[colnames(pointsObj)[1]]],
                                .data[[colnames(pointsObj)[2]]]),
                            data=pointsObj,
                            alpha=alpha,
                            shape=pointShape,
                            size=pointSize,
                            color=pointColor)
        if(!is.null(labelSize))
            p <- p + geom_text_repel(aes(.data[[colnames(pointsObj)[1]]],
                                         .data[[colnames(pointsObj)[2]]]),
                                     data=pointsObj,
                                     label=rownames(pointsObj),
                                     size=labelSize,
                                     max.overlaps=maxOverlaps)
    }
    p <- centerTitle(p, title)
    return(p)
}

#' Create a single-cell dimensionality reduction plot with added labeled points
#' for genes
#'
#' This function creates a single-cell dimensionality reduction plot with
#' added labeled points for genes.
#'
#' @inheritParams geneCenters
#' @param genes Genes whose centers of mass will be plotted.
#' @param ... Additional parameters passed to \code{pointsDimPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' genesDimPlot(sceObj, c('Gene_0364', 'Gene_0388', 'Gene_0477'))
#'
#' @export
#'
genesDimPlot <- function(scObj, genes, ...){
    centersDF <- geneCenters(scObj, genes)
    return(pointsDimPlot(scObj, pointsObj=centersDF, ...))
}

#' Create a single-cell dimensionality reduction plot with added labeled points
#' for numeric columns
#'
#' This function creates a single-cell dimensionality reduction plot with
#' added labeled points for metadata numeric columns.
#'
#' @inheritParams geneCenters
#' @param cols Columns whose centers of mass will be plotted.
#' @param ... Additional parameters passed to \code{pointsDimPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' colsDimPlot(sceObj, c('sizeFactor', 'silhouette'))
#'
#' @export
#'
colsDimPlot <- function(scObj, cols, ...){
    centersDF <- colCenters(scObj, cols)
    return(pointsDimPlot(scObj, pointsObj=centersDF, ...))
}
