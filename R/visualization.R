#' @importFrom henna centerTitle riverPlot
#'
NULL

#' Plot the distribution of cells across two columns
#'
#' This function plots the distribution of cells across two columns.
#'
#' @inheritParams metadataDF
#' @param plotTitle Plot title.
#' @param col1 Column 1
#' @param col2 Column 2.
#' @param xLab x axis label.
#' @param yLab y axis label.
#' @param legendLab Legend label.
#' @param palette Color palette.
#'
#' @export
#'
distributionPlot <- function(scObj,
                             plotTitle = 'Distribution plot',
                             col1 = seurat_clusters,
                             col2 = orig.ident,
                             xLab = 'Column 1',
                             yLab = 'Count',
                             legendLab = 'Column 2',
                             palette = 'Spectral'){
    df <- dplyr::count(metadataDF(scObj), {{col1}}, {{col2}})
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
#' @export
#'
pvalRiverPlot <- function(df, weightExp = 1/2, ...){
    resDF <- prepAlluvial(df)
    p <- riverPlot(resDF, weightExp, ...)
    return(p)
}


