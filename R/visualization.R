#' @importFrom henna riverPlot
#'
NULL

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
