#' @importFrom henna riverPlot
#'
NULL

#' Plot representation data frame
#'
#' This function plots representation data frame as an alluvial plot.
#'
#' @param df A representation data frame.
#' @param ... Additional parameters passed to \code{henna::riverPlot}
#'
#' @return A ggplot object
#'
#' @export
#'
repPlot <- function(df, ...){
    resDF <- prepAlluvial(df)
    p <- riverPlot(resDF, ...)
    return(p)
}
