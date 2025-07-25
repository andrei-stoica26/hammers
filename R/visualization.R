#' @importFrom grDevices dev.new dev.off
#' @importFrom henna riverPlot
#'
NULL

#' @rdname devPlot
#' @export
#'
devPlot.default <- function(plotObject, ...)
    stop('Unrecognized input type: plotObject must be a function,',
         ' a ggplot object or a list of ggplot objects')

#' @rdname devPlot
#' @export
#'
devPlot.function <- function(plotObject, ...){
    dev.new(noRStudioGD = TRUE)
    print(plotObject(...))
    dev.off()
}

#' @rdname devPlot
#' @export
#'
devPlot.ggplot <- function(plotObject, ...)
    devPlot.function(identity, plotObject)

#' @rdname devPlot
#' @export
#'
devPlot.list <- function(plotObject, ...)
    invisible(lapply(plotObject, devPlot.ggplot))

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
pvalRiverPlot <- function(df, ...){
    resDF <- prepAlluvial(df)
    p <- riverPlot(resDF, ...)
    return(p)
}
