#' @importFrom grDevices dev.new dev.off hcl.colors
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
