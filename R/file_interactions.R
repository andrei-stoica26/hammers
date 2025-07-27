#' @importFrom grDevices dev.new dev.off
#' @importFrom qs qread
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

#' Read and delete a .qs file
#'
#' This functions reads a .qs file, deletes it, and returns its content.
#'
#' @param qsFile Name of .qs file with path.
#'
#' @return The content of the .qs file.
#'
#' @examples
#' library(qs)
#' qsave(c(1, 2, 3), 'temp.qs')
#' qGrab('temp.qs')
#'
#' @export
#'
qGrab <- function(qsFile){
    res <- qread(qsFile)
    file.remove(qsFile)
    return(res)
}
