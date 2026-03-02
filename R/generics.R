#' Saves plot or list of plots
#'
#' This function saves a plot or list of plots as a pdf. Can also take as input
#' a function that returns a ggplot object together with its arguments.
#'
#' @param plotObject A function, ggplot object, or list of ggplot objects.
#' @param ... Additional arguments.
#'
#' @return No value. This function is called for its side effect.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x = c(1, 2), y = c(3, 5))
#' p <- ggplot(df) + geom_point(aes(x, y))
#' devPlot(p)
#'
#' simplePlot <- function(df, title)
#'     return(ggplot(df) + geom_point(aes(x, y)) + ggtitle(title))
#'
#' devPlot(simplePlot, df, 'Plot title')
#'
#' @export
#'
#'
devPlot <- function(plotObject, ...)
    UseMethod(generic='devPlot', object=plotObject)
