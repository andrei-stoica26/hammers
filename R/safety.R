#' Perform min-max normalization when possible; otherwise return a single-value
#' vector.
#'
#' This function min-max-normalizes a vector when possible, and otherwise returns
#' a single-value vector.
#'
#' @param scores Numeric vector.
#' @param safeVal Value to replace all values with when all values in the vector
#' are the same.
#'
#' @return Min-max-normalized scores or a single-value vector.
#' @examples
#' safeMinmax(c(0, 3, 2, 1, 4, 5.5, 6.32, 8, 1.1))
#'
#' @export
#'
safeMinmax <- function(scores, safeVal = 0){
    if(length(unique(scores)) < 2)
        return(rep(safeVal, length(scores)))
    return(liver::minmax(scores))
}

#' Message an input if verbose is set to TRUE
#'
#' This function messages an input if \code{verbose} is set to \code{TRUE}.
#'
#' @param msg Message
#' @param verbose Whether the message should be displayed.
#'
#' @return No return value. This function is called for its side effect
#' (messaging the input if \code{verbose} is set to \code{TRUE}).
#'
#' @examples
#' safeMessage('message')
#'
#' @export
#'
safeMessage <- function(msg, verbose = TRUE)
    if(verbose)
        message(msg)
