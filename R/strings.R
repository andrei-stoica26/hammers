#' Join all combinations of elements from character vectors
#'
#' This function joins all combinations of elements from character vectors
#' with a separating character.
#'
#' @param ... Vectors passed to \code{expandGrid}.
#' @param joinChar Character used to join combinations.
#'
#' @return A character vector.
#'
#' @examples
#' joinCharCombs(c('a', 'b', 'c', 'd'), c('eee', 'ff'), c(1, 2, 3))
#'
#' @export
#'
joinCharCombs <- function(..., joinChar = '_')
    return(sort(apply(expand.grid(...), 1,
               function(x) paste0(x, collapse = joinChar))))
