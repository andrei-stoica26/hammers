#' Create a hash map
#'
#' This function creates a hash map.
#'
#' @param hashKeys A list of hash keys. If vectors are part of the hash keys,
#' each of their elements will be assigned the corresponding value.
#' @param hashValues A vector of hash values.
#' Must have the same length as \code{hashKeys}.
#'
#' @return A named vector.
#'
#' @examples
#' hashMap(list(2, c(3, 4, 5), 6, 8), c('a', 'b', 'c', 'd'))
#'
#' @export
#'
hashMap <- function(hashKeys, hashValues){
    if (length(hashKeys) != length(hashValues))
        stop('`hashKeys` and `hashValues` must have the same length.')
    if (max(table(unlist(hashKeys))) > 1)
        stop('All values in `hashKeys` must be unique. ',
             names(which.max(table(unlist(hashKeys)))), ' is repeated.')
    hashNames <- unlist(hashKeys)
    hash <- unlist(mapply(function(x, y) rep(y, length(x)), hashKeys,
                          hashValues))
    names(hash) <- hashNames
    return(hash)
}

#' Add a categorical column to a data frame based on another column
#'
#' This function adds a categorical column to a data frame based on another
#' column.
#'
#' @param df A data frame.
#' @param col Column whose values will be used for creating the new column.
#' @param newCol Column to be added.
#' @inheritParams hashMap
#'
#' @return A data frame with a new categorical column.
#'
#' @examples
#' df <- data.frame(fruit = c('apple', 'banana', 'cherry', 'grape'))
#' df <- addCategory(df,
#'                 'fruit',
#'                 'color',
#'                 list(c('apple', 'cherry'),
#'                 'banana',
#'                 'grape'),
#'                 c('red', 'yellow', 'purple'))
#'
#' @export
#'
addCategory <- function(df, col, newCol, hashKeys, hashValues){
    remainder <- setdiff(unique(df[[col]]), unlist(hashKeys))
    if (length(hashKeys) + 1 == length(hashValues))
        hashKeys[[length(hashKeys) + 1]] <- setdiff(unique(df[[col]]),
                                                    unlist(hashKeys))
    else if (length(remainder) > 0){
        hashKeys <- c(hashKeys, remainder)
        hashValues <- c(hashValues, remainder)
    }
    hash <- hashMap(hashKeys, hashValues)
    df[[newCol]] <- hash[as.character(df[[col]])]
    return(df)
}


#' Add a categorical column to a Seurat metadata or SingleCellExperiment
#' coldata
#'
#' @inheritParams metadataDF
#' @inheritParams addCategory
#' @param newCol2 A second column to be added based on the same hash keys.
#' Default is \code{NULL} (no second column will be added).
#' @param hashValues2 A vector of hash values corresponding to the second
#' column. Default is \code{NULL} (no second column will be added).
#'
#' @return A Seurat or SingleCellExpression object with one or two new
#' categorical column(s) in the metadata/coldata.
#'
#' @examples
#' df <- data.frame(fruit = c('apple', 'banana', 'cherry', 'grape'))
#' df <- addCategory(df,
#'                 'fruit',
#'                 'color',
#'                 list(c('apple', 'cherry'),
#'                 'banana',
#'                 'grape'),
#'                 c('red', 'yellow', 'purple'))
#'
#' @export
#'
addMetadataCategory <- function(scObj,
                                col,
                                newCol,
                                hashKeys,
                                hashValues,
                                newCol2 = NULL,
                                hashValues2 = NULL){
    metadataDF(scObj) <- addCategory(metadataDF(scObj), col, newCol,
                                 hashKeys, hashValues)
    if (!is.null(newCol2))
    {
        hashValues2 <- unlist(lapply(hashValues2, function(x)
            paste0(x, collapse = '/')))
        metadataDF(scObj) <- addCategory(metadataDF(scObj), col, newCol2,
                                     hashKeys, hashValues2)
    }
    return(scObj)
}
