#' Create a map from keys to values
#'
#' This function creates a map from keys to values.
#'
#' @param keys A list of keys. If vectors are part of the keys,
#' each of their elements will be assigned the corresponding value.
#' @param values A vector of values.
#' Must have the same length as \code{keys}.
#'
#' @return A named vector.
#'
#' @examples
#' keyvalMap(list(2, c(3, 4, 5), 6, 8), c('a', 'b', 'c', 'd'))
#'
#' @export
#'
keyvalMap <- function(keys, values){
    if (length(keys) != length(values))
        stop('`keys` and `values` must have the same length.')
    if (max(table(unlist(keys))) > 1)
        stop('All values in `keys` must be unique. ',
             names(which.max(table(unlist(keys)))), ' is repeated.')
    mapNames <- unlist(keys)
    map <- unlist(mapply(function(x, y) rep(y, length(x)), keys,
                          values))
    names(map) <- mapNames
    return(map)
}

#' Add a categorical column to a data frame based on another column
#'
#' This function adds a categorical column to a data frame based on another
#' column.
#'
#' @param df A data frame.
#' @param col Column whose values will be used for creating the new column.
#' @param newCol Column to be added.
#' @inheritParams keyvalMap
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
#' @export
#'
addCategory <- function(df, col, newCol, keys, values){
    remainder <- setdiff(unique(df[[col]]), unlist(keys))
    if (length(keys) + 1 == length(values))
        keys[[length(keys) + 1]] <- setdiff(unique(df[[col]]),
                                                    unlist(keys))
    else if (length(remainder) > 0){
        keys <- c(keys, remainder)
        values <- c(values, remainder)
    }
    map <- keyvalMap(keys, values)
    df[[newCol]] <- map[as.character(df[[col]])]
    return(df)
}

#' Add a categorical column to a Seurat metadata or SingleCellExperiment
#' coldata
#'
#' @inheritParams geneCenters
#' @inheritParams addCategory
#' @param newCol2 A second column to be added based on the same keys.
#' Default is \code{NULL} (no second column will be added).
#' @param values2 A vector of values corresponding to the second
#' column. Default is \code{NULL} (no second column will be added).
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object with one or
#' two new categorical column(s) in the metadata/coldata.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='hammers')
#' sceObj <- qs2::qs_read(scePath)
#' sceObj <- addMetadataCategory(sceObj, 'Cell_Cycle', 'Type',
#' list(c('G0', 'G1'), 'G2M', 'S'), c(2, 3, 1))
#'
#' @export
#'
addMetadataCategory <- function(scObj,
                                col,
                                newCol,
                                keys,
                                values,
                                newCol2 = NULL,
                                values2 = NULL){
    df <- metadataDF(scObj)
    df <- addCategory(df, col, newCol, keys, values)
    scCol(scObj, newCol) <- df[[newCol]]
    if (!is.null(newCol2))
    {
        values2 <- unlist(lapply(values2, function(x)
            paste0(x, collapse = '/')))
        df <- addCategory(df, col, newCol2, keys, values2)
        scCol(scObj, newCol2) <- df[[newCol2]]
    }
    return(scObj)
}
