#' @importFrom sgof BH BY
#'
NULL

#' Perform the Bonferroni multiple testing correction on a vector of p-values
#'
#' This function performs the Bonferroni multiple testing correction on a
#' vector of p-values.
#'
#' @param pvals A vector of p-values.
#' @param nTests Number of tests. Defaults to the number of p-values.
#' @param replaceHigh Whether to replace adjusted p-values greater than 1
#' with 1.
#'
#' @return A vector of Bonferroni-corrected p-values.
#'
#' @keywords internal
#'
bfCorrectV <- function(pvals, nTests = length(pvals), replaceHigh = TRUE){
    pvalsAdj <- pvals * nTests
    if (replaceHigh)
        pvalsAdj[pvalsAdj > 1] <- 1
    return(pvalsAdj)
}

#' Perform the Benjamini-Hochberg multiple testing correction on a vector of
#' p-values
#'
#' This function performs the Benjamini-Hochberg multiple testing correction on
#' a vector of p-values.
#'
#' @inheritParams bfCorrectV
#' @param keepInputOrder Whether to keep the input order of p-values as opposed
#' to returning them as sorted.
#'
#' @return A vector of Benjamini-Hochberg-corrected p-values.
#'
#' @keywords internal
#'
bhCorrectV <- function(pvals, keepInputOrder = TRUE){
    pvalObj <- BH(pvals)
    pvalsAdj <- pvalObj$Adjusted.pvalues
    if (keepInputOrder)
        return(unname(setNames(pvalsAdj,
                               pvalObj$data)[as.character(pvals)]))
    return(pvalsAdj)
}

#' Perform the Benjamini-Yekutieli multiple testing correction on a vector of
#' p-values
#'
#' This function performs the Benjamini-Yekutieli multiple testing correction
#' on a vector of p-values.
#'
#' @inheritParams bfCorrectV
#' @param keepInputOrder Whether to keep the input order of p-values as opposed
#' to returning them as sorted.
#'
#' @return A vector of Benjamini-Yekutieli-corrected p-values.
#'
#' @keywords internal
#'
byCorrectV <- function(pvals, keepInputOrder = TRUE){
    pvalObj <- BY(pvals)
    pvalsAdj <- pvalObj$Adjusted.pvalues
    if (keepInputOrder)
        return(unname(setNames(pvalsAdj,
                               pvalObj$data)[as.character(pvals)]))
    return(pvalsAdj)
}

#' Perform multiple testing correction with Benjamini-Yekutieli on a vector of
#' p-values
#'
#' This function perform multiple testing correction with Benjamini-Yekutieli
#' on a vector of p-values.
#'
#' @inheritParams bfCorrectV
#' @param mtMethod Multiple testing correction method. Choices are
#' Bonferroni ('bf'), Benjamini-Hochberg('bh'), and Benjamini-Yekutieli ('by').
#' @param mtStat A statistics to be optionally computed. Choices are 'identity'
#' (no statistics will be computed and the adjusted p-values will be returned
#' as such), 'median', 'mean', 'max' and 'min'.
#' @param ... Additional parameters passed to the multiple testing function.
#'
#' @return A vector of p-values corrected for multiple testing.
#'
#' @examples
#' pvals <- c(0.032, 0.001, 0.0045, 0.051, 0.048)
#' mtCorrectV(pvals)
#'
#' @export
#'
mtCorrectV <- function(pvals,
                       mtMethod = c('bf', 'bh', 'by'),
                       mtStat = c('identity', 'median', 'mean', 'max', 'min'),
                       ...){

    mtMethod <- match.arg(mtMethod, c('bf', 'bh', 'by'))
    mtStat <- match.arg(mtStat, c('identity', 'median', 'mean', 'max', 'min'))

    methodFun <- eval(as.name(paste0(mtMethod, 'CorrectV')))
    statFun <- eval(as.name(mtStat))

    if(mtMethod %in% c('bh', 'by') & mtStat == 'identity')
        return(statFun(methodFun(pvals, FALSE, ...))) else
            return(statFun(methodFun(pvals, ...)))
}

#' Perform multiple testing correction on a data frame column
#'
#' This function orders a data frame based on a column of p-values, performs
#' multiple testing on the column, and filters the data-frame based on it.
#'
#' @param df A data frame with a p-values columnn.
#' @inheritParams mtCorrectV
#' @param pvalThr p-value threshold.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
#' @param ... Additional arguments passed to the multiple testing correction
#' method.
#'
#' @examples
#' df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
#' pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
#' mtCorrectDF(df)
#'
#' @export
#'
mtCorrectDF <- function(df,
                        mtMethod = c('bf', 'bh', 'by'),
                        pvalThr = 0.05,
                        colStr = 'pval',
                        newColStr = 'pvalAdj',
                        ...){
    mtMethod <- match.arg(mtMethod, c('bf', 'bh', 'by'))
    fun <- eval(as.name(paste0(mtMethod, 'CorrectV')))
    df <- df[order(df[[colStr]]), ]
    df[[newColStr]] <- fun(df[[colStr]], ...)
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}
