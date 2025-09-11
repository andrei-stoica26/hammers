#' Perform multiple testing correction and filtering with Bonferroni
#'
#' This function performs the Bonferroni correction for multiple
#' testing in a dataframe column of p-values and filters the data-frame
#' based on p-values.
#'
#' @param df A dataframe with a column of p-values.
#' @param nTests Number of tests.
#' @param pvalThr p-value threshold.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
#'
#' @return The data frame with Benjamini-Yekutieli-corrected p-values.
#'
#' @examples
#' df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
#' pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
#' bfCorrectDF(df, 5)
#'
#' @export
#'
bfCorrectDF <- function(df,
                        nTests,
                        pvalThr = 0.05,
                        colStr = 'pval',
                        newColStr = 'pvalAdj'){
    df[[newColStr]] <- df[[colStr]] * nTests
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}

#' Perform multiple testing correction and filtering with Benjamini-Yekutieli
#'
#' This function performs the Benjamini-Yekutieli correction for multiple
#' testing in a dataframe column of p-values and filters the data-frame
#' based on p-values.
#'
#' @inheritParams bfCorrectDF
#'
#' @return The data frame with Benjamini-Yekutieli-corrected p-values.
#'
#' @examples
#' df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
#' pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
#' byCorrectDF(df)
#'
#' @export
#'
byCorrectDF <- function(df,
                        pvalThr = 0.05,
                        colStr = 'pval',
                        newColStr = 'pvalAdj'){
    df <- df[order(df[[colStr]]), ]
    df[[newColStr]] <- BY(df[[colStr]], pvalThr)$Adjusted.pvalues
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}
