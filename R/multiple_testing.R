#' Perform multiple testing correction and filtering with Bonferroni
#'
#' This function performs the Bonferroni correction for multiple
#' testing in a dataframe column of p-values and filters the data-frame
#' based on p-values.
#'
#' @param df A dataframe with a column of p-values.
#' @param pvalThr p-value threshold.
#' @param colStr Name of the column of p-values.
#'
#' @return The data frame with Benjamini-Yekutieli-corrected p-values.
#'
#' @examples
#' df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
#' pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
#' bfCorrectDF(df)
#'
#' @export
#'
bfCorrectDF <- function(df, nTests, pvalThr = 0.05, colStr = 'pval'){
    df$pvalAdj <- df[[colStr]] * nTests
    df <- subset(df, pvalAdj < pvalThr)
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
#'
byCorrectDF <- function(df, pvalThr = 0.05, colStr = 'pval'){
    df <- df[order(df[[colStr]]), ]
    df$pvalAdj <- BY(df[[colStr]], pvalThr)$Adjusted.pvalues
    df <- subset(df, pvalAdj < pvalThr)
    return(df)
}
