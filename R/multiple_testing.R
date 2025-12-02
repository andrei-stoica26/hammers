#' @importFrom stats p.adjust
#'
NULL

#' Perform multiple testing correction on a data frame column
#'
#' This function orders a data frame based on a column of p-values, performs
#' multiple testing on the column, and filters the data-frame based on it.
#'
#' @param df A data frame with a p-values columnn.
#' @param mtMethod Multiple testing correction method. Choices are 'holm',
#' 'hochberg', hommel', 'bonferroni', 'BH', 'BY', 'fdr' and 'none'. Default is
#' 'holm'.
#' @param doOrder Whether to increasingly order the data frame based on the
#' adjusted p-values.
#' @param pvalThr p-value threshold.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
#' @param ... Additional arguments passed to the multiple testing correction
#' method.
#'
#' @return A data frame with the p-value column corrected for multiple testing.
#'
#' @examples
#' df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
#' pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
#' mtCorrectDF(df)
#'
#' @export
#'
mtCorrectDF <- function(df,
                        mtMethod = c('holm', 'hochberg', 'hommel',
                                     'bonferroni', 'BH', 'BY',
                                     'fdr', 'none'),
                        doOrder = TRUE,
                        pvalThr = 0.05,
                        colStr = 'pval',
                        newColStr = 'pvalAdj',
                        ...){
    mtMethod <- match.arg(mtMethod, c('holm', 'hochberg', 'hommel',
                                      'bonferroni', 'BH', 'BY',
                                      'fdr', 'none'))
    df[[newColStr]] <- p.adjust(df[[colStr]], mtMethod, ...)
    if (doOrder)
        df <- df[order(df[[newColStr]]), ]
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}
