#' @importFrom qs qread
#'
NULL

#' Create and save test data
#'
#' This function creates and saves test data.
#'
#' @return Nothing. This function is included for documentation purposes.
#'
#' @noRd
#'
createTestData <- function(){
    if (requireNamespace(c('qs', 'scater', 'scuttle', 'Seurat', 'withr'),
                         quietly=TRUE)){
        sceObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=500))
        sceObj <- scuttle::logNormCounts(sceObj)
        sceObj <- scater::runPCA(sceObj)
        sceObj <- withr::with_seed(1, scater::runUMAP(sceObj))
        qs::qsave(sceObj, 'inst/testdata/sceObj.qs')
        seuratObj <- suppressWarnings(Seurat::as.Seurat(sceObj, data=NULL))
        qs::qsave(seuratObj, 'inst/testdata/seuratObj.qs')
    }
}


