# A self-contained script to generate the data in extdata
#'
createExtData <- function(){
    if (requireNamespace(c('qs', 'scater', 'scuttle', 'Seurat', 'withr'),
                         quietly=TRUE)){
        sceObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=500))
        sceObj <- scuttle::logNormCounts(sceObj)
        sceObj <- scater::runPCA(sceObj)
        sceObj <- withr::with_seed(1, scater::runUMAP(sceObj))
        scCol(sceObj, 'Cluster') <- withr::with_seed(1,
                                                     sample(paste0('Cluster',
                                                                   seq(5)),
                                                            dim(sceObj)[2],
                                                            replace=TRUE))
        scCol(sceObj, 'Donor') <- rep('Donor1', dim(sceObj)[2])
        for (i in seq(5)){
            scCol(sceObj, 'Donor')[withr::with_seed(1,
                                                    sample(which(scCol(sceObj, 'Cluster') == paste0('Cluster', i))
                                                           , 15))]<- paste0('Donor', i)
            scCol(sceObj, 'Donor')[withr::with_seed(1,
                                                    sample(which(scCol(sceObj,
                                                                       'Cluster') ==
                                                                     paste0('Cluster', i))
                                                           , 15))]<- paste0('Donor', i + 1)
        }
        sceObj <- computeSilhouette(sceObj, 'Cluster')
        qs::qsave(sceObj, 'inst/extdata/sceObj.qs')
        seuratObj <- suppressWarnings(Seurat::as.Seurat(sceObj, data=NULL))
        qs::qsave(seuratObj, 'inst/extdata/seuratObj.qs')
    }
}
