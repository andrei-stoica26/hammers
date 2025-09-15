test_that("centerOfMass works", {
    dimMat <- matrix(data=c(2, 3, 1, 3, 6, 8), nrow=3, ncol=2)
    weights <- c(0.8, 6, 16)
    res <- centerOfMass(dimMat, weights)
    expect_equal(res, c(1.561404, 7.298246), tolerance=0.001)
})

test_that("geneCenters and colCenters work", {
    df <- geneCenters(seuratObj, c('AURKA', 'MKI67', 'TOP2A'))
    expect_equal(df$umap_1, c(1.563351, 6.073028, 5.656954),
                 tolerance=0.001)
    expect_equal(df$umap_2, c(-2.3658094, -1.1434781, -0.5404839),
                 tolerance=0.001)

    df <- colCenters(seuratObj, c('nCount_originalexp',
                              'nFeature_originalexp',
                              'sizeFactor'))
    expect_equal(df$umap_1, c(0.1617853, -0.1486735, 0.1617853),
                 tolerance=0.001)
    expect_equal(df$umap_2, c(-1.197215, -0.509185, -1.197215),
                 tolerance=0.001)
})

test_that("compatibility functions and checks work", {
    expect_null(checkGenes(sceObj, c('AURKA', 'TOP2A', 'MKI67')))
    expect_error(checkGenes(sceObj, c('AURKA', 'TOP2A', 'MKI67', 'DSFDGDG')))

    expect_equal(colnames(metadataDF(sceObj)), c('donor', 'label', 'sizeFactor',
                                                 'silhouette'))
    expect_equal(metadataNames(sceObj), c('donor', 'label', 'sizeFactor',
                                          'silhouette'))

    expect_equal(colnames(metadataDF(seuratObj)), c('orig.ident',
                                                    'nCount_originalexp',
                                                    'nFeature_originalexp',
                                                    'donor',
                                                    'label',
                                                    'sizeFactor'))
    expect_equal(metadataNames(seuratObj), c('orig.ident',
                                             'nCount_originalexp',
                                             'nFeature_originalexp',
                                             'donor',
                                             'label',
                                             'sizeFactor'))

    expect_error(metadataDF(c(1, 2, 3)))
    expect_error(metadataNames(c(1, 2, 3)))

    expect_equal(scCol(sceObj, 'label'), scCol(seuratObj, 'label'))
    expect_equal(length(scCol(sceObj, 'label')), 8569)

    expect_equal(scColCounts(sceObj, 'label')['acinar'],
                 setNames(958, 'acinar'))
    expect_equal(scColPairCounts(sceObj, 'label', 'donor')[1, 3], 110)

    expect_equal(scGeneExp(sceObj, 'AURKA'), scGeneExp(seuratObj, 'AURKA'))
    expect_equal(length(scGeneExp(sceObj, 'AURKA')), 8569)

    res <- scExpMat(sceObj)
    expect_equal(res, scExpMat(seuratObj))
    expect_equal(dim(res), c(20125, 8569))

    convSeuratObj <- suppressWarnings(as.Seurat(sceObj))

    v <- scPCAMat(sceObj)
    w <- scPCAMat(convSeuratObj)
    colnames(v) <- paste0('PC_', seq(50))
    expect_equal(v, w)

    v <- scUMAPMat(sceObj)
    w <- scUMAPMat(convSeuratObj)
    colnames(v) <- paste0('UMAP_', seq(2))
    expect_equal(v, w)
})

test_that("repAnalysis and pvalRiverPlot work", {
    df <- repAnalysis(sceObj, 'donor', 'label')
    expect_equal(ncol(df), 9)
    expect_equal(mean(df$pvalAdj), 0.0005422197, tolerance=0.001)
    expect_equal(is(pvalRiverPlot(df)), 'gg')
})

test_that("gene enrichment functions work", {
    m <- genesER(c('AURKA', 'TOP2A', 'CENPF', 'PTTG2', 'MKI67', 'BIRC5', 'RRM2'),
                 'human')
    expect_true('chromosome segregation' %in% m@result$Description)
    expect_equal(termGenes(m, 'chromosome segregation', 'meiosis I'),
                 c('BIRC5', 'CENPF', 'MKI67'))
})

test_that("addMetadataCategory works", {
    scObjExt <- addMetadataCategory(scObj,
                                    'label',
                                    'initial',
                                    list(c('acinar',
                                           'activated_stellate',
                                           'alpha'),
                                         'beta',
                                         c('delta', 'ductal'),
                                         c('epsilon', 'endothelial'),
                                         'gamma',
                                         'quiescent_stellate',
                                         c('macrophage', 'mast'),
                                         'schwann',
                                         't_cell'),
                                    c('a', 'b', 'd', 'e', 'g', 'q', 'm',
                                      's', 't'))
    expect_equal(unique(metadataDF(scObjExt)[['initial']]),
                 c('a', 'b', 'd', 'e', 'g', 'q', 'm', 's', 't'))
})

test_that("multiple testing functions work", {
    df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
                     pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
    res <- bfCorrectDF(df, 5)
    expect_equal(res$pvalAdj, c(0.0050, 0.0225), tolerance=0.001)
    res <- bhCorrectDF(df)
    expect_equal(res$pvalAdj, c(0.00500, 0.01125), tolerance=0.001)
    res <- byCorrectDF(df)
    expect_equal(res$pvalAdj, c(0.01141667, 0.02568750), tolerance=0.001)
})

test_that("silhouette functions work", {
    expect_equal(mean(scCol(sceObj, 'silhouette')), 0.3981082, tolerance=0.001)
})

test_that("joinCharCombs works", {
    res <- joinCharCombs(c('a', 'b', 'c', 'd'), c('eee', 'ff'), c(1, 2, 3))
    expected <- c('a_eee_1', 'a_eee_2', 'a_eee_3', 'a_ff_1', 'a_ff_2', 'a_ff_3',
                  'b_eee_1', 'b_eee_2', 'b_eee_3', 'b_ff_1', 'b_ff_2', 'b_ff_3',
                  'c_eee_1', 'c_eee_2', 'c_eee_3', 'c_ff_1', 'c_ff_2', 'c_ff_3',
                  'd_eee_1', 'd_eee_2', 'd_eee_3', 'd_ff_1', 'd_ff_2', 'd_ff_3')
    expect_equal(res, expected)
})

test_that("nearestNeighbors works", {
    df <- data.frame(v = c(1, 2, 4, 5, 6),
                     w = c(2, 3, 1, 5, 8),
                     x = c(2, 8, 7, 1, 1),
                     y = c(2, 3, 2, 2, 4),
                     z = c(1, 9, 9, 7, 6))
    distMat <- as.matrix(stats::dist(df))
    rownames(distMat) <- c('v', 'w', 'x', 'y', 'z')
    colnames(distMat) <- c('v', 'w', 'x', 'y', 'z')
    res <- nearestNeighbors(distMat)
    expected <- setNames(c('y', 'x', 'w', 'z', 'y'), rownames(distMat))
    expect_equal(res, expected)
})

test_that("proximity works", {
    expect_equal(proximity(2, 3, 6), 0.8333333, tolerance=0.001)
})

test_that("safeMinmax works", {
    expect_equal(safeMinmax(c(2.1, 3.2, 2.8)), c(0, 1, 0.6363636),
                 tolerance=0.001)
    expect_equal(safeMinmax(c(8, 8, 8)), c(0, 0, 0))
})

test_that("safeMessage works", {
    expect_message(safeMessage('message'), 'message')
    expect_null(safeMessage('message'))
    expect_message(safeMessage('message', FALSE), NA)
})

test_that("shuffleGenes works", {
    genes <- c('TOP2A', 'BIRC5', 'MKI67', 'RRM2', 'CENPF', 'PTTG2', 'CLSPN')
    newGenes <- shuffleGenes(sceObj, genes, 0.3, 0.9)
    expect_equal(length(intersect(genes, newGenes)), 5)
    expect_equal(length(newGenes), 50)
})

test_that("tabulateVector works", {
    v <- c(2, 3, 4, 19, 15, 25, 32, 8)
    res <- tabulateVector(v, paste0('r', seq(4)), paste0('c', seq(2)))
    df <- data.frame(c1 = c(2, 3, 4, 19),
                     c2 = c(15, 25, 32, 8),
                     row.names = paste0('r', seq(4)))
    expect_equal(res, df)
})

test_that("timeFun works", {
    expect_output(res <- timeFun(sum, 2, 3, 4))
    expect_equal(res, 9)
})

test_that("distributionPlot works", {
    expect_equal(is(distributionPlot(sceObj,
                                     col1='donor', col2='label')), 'gg')
})

test_that("dimPlot functions work", {
    pointsDF <- data.frame(x = c(2, 3),
                           y = c(1, 6),
                           row.names = c('P1', 'P2'))
    pointsDimPlot(seuratObj, pointsDF=pointsDF)
    expect_equal(is(pointsDimPlot(seuratObj, pointsDF=pointsDF)), 'patchwork')
    expect_equal(is(genesDimPlot(seuratObj, c('AURKA', 'TOP2A', 'MKI67'))),
                    'patchwork')
    expect_equal(is(colsDimPlot(seuratObj,
                                c('nCount_originalexp',
                                  'nFeature_originalexp'))),
                    'patchwork')
})
