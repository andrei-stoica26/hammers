test_that("centerOfMass works", {
    dimMat <- matrix(data=c(2, 3, 1, 3, 6, 8), nrow=3, ncol=2)
    weights <- c(0.8, 6, 16)
    res <- centerOfMass(dimMat, weights)
    expect_equal(res, c(1.561404, 7.298246), tolerance=0.001)
})

test_that("geneCenters and colCenters work", {
    df <- geneCenters(scObj, c('Gene_0980', 'Gene_0981', 'Gene_0982'))
    expect_equal(df$UMAP1, c(0.022384739, -0.003861012, -0.003507070),
                 tolerance=0.001)
    expect_equal(df$UMAP2, c(0.0111317786, -0.0004384431, -0.0139794494),
                 tolerance=0.001)

    df <- colCenters(scObj, c('sizeFactor'))
    expect_equal(df$UMAP1, c(0.0008518303), tolerance=0.001)
    expect_equal(df$UMAP1, c(0.0006236511), tolerance=0.001)
})

test_that("compatibility functions and checks work", {
    expect_null(checkGenes(scObj, c('Gene_0980', 'Gene_0981', 'Gene_0982')))
    expect_error(checkGenes(scObj, c('Gene_0980', 'Gene_0981', 'Gene_0982',
                                     'DSFDGDG')))

    expect_equal(metadataNames(scObj), c('Mutation_Status',
                                         'Cell_Cycle',
                                         'Treatment',
                                         'sizeFactor'))

    expect_error(metadataNames(c(1, 2, 3)))
    expect_equal(length(scCol(scObj, 'Cell_Cycle')), 200)

    expect_equal(scColCounts(scObj, 'Cell_Cycle')['G2M'],
                 setNames(47, 'G2M'))
    expect_equal(scColPairCounts(scObj, 'Mutation_Status', 'Cell_Cycle')[1, 3],
                 28)

    expect_equal(length(scGeneExp(scObj, 'Gene_0980')), 200)
    expect_equal(dim(scExpMat(scObj)), c(20000, 200))

    v <- scPCAMat(scObj)
    w <- scPCAMat(seuratObj)
    colnames(v) <- paste0('PC_', seq(50))
    expect_equal(v, w)

    v <- scUMAPMat(scObj)
    w <- scUMAPMat(seuratObj)
    colnames(v) <- paste0('UMAP_', seq(2))
    expect_equal(v, w)
})

test_that("genePresence works", {
    scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=200))
    df <- genePresence(scObj)
    expect_identical(mean(df[, 2]), 144.775, tolerance=0.001)
})

test_that("repAnalysis and pvalRiverPlot work", {
    scObj <- withr::with_seed(1, scuttle::mockSCE(ngenes=20000))
    scCol(scObj, 'Cluster') <- withr::with_seed(1,
                                                sample(paste0('Cluster',
                                                              seq(5)),
                                                       dim(scObj)[2],
                                                       replace=TRUE))
    scCol(scObj, 'Donor') <- rep('Donor1', dim(scObj)[2])
    for (i in seq(5)){
        scCol(scObj, 'Donor')[withr::with_seed(1,
                                               sample(which(scCol(scObj,
                                                                  'Cluster') ==
                                               paste0('Cluster', i))
                                     , 15))]<- paste0('Donor', i)
        scCol(scObj, 'Donor')[withr::with_seed(1,
                                               sample(which(scCol(scObj,
                                                                  'Cluster') ==
                                               paste0('Cluster', i))
                                     , 15))]<- paste0('Donor', i + 1)
    }
    df <- repAnalysis(scObj, 'Cluster', 'Donor')
    expect_equal(ncol(df), 9)
    expect_equal(mean(df$pvalAdj), 7.30162e-11, tolerance=0.001)
    p <- pvalRiverPlot(df)
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("addMetadataCategory works", {
    scObj <- addMetadataCategory(scObj, 'Cell_Cycle', 'Type',
                                 list(c('G0', 'G1'), 'G2M', 'S'), c(2, 3, 1))
    expect_equal(unique(metadataDF(scObj)[['Type']]),
                 c(2, 3, 1))
})

test_that("multiple testing functions work", {
    df <- data.frame(elem = c('A', 'B', 'C', 'D', 'E'),
                     pval = c(0.032, 0.001, 0.0045, 0.051, 0.048))
    res <- mtCorrectDF(df, 'bf', nTests=5)
    expect_equal(res$pvalAdj, c(0.0050, 0.0225), tolerance=0.001)
    res <- mtCorrectDF(df, 'bh')
    expect_equal(res$pvalAdj, c(0.00500, 0.01125), tolerance=0.001)
    res <- mtCorrectDF(df, 'by')
    expect_equal(res$pvalAdj, c(0.01141667, 0.02568750), tolerance=0.001)
})

test_that("silhouette functions and scCol work", {
    scObj <- computeSilhouette(scObj, 'Cell_Cycle')
    df <- normalizeSilhouette(scObj, 'Cell_Cycle')
    expect_equal(sum(df), 112.3134, tolerance=0.001)
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
    genes <- c('Gene_0826', 'Gene_0610', 'Gene_0380', 'Gene_0602',
               'Gene_0613', 'Gene_0201', 'Gene_0295')
    newGenes <- shuffleGenes(scObj, genes, 0.3, 0.9)
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
    p <- distributionPlot(scObj, col1='Mutation_Status', col2='Cell_Cycle')
    expect_equal(length(intersect(is(p), c('gg', 'ggplot2::ggplot'))), 1)
})

test_that("dimPlot functions work", {
    pointsDF <- data.frame(x = c(2, 3),
                           y = c(1, 0.5),
                           row.names = c('P1', 'P2'))
    pointsDimPlot(seuratObj, pointsDF=pointsDF)
    expect_equal(is(pointsDimPlot(seuratObj, pointsDF=pointsDF)), 'patchwork')
    expect_equal(is(genesDimPlot(seuratObj, c('Spike-0001', 'Spike-0074',
                                              'Spike-0067'))),
                    'patchwork')
    expect_equal(is(colsDimPlot(seuratObj,
                                c('nCount_originalexp',
                                  'nFeature_originalexp'))),
                    'patchwork')
})
