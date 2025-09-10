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



