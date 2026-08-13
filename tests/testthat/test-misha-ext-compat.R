test_that("intervs_to_mat keeps all value columns and iceqream-style rownames", {
    intervs <- data.frame(
        chrom = c("chr1", "chr1", "chr2"),
        start = c(100, 500, 100),
        end = c(400, 800, 400),
        ct1 = c(1, 2, 3),
        ct2 = c(4, 5, 6)
    )

    mat <- intervs_to_mat(intervs)

    expect_true(is.matrix(mat))
    expect_null(attr(mat, "class"))
    expect_equal(colnames(mat), c("ct1", "ct2"))
    expect_equal(rownames(mat), c("chr1_100_400", "chr1_500_800", "chr2_100_400"))
    expect_equal(unname(mat[, "ct1"]), c(1, 2, 3))
    # single-column selection must drop, like a plain matrix
    expect_null(dim(mat[, "ct1"]))
})

test_that("intervs_to_mat drops intervalID from the values", {
    intervs <- data.frame(
        chrom = "chr1", start = 100, end = 400, intervalID = 0, ct1 = 7
    )
    expect_equal(colnames(intervs_to_mat(intervs)), "ct1")
})

test_that("mat_to_intervs round-trips intervs_to_mat", {
    intervs <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 500),
        end = c(400, 800),
        ct1 = c(1, 2)
    )
    expect_equal(mat_to_intervs(intervs_to_mat(intervs)), intervs)
})
