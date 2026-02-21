# Tests for validation functions in trajectory-regression.R

# ---- resolve_bin_index ----

test_that("resolve_bin_index returns integer passthrough for integer input", {
    mat <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(mat) <- c("col1", "col2", "col3", "col4")
    result <- iceqream:::resolve_bin_index(3L, mat, "test")
    expect_identical(result, 3L)
})

test_that("resolve_bin_index coerces numeric to integer", {
    mat <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(mat) <- c("col1", "col2", "col3", "col4")
    result <- iceqream:::resolve_bin_index(2.7, mat, "test")
    expect_identical(result, 2L)
    expect_type(result, "integer")
})

test_that("resolve_bin_index resolves character column name to correct index", {
    mat <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(mat) <- c("col1", "col2", "col3", "col4")
    result <- iceqream:::resolve_bin_index("col2", mat, "test")
    expect_equal(result, 2L)
})

test_that("resolve_bin_index resolves first column name", {
    mat <- matrix(1:9, nrow = 3, ncol = 3)
    colnames(mat) <- c("alpha", "beta", "gamma")
    expect_equal(iceqream:::resolve_bin_index("alpha", mat, "test"), 1L)
})

test_that("resolve_bin_index resolves last column name", {
    mat <- matrix(1:9, nrow = 3, ncol = 3)
    colnames(mat) <- c("alpha", "beta", "gamma")
    expect_equal(iceqream:::resolve_bin_index("gamma", mat, "test"), 3L)
})

test_that("resolve_bin_index errors on missing character column name", {
    mat <- matrix(1:12, nrow = 3, ncol = 4)
    colnames(mat) <- c("col1", "col2", "col3", "col4")
    expect_error(
        iceqream:::resolve_bin_index("nonexistent", mat, "bin_start"),
        "not found"
    )
})

test_that("resolve_bin_index errors when matrix has no column names and character is given", {
    mat <- matrix(1:12, nrow = 3, ncol = 4)
    expect_error(
        iceqream:::resolve_bin_index("col1", mat, "bin_start"),
        "not found"
    )
})

test_that("resolve_bin_index passes through 1 as integer", {
    mat <- matrix(1:6, nrow = 2, ncol = 3)
    result <- iceqream:::resolve_bin_index(1, mat, "test")
    expect_identical(result, 1L)
})

# ---- validate_atac_scores ----

test_that("validate_atac_scores returns correct indices for valid integer inputs", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    result <- iceqream:::validate_atac_scores(mat, 1, 3)
    expect_equal(result$bin_start, 1L)
    expect_equal(result$bin_end, 3L)
})

test_that("validate_atac_scores resolves valid character column names", {
    mat <- matrix(rnorm(40), nrow = 10, ncol = 4)
    colnames(mat) <- c("bin_A", "bin_B", "bin_C", "bin_D")
    result <- iceqream:::validate_atac_scores(mat, "bin_A", "bin_D")
    expect_equal(result$bin_start, 1L)
    expect_equal(result$bin_end, 4L)
})

test_that("validate_atac_scores resolves mixed character and integer inputs", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(mat) <- c("first", "middle", "last")
    result <- iceqream:::validate_atac_scores(mat, "first", 3)
    expect_equal(result$bin_start, 1L)
    expect_equal(result$bin_end, 3L)
})

test_that("validate_atac_scores errors when bin_start is 0", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    expect_error(
        iceqream:::validate_atac_scores(mat, 0, 3),
        "bin_start"
    )
})

test_that("validate_atac_scores errors when bin_start exceeds ncol", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    expect_error(
        iceqream:::validate_atac_scores(mat, 5, 3),
        "bin_start"
    )
})

test_that("validate_atac_scores errors when bin_end is 0", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    expect_error(
        iceqream:::validate_atac_scores(mat, 1, 0),
        "bin_end"
    )
})

test_that("validate_atac_scores errors when bin_end exceeds ncol", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    expect_error(
        iceqream:::validate_atac_scores(mat, 1, 10),
        "bin_end"
    )
})

test_that("validate_atac_scores errors when bin_start equals bin_end", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    expect_error(
        iceqream:::validate_atac_scores(mat, 2, 2),
        "different columns"
    )
})

test_that("validate_atac_scores allows bin_start > bin_end (reverse trajectory)", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    result <- iceqream:::validate_atac_scores(mat, 3, 1)
    expect_equal(result$bin_start, 3L)
    expect_equal(result$bin_end, 1L)
})

test_that("validate_atac_scores errors on character bin_start not found", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(mat) <- c("a", "b", "c")
    expect_error(
        iceqream:::validate_atac_scores(mat, "missing_col", "c"),
        "not found"
    )
})

test_that("validate_atac_scores works with two-column matrix", {
    mat <- matrix(rnorm(20), nrow = 10, ncol = 2)
    result <- iceqream:::validate_atac_scores(mat, 1, 2)
    expect_equal(result$bin_start, 1L)
    expect_equal(result$bin_end, 2L)
})

test_that("validate_atac_scores returns a list with bin_start and bin_end names", {
    mat <- matrix(rnorm(30), nrow = 10, ncol = 3)
    result <- iceqream:::validate_atac_scores(mat, 1, 3)
    expect_named(result, c("bin_start", "bin_end"))
})

# ---- validate_peak_intervals ----

test_that("validate_peak_intervals accepts valid data frame with equal-size peaks", {
    peaks <- data.frame(
        chrom = rep("chr1", 5),
        start = seq(1000, 5000, by = 1000),
        end = seq(1500, 5500, by = 1000)
    )
    expect_no_error(iceqream:::validate_peak_intervals(peaks))
})

test_that("validate_peak_intervals errors when chrom column is missing", {
    peaks <- data.frame(
        start = c(100, 200),
        end = c(600, 700)
    )
    expect_error(
        iceqream:::validate_peak_intervals(peaks),
        "peak_intervals"
    )
})

test_that("validate_peak_intervals errors when start column is missing", {
    peaks <- data.frame(
        chrom = c("chr1", "chr2"),
        end = c(600, 700)
    )
    expect_error(
        iceqream:::validate_peak_intervals(peaks),
        "peak_intervals"
    )
})

test_that("validate_peak_intervals errors when end column is missing", {
    peaks <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 200)
    )
    expect_error(
        iceqream:::validate_peak_intervals(peaks),
        "peak_intervals"
    )
})

test_that("validate_peak_intervals errors on unequal peak sizes", {
    peaks <- data.frame(
        chrom = c("chr1", "chr1", "chr2"),
        start = c(100, 200, 300),
        end = c(600, 800, 900)
    )
    expect_error(
        iceqream:::validate_peak_intervals(peaks),
        "same size"
    )
})

test_that("validate_peak_intervals warns on peak sizes smaller than 100bp", {
    peaks <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 200),
        end = c(150, 250) # 50bp peaks
    )
    expect_warning(
        iceqream:::validate_peak_intervals(peaks),
        "smaller than 100bp"
    )
})

test_that("validate_peak_intervals passes with exactly 100bp peaks", {
    peaks <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 200),
        end = c(200, 300) # 100bp peaks
    )
    expect_no_warning(iceqream:::validate_peak_intervals(peaks))
    expect_no_error(iceqream:::validate_peak_intervals(peaks))
})

test_that("validate_peak_intervals passes with large peaks (500bp)", {
    peaks <- data.frame(
        chrom = c("chr1", "chr2", "chr3"),
        start = c(100, 200, 300),
        end = c(600, 700, 800) # 500bp peaks
    )
    expect_no_error(iceqream:::validate_peak_intervals(peaks))
})

test_that("validate_peak_intervals works with single peak", {
    peaks <- data.frame(
        chrom = "chr1",
        start = 1000,
        end = 1500
    )
    expect_no_error(iceqream:::validate_peak_intervals(peaks))
})

test_that("validate_peak_intervals respects custom columns parameter", {
    peaks <- data.frame(
        chrom = "chr1",
        start = 100,
        end = 600,
        extra_col = "x"
    )
    # Default columns work fine
    expect_no_error(iceqream:::validate_peak_intervals(peaks))
    # Custom columns can require additional columns
    expect_error(
        iceqream:::validate_peak_intervals(peaks, columns = c("chrom", "start", "end", "score")),
        "peak_intervals"
    )
})

# ---- validate_additional_features ----

test_that("validate_additional_features passes with NULL additional_features", {
    peaks <- data.frame(chrom = "chr1", start = 100, end = 600)
    expect_no_error(iceqream:::validate_additional_features(NULL, peaks))
})

test_that("validate_additional_features passes with matching row count", {
    peaks <- data.frame(
        chrom = rep("chr1", 5),
        start = seq(100, 500, 100),
        end = seq(600, 1000, 100)
    )
    feats <- data.frame(
        gc_content = runif(5),
        dist_to_tss = runif(5, 1000, 50000)
    )
    expect_no_error(iceqream:::validate_additional_features(feats, peaks))
})

test_that("validate_additional_features errors on mismatched row count", {
    peaks <- data.frame(
        chrom = rep("chr1", 5),
        start = seq(100, 500, 100),
        end = seq(600, 1000, 100)
    )
    feats <- data.frame(
        gc_content = runif(3),
        dist_to_tss = runif(3, 1000, 50000)
    )
    expect_error(
        iceqream:::validate_additional_features(feats, peaks),
        "same number of rows"
    )
})

test_that("validate_additional_features warns when additional_features contains NAs", {
    peaks <- data.frame(
        chrom = rep("chr1", 3),
        start = c(100, 200, 300),
        end = c(600, 700, 800)
    )
    feats <- data.frame(
        gc_content = c(0.5, NA, 0.3),
        dist_to_tss = c(1000, 2000, NA)
    )
    expect_warning(
        iceqream:::validate_additional_features(feats, peaks),
        "NA"
    )
})

test_that("validate_additional_features does not warn when no NAs present", {
    peaks <- data.frame(
        chrom = rep("chr1", 3),
        start = c(100, 200, 300),
        end = c(600, 700, 800)
    )
    feats <- data.frame(
        gc_content = c(0.5, 0.6, 0.3),
        dist_to_tss = c(1000, 2000, 3000)
    )
    expect_no_warning(iceqream:::validate_additional_features(feats, peaks))
})

test_that("validate_additional_features errors when features have more rows than peaks", {
    peaks <- data.frame(
        chrom = rep("chr1", 2),
        start = c(100, 200),
        end = c(600, 700)
    )
    feats <- data.frame(
        gc_content = runif(5)
    )
    expect_error(
        iceqream:::validate_additional_features(feats, peaks),
        "same number of rows"
    )
})

test_that("validate_additional_features errors when features have fewer rows than peaks", {
    peaks <- data.frame(
        chrom = rep("chr1", 5),
        start = seq(100, 500, 100),
        end = seq(600, 1000, 100)
    )
    feats <- data.frame(
        gc_content = runif(2)
    )
    expect_error(
        iceqream:::validate_additional_features(feats, peaks),
        "same number of rows"
    )
})

test_that("validate_additional_features passes with single-column additional features", {
    peaks <- data.frame(
        chrom = rep("chr1", 4),
        start = c(100, 200, 300, 400),
        end = c(600, 700, 800, 900)
    )
    feats <- data.frame(
        gc_content = runif(4)
    )
    expect_no_error(iceqream:::validate_additional_features(feats, peaks))
})
