# Tests for preprocess-related functions (no misha genome database required)
# Functions tested: psum, normalize_feature, load_peaks (data frame / file input only)

# --------------------------------------------------------------------------
# psum tests
# --------------------------------------------------------------------------

test_that("psum: two equal-length vectors give element-wise sum", {
    result <- iceqream:::psum(c(1, 2, 3), c(10, 20, 30))
    expect_equal(result, c(11, 22, 33))
})

test_that("psum: three vectors give element-wise sum", {
    result <- iceqream:::psum(c(1, 2, 3), c(10, 20, 30), c(100, 200, 300))
    expect_equal(result, c(111, 222, 333))
})

test_that("psum: single vector returns the same values (as doubles)", {
    result <- iceqream:::psum(c(4, 5, 6))
    expect_equal(result, c(4, 5, 6))
})

test_that("psum: na.rm = FALSE keeps NAs in result", {
    result <- iceqream:::psum(c(1, NA, 3), c(10, 20, 30), na.rm = FALSE)
    expect_equal(result[1], 11)
    expect_true(is.na(result[2]))
    expect_equal(result[3], 33)
})

test_that("psum: na.rm = TRUE sums over non-NA values", {
    result <- iceqream:::psum(c(1, NA, 3), c(10, 20, 30), na.rm = TRUE)
    expect_equal(result, c(11, 20, 33))
})

test_that("psum: all-NA row returns NA even with na.rm = TRUE", {
    result <- iceqream:::psum(c(NA_real_, 1), c(NA_real_, 2), na.rm = TRUE)
    expect_true(is.na(result[1]))
    expect_equal(result[2], 3)
})

test_that("psum: no arguments returns numeric(0)", {
    result <- iceqream:::psum()
    expect_equal(result, numeric(0))
})

test_that("psum: recycling shorter vectors to match longest", {
    # c(1) is recycled to c(1, 1, 1)
    result <- iceqream:::psum(c(1), c(10, 20, 30))
    expect_equal(result, c(11, 21, 31))
})

test_that("psum: zero-length vectors return numeric(0)", {
    result <- iceqream:::psum(numeric(0), numeric(0))
    expect_equal(result, numeric(0))
})

test_that("psum: integer input is coerced to double", {
    result <- iceqream:::psum(1L, 2L)
    expect_type(result, "double")
    expect_equal(result, 3)
})

test_that("psum: works with negative values", {
    result <- iceqream:::psum(c(-5, 10), c(3, -7))
    expect_equal(result, c(-2, 3))
})

# --------------------------------------------------------------------------
# normalize_feature tests
# --------------------------------------------------------------------------

test_that("normalize_feature: result is in [0, 10] range", {
    set.seed(42)
    x <- rnorm(200, mean = 50, sd = 10)
    result <- iceqream:::normalize_feature(x, quant = 0.05)
    expect_true(all(result >= 0))
    expect_true(all(result <= 10))
})

test_that("normalize_feature: values beyond quantiles are capped", {
    # Create a vector with clear spread and known outliers
    set.seed(123)
    x <- c(seq(0, 100, length.out = 100), -500, 500)
    result <- iceqream:::normalize_feature(x, quant = 0.05)
    # After capping at quantile boundaries and normalizing, the result
    # should span [0, 10] and the outlier values should be capped to the
    # same values as the quantile boundaries
    expect_equal(min(result), 0)
    expect_equal(max(result), 10)
    # The two outlier values at positions 101 and 102 should equal
    # the min and max respectively (they were capped)
    expect_equal(result[101], 0)
    expect_equal(result[102], 10)
})

test_that("normalize_feature: all same values gives constant output", {
    x <- rep(7, 50)
    result <- iceqream:::normalize_feature(x, quant = 0.05)
    # norm01 of a constant vector returns all zeros, times 10 => all zeros
    expect_true(all(result == 0))
})

test_that("normalize_feature: endpoints are exactly 0 and 10 for spread data", {
    # With enough spread, quantile capping won't collapse everything
    x <- seq(1, 100, length.out = 200)
    result <- iceqream:::normalize_feature(x, quant = 0.05)
    expect_equal(min(result), 0)
    expect_equal(max(result), 10)
})

test_that("normalize_feature: quant = 0 does no capping (just norm01 * 10)", {
    x <- c(0, 5, 10)
    result <- iceqream:::normalize_feature(x, quant = 0)
    expected <- iceqream:::norm01(x) * 10
    expect_equal(result, expected)
})

# --------------------------------------------------------------------------
# load_peaks tests (data frame path only - no misha)
# --------------------------------------------------------------------------

test_that("load_peaks: valid data frame with chrom, start, end is accepted", {
    peaks_df <- data.frame(
        chrom = c("chr1", "chr1", "chr2"),
        start = c(100, 500, 200),
        end   = c(200, 600, 300),
        stringsAsFactors = FALSE
    )
    # load_peaks calls gintervals.canonic which needs misha, so we mock it
    # We need to stub out the gintervals.canonic call
    # Since load_peaks is not exported and calls misha internals for overlap
    # removal, we use mockery/local_mocked_bindings to skip those parts
    local_mocked_bindings(
        gintervals.canonic = function(x) x,
        .package = "misha"
    )
    result <- iceqream:::load_peaks(peaks_df)
    expect_s3_class(result, "data.frame")
    expect_true(all(c("chrom", "start", "end") %in% colnames(result)))
    expect_equal(nrow(result), 3)
})

test_that("load_peaks: missing required columns throws error", {
    bad_df <- data.frame(
        chromosome = c("chr1"),
        begin = c(100),
        finish = c(200)
    )
    expect_error(
        iceqream:::load_peaks(bad_df),
        "chrom.*start.*end"
    )
})

test_that("load_peaks: non-data-frame, non-string input throws error", {
    expect_error(
        iceqream:::load_peaks(42),
        "data frame or a path"
    )
    expect_error(
        iceqream:::load_peaks(list(a = 1)),
        "data frame or a path"
    )
})

test_that("load_peaks: reads TSV file with header correctly", {
    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(
        c("chrom\tstart\tend", "chr1\t100\t200", "chr1\t500\t600"),
        tmp
    )
    local_mocked_bindings(
        gintervals.canonic = function(x) x,
        .package = "misha"
    )
    result <- iceqream:::load_peaks(tmp)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_true(all(c("chrom", "start", "end") %in% colnames(result)))
})

test_that("load_peaks: reads BED file (no header) correctly", {
    tmp <- tempfile(fileext = ".bed")
    on.exit(unlink(tmp), add = TRUE)
    writeLines(
        c("chr1\t100\t200", "chr2\t300\t400"),
        tmp
    )
    local_mocked_bindings(
        gintervals.canonic = function(x) x,
        .package = "misha"
    )
    result <- iceqream:::load_peaks(tmp)
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 2)
    expect_equal(result$chrom[1], "chr1")
    expect_equal(result$start[1], 100)
    expect_equal(result$end[1], 200)
})

test_that("load_peaks: data frame with extra columns preserves them", {
    peaks_df <- data.frame(
        chrom = c("chr1", "chr2"),
        start = c(100, 200),
        end   = c(200, 300),
        score = c(5.5, 3.2),
        stringsAsFactors = FALSE
    )
    local_mocked_bindings(
        gintervals.canonic = function(x) x,
        .package = "misha"
    )
    result <- iceqream:::load_peaks(peaks_df)
    expect_true("score" %in% colnames(result))
    expect_equal(result$score, c(5.5, 3.2))
})
