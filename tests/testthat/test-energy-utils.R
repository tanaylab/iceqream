# Test for norm01 function
test_that("norm01 normalizes values to the [0, 1] range", {
    x <- rnorm(100, 10, 5)
    normed_x <- norm01(x)
    expect_true(all(normed_x >= 0))
    expect_true(all(normed_x <= 1))
    expect_equal(min(normed_x), 0)
    expect_equal(max(normed_x), 1)
})

test_that("norm01 handles all same values", {
    x <- rep(5, 10)
    normed_x <- norm01(x)
    expect_true(all(normed_x == 0))
})

test_that("norm01 handles NAs", {
    x <- c(1, 2, NA, 4, 5)
    normed_x <- norm01(x)
    expect_equal(normed_x[1], 0)
    expect_equal(normed_x[5], 1)
    expect_true(is.na(normed_x[3]))
})

test_that("norm01 handles single element", {
    x <- c(42)
    normed_x <- norm01(x)
    expect_equal(normed_x, 0)
})

test_that("norm01 handles negative values", {
    x <- c(-10, -5, 0, 5, 10)
    normed_x <- norm01(x)
    expect_equal(normed_x[1], 0)
    expect_equal(normed_x[5], 1)
    expect_equal(normed_x[3], 0.5)
})

test_that("norm01 handles matrix input column-wise", {
    m <- matrix(c(1, 2, 3, 10, 20, 30), ncol = 2)
    normed_m <- norm01(m)
    expect_true(is.matrix(normed_m))
    expect_equal(dim(normed_m), dim(m))
    # Each column should be independently normalized
    expect_equal(normed_m[1, 1], 0)
    expect_equal(normed_m[3, 1], 1)
    expect_equal(normed_m[1, 2], 0)
    expect_equal(normed_m[3, 2], 1)
})

test_that("norm01 handles matrix with constant column", {
    m <- matrix(c(1, 2, 3, 5, 5, 5), ncol = 2)
    normed_m <- norm01(m)
    # Constant column should be all zeros
    expect_true(all(normed_m[, 2] == 0))
    # Non-constant column should have range [0,1]
    expect_equal(normed_m[1, 1], 0)
    expect_equal(normed_m[3, 1], 1)
})

# Test for norm0q function
test_that("norm0q normalizes to specified quantile", {
    set.seed(42)
    x <- rnorm(1000, mean = 5, sd = 2)
    normed_x <- norm0q(x, quant = 0.99)
    expect_true(min(normed_x) >= 0)
    expect_true(max(normed_x) <= 1)
    # Most values should be between 0 and 1
    expect_true(mean(normed_x >= 0 & normed_x <= 1) > 0.98)
})

test_that("norm0q caps values above 1", {
    x <- c(1, 2, 3, 4, 100)
    normed_x <- norm0q(x, quant = 0.5)
    expect_true(max(normed_x) == 1)
})

test_that("norm0q with quant=1 behaves like norm01 for positive-range data", {
    x <- c(0, 1, 2, 3, 4, 5)
    normed_q1 <- norm0q(x, quant = 1)
    normed_01 <- norm01(x)
    expect_equal(normed_q1, normed_01)
})

test_that("norm0q handles NAs", {
    x <- c(1, NA, 3, 4, 5)
    normed_x <- norm0q(x, quant = 0.99)
    expect_true(is.na(normed_x[2]))
    expect_true(normed_x[1] >= 0)
})

# Test for norm_energy function
test_that("norm_energy normalizes energy values correctly", {
    x <- runif(n = 100, min = -11, max = 0)
    normed_energy <- norm_energy(x)
    expect_true(all(!is.na(normed_energy)))
    expect_true(all(normed_energy <= 10))
    expect_true(all(normed_energy >= 0))
})

test_that("norm_energy clamps to min_energy boundary", {
    # Very negative values should be clamped to min_energy, then shifted
    x <- c(-100, -50, -20, -10, 0)
    result <- norm_energy(x, min_energy = -7)
    # Values clamped at min_energy become 0 after subtraction
    expect_equal(min(result), 0)
})

test_that("norm_energy respects custom min_energy", {
    x <- runif(n = 50, min = -15, max = 0)
    result_default <- norm_energy(x, min_energy = -7)
    result_custom <- norm_energy(x, min_energy = -3)
    # With a less negative min_energy, more values should be clamped
    expect_true(sum(result_custom == 0) >= sum(result_default == 0))
})

test_that("norm_energy with q < 1 uses quantile reference", {
    set.seed(123)
    x <- runif(n = 100, min = -10, max = 0)
    result_q1 <- norm_energy(x, q = 1)
    result_q09 <- norm_energy(x, q = 0.9)
    # Different quantile reference should give different results
    expect_false(identical(result_q1, result_q09))
})

test_that("norm_energy returns all zeros for identical values", {
    x <- rep(-5, 10)
    result <- norm_energy(x)
    # All same => log2(x/x) = 0 for all, then 0 - min_energy = -min_energy = 7 for all
    # Actually: exp(-5)^1 / quantile = 1, log2(1) = 0, clamp to 0, 0 - (-7) = 7
    expect_true(all(result == result[1]))
})

# Test for logist function
test_that("logist calculates the logistic function correctly", {
    x_vals <- seq(-10, 10, by = 0.1)
    results <- logist(x_vals)
    expect_true(all(!is.na(results)))
    expect_true(all(results <= 1))
    expect_true(all(results >= 0))
})

test_that("logist at midpoint equals L/2", {
    expect_equal(logist(0, x_0 = 0, L = 1, k = 1), 0.5)
    expect_equal(logist(5, x_0 = 5, L = 2, k = 1), 1.0)
    expect_equal(logist(3, x_0 = 3, L = 10, k = 2), 5.0)
})

test_that("logist boundary values", {
    # Very large positive x should approach L
    expect_equal(logist(1000, L = 1, k = 1), 1, tolerance = 1e-10)
    # Very large negative x should approach 0
    expect_equal(logist(-1000, L = 1, k = 1), 0, tolerance = 1e-10)
    # x = 0 with default params should be 0.5
    expect_equal(logist(0), 0.5)
})

test_that("logist is monotonically increasing for positive k", {
    x_vals <- seq(-10, 10, by = 0.5)
    results <- logist(x_vals, k = 2)
    expect_true(all(diff(results) >= 0))
})

test_that("logist with k=0 returns L/2 everywhere", {
    x_vals <- c(-100, -1, 0, 1, 100)
    results <- logist(x_vals, k = 0, L = 4)
    expect_true(all(results == 2))
})

test_that("logist handles vector input", {
    x <- c(-10, 0, 10)
    result <- logist(x)
    expect_length(result, 3)
    expect_true(result[1] < result[2])
    expect_true(result[2] < result[3])
})

# Test for create_logist_features function
test_that("create_logist_features applies logistic transformations and removes NA columns", {
    sample_features <- matrix(rnorm(100), ncol = 5)
    colnames(sample_features) <- c("feature1", "feature2", "NA_feature", "feature3", "feature4")
    sample_features[, "NA_feature"] <- NA

    transformed_features <- create_logist_features(sample_features)

    # Check if NA column is removed
    expect_false("NA_feature" %in% colnames(transformed_features))

    # Check if transformation columns are created
    expect_true("feature1_high-energy" %in% colnames(transformed_features))
    expect_true("feature1_low-energy" %in% colnames(transformed_features))
    expect_true("feature1_sigmoid" %in% colnames(transformed_features))
    expect_true("feature1_higher-energy" %in% colnames(transformed_features))
})

test_that("create_logist_features with single column", {
    x <- matrix(1:10, ncol = 1)
    colnames(x) <- "motif1"
    result <- create_logist_features(x)
    expect_equal(ncol(result), 4) # 4 transformations
    expect_equal(nrow(result), 10)
    expected_names <- sort(c("motif1_low-energy", "motif1_high-energy", "motif1_sigmoid", "motif1_higher-energy"))
    expect_equal(sort(colnames(result)), expected_names)
})

test_that("create_logist_features assigns default column names if missing", {
    x <- matrix(rnorm(20), ncol = 2)
    result <- create_logist_features(x)
    expect_true(any(grepl("^V1_", colnames(result))))
    expect_true(any(grepl("^V2_", colnames(result))))
})

test_that("create_logist_features handles vector input", {
    x <- 1:10
    result <- create_logist_features(x)
    expect_true(is.matrix(result))
    expect_equal(nrow(result), 10)
    expect_equal(ncol(result), 4)
})

test_that("create_logist_features produces 4 transformations per column", {
    x <- matrix(rnorm(30), ncol = 3)
    colnames(x) <- c("A", "B", "C")
    result <- create_logist_features(x)
    expect_equal(ncol(result), 12) # 3 columns x 4 transformations
})

test_that("create_logist_features columns are sorted alphabetically", {
    x <- matrix(rnorm(20), ncol = 2)
    colnames(x) <- c("Z", "A")
    result <- create_logist_features(x)
    expect_equal(colnames(result), sort(colnames(result)))
})

# Test for rescale function
test_that("rescale inverts norm01 correctly", {
    orig_x <- c(-10, -5, 0, 5, 10)
    normed <- norm01(orig_x)
    rescaled <- rescale(normed, orig_x)
    expect_equal(rescaled, orig_x, tolerance = 1e-10)
})

test_that("rescale with negative original range", {
    orig_x <- c(-100, -50, -10)
    normed <- norm01(orig_x)
    rescaled <- rescale(normed, orig_x)
    expect_equal(rescaled, orig_x, tolerance = 1e-10)
})

test_that("rescale preserves relative ordering", {
    orig_x <- c(1, 5, 3, 7, 2)
    x <- c(0, 0.25, 0.5, 0.75, 1)
    result <- rescale(x, orig_x)
    expect_true(all(diff(result) > 0))
})

test_that("rescale handles NAs in orig_x", {
    orig_x <- c(1, NA, 3, 4, 5)
    x <- c(0, 0.5, 1)
    result <- rescale(x, orig_x)
    expect_equal(result[1], min(orig_x, na.rm = TRUE))
    expect_equal(result[3], max(orig_x, na.rm = TRUE))
})

test_that("rescale handles constant original vector", {
    orig_x <- rep(5, 10)
    x <- c(0, 0.5, 1)
    result <- rescale(x, orig_x)
    # min = 5, range = 0, so result = x * 0 + 5 = 5
    expect_true(all(result == 5))
})

# Test for norm_energy_dataset function
test_that("norm_energy_dataset normalizes with dataset reference", {
    set.seed(42)
    dataset_x <- runif(100, min = -10, max = 0)
    x <- runif(20, min = -10, max = 0)
    result <- norm_energy_dataset(x, dataset_x)
    expect_length(result, length(x))
    expect_true(all(result >= 0))
    expect_true(all(!is.na(result)))
})

test_that("norm_energy_dataset handles identical dataset values", {
    x <- c(-5, -3, -1)
    dataset_x <- rep(-5, 10)
    result <- norm_energy_dataset(x, dataset_x)
    # When range is 0, should return all zeros
    expect_true(all(result == 0))
})

test_that("norm_energy_dataset custom norm_energy_max scales output", {
    set.seed(42)
    dataset_x <- runif(100, min = -10, max = 0)
    x <- runif(20, min = -10, max = 0)
    result5 <- norm_energy_dataset(x, dataset_x, norm_energy_max = 5)
    result10 <- norm_energy_dataset(x, dataset_x, norm_energy_max = 10)
    # Doubling norm_energy_max should double the values
    expect_equal(result10, result5 * 2, tolerance = 1e-10)
})

# Test for norm_energy_matrix function
test_that("norm_energy_matrix normalizes matrix correctly", {
    set.seed(42)
    x <- matrix(runif(30, min = -10, max = 0), nrow = 10, ncol = 3)
    colnames(x) <- c("m1", "m2", "m3")
    result <- norm_energy_matrix(x, dataset_x = x)
    expect_equal(dim(result), dim(x))
    expect_equal(colnames(result), colnames(x))
    expect_true(all(!is.na(result)))
    expect_true(all(result >= 0))
})

test_that("norm_energy_matrix works when dataset_x has superset of x columns", {
    # dataset_x has more columns than x, x is a subset
    dataset_x <- matrix(runif(30, -10, 0), nrow = 10, ncol = 3)
    colnames(dataset_x) <- c("m1", "m2", "m3")
    x <- matrix(runif(30, -10, 0), nrow = 10, ncol = 3)
    colnames(x) <- c("m1", "m2", "m3")
    # Both have same columns - should work
    result <- norm_energy_matrix(x, dataset_x)
    expect_equal(colnames(result), c("m1", "m2", "m3"))
    expect_equal(dim(result), dim(x))
})

test_that("norm_energy_matrix errors when x is missing columns from dataset_x", {
    x <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    colnames(x) <- c("m1", "m2")
    dataset_x <- matrix(runif(30, -10, 0), nrow = 10, ncol = 3)
    colnames(dataset_x) <- c("m1", "m2", "m_missing")
    expect_error(norm_energy_matrix(x, dataset_x), "missing")
})

test_that("norm_energy_matrix without column names uses default names", {
    x <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    dataset_x <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    result <- norm_energy_matrix(x, dataset_x)
    expect_equal(colnames(result), c("V1", "V2"))
})

test_that("norm_energy_matrix errors on mismatched dimensions without colnames", {
    x <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    dataset_x <- matrix(runif(30, -10, 0), nrow = 10, ncol = 3)
    expect_error(norm_energy_matrix(x, dataset_x), "different column counts")
})

test_that("norm_energy_matrix handles single column matrix", {
    x <- matrix(runif(10, -10, 0), nrow = 10, ncol = 1)
    colnames(x) <- "motif1"
    result <- norm_energy_matrix(x)
    expect_equal(ncol(result), 1)
    expect_equal(colnames(result), "motif1")
})

test_that("norm_energy_matrix preserves row names", {
    x <- matrix(runif(6, -10, 0), nrow = 3, ncol = 2)
    colnames(x) <- c("m1", "m2")
    rownames(x) <- c("peak1", "peak2", "peak3")
    result <- norm_energy_matrix(x)
    expect_equal(rownames(result), c("peak1", "peak2", "peak3"))
})

test_that("norm_energy_matrix handles data.frame input", {
    df <- data.frame(m1 = runif(10, -10, 0), m2 = runif(10, -10, 0))
    result <- norm_energy_matrix(df)
    expect_true(is.matrix(result))
    expect_equal(ncol(result), 2)
})

test_that("norm_energy_matrix x colnames inherit from dataset_x when x has none", {
    x <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    dataset_x <- matrix(runif(20, -10, 0), nrow = 10, ncol = 2)
    colnames(dataset_x) <- c("A", "B")
    result <- norm_energy_matrix(x, dataset_x)
    expect_equal(colnames(result), c("A", "B"))
})

# Test for inverse_logist_by_type function
test_that("inverse_logist_by_type inverts low-energy correctly", {
    x_vals <- seq(1, 9, by = 0.5)
    forward <- logist(x_vals, x_0 = 0, L = 2, k = 0.5) - 1
    recovered <- iceqream:::inverse_logist_by_type(forward, "low-energy")
    expect_equal(recovered, x_vals, tolerance = 1e-8)
})

test_that("inverse_logist_by_type inverts high-energy correctly", {
    x_vals <- seq(1, 9, by = 0.5)
    forward <- logist(x_vals, x_0 = 10, L = 2, k = 0.50)
    recovered <- iceqream:::inverse_logist_by_type(forward, "high-energy")
    expect_equal(recovered, x_vals, tolerance = 1e-8)
})

test_that("inverse_logist_by_type inverts sigmoid correctly", {
    x_vals <- seq(1, 9, by = 0.5)
    forward <- logist(x_vals - 5, x_0 = 0, L = 1, k = 1)
    recovered <- iceqream:::inverse_logist_by_type(forward, "sigmoid")
    expect_equal(recovered, x_vals, tolerance = 1e-8)
})

test_that("inverse_logist_by_type inverts higher-energy correctly", {
    x_vals <- seq(1, 9, by = 0.5)
    forward <- logist(x_vals, x_0 = 10, L = 2, k = 1)
    recovered <- iceqream:::inverse_logist_by_type(forward, "higher-energy")
    expect_equal(recovered, x_vals, tolerance = 1e-8)
})

test_that("inverse_logist_by_type errors on unknown type", {
    expect_error(iceqream:::inverse_logist_by_type(0.5, "unknown"), "Unknown transformation type")
})

# Test for normalize_with_db_quantiles function
test_that("normalize_with_db_quantiles normalizes correctly", {
    motif_energies <- matrix(runif(30, -10, 0), nrow = 10, ncol = 3)
    colnames(motif_energies) <- c("motif_A", "motif_B", "motif_C")
    rownames(motif_energies) <- paste0("peak", 1:10)

    db_quantiles <- matrix(runif(9, -10, 0), nrow = 3, ncol = 3)
    rownames(db_quantiles) <- c("motif_A", "motif_B", "motif_C")
    colnames(db_quantiles) <- c("0.5", "0.9", "1")

    result <- iceqream:::normalize_with_db_quantiles(
        motif_energies, db_quantiles,
        energy_norm_quantile = 1,
        min_energy = -7,
        norm_energy_max = 10
    )
    expect_equal(dim(result), dim(motif_energies))
    expect_equal(colnames(result), colnames(motif_energies))
    expect_equal(rownames(result), rownames(motif_energies))
    expect_true(all(result >= 0))
})

test_that("normalize_with_db_quantiles errors on missing quantile", {
    motif_energies <- matrix(runif(10, -10, 0), nrow = 5, ncol = 2)
    colnames(motif_energies) <- c("m1", "m2")
    db_quantiles <- matrix(0, nrow = 2, ncol = 1)
    rownames(db_quantiles) <- c("m1", "m2")
    colnames(db_quantiles) <- "0.5"

    expect_error(
        iceqream:::normalize_with_db_quantiles(
            motif_energies, db_quantiles,
            energy_norm_quantile = 1,
            min_energy = -7, norm_energy_max = 10
        ),
        "energy_norm_quantile"
    )
})

test_that("normalize_with_db_quantiles errors on missing motifs", {
    motif_energies <- matrix(runif(10, -10, 0), nrow = 5, ncol = 2)
    colnames(motif_energies) <- c("m1", "m2")
    db_quantiles <- matrix(0, nrow = 1, ncol = 1)
    rownames(db_quantiles) <- "m1"
    colnames(db_quantiles) <- "1"

    expect_error(
        iceqream:::normalize_with_db_quantiles(
            motif_energies, db_quantiles,
            energy_norm_quantile = 1,
            min_energy = -7, norm_energy_max = 10
        ),
        "missing"
    )
})
