# Test for norm01 function
test_that("norm01 normalizes values to the [0, 1] range", {
    x <- rnorm(100, 10, 5)
    normed_x <- norm01(x)
    expect_true(all(normed_x >= 0))
    expect_true(all(normed_x <= 1))
    expect_equal(min(normed_x), 0)
    expect_equal(max(normed_x), 1)
})

# Test for norm_energy function
test_that("norm_energy normalizes energy values correctly", {
    x <- runif(n = 100, min = -11, max = 0)
    normed_energy <- norm_energy(x)
    expect_true(all(!is.na(normed_energy)))
    expect_true(all(normed_energy <= 10))
    expect_true(all(normed_energy >= 0))
})

# Test for logist function
test_that("logist calculates the logistic function correctly", {
    x_vals <- seq(-10, 10, by = 0.1)
    results <- logist(x_vals)
    expect_true(all(!is.na(results)))
    expect_true(all(results <= 1))
    expect_true(all(results >= 0))
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
