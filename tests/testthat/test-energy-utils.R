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

# # Test for norm_energy_dataset function
# test_that("norm_energy_dataset normalizes energy values using an external dataset", {
#     x <- runif(n = 50, min = -11, max = 0)
#     ds_x <- runif(n = 50, min = -11, max = 0)
#     y <- norm_energy_dataset(x, ds_x)
#     expect_equal(length(y), length(x))
#     expect_true(all(!is.na(y)))
#     expect_true(all(y <= 10))
#     expect_true(all(y >= 0))
# })

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
    expect_true("feature1_early" %in% colnames(transformed_features))
    expect_true("feature1_linear" %in% colnames(transformed_features))
    expect_true("feature1_late" %in% colnames(transformed_features))

    # Check if transformations are correctly bounded
    early_cols <- grep("_early$", colnames(transformed_features), value = TRUE)
    linear_cols <- grep("_linear$", colnames(transformed_features), value = TRUE)
    late_cols <- grep("_late$", colnames(transformed_features), value = TRUE)

    expect_true(all(transformed_features[, early_cols] <= 1))
    expect_true(all(transformed_features[, early_cols] >= -1))
    expect_true(all(transformed_features[, linear_cols] <= 1))
    expect_true(all(transformed_features[, linear_cols] >= -1))
    expect_true(all(transformed_features[, late_cols] <= 1))
    expect_true(all(transformed_features[, late_cols] >= -1))
})
